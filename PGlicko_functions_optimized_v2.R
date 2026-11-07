# ============================================================
# CONSTANTS (computed once globally)
# ============================================================
.Q  <- log(10) / 400
.Q2 <- .Q^2
.K  <- 3 * .Q2 / pi^2

# ============================================================
# CACHE ENVIRONMENT
# ============================================================
.cache <- new.env(hash = TRUE, parent = emptyenv())

col_rename <- function(x) {
  setNames(x, c("Season", "Player1", "Player2", "Game_Result"))
}

periods <- function(x) {
  sort(unique(x$Season))
}

players <- function(x) {
  sort(unique(c(x$Player1, x$Player2)))
}

build_cache <- function(dataset) {
  
  by_season <- split(dataset, dataset$Season)
  .cache$by_season <- by_season
  
  # Per season: split by player
  .cache$by_season_player <- lapply(by_season, function(d) {
    all_p <- unique(c(d$Player1, d$Player2))
    setNames(
      lapply(all_p, function(pl)
        d[d$Player1 == pl | d$Player2 == pl, , drop = FALSE]),
      all_p
    )
  })
  
  # Per season: split by pair key
  .cache$by_season_pair <- lapply(by_season, function(d) {
    keys <- paste(d$Player1, d$Player2, sep = "_")
    split(d, keys)
  })
  
  # Pre-compute pair keys for every row
  .cache$pair_keys <- paste(dataset$Player1, dataset$Player2, sep = "_")
  
  # FIX 1: cache stable values — periods, players, pair list
  # so get_like never recomputes these across ~175 optim() calls
  .cache$p_x      <- periods(dataset)
  .cache$Players  <- players(dataset)
  .cache$Pairwise_List_Vector <- Pairwise_initialize(.cache$Players)[[2]]
  
  # FIX 5: pre-build colname lookup set for winning_prob_vectorized
  # avoids repeated %in% colnames(PS_ts) inside the row loop
  .cache$ps_col_set <- NULL  # populated after first get_like builds PS_ts colnames
  # (set lazily on first call since PS_ts colnames = Pairwise_List_Vector)
  .cache$ps_col_idx <- setNames(
    seq_along(.cache$Pairwise_List_Vector),
    .cache$Pairwise_List_Vector
  )
  
  invisible(NULL)
}

# ============================================================
# FAST LOOKUP REPLACEMENTS
# ============================================================
get_player_games <- function(p, name) {
  sp <- .cache$by_season_player[[as.character(p)]]
  if (is.null(sp)) return(NULL)
  sp[[name]]
}

get_pair_games <- function(p, name1, name2) {
  .cache$by_season_pair[[as.character(p)]][[paste(name1, name2, sep = "_")]]
}

# ============================================================
# MATH FUNCTIONS
# ============================================================
q_function <- function(k = 400) log(10) / k

g_function <- function(A, k = 400) {
  q <- if (k == 400) .Q else log(10) / k
  1 / sqrt(1 + 3 * q^2 * A^2 / pi^2)
}

Gg_function <- function(A, B, k = 400) {
  q <- if (k == 400) .Q else log(10) / k
  1 / sqrt(1 + 3 * q^2 * A^2 / pi^2 + 3 * q^2 * B^2 / pi^2)
}

Es_ABCDE <- function(A, B, C, D, E, k = 400) {
  1 / (1 + 10^(-Gg_function(C, E, k) * (A + D - B) / k))
}

# ============================================================
# DELTA FUNCTIONS
# ============================================================
delta_player_rating <- function(p, name, x, rating_tbl, pairwise_r_tbl, pairwise_sd_tbl, k = 400) {
  
  this_period <- get_player_games(p, name)
  if (is.null(this_period) || nrow(this_period) == 0) return(NA_real_)
  
  opponents <- sort(unique(c(this_period$Player1, this_period$Player2)))
  opponents <- opponents[opponents != name]
  if (length(opponents) == 0) return(NA_real_)
  
  opp_ratings <- rating_tbl$Rating[opponents]
  opp_sigmas  <- rating_tbl$Sigma[opponents]
  pair_adv    <- pairwise_r_tbl[name, opponents]
  pair_adv_sd <- pairwise_sd_tbl[name, opponents]
  player_rating <- rating_tbl$Rating[name]
  
  all_players <- c(this_period$Player1, this_period$Player2)
  game_cnt    <- tabulate(match(all_players, opponents), nbins = length(opponents))
  
  G_vals  <- Gg_function(opp_sigmas, pair_adv_sd, k)
  Es_vals <- Es_ABCDE(player_rating, opp_ratings, opp_sigmas, pair_adv, pair_adv_sd, k)
  
  sqrt((.Q2 * sum(game_cnt * G_vals^2 * Es_vals * (1 - Es_vals)))^-1)
}

delta_pairwise_strength <- function(p, name1, name2, x, rating_tbl, pairwise_r_tbl, pairwise_sd_tbl, k = 400) {
  
  this_period <- get_pair_games(p, name1, name2)
  if (is.null(this_period) || nrow(this_period) == 0) stop("No games found for pair.")
  
  opp          <- name2
  player_sigma <- rating_tbl$Sigma[name1]
  opp_rating   <- rating_tbl$Rating[opp]
  opp_sigma    <- rating_tbl$Sigma[opp]
  pair_adv     <- pairwise_r_tbl[name1, opp]
  pair_adv_sd  <- pairwise_sd_tbl[name1, opp]
  game_cnt     <- nrow(this_period)
  
  G_val  <- Gg_function(player_sigma, opp_sigma, k)
  Es_val <- Es_ABCDE(pair_adv, opp_rating, opp_sigma, rating_tbl$Rating[name1], player_sigma, k)
  
  sqrt((.Q2 * game_cnt * G_val^2 * Es_val * (1 - Es_val))^-1)
}

# ============================================================
# UPDATE FUNCTIONS
# ============================================================
sigma_update <- function(A, B) sqrt((1/A^2 + 1/B^2)^-1)

player_rating_update <- function(p, name, x, rating_tbl, pairwise_r_tbl, pairwise_sd_tbl, new_sigma, k = 400) {
  
  player_rating <- rating_tbl$Rating[name]
  record        <- get_player_games(p, name)
  if (is.null(record) || nrow(record) == 0) return(player_rating)
  
  opp         <- ifelse(record$Player1 == name, record$Player2, record$Player1)
  opp_rating  <- rating_tbl$Rating[opp]
  opp_sigma   <- rating_tbl$Sigma[opp]
  pair_adv    <- pairwise_r_tbl[name, opp]
  pair_adv_sd <- pairwise_sd_tbl[name, opp]
  result      <- ifelse(record$Player2 == opp, record$Game_Result, 1 - record$Game_Result)
  
  Es_vals <- Es_ABCDE(player_rating, opp_rating, opp_sigma, pair_adv, pair_adv_sd, k)
  G_vals  <- Gg_function(opp_sigma, pair_adv_sd, k)
  
  player_rating + .Q * new_sigma^2 * sum(G_vals * (result - Es_vals))
}

pairwise_strength_update <- function(p, name1, name2, x, rating_tbl, pairwise_r_tbl, pairwise_sd_tbl, new_sigma, k = 400) {
  
  record <- get_pair_games(p, name1, name2)
  if (is.null(record) || nrow(record) == 0) return(pairwise_r_tbl[name1, name2])
  
  opp         <- name2
  opp_rating  <- rating_tbl$Rating[opp]
  opp_sigma   <- rating_tbl$Sigma[opp]
  pair_adv    <- pairwise_r_tbl[name1, opp]
  pair_adv_sd <- pairwise_sd_tbl[name1, opp]
  result      <- ifelse(record$Player2 == opp, record$Game_Result, 1 - record$Game_Result)
  
  Es_vals <- Es_ABCDE(pair_adv, opp_rating, opp_sigma, rating_tbl$Rating[name1], rating_tbl$Sigma[name1], k)
  G_vals  <- Gg_function(rating_tbl$Sigma[name1], opp_sigma, k)
  
  pairwise_r_tbl[name1, opp] + .Q * new_sigma^2 * sum(G_vals * (result - Es_vals))
}

# ============================================================
# INITIALIZERS
# ============================================================
Rating_Table_initialize <- function(Players, initial_rating, initial_sd) {
  n <- length(Players)
  list(
    Rating  = setNames(rep(initial_rating, n), Players),
    Sigma   = setNames(rep(initial_sd,     n), Players),
    Games   = setNames(rep(0L,             n), Players),
    Win     = setNames(rep(0L,             n), Players),
    Loss    = setNames(rep(0L,             n), Players),
    Lag     = setNames(rep(0L,             n), Players),
    debut   = setNames(rep("No",           n), Players),
    Players = Players
  )
}

rating_tbl_to_df <- function(rt) {
  data.frame(
    Players = rt$Players,
    Rating  = rt$Rating,
    Sigma   = rt$Sigma,
    Games   = rt$Games,
    Win     = rt$Win,
    Loss    = rt$Loss,
    Lag     = rt$Lag,
    debut   = rt$debut,
    stringsAsFactors = FALSE
  )
}

Pairwise_initialize <- function(Players) {
  list(0, combn(Players, 2, FUN = function(x) paste(x, collapse = "_")))
}

init_square_matrix <- function(Players, fill_value = 0) {
  n   <- length(Players)
  mat <- matrix(fill_value, nrow = n, ncol = n, dimnames = list(Players, Players))
  diag(mat) <- NA
  mat
}

Pairwise_Strength_initialize    <- function(Players)                init_square_matrix(Players, 0)
Pairwise_Strength_sd_initialize <- function(Players, init_pw_sigma) init_square_matrix(Players, init_pw_sigma)

# ============================================================
# WINNING PROBABILITY
# ============================================================
winning_prob <- function(mu_i, sig_i, mu_j, sig_j, pairwise_strength_ij, pairwise_sd_ij) {
  1 / (1 + 10^(-g_function(sqrt(sig_i^2 + sig_j^2 + pairwise_sd_ij^2)) *
                 (mu_i + pairwise_strength_ij - mu_j) / 400))
}

winning_prob_vectorized <- function(dataset, ratings_ts, sigma_ts, PS_ts, PS_sd_ts) {
  
  seasons  <- dataset$Season
  p1_names <- dataset$Player1
  p2_names <- dataset$Player2
  n        <- nrow(dataset)
  wp       <- numeric(n)
  
  early     <- seasons < 7
  wp[early] <- 0.5
  
  # FIX 5: pre-build named index once — O(1) lookup replaces %in% colnames() each row
  ps_col_idx <- .cache$ps_col_idx
  
  if (any(!early)) {
    idx <- which(!early)
    for (i in idx) {
      s   <- seasons[i];  p1 <- p1_names[i];  p2 <- p2_names[i]
      prev  <- (s - 6):(s - 1)
      mu_i  <- mean(ratings_ts[prev, p1])
      mu_j  <- mean(ratings_ts[prev, p2])
      sig_i <- sigma_ts[s - 1, p1]
      sig_j <- sigma_ts[s - 1, p2]
      
      # FIX 5: named index lookup instead of %in% colnames()
      col_nm  <- paste0(p1, "_", p2)
      flipped <- is.na(ps_col_idx[col_nm])
      if (flipped) col_nm <- paste0(p2, "_", p1)
      
      ps_val  <- PS_ts[s - 1, col_nm]
      psd_val <- PS_sd_ts[s - 1, col_nm]
      if (flipped) ps_val <- -ps_val
      ps_val  <- if (is.na(ps_val))  0  else ps_val
      psd_val <- if (is.na(psd_val)) 40 else psd_val
      
      wp[i] <- winning_prob(mu_i, sig_i, mu_j, sig_j, ps_val, psd_val)
    }
  }
  wp
}

prediction <- function(x) as.integer(x >= 0.5)

# ============================================================
# MAIN LIKELIHOOD FUNCTION
# ============================================================
get_like <- function(para, data, rdmax = 350, metric = "likelihood") {
  
  # Local alias — all internal references to 'dataset' now use this
  dataset <- data
  
  tic <- Sys.time()
  if (any(para <= 0)) return(-Inf)
  
  init_sigm    <- para[1]
  init_pw_sigm <- para[2]
  
  # FIX 1: pull stable values from cache — never recomputed across optim() calls
  p_x                  <- .cache$p_x
  Players              <- .cache$Players
  Pairwise_List_Vector <- .cache$Pairwise_List_Vector
  n_p                  <- length(p_x)
  n_pl                 <- length(Players)
  period_labels        <- paste0("p", p_x)
  
  # FIX 2: pre-compute period index map — O(1) lookup replaces which(p_x == p) each iteration
  p_idx <- setNames(seq_along(p_x), as.character(p_x))
  
  Ratings                     <- Rating_Table_initialize(Players, 1500, init_sigm)
  Pairwise_Strength_Matrix    <- Pairwise_Strength_initialize(Players)
  Pairwise_Strength_sd_Matrix <- Pairwise_Strength_sd_initialize(Players, init_pw_sigm)
  
  attendance <- matrix(NA_character_, nrow = n_p, ncol = n_pl,
                       dimnames = list(period_labels, Players))
  ratings_ts <- matrix(1500,          nrow = n_p, ncol = n_pl,
                       dimnames = list(period_labels, Players))
  sigma_ts   <- matrix(init_sigm,     nrow = n_p, ncol = n_pl,
                       dimnames = list(period_labels, Players))
  PS_ts      <- matrix(0,             nrow = n_p, ncol = length(Pairwise_List_Vector),
                       dimnames = list(period_labels, Pairwise_List_Vector))
  PS_sd_ts   <- matrix(init_pw_sigm,  nrow = n_p, ncol = length(Pairwise_List_Vector),
                       dimnames = list(period_labels, Pairwise_List_Vector))
  
  for (p in p_x) {
    
    period_idx  <- paste0("p", p)
    period_data <- .cache$by_season[[as.character(p)]]
    players_period <- sort(unique(c(period_data$Player1, period_data$Player2)))
    
    # Vectorized game stats
    p1_tab <- table(factor(period_data$Player1, levels = Players))
    p2_tab <- table(factor(period_data$Player2, levels = Players))
    Ratings$Games[Players]        <- Ratings$Games[Players] + as.integer(p1_tab) + as.integer(p2_tab)
    Ratings$Win[Players]          <- Ratings$Win[Players]   + as.integer(p1_tab)
    Ratings$Loss[Players]         <- Ratings$Loss[Players]  + as.integer(p2_tab)
    Ratings$Lag[players_period]   <- 0L
    Ratings$debut[players_period] <- "Yes"
    
    debuted      <- Players[Ratings$debut == "Yes"]
    idle_players <- setdiff(debuted, players_period)
    
    if (length(idle_players) > 0) {
      attendance[period_idx, idle_players] <- "No"
      ratings_ts[period_idx, idle_players] <- Ratings$Rating[idle_players]
      sigma_ts[period_idx,   idle_players] <- Ratings$Sigma[idle_players]
      Ratings$Lag[idle_players] <- Ratings$Lag[idle_players] + 1L
    }
    
    old_ratings <- Ratings
    new_ratings <- Ratings
    old_pw      <- Pairwise_Strength_Matrix
    old_pw_sd   <- Pairwise_Strength_sd_Matrix
    new_pw      <- old_pw
    new_pw_sd   <- old_pw_sd
    
    # Player rating updates
    for (i in players_period) {
      delta_i     <- delta_player_rating(p, i, dataset, old_ratings, old_pw, old_pw_sd)
      new_sigma_i <- min(rdmax, sigma_update(old_ratings$Sigma[i], delta_i))
      new_ratings$Sigma[i]  <- new_sigma_i
      new_ratings$Rating[i] <- player_rating_update(p, i, dataset, old_ratings, old_pw, old_pw_sd,
                                                    new_sigma = new_sigma_i)
      attendance[period_idx, i] <- "Yes"
    }
    
    # FIX 3: parse pairs once as matrix — no strsplit, no repeated apply()
    pairs_mat <- unique(as.matrix(period_data[, c("Player1", "Player2")]))
    
    # Carry forward previous period pairwise timeseries BEFORE the pairs loop
    # (must be outside — moving it inside the loop would overwrite per-pair writes)
    cur_i  <- p_idx[as.character(p)]   # FIX 2: O(1) lookup
    prev_i <- max(1L, cur_i - 1L)
    PS_ts[period_idx, ]    <- PS_ts[prev_i, ]
    PS_sd_ts[period_idx, ] <- PS_sd_ts[prev_i, ]
    
    # FIX 4: single pass over pairs for both pairwise updates AND timeseries assignment
    for (pi in seq_len(nrow(pairs_mat))) {
      p1 <- pairs_mat[pi, 1];  p2 <- pairs_mat[pi, 2]
      
      # Pairwise update
      new_pw_sd[p1, p2] <- min(rdmax, sigma_update(old_pw_sd[p1, p2],
                                                   delta_pairwise_strength(p, p1, p2, dataset, old_ratings, old_pw, old_pw_sd)))
      new_pw[p1, p2]    <- pairwise_strength_update(p, p1, p2, dataset, old_ratings, old_pw, old_pw_sd,
                                                    new_sigma = new_pw_sd[p1, p2])
      new_pw_sd[p2, p1] <- new_pw_sd[p1, p2]
      new_pw[p2, p1]    <- -new_pw[p1, p2]
      
      # Timeseries assignment in same pass (FIX 4)
      col_nm <- if (p1 < p2) paste0(p1, "_", p2) else paste0(p2, "_", p1)
      if (p1 < p2) {
        PS_ts[period_idx,    col_nm] <- new_pw[p1, p2]
        PS_sd_ts[period_idx, col_nm] <- new_pw_sd[p1, p2]
      } else {
        PS_ts[period_idx,    col_nm] <- new_pw[p2, p1]
        PS_sd_ts[period_idx, col_nm] <- new_pw_sd[p2, p1]
      }
      Pairwise_Strength_Matrix[p1, p2]    <- new_pw[p1, p2]
      Pairwise_Strength_sd_Matrix[p1, p2] <- new_pw_sd[p1, p2]
    }
    
    # Batch commit active players
    Ratings$Rating[players_period] <- new_ratings$Rating[players_period]
    Ratings$Sigma[players_period]  <- new_ratings$Sigma[players_period]
    
    # Batch timeseries update
    all_updated <- c(players_period, idle_players)
    ratings_ts[period_idx, all_updated] <- Ratings$Rating[all_updated]
    sigma_ts[period_idx,   all_updated] <- Ratings$Sigma[all_updated]
  }
  
  # Vectorized winning probability
  wp          <- winning_prob_vectorized(dataset, ratings_ts, sigma_ts, PS_ts, PS_sd_ts)
  gr          <- dataset$Game_Result
  result_prob <- wp^gr * (1 - wp)^(1 - gr)
  
  n         <- length(gr)
  train_end <- floor(0.8 * n)
  test_idx  <- (train_end + 1L):n
  
  likelihood <- sum(log(result_prob[seq_len(train_end)]))
  roc_curve  <- roc(gr[test_idx], wp[test_idx])
  
  auc        <- as.numeric(roc_curve$auc)
  acc        <- mean((wp[test_idx] > 0.5) == gr[test_idx])
  
  toc <- Sys.time()
  print(paste(init_sigm, init_pw_sigm, "time(mins):", round(as.numeric(toc - tic, units = "mins"), 4),
              "(auc, acc, likelihood):", paste(round(c(auc, acc, likelihood), 4), collapse = ", ")))
  
  record <<- rbind(record, c(para[1], para[2], auc, acc, likelihood))
  
  switch(metric,
         all        = c(auc, acc, likelihood),
         auc        = auc,
         acc        = acc,
         likelihood = likelihood,
         stop("choose a metric"))
}