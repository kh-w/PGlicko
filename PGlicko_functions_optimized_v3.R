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
  
  # Per season: split by UNORDERED pair key, so all n_j games between two
  # players are pooled regardless of who is listed as Player1 (eqs 30, 32)
  .cache$by_season_pair <- lapply(by_season, function(d) {
    keys <- paste(pmin(d$Player1, d$Player2), pmax(d$Player1, d$Player2), sep = "_")
    split(d, keys)
  })

  # Remember which dataset this cache was built from (checked in get_like)
  .cache$dataset <- dataset
  
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
  key <- paste(pmin(name1, name2), pmax(name1, name2), sep = "_")
  .cache$by_season_pair[[as.character(p)]][[key]]
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

  # Fully vectorised: mean rating over the previous 6 periods, RD and pairwise
  # values at the previous period; first 6 periods fixed at 0.5.
  # Period position (not raw Season value) indexes the time-series rows.
  seasons <- match(paste0("p", dataset$Season), rownames(ratings_ts))
  wp      <- rep(0.5, nrow(dataset))
  ok      <- which(seasons >= 7)

  if (length(ok) > 0) {
    t  <- seasons[ok]
    p1 <- dataset$Player1[ok];  p2 <- dataset$Player2[ok]
    i1 <- match(p1, colnames(ratings_ts));  i2 <- match(p2, colnames(ratings_ts))

    cs   <- rbind(0, apply(ratings_ts, 2, cumsum))       # cs[k + 1, ] = sum of rows 1..k
    mu_i <- (cs[cbind(t, i1)] - cs[cbind(t - 6, i1)]) / 6
    mu_j <- (cs[cbind(t, i2)] - cs[cbind(t - 6, i2)]) / 6

    col     <- match(paste0(p1, "_", p2), colnames(PS_ts))
    flipped <- is.na(col)
    col[flipped] <- match(paste0(p2, "_", p1)[flipped], colnames(PS_ts))
    ps_val  <- PS_ts[cbind(t - 1, col)]
    psd_val <- PS_sd_ts[cbind(t - 1, col)]
    ps_val[flipped]       <- -ps_val[flipped]
    ps_val[is.na(ps_val)]   <- 0
    psd_val[is.na(psd_val)] <- 40

    wp[ok] <- winning_prob(mu_i, sigma_ts[cbind(t - 1, i1)], mu_j, sigma_ts[cbind(t - 1, i2)],
                           ps_val, psd_val)
  }
  wp
}

prediction <- function(x) as.integer(x >= 0.5)

# ============================================================
# VECTORISED ENGINE
# Runs P-Glicko over every period of `dataset`. Within a period all games are
# processed at once with integer player indices and plain matrices:
#   step 1  RD growth for debuted players: sigma <- min(sqrt(sigma^2 + c^2), rdmax)
#   player updates   eqs (16)-(21), simultaneous, using pre-period values
#   pairwise updates eqs (28)-(32), all games of an unordered pair pooled,
#                    beta_ji = -beta_ij, b_ji = b_ij
# Predictions (eq 34) use the mean rating over the previous 6 periods and the
# RD / pairwise values at the end of the previous period.
# ============================================================
.pglicko_engine <- function(para, dataset, rdmax = 350, c_rd = 63.2) {

  init_sigm    <- para[1]
  init_pw_sigm <- para[2]

  p_x     <- periods(dataset)
  Players <- players(dataset)
  np      <- length(Players)
  TT      <- length(p_x)
  period_labels <- paste0("p", p_x)

  a_all <- match(dataset$Player1, Players)
  b_all <- match(dataset$Player2, Players)
  s_all <- dataset$Game_Result
  t_all <- match(dataset$Season, p_x)
  rows  <- split(seq_along(t_all), factor(t_all, levels = seq_len(TT)))

  mu    <- rep(1500, np)
  sg    <- rep(init_sigm, np)
  debut <- logical(np)
  games <- win <- loss <- lag <- integer(np)
  B     <- matrix(0, np, np)              # pairwise advantage beta[i, j]
  Bs    <- matrix(init_pw_sigm, np, np)   # pairwise deviation b[i, j]

  ratings_ts <- matrix(1500,      TT, np, dimnames = list(period_labels, Players))
  sigma_ts   <- matrix(init_sigm, TT, np, dimnames = list(period_labels, Players))
  attendance <- matrix(NA_character_, TT, np, dimnames = list(period_labels, Players))

  # Per-game values at the end of the previous period (for prediction)
  g_beta <- g_b <- numeric(length(s_all))
  q <- .Q; K <- .K

  for (t in seq_len(TT)) {
    r  <- rows[[t]]
    a  <- a_all[r];  b <- b_all[r];  s <- s_all[r]
    ab <- cbind(a, b)
    bet <- B[ab];  bsd <- Bs[ab]
    g_beta[r] <- bet;  g_b[r] <- bsd

    # Step 1: RD growth for players who have already debuted
    sg[debut] <- pmin(sqrt(sg[debut]^2 + c_rd^2), rdmax)
    mu0 <- mu;  sg0 <- sg

    # Bookkeeping (Game_Result = 1 means Player1 won)
    who   <- c(a, b)
    act   <- sort(unique(who))
    games <- games + tabulate(who, np)
    win   <- win   + tabulate(c(a[s == 1], b[s == 0]), np)
    loss  <- loss  + tabulate(c(a[s == 0], b[s == 1]), np)
    debut[act] <- TRUE
    idle  <- setdiff(which(debut), act)
    lag[idle] <- lag[idle] + 1L
    lag[act]  <- 0L
    attendance[t, idle] <- "No"
    attendance[t, act]  <- "Yes"

    # Player updates, eqs (16)-(21): Player1's side and Player2's side of every game
    G1 <- 1 / sqrt(1 + K * (sg0[b]^2 + bsd^2))
    E1 <- 1 / (1 + 10^(-G1 * (mu0[a] + bet - mu0[b]) / 400))
    G2 <- 1 / sqrt(1 + K * (sg0[a]^2 + bsd^2))
    E2 <- 1 / (1 + 10^(-G2 * (mu0[b] - bet - mu0[a]) / 400))
    v  <- rowsum(cbind(c(G1^2 * E1 * (1 - E1), G2^2 * E2 * (1 - E2)),   # 1/(q^2 delta^2)
                       c(G1 * (s - E1),        G2 * (1 - s - E2))),      # score residuals
                 who)                                                    # rows in order of `act`
    sg[act] <- pmin(rdmax, sqrt(1 / (1 / sg0[act]^2 + q^2 * v[, 1])))
    mu[act] <- mu0[act] + q * sg[act]^2 * v[, 2]

    # Pairwise updates, eqs (28)-(32): unordered pairs, all games pooled
    lo  <- pmin(a, b);  hi <- pmax(a, b)
    s_lo <- ifelse(a == lo, s, 1 - s)                     # result from `lo`'s side
    key <- lo * (np + 1) + hi
    uk  <- sort(unique(key))
    w   <- rowsum(cbind(1, s_lo), match(key, uk))         # n_j and sum of scores per pair
    pl  <- uk %/% (np + 1);  ph <- uk %% (np + 1)
    lh  <- cbind(pl, ph);    hl <- cbind(ph, pl)
    G   <- 1 / sqrt(1 + K * (sg0[pl]^2 + sg0[ph]^2))
    E   <- 1 / (1 + 10^(-G * (mu0[pl] + B[lh] - mu0[ph]) / 400))
    b_new    <- pmin(rdmax, sqrt(1 / (1 / Bs[lh]^2 + q^2 * w[, 1] * G^2 * E * (1 - E))))
    beta_new <- B[lh] + q * b_new^2 * G * (w[, 2] - w[, 1] * E)
    B[lh]  <- beta_new;  B[hl]  <- -beta_new
    Bs[lh] <- b_new;     Bs[hl] <- b_new

    ratings_ts[t, ] <- mu
    sigma_ts[t, ]   <- sg
  }

  # Predictions, eq (34)
  wp <- rep(0.5, length(s_all))
  ok <- which(t_all >= 7)
  if (length(ok) > 0) {
    tt <- t_all[ok];  i1 <- a_all[ok];  i2 <- b_all[ok]
    cs <- rbind(0, apply(ratings_ts, 2, cumsum))
    mu_i <- (cs[cbind(tt, i1)] - cs[cbind(tt - 6, i1)]) / 6
    mu_j <- (cs[cbind(tt, i2)] - cs[cbind(tt - 6, i2)]) / 6
    wp[ok] <- winning_prob(mu_i, sigma_ts[cbind(tt - 1, i1)], mu_j, sigma_ts[cbind(tt - 1, i2)],
                           g_beta[ok], g_b[ok])
  }

  diag(B) <- NA;  diag(Bs) <- NA
  dimnames(B) <- dimnames(Bs) <- list(Players, Players)
  nm <- function(x) setNames(x, Players)
  list(Ratings = list(Rating = nm(mu), Sigma = nm(sg), Games = nm(games), Win = nm(win),
                      Loss = nm(loss), Lag = nm(lag), debut = nm(ifelse(debut, "Yes", "No")),
                      Players = Players),
       Pairwise_Strength_Matrix    = B,
       Pairwise_Strength_sd_Matrix = Bs,
       ratings_ts = ratings_ts,
       sigma_ts   = sigma_ts,
       attendance = attendance,
       wp         = wp)
}

# ============================================================
# MAIN LIKELIHOOD FUNCTION
# ============================================================
get_like <- function(para, data, rdmax = 350, metric = "likelihood",
                     c_rd = 63.2, train_frac = 0.8) {

  dataset <- data

  tic <- Sys.time()
  if (any(para <= 0)) return(-Inf)
  if (!(train_frac > 0 && train_frac <= 1)) stop("train_frac must be in (0, 1]")

  init_sigm    <- para[1]
  init_pw_sigm <- para[2]

  # Whole rating run and predictions (no dependence on build_cache())
  wp          <- .pglicko_engine(para, dataset, rdmax = rdmax, c_rd = c_rd)$wp
  gr          <- dataset$Game_Result
  result_prob <- wp^gr * (1 - wp)^(1 - gr)

  # Fit on the first train_frac of games, evaluate on the rest
  n         <- length(gr)
  train_end <- floor(train_frac * n)

  likelihood <- sum(log(result_prob[seq_len(train_end)]))
  if (train_end < n) {
    test_idx <- (train_end + 1L):n
    auc <- as.numeric(pROC::roc(gr[test_idx], wp[test_idx], quiet = TRUE)$auc)
    acc <- mean((wp[test_idx] > 0.5) == gr[test_idx])
  } else {
    auc <- NA_real_
    acc <- NA_real_
  }

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
