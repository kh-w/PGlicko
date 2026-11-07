library(compiler)

# ---------------------------------------------------------------------------
# Constants (computed once at load time instead of recomputed every call)
# ---------------------------------------------------------------------------
.q  <- log(10) / 400          # q constant used throughout Glicko
.q2 <- .q^2
.k  <- 3 * .q2 / pi^2         # 3*q^2/pi^2 used in g_function

# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------
col_rename <- function(x) {
  setNames(x, c("Season", "Player1", "Player2", "Game_Result"))
}

periods <- function(x) {
  sort(unique(x$Season))
}

players <- function(x) {
  sort(unique(c(x$Player1, x$Player2)))
}

# Games played by a player (=name) in a period (=p)
games_in_p <- function(p, name, x) {
  x[x$Season == p & (x$Player1 == name | x$Player2 == name), ]
}

# ---------------------------------------------------------------------------
# Core Glicko math
# ---------------------------------------------------------------------------
g_function <- function(sigma) {
  1 / sqrt(1 + .k * sigma^2)
}

Es_uusj2 <- function(u, uj, sj) {
  1 / (1 + 10^(-(g_function(sj) * (u - uj) / 400)))
}

delta <- cmpfun(function(p, name, x, rating_tbl) {
  this_period <- .games_in_p_cached(p, name, x)
  
  opponents <- sort(unique(c(this_period$Player1, this_period$Player2)))
  opponents <- opponents[opponents != name]
  
  opp_rows     <- match(opponents, rating_tbl$Players)
  sigma_vals   <- rating_tbl$Sigma[opp_rows]
  rating_vals  <- rating_tbl$Rating[opp_rows]
  player_rating <- rating_tbl$Rating[rating_tbl$Players == name]
  
  # Count games per opponent with a single tabulation pass
  all_opps    <- c(this_period$Player1, this_period$Player2)
  game_counts <- tabulate(match(all_opps[all_opps %in% opponents], opponents),
                          nbins = length(opponents))
  
  g_vals  <- g_function(sigma_vals)
  Es_vals <- Es_uusj2(player_rating, rating_vals, sigma_vals)
  
  sqrt(1 / (.q2 * sum(game_counts * g_vals^2 * Es_vals * (1 - Es_vals))))
})

player_sigma_update <- function(sigma, delta) {
  sqrt((1 / sigma^2 + 1 / delta^2)^-1)
}

player_rating_update <- cmpfun(function(p, name, x, Ratings_tbl, new_sigma) {
  player_idx    <- match(name, Ratings_tbl$Players)
  player_rating <- Ratings_tbl$Rating[player_idx]
  
  record <- .games_in_p_cached(p, name, x)
  
  opp        <- ifelse(record$Player1 == name, record$Player2, record$Player1)
  opp_idx    <- match(opp, Ratings_tbl$Players)
  opp_rating <- Ratings_tbl$Rating[opp_idx]
  opp_sigma  <- Ratings_tbl$Sigma[opp_idx]
  result     <- ifelse(record$Player2 == opp, record$Game_Result, 1 - record$Game_Result)
  
  g_vals  <- g_function(opp_sigma)
  Es_vals <- Es_uusj2(player_rating, opp_rating, opp_sigma)
  
  player_rating + .q * new_sigma^2 * sum(g_vals * (result - Es_vals))
})

Rating_Table_initialize <- function(Players, initial_rating, initial_sd) {
  data.frame(
    Players = Players,
    Rating  = initial_rating,
    Sigma   = initial_sd,
    Games   = 0L,
    Win     = 0L,
    Loss    = 0L,
    Lag     = 0L,
    debut   = "No",
    stringsAsFactors = FALSE
  )
}

winning_prob <- function(mu_i, sig_i, mu_j, sig_j) {
  1 / (1 + 10^(-g_function(sqrt(sig_i^2 + sig_j^2)) * (mu_i - mu_j) / 400))
}

get_winning_prob <- cmpfun(function(row, ratings_ts, sigma_ts) {
  season <- as.integer(row[1])
  if (season < 7L) return(0.5)
  
  p1 <- row[2]
  p2 <- row[3]
  idx <- (season - 6L):(season - 1L)
  
  mu_i  <- mean(ratings_ts[idx, p1], na.rm = TRUE)
  mu_j  <- mean(ratings_ts[idx, p2], na.rm = TRUE)
  sig_i <- sigma_ts[season - 1L, p1]
  sig_j <- sigma_ts[season - 1L, p2]
  
  winning_prob(mu_i, sig_i, mu_j, sig_j)
})

prediction <- function(x) {
  as.integer(x >= 0.5)
}

# ---------------------------------------------------------------------------
# Cache environment – avoids re-indexing the same dataset in tight loops
# ---------------------------------------------------------------------------
.glicko_cache <- new.env(parent = emptyenv())

.build_period_cache <- function(dataset, p_x) {
  # Pre-split dataset by period (single pass), stored in the cache env
  by_p <- split(dataset, factor(dataset$Season, levels = p_x))
  for (k in seq_along(p_x)) {
    key <- paste0("period_", p_x[k])
    if (!exists(key, envir = .glicko_cache)) {
      .glicko_cache[[key]] <- by_p[[k]]
    }
  }
}

.clear_period_cache <- function() {
  rm(list = ls(.glicko_cache), envir = .glicko_cache)
}

# games_in_p variant that uses the cache (same signature for internal use)
.games_in_p_cached <- function(p, name, x) {
  key  <- paste0("period_", p)
  sub  <- if (exists(key, envir = .glicko_cache)) .glicko_cache[[key]] else
    x[x$Season == p, ]
  sub[sub$Player1 == name | sub$Player2 == name, ]
}

# ---------------------------------------------------------------------------
# Period update (uses cache internally)
# ---------------------------------------------------------------------------
update_ratings_for_period <- cmpfun(function(period, row_label, dataset, Ratings, ratings_ts, sigma_ts, attendance,
                                             rdmax = 350, c_rd = 63.2) {
  key <- paste0("period_", period)
  period_data <- if (exists(key, envir = .glicko_cache)) .glicko_cache[[key]] else
    dataset[dataset$Season == period, ]

  # Glicko step 1: RD growth over one rating period for every player who has
  # already debuted: sigma <- min(sqrt(sigma^2 + c^2), rdmax)
  prev_debut <- Ratings$debut == "Yes"
  if (any(prev_debut)) {
    Ratings$Sigma[prev_debut] <- pmin(sqrt(Ratings$Sigma[prev_debut]^2 + c_rd^2), rdmax)
  }

  # Integer indices of both players in every game of the period
  # (Game_Result = 1 means Player1 won)
  n_all <- nrow(Ratings)
  a <- match(period_data$Player1, Ratings$Players)
  b <- match(period_data$Player2, Ratings$Players)
  s <- period_data$Game_Result
  act <- sort(unique(c(a, b)))

  # Game / win / loss bookkeeping, debut and lag
  Ratings$Games <- Ratings$Games + tabulate(c(a, b), nbins = n_all)
  Ratings$Win   <- Ratings$Win   + tabulate(c(a[s == 1], b[s == 0]), nbins = n_all)
  Ratings$Loss  <- Ratings$Loss  + tabulate(c(a[s == 0], b[s == 1]), nbins = n_all)
  Ratings$debut[act] <- "Yes"

  idle_idx <- setdiff(which(Ratings$debut == "Yes"), act)
  if (length(idle_idx) > 0L) {
    idle_players <- Ratings$Players[idle_idx]
    Ratings$Lag[idle_idx]               <- Ratings$Lag[idle_idx] + 1L
    ratings_ts[row_label, idle_players] <- Ratings$Rating[idle_idx]
    sigma_ts[row_label,   idle_players] <- Ratings$Sigma[idle_idx]
    attendance[row_label, idle_players] <- "No"
  }

  # Active players: eqs (4)-(9), all games of the period at once, using the
  # pre-period values of every player (simultaneous update)
  if (length(act) > 0L) {
    mu0 <- Ratings$Rating
    sg0 <- Ratings$Sigma

    # Player1's side and Player2's side of every game
    G1 <- g_function(sg0[b]); E1 <- Es_uusj2(mu0[a], mu0[b], sg0[b])
    G2 <- g_function(sg0[a]); E2 <- Es_uusj2(mu0[b], mu0[a], sg0[a])
    v  <- rowsum(cbind(c(G1^2 * E1 * (1 - E1), G2^2 * E2 * (1 - E2)),   # 1/(q^2 delta^2) terms
                       c(G1 * (s - E1),        G2 * (1 - s - E2))),      # score residual terms
                 c(a, b))
    id <- as.integer(rownames(v))

    new_sigma <- pmin(rdmax, sqrt(1 / (1 / sg0[id]^2 + .q2 * v[, 1])))
    new_rating <- mu0[id] + .q * new_sigma^2 * v[, 2]

    Ratings$Sigma[id]  <- new_sigma
    Ratings$Rating[id] <- new_rating
    Ratings$Lag[id]    <- 0L

    act_players <- Ratings$Players[id]
    ratings_ts[row_label, act_players] <- new_rating
    sigma_ts[row_label,   act_players] <- new_sigma
    attendance[row_label, act_players] <- "Yes"
  }

  list(Ratings    = Ratings,
       ratings_ts = ratings_ts,
       sigma_ts   = sigma_ts,
       attendance = attendance)
})

# ---------------------------------------------------------------------------
# Likelihood (model fitting)
# ---------------------------------------------------------------------------
get_like <- cmpfun(function(para, train_data, rdmax = 350, c_rd = 63.2,
                            train_frac = 0.8, metric = "likelihood") {
  tic <- Sys.time()
  if (!(train_frac > 0 && train_frac <= 1)) stop("train_frac must be in (0, 1]")
  
  init_sigm <- para[1]
  dataset   <- train_data
  p_x       <- periods(dataset)
  
  Players  <- sort(unique(c(dataset$Player1, dataset$Player2)))
  nPlayers <- length(Players)
  nPeriods <- length(p_x)
  
  Ratings    <- Rating_Table_initialize(Players, 1500, init_sigm)
  row_labels <- paste0("p", p_x)
  
  ratings_ts <- matrix(1500,     nrow = nPeriods, ncol = nPlayers, dimnames = list(row_labels, Players))
  sigma_ts   <- matrix(init_sigm, nrow = nPeriods, ncol = nPlayers, dimnames = list(row_labels, Players))
  attendance <- matrix(NA,        nrow = nPeriods, ncol = nPlayers, dimnames = list(row_labels, Players))
  
  # Build period cache once per call
  .clear_period_cache()
  .build_period_cache(dataset, p_x)
  
  for (p in seq_along(p_x)) {
    res        <- update_ratings_for_period(p_x[p], row_labels[p], dataset,
                                            Ratings, ratings_ts, sigma_ts, attendance,
                                            rdmax = rdmax, c_rd = c_rd)
    Ratings    <- res$Ratings
    ratings_ts <- res$ratings_ts
    sigma_ts   <- res$sigma_ts
    attendance <- res$attendance
  }
  
  .clear_period_cache()
  
  # Vectorised winning probability over all rows
  # (period position, not raw Season value, indexes the time-series rows)
  # Vectorised equivalent of get_winning_prob(): mean rating over the previous
  # 6 periods and RD at the previous period; first 6 periods fixed at 0.5
  season_vec <- match(dataset$Season, p_x)
  wp <- rep(0.5, nrow(dataset))
  ok <- season_vec >= 7L
  if (any(ok)) {
    t  <- season_vec[ok]
    i1 <- match(dataset$Player1[ok], Players)
    i2 <- match(dataset$Player2[ok], Players)
    cs <- rbind(0, apply(ratings_ts, 2, cumsum))            # cs[k + 1, ] = sum of rows 1..k
    mu_i  <- (cs[cbind(t, i1)] - cs[cbind(t - 6L, i1)]) / 6
    mu_j  <- (cs[cbind(t, i2)] - cs[cbind(t - 6L, i2)]) / 6
    wp[ok] <- winning_prob(mu_i, sigma_ts[cbind(t - 1L, i1)], mu_j, sigma_ts[cbind(t - 1L, i2)])
  }

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
  print(paste(init_sigm, "time(mins):", round(as.numeric(toc - tic, units = "mins"), 4),
              "(auc, acc, likelihood):", paste(round(c(auc, acc, likelihood), 4), collapse = ", ")))

  record <<- rbind(record, c(init_sigm, auc, acc, likelihood))

  switch(metric,
         all        = c(auc, acc, likelihood),
         auc        = auc,
         acc        = acc,
         likelihood = likelihood,
         stop("choose a metric"))
})

# ---------------------------------------------------------------------------
# Core processor: handles both initial and incremental updates
# ---------------------------------------------------------------------------
process_standing <- cmpfun(function(para, dataset, prev_standings = NULL, rdmax = 350, c_rd = 63.2) {
  init_sigm <- para[1]
  p_x       <- periods(dataset)
  
  incoming_players <- sort(unique(c(dataset$Player1, dataset$Player2)))
  
  if (is.null(prev_standings)) {
    Players <- incoming_players
  } else {
    new_players <- setdiff(incoming_players, prev_standings$Player_Names)
    Players     <- sort(c(new_players, prev_standings$Player_Names))
  }
  
  extend_mat <- function(mat, new_row_labels, Players, new_val) {
    mat          <- as.matrix(mat)
    missing_cols <- setdiff(Players, colnames(mat))
    if (length(missing_cols) > 0L) {
      extra <- matrix(new_val, nrow = nrow(mat), ncol = length(missing_cols),
                      dimnames = list(rownames(mat), missing_cols))
      mat   <- cbind(mat, extra)
    }
    mat        <- mat[, Players, drop = FALSE]
    extra_rows <- matrix(new_val, nrow = length(new_row_labels), ncol = ncol(mat),
                         dimnames = list(new_row_labels, Players))
    rbind(mat, extra_rows)
  }
  
  if (is.null(prev_standings)) {
    Ratings    <- Rating_Table_initialize(Players, 1500, init_sigm)
    row_labels <- paste0("p", p_x)
    n          <- length(p_x)
    np         <- length(Players)
    ratings_ts <- matrix(1500,      nrow = n, ncol = np, dimnames = list(row_labels, Players))
    sigma_ts   <- matrix(init_sigm, nrow = n, ncol = np, dimnames = list(row_labels, Players))
    attendance <- matrix(NA,        nrow = n, ncol = np, dimnames = list(row_labels, Players))
  } else {
    Ratings     <- prev_standings$Current_Standings
    new_players <- setdiff(incoming_players, prev_standings$Player_Names)
    if (length(new_players) > 0L) {
      Ratings <- rbind(Ratings, Rating_Table_initialize(new_players, 1500, init_sigm))
    }
    last_num       <- max(as.numeric(sub("^p", "", rownames(prev_standings$Past_Ratings_time_series))))
    new_row_labels <- paste0("p", seq(last_num + 1L, length.out = length(p_x)))
    ratings_ts     <- extend_mat(prev_standings$Past_Ratings_time_series,              new_row_labels, Players, 1500)
    sigma_ts       <- extend_mat(prev_standings$Past_Rating_Deviations_time_series,    new_row_labels, Players, init_sigm)
    attendance     <- extend_mat(prev_standings$Attendance,                            new_row_labels, Players, NA)
    row_labels     <- new_row_labels
  }
  
  .clear_period_cache()
  .build_period_cache(dataset, p_x)
  
  for (pi in seq_along(p_x)) {
    res        <- update_ratings_for_period(p_x[pi], row_labels[pi], dataset,
                                            Ratings, ratings_ts, sigma_ts, attendance,
                                            rdmax = rdmax, c_rd = c_rd)
    Ratings    <- res$Ratings
    ratings_ts <- res$ratings_ts
    sigma_ts   <- res$sigma_ts
    attendance <- res$attendance
  }
  
  .clear_period_cache()
  
  list(Player_Names                         = Players,
       Current_Standings                    = Ratings,
       Past_Ratings_time_series             = ratings_ts,
       Past_Rating_Deviations_time_series   = sigma_ts,
       Attendance                           = attendance)
})

# ---------------------------------------------------------------------------
# User-facing wrappers (identical signatures to originals)
# ---------------------------------------------------------------------------
current_standing <- function(para, data, rdmax = 350, c_rd = 63.2) {
  process_standing(para, data, prev_standings = NULL, rdmax = rdmax, c_rd = c_rd)
}

further_standing <- function(current_standings_obj, new_data, para, rdmax = 350, c_rd = 63.2) {
  process_standing(para, new_data, prev_standings = current_standings_obj, rdmax = rdmax, c_rd = c_rd)
}