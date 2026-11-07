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
  this_period <- games_in_p(p, name, x)
  
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
  
  record <- games_in_p(p, name, x)
  
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
  # Pre-split dataset by period, stored in the cache env
  for (p in p_x) {
    key <- paste0("period_", p)
    if (!exists(key, envir = .glicko_cache)) {
      .glicko_cache[[key]] <- dataset[dataset$Season == p, ]
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
update_ratings_for_period <- cmpfun(function(period, row_label, dataset, Ratings, ratings_ts, sigma_ts, attendance) {
  key <- paste0("period_", period)
  period_data <- if (exists(key, envir = .glicko_cache)) .glicko_cache[[key]] else
    dataset[dataset$Season == period, ]
  
  players_period <- sort(unique(c(period_data$Player1, period_data$Player2)))
  
  # Debut flags
  debut_idx <- match(players_period, Ratings$Players)
  debut_idx <- debut_idx[!is.na(debut_idx)]
  Ratings$debut[debut_idx] <- "Yes"
  
  # Idle players
  idle_mask    <- Ratings$debut == "Yes" & !(Ratings$Players %in% players_period)
  idle_players <- Ratings$Players[idle_mask]
  
  if (length(idle_players) > 0L) {
    idle_idx <- which(idle_mask)
    Ratings$Lag[idle_idx]        <- Ratings$Lag[idle_idx] + 1L
    ratings_ts[row_label, idle_players] <- Ratings$Rating[idle_idx]
    sigma_ts[row_label,   idle_players] <- Ratings$Sigma[idle_idx]
    attendance[row_label, idle_players] <- "No"
  }
  
  # Active players
  if (length(players_period) > 0L) {
    active_mask  <- Ratings$Players %in% players_period
    old_ratings  <- Ratings[active_mask, , drop = FALSE]
    new_ratings  <- old_ratings
    
    for (i in players_period) {
      idx  <- which(Ratings$Players == i)
      nidx <- which(new_ratings$Players == i)
      oidx <- which(old_ratings$Players == i)
      
      gip <- .games_in_p_cached(period, i, dataset)
      
      # Vectorised win/loss counts
      p1_i    <- gip$Player1 == i
      win_c   <- sum(( p1_i & gip$Game_Result == 1L) | (!p1_i & gip$Game_Result == 0L))
      loss_c  <- sum(( p1_i & gip$Game_Result == 0L) | (!p1_i & gip$Game_Result == 1L))
      
      Ratings$Games[idx] <- Ratings$Games[idx] + nrow(gip)
      Ratings$Win[idx]   <- Ratings$Win[idx]   + win_c
      Ratings$Loss[idx]  <- Ratings$Loss[idx]  + loss_c
      Ratings$Lag[idx]   <- 0L
      Ratings$debut[idx] <- "Yes"
      
      old_sigma_i <- old_ratings$Sigma[oidx]
      new_sigma   <- min(rdmax,
                         player_sigma_update(old_sigma_i,
                                             delta(period, i, dataset, old_ratings)))
      new_rating  <- player_rating_update(period, i, dataset, old_ratings, new_sigma)
      
      new_ratings$Sigma[nidx]  <- new_sigma
      new_ratings$Rating[nidx] <- new_rating
      
      ratings_ts[row_label, i] <- new_rating
      sigma_ts[row_label,   i] <- new_sigma
      attendance[row_label, i] <- "Yes"
    }
    
    Ratings[active_mask, c("Rating", "Sigma")] <- new_ratings[, c("Rating", "Sigma")]
  }
  
  list(Ratings    = Ratings,
       ratings_ts = ratings_ts,
       sigma_ts   = sigma_ts,
       attendance = attendance)
})

# ---------------------------------------------------------------------------
# Likelihood (model fitting)
# ---------------------------------------------------------------------------
get_like <- cmpfun(function(para, train_data) {
  tic <- Sys.time()
  
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
                                            Ratings, ratings_ts, sigma_ts, attendance)
    Ratings    <- res$Ratings
    ratings_ts <- res$ratings_ts
    sigma_ts   <- res$sigma_ts
    attendance <- res$attendance
  }
  
  .clear_period_cache()
  
  # Vectorised winning probability over all rows
  season_vec <- dataset$Season
  wp <- vapply(seq_len(nrow(dataset)), function(ri) {
    get_winning_prob(c(season_vec[ri], dataset$Player1[ri], dataset$Player2[ri]),
                     ratings_ts, sigma_ts)
  }, numeric(1))
  
  gr          <- dataset$Game_Result
  result_prob <- wp^gr * (1 - wp)^(1 - gr)
  likelihood  <- sum(log(result_prob))
  
  toc <- Sys.time()
  print(paste(init_sigm, "time(s):", toc - tic, "likelihood:", likelihood))
  
  record <<- rbind(record, c(init_sigm, likelihood))
  
  likelihood
})

# ---------------------------------------------------------------------------
# Core processor: handles both initial and incremental updates
# ---------------------------------------------------------------------------
process_standing <- cmpfun(function(para, dataset, prev_standings = NULL) {
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
                                            Ratings, ratings_ts, sigma_ts, attendance)
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
current_standing <- function(para, data) {
  process_standing(para, data, prev_standings = NULL)
}

further_standing <- function(current_standings_obj, new_data, para) {
  process_standing(para, new_data, prev_standings = current_standings_obj)
}