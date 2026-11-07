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

  # The per-player and per-pair splits of every season are only needed by the
  # per-player helper functions (delta_player_rating() etc.), not by get_like()
  # or pglicko_phase(). They hold one small data frame per player and per pair
  # per season, so they are built lazily, on first use (see .ensure_splits()).
  .cache$by_season_player <- NULL
  .cache$by_season_pair   <- NULL

  # Remember which dataset this cache was built from (checked in get_like)
  .cache$dataset <- dataset
  
  # Pre-compute pair keys for every row
  .cache$pair_keys <- paste(dataset$Player1, dataset$Player2, sep = "_")
  
  # FIX 1: cache stable values — periods, players, pair list
  # so get_like never recomputes these across ~175 optim() calls
  .cache$p_x      <- periods(dataset)
  .cache$Players  <- players(dataset)
  # Pairs that actually meet in the data ("p1_p2", p1 first in sort order).
  # Listing every possible pair (combn) would need players^2 / 2 strings.
  lo <- pmin(dataset$Player1, dataset$Player2); hi <- pmax(dataset$Player1, dataset$Player2)
  .cache$Pairwise_List_Vector <- sort(unique(paste(lo, hi, sep = "_")))
  .cache$ps_col_idx <- setNames(
    seq_along(.cache$Pairwise_List_Vector),
    .cache$Pairwise_List_Vector
  )
  
  invisible(NULL)
}

# ============================================================
# FAST LOOKUP REPLACEMENTS
# ============================================================
.ensure_splits <- function() {
  if (is.null(.cache$by_season_player)) {
    # Per season: split by player
    .cache$by_season_player <- lapply(.cache$by_season, function(d) {
      all_p <- unique(c(d$Player1, d$Player2))
      setNames(
        lapply(all_p, function(pl)
          d[d$Player1 == pl | d$Player2 == pl, , drop = FALSE]),
        all_p
      )
    })
    # Per season: split by UNORDERED pair key, so all n_j games between two
    # players are pooled regardless of who is listed as Player1 (eqs 30, 32)
    .cache$by_season_pair <- lapply(.cache$by_season, function(d) {
      keys <- paste(pmin(d$Player1, d$Player2), pmax(d$Player1, d$Player2), sep = "_")
      split(d, keys)
    })
  }
  invisible(NULL)
}

get_player_games <- function(p, name) {
  .ensure_splits()
  sp <- .cache$by_season_player[[as.character(p)]]
  if (is.null(sp)) return(NULL)
  sp[[name]]
}

get_pair_games <- function(p, name1, name2) {
  .ensure_splits()
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
# processed at once with integer player indices and plain vectors:
#   step 1  RD growth for debuted players: sigma <- min(sqrt(sigma^2 + c^2), rdmax)
#   player updates   eqs (16)-(21), simultaneous, using pre-period values
#   pairwise updates eqs (28)-(32), all games of an unordered pair pooled,
#                    beta_ji = -beta_ij, b_ji = b_ij
#
# Pairwise strengths and their deviations are stored once per pair of players
# that can meet (a compact pair table), not as players x players matrices, so
# memory grows with the number of pairs that play, not with players^2.
#
# predict = "avg6":    eq (34) with the mean rating over the previous 6 periods
#                      and the RD / pairwise values at the end of the previous
#                      period (used by get_like); the first 6 periods get 0.5.
# predict = "current": eq (34) with each player's rating and RD at the start of
#                      the period (used by pglicko_phase).
#
# state (optional) continues from an earlier run: list(Ratings, pairs) as
# returned in $state. New players start at init_rating / para[1]; new pairs at
# strength 0 and deviation para[2].
# ============================================================
.pglicko_engine <- function(para, dataset, rdmax = 350, c_rd = 63.2, state = NULL,
                            predict = c("avg6", "current"), init_rating = 1500) {

  predict <- match.arg(predict)
  init_sigm    <- para[1]
  init_pw_sigm <- para[2]

  p_x <- periods(dataset)
  TT  <- length(p_x)
  period_labels <- paste0("p", p_x)

  # Players: sorted when starting fresh; otherwise the state's players followed
  # by new players in order of first appearance (period by period)
  if (is.null(state)) {
    Players <- players(dataset)
    np <- length(Players)
    mu <- rep(init_rating, np); sg <- rep(init_sigm, np); debut <- logical(np)
    games <- win <- loss <- lag <- integer(np)
  } else {
    R0 <- state$Ratings
    first_seen <- unique(unlist(lapply(split(dataset, factor(dataset$Season, levels = p_x)),
                                       function(d) unique(c(d$Player1, d$Player2)))))
    new_pl  <- setdiff(first_seen, R0$Players)
    Players <- c(R0$Players, new_pl)
    np <- length(Players); nn <- length(new_pl)
    mu <- c(unname(R0$Rating[R0$Players]), rep(init_rating, nn))
    sg <- c(unname(R0$Sigma[R0$Players]),  rep(init_sigm, nn))
    debut <- c(unname(R0$debut[R0$Players]) == "Yes", rep(FALSE, nn))
    games <- c(unname(R0$Games[R0$Players]), integer(nn)); win <- c(unname(R0$Win[R0$Players]), integer(nn))
    loss  <- c(unname(R0$Loss[R0$Players]),  integer(nn)); lag <- c(unname(R0$Lag[R0$Players]), integer(nn))
  }

  a_all <- match(dataset$Player1, Players)
  b_all <- match(dataset$Player2, Players)
  s_all <- dataset$Game_Result
  t_all <- match(dataset$Season, p_x)
  rows  <- split(seq_along(t_all), factor(t_all, levels = seq_len(TT)))

  # ---- pair table: one strength (beta, from the lower-index player's side) and
  #      one deviation per unordered pair; pairs from `state` are carried over
  lo_all <- pmin(a_all, b_all); hi_all <- pmax(a_all, b_all)
  key_all <- lo_all * (np + 1) + hi_all
  pkey <- unique(key_all)
  beta <- numeric(length(pkey)); bsd <- rep(init_pw_sigm, length(pkey))
  if (!is.null(state) && nrow(state$pairs) > 0) {
    i1 <- match(state$pairs$Player1, Players); i2 <- match(state$pairs$Player2, Players)
    sk <- pmin(i1, i2) * (np + 1) + pmax(i1, i2)
    sv <- ifelse(i1 <= i2, state$pairs$PairAdv, -state$pairs$PairAdv)
    extra <- !(sk %in% pkey)
    pkey <- c(pkey, sk[extra]); beta <- c(beta, numeric(sum(extra))); bsd <- c(bsd, numeric(sum(extra)))
    m <- match(sk, pkey); beta[m] <- sv; bsd[m] <- state$pairs$PairSD
  }
  pid_all <- match(key_all, pkey)
  sgn_all <- ifelse(a_all < b_all, 1, ifelse(a_all > b_all, -1, 0))

  keep_ts <- predict == "avg6"
  if (keep_ts) {
    ratings_ts <- matrix(init_rating, TT, np, dimnames = list(period_labels, Players))
    sigma_ts   <- matrix(init_sigm,   TT, np, dimnames = list(period_labels, Players))
  }

  # Per-game values at the start of the period (for prediction)
  g_beta <- g_b <- numeric(length(s_all))
  if (!keep_ts) g_mua <- g_mub <- g_sga <- g_sgb <- numeric(length(s_all))
  q <- .Q; K <- .K

  for (t in seq_len(TT)) {
    r  <- rows[[t]]
    a  <- a_all[r];  b <- b_all[r];  s <- s_all[r]
    k  <- pid_all[r]
    bet <- sgn_all[r] * beta[k];  bsdv <- bsd[k]
    g_beta[r] <- bet;  g_b[r] <- bsdv
    if (!keep_ts) { g_mua[r] <- mu[a]; g_mub[r] <- mu[b]; g_sga[r] <- sg[a]; g_sgb[r] <- sg[b] }

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

    # Player updates, eqs (16)-(21): Player1's side and Player2's side of every game
    G1 <- 1 / sqrt(1 + K * (sg0[b]^2 + bsdv^2))
    E1 <- 1 / (1 + 10^(-G1 * (mu0[a] + bet - mu0[b]) / 400))
    G2 <- 1 / sqrt(1 + K * (sg0[a]^2 + bsdv^2))
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
    kk  <- match(uk, pkey)
    G   <- 1 / sqrt(1 + K * (sg0[pl]^2 + sg0[ph]^2))
    E   <- 1 / (1 + 10^(-G * (mu0[pl] + beta[kk] - mu0[ph]) / 400))
    b_new    <- pmin(rdmax, sqrt(1 / (1 / bsd[kk]^2 + q^2 * w[, 1] * G^2 * E * (1 - E))))
    beta[kk] <- beta[kk] + q * b_new^2 * G * (w[, 2] - w[, 1] * E)
    bsd[kk]  <- b_new

    if (keep_ts) { ratings_ts[t, ] <- mu; sigma_ts[t, ] <- sg }
  }

  # Predictions, eq (34)
  if (keep_ts) {
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
  } else {
    wp <- winning_prob(g_mua, g_sga, g_mub, g_sgb, g_beta, g_b)
  }

  nm <- function(x) setNames(x, Players)
  pk_lo <- pkey %/% (np + 1); pk_hi <- pkey %% (np + 1)
  st <- list(
    Ratings = list(Rating = nm(mu), Sigma = nm(sg), Games = nm(games), Win = nm(win),
                   Loss = nm(loss), Lag = nm(lag), debut = nm(ifelse(debut, "Yes", "No")),
                   Players = Players),
    pairs = data.frame(Player1 = Players[pk_lo], Player2 = Players[pk_hi],
                       PairAdv = beta, PairSD = bsd, stringsAsFactors = FALSE))
  out <- list(state = st, wp = wp)
  if (keep_ts) { out$ratings_ts <- ratings_ts; out$sigma_ts <- sigma_ts }
  out
}

# ============================================================
# WALK-FORWARD PHASE (predict, then update), memory-efficient
# ============================================================
# Runs P-Glicko over `data` period by period: each period's games are first
# predicted from the ratings, RDs and pairwise values at the start of the
# period (eq 34), then all of the period's updates are applied together.
#
#   para   c(initial sigma, initial pairwise sigma)
#   data   data frame with Season, Player1, Player2, Game_Result
#   state  NULL to start fresh (all players at init_rating), or the $state of
#          an earlier call, to continue (e.g. from training into testing).
#          New players and pairs are added automatically.
#
# Returns list(preds, state). state$pairs is a compact pair table with columns
# Player1, Player2, PairAdv (Player1's pairwise strength over Player2) and
# PairSD, one row per pair that has played; state$Ratings is a rating list as
# from Rating_Table_initialize().
pglicko_phase <- function(para, data, state = NULL, rdmax = 350, c_rd = 63.2, init_rating = 1500) {
  e <- .pglicko_engine(para, data, rdmax = rdmax, c_rd = c_rd, state = state,
                       predict = "current", init_rating = init_rating)
  list(preds = e$wp, state = e$state)
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
