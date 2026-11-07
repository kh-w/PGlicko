# ==============================================================================
# optimize_rating_systems.R
#
# Optimizes hyperparameters for four rating systems on the first 80% of data
# (by unique season), then evaluates all models on the held-out 20%.
#
# Models & parameters optimized
# ──────────────────────────────
#   Elo    (PlayerRatings::elo)   : K
#   PElo   (PELO::run_pelo)       : K, K2
#   Glicko (PlayerRatings::glicko): init_sig
#   PGlicko (custom PGlicko fns)  : init_sig, init_pw_sig
#
# Dataset requirements
# ────────────────────
#   CSV with four columns (names flexible – see col_rename() below):
#     Season      – integer or factor identifying the time period
#     Player1     – winner identifier
#     Player2     – loser identifier
#     Game_Result – 1 if Player1 won, 0 if Player2 won
#
# Dependencies
# ────────────
#   install.packages(c("PlayerRatings", "dplyr", "pROC", "PELO"))
#   # PELO from GitHub: remotes::install_github("kh-w/PELO")
#   # PGlicko helper functions sourced below
# ==============================================================================


# ==============================================================================
# 0. LIBRARIES & SOURCES
# ==============================================================================
library(PlayerRatings)
library(dplyr)
library(pROC)
remotes::install_github("kinhang8094-claude/PELO", ref = "speedup-loops")
library(PELO)
setwd("C:/Users/woodr/OneDrive - Wayne State University/PhD Wayne State/Research/Mathematics/Shinki/Bayesian Pairwise Strengths/")

PGLICKO_FUNCTIONS_PATH <- "PGlicko_functions_optimized_v3.R"
source(PGLICKO_FUNCTIONS_PATH)

# Glicko functions go in their own environment: several names clash with PGlicko's
GLICKO_FUNCTIONS_PATH <- "Glicko_functions_optimized_v4.R"
glicko_env <- new.env()
sys.source(GLICKO_FUNCTIONS_PATH, envir = glicko_env)


# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
SCENARIO_NAME <- "Sumo"
DATASET_PATH  <- paste0("datasets/simpledata_pairwise_glicko_compatible_", SCENARIO_NAME, ".csv")
INIT_RATING   <- 1500
RDMAX         <- 350
CVAL_GLICKO   <- 63.2
CVAL_PGLICKO  <- 63.2
PGLICKO_TRAIN_FRAC <- 0.8   # share of train_df used for the P-Glicko log-likelihood in get_like()
OUTPUT_DIR    <- paste0("outputs_", SCENARIO_NAME, "_v2")

ELO_K_INIT       <- 32
ELO_K_LOWER      <- 1
ELO_K_UPPER      <- 200

PELO_KK2_INIT    <- c(32, 16)
PELO_K_LOWER     <- 1;  PELO_K_UPPER  <- 200
PELO_K2_LOWER    <- 1;  PELO_K2_UPPER <- 200

GLICKO_SIG_INIT  <- 300
GLICKO_SIG_LOWER <- 50
GLICKO_SIG_UPPER <- 500

PGLICKO_INIT        <- c(300, 300)
PGLICKO_SIG_LOWER   <- 50;  PGLICKO_SIG_UPPER   <- 500
PGLICKO_PWSIG_LOWER <- 50;  PGLICKO_PWSIG_UPPER <- 500


# ==============================================================================
# 2. HELPER UTILITIES
# ==============================================================================

col_rename <- function(df) {
  if (ncol(df) != 4) { stop("Make sure your dataset is using format: c('Time','A','B','match_result')") }
  if (all(sort(unique(df[, 4])) != c(0, 1))) { stop("Make sure match result is binary") }
  colnames(df) <- c("Season", "Player1", "Player2", "Game_Result")
  return(df)
}

log_loss <- function(actual, prob) {
  prob <- pmax(pmin(prob, 1 - 1e-15), 1e-15)
  -mean(actual * log(prob) + (1 - actual) * log(1 - prob))
}

calc_auc <- function(actual, prob) {
  suppressWarnings(as.numeric(roc(actual, prob)$auc))
}

calc_acc <- function(actual, prob) mean((prob >= 0.5) == actual)


# ==============================================================================
# 3. LOAD & PREPARE DATASET
# ==============================================================================
cat("Loading dataset:", DATASET_PATH, "\n")
dataset <- read.csv(DATASET_PATH)
dataset <- col_rename(dataset)

unique_seasons <- sort(unique(dataset$Season))
dataset$Season <- match(dataset$Season, unique_seasons)

n_seasons    <- length(unique(dataset$Season))
train_cutoff <- floor(0.8 * n_seasons)

train_df <- dataset[dataset$Season <= train_cutoff, ]
test_df  <- dataset[dataset$Season >  train_cutoff, ]

cat(sprintf("Seasons total: %d  |  train: 1–%d  |  test: %d–%d\n",
            n_seasons, train_cutoff, train_cutoff + 1, n_seasons))
cat(sprintf("Train rows: %d  |  Test rows: %d\n",
            nrow(train_df), nrow(test_df)))

all_players   <- sort(unique(c(dataset$Player1, dataset$Player2)))
train_players <- sort(unique(c(train_df$Player1, train_df$Player2)))


# ==============================================================================
# 4. ELO  (PlayerRatings::elo)
# ==============================================================================
cat("\n── Optimising Elo ──────────────────────────────────────────────────────\n")

elo_walk_forward <- function(K, data, start_status, init_r = 1500) {
  seasons    <- sort(unique(data$Season))
  preds      <- numeric(nrow(data))
  actual     <- data$Game_Result
  cur_status <- start_status
  
  for (s in seasons) {
    idx   <- which(data$Season == s)
    chunk <- data[idx, ]
    
    new_pl <- setdiff(unique(c(chunk$Player1, chunk$Player2)), cur_status$Player)
    if (length(new_pl) > 0) {
      cur_status <- rbind(cur_status,
                          data.frame(Player = new_pl, Rating = init_r,
                                     Games = 0L, Win = 0L, Draw = 0L,
                                     Loss  = 0L, Lag = 0L,
                                     stringsAsFactors = FALSE))
    }
    
    r_lookup <- setNames(cur_status$Rating, cur_status$Player)
    for (g in seq_len(nrow(chunk))) {
      ri <- r_lookup[chunk$Player1[g]]
      rj <- r_lookup[chunk$Player2[g]]
      preds[idx[g]] <- 1 / (1 + 10^((rj - ri) / 400))
    }
    
    cur_status <- elo(x = chunk, status = cur_status,
                      kfac = K, history = FALSE)$ratings
  }
  
  list(preds = preds, actual = actual, end_status = cur_status)
}

# Runs walk-forward over `data` starting from `start_status`.
# Returns preds, metrics, and end_status for chaining train → test.
run_elo_playerratings <- function(K, data, start_status = NULL, init_r = INIT_RATING) {
  if (is.null(start_status)) {
    players      <- sort(unique(c(data$Player1, data$Player2)))
    start_status <- data.frame(Player = players, Rating = init_r,
                               Games  = 0L, Win = 0L, Draw = 0L,
                               Loss   = 0L, Lag = 0L,
                               stringsAsFactors = FALSE)
  }
  
  res <- elo_walk_forward(K, data, start_status, init_r)
  
  list(preds      = res$preds,
       actual     = res$actual,
       auc        = calc_auc(res$actual, res$preds),
       acc        = calc_acc(res$actual, res$preds),
       ll         = log_loss(res$actual, res$preds),
       end_status = res$end_status)
}

elo_objective <- function(K, train) {
  if (K <= 0) return(-Inf)
  run_elo_playerratings(K, train)$auc
}

elo_opt <- optimise(f        = function(K) elo_objective(K, train_df),
                    interval = c(ELO_K_LOWER, ELO_K_UPPER),
                    maximum  = TRUE)
best_K_elo <- elo_opt$maximum
cat(sprintf("  Best K = %.4f  (train AUC = %.4f)\n", best_K_elo, elo_opt$objective))


# ==============================================================================
# 5. P-ELO  (PELO package)
# ==============================================================================
cat("\n── Optimising P-Elo ────────────────────────────────────────────────────\n")

to_pelo_fmt <- function(df) {
  data.frame(Time         = df$Season,
             A            = df$Player1,
             B            = df$Player2,
             match_result = df$Game_Result,
             stringsAsFactors = FALSE)
}

train_pelo <- to_pelo_fmt(train_df)
test_pelo  <- to_pelo_fmt(test_df)

players_train     <- sort(unique(c(train_pelo$A, train_pelo$B)))
init_ratings_pelo <- initialize_ratings(players_train, INIT_RATING)
init_pair_adv     <- initialize_pair_adv_mat(players_train, 0)

best_KK2 <- optimize_pelo_kk2(train_pelo, init_ratings_pelo, init_pair_adv)
cat(sprintf("  Best K = %.4f, K2 = %.4f\n", best_KK2[1], best_KK2[2]))

players_id_pelo <- get_player_id_mapping(train_pelo)

# Train phase: fresh state → end-of-train state
pelo_train_obj <- run_pelo(best_KK2, train_pelo,
                           init_ratings_pelo, init_pair_adv,
                           players_id_pelo,
                           new_player_initial_rating = INIT_RATING)

# Test phase: continues from end-of-train ratings and pair_adv
pelo_test_obj <- run_pelo(best_KK2, test_pelo,
                          pelo_train_obj$ratings, pelo_train_obj$pair_adv,
                          players_id_pelo,
                          new_player_initial_rating = INIT_RATING)


# ==============================================================================
# 6. GLICKO  (PlayerRatings::glicko)
# ==============================================================================
cat("\n── Optimising Glicko ───────────────────────────────────────────────────\n")

glicko_walk_forward <- function(init_sig, data, start_status, init_r = 1500,
                                cval = 0, rdmax = 350) {
  seasons    <- sort(unique(data$Season))
  preds      <- numeric(nrow(data))
  actual     <- data$Game_Result
  cur_status <- start_status
  
  for (s in seasons) {
    idx   <- which(data$Season == s)
    chunk <- data[idx, ]

    # cur_status is a Glicko v4 standings object (NULL before the first period);
    # players not yet rated start at init_r / init_sig
    all_pl    <- unique(c(chunk$Player1, chunk$Player2))
    r_lookup  <- setNames(rep(init_r,   length(all_pl)), all_pl)
    rd_lookup <- setNames(rep(init_sig, length(all_pl)), all_pl)
    if (!is.null(cur_status)) {
      cs <- cur_status$Current_Standings
      r_lookup  <- c(setNames(cs$Rating, cs$Players), r_lookup)
      rd_lookup <- c(setNames(cs$Sigma,  cs$Players), rd_lookup)
    }
    for (g in seq_len(nrow(chunk))) {
      ri  <- r_lookup[chunk$Player1[g]];  rj  <- r_lookup[chunk$Player2[g]]
      rdi <- rd_lookup[chunk$Player1[g]]; rdj <- rd_lookup[chunk$Player2[g]]
      q    <- log(10) / 400
      g_rd <- 1 / sqrt(1 + 3 * q^2 * (rdi^2 + rdj^2) / pi^2)
      preds[idx[g]] <- 1 / (1 + 10^(-g_rd * (ri - rj) / 400))
    }

    cur_status <- if (is.null(cur_status)) {
      glicko_env$current_standing(init_sig, chunk, rdmax = rdmax, c_rd = cval)
    } else {
      glicko_env$further_standing(cur_status, chunk, init_sig, rdmax = rdmax, c_rd = cval)
    }
  }
  
  list(preds = preds, actual = actual, end_status = cur_status)
}

# Runs walk-forward over `data` starting from `start_status`.
# Returns preds, metrics, and end_status for chaining train → test.
run_glicko_playerratings <- function(init_sig, data, start_status = NULL,
                                     init_r = INIT_RATING, cval = 0, rdmax = 350) {
  # start_status = NULL means a fresh start (handled inside glicko_walk_forward)
  res <- glicko_walk_forward(init_sig, data, start_status, init_r, cval, rdmax)
  
  list(preds      = res$preds,
       actual     = res$actual,
       auc        = calc_auc(res$actual, res$preds),
       acc        = calc_acc(res$actual, res$preds),
       ll         = log_loss(res$actual, res$preds),
       end_status = res$end_status)
}

glicko_objective <- function(init_sig, train) {
  if (init_sig <= 0) return(-Inf)
  run_glicko_playerratings(init_sig, train,
                           init_r = INIT_RATING, cval = CVAL_GLICKO, rdmax = RDMAX)$auc
}

glicko_opt <- optimise(f        = function(s) glicko_objective(s, train_df),
                       interval = c(GLICKO_SIG_LOWER, GLICKO_SIG_UPPER),
                       maximum  = TRUE)
best_sig_glicko <- glicko_opt$maximum
cat(sprintf("  Best init_sig = %.4f  (train AUC = %.4f)\n",
            best_sig_glicko, glicko_opt$objective))


# ==============================================================================
# 7. P-GLICKO  (custom PGlicko functions)
# ==============================================================================
cat("\n── Optimising PGlicko ──────────────────────────────────────────────────\n")

build_cache(train_df)
p_x_train      <- .cache$p_x
Players_train  <- .cache$Players
PairList_train <- .cache$Pairwise_List_Vector

record <- data.frame(init_sigm    = numeric(),
                     init_pw_sigm = numeric(),
                     auc          = numeric(),
                     acc          = numeric(),
                     likelihood   = numeric())

pglicko_objective <- function(para) {
  if (any(para <= 0)) return(-Inf)
  get_like(para, train_df, rdmax = RDMAX, c_rd = CVAL_PGLICKO, train_frac = PGLICKO_TRAIN_FRAC)
}

pglicko_opt <- optim(par     = PGLICKO_INIT,
                     fn      = pglicko_objective,
                     method  = "Nelder-Mead",
                     control = list(fnscale = -1,
                                    maxit   = 500))
# stop("stop")
best_sig_pglicko   <- pglicko_opt$par[1]
best_pwsig_pglicko <- pglicko_opt$par[2]

cat(sprintf("  Best init_sig = %.4f, init_pw_sig = %.4f  (train log-lik = %.4f)\n",
            best_sig_pglicko, best_pwsig_pglicko, pglicko_opt$value))


# ==============================================================================
# 8. WALK-FORWARD ENGINE & WRAPPERS
# ==============================================================================

# Walk-forward for PGlicko: predict each period's games from the ratings at the
# start of the period, then update (pglicko_phase() in the PGlicko functions).
# The state - ratings plus a compact table of pairwise strengths and deviations
# for the pairs that have played - is returned for chaining train -> test.
run_pglicko_phase <- function(data, cur_state, init_sig, init_pw_sig,
                              rdmax = RDMAX, c_rd = CVAL_PGLICKO) {
  out    <- pglicko_phase(c(init_sig, init_pw_sig), data, state = cur_state,
                          rdmax = rdmax, c_rd = c_rd, init_rating = INIT_RATING)
  preds  <- out$preds
  actual <- data$Game_Result
  list(preds      = preds,
       pred_class = as.integer(preds >= 0.5),
       actual     = actual,
       auc        = calc_auc(actual, preds),
       acc        = calc_acc(actual, preds),
       ll         = log_loss(actual, preds),
       end_state  = out$state)
}

# Wrapper: fresh start (every player at INIT_RATING) and walk forward over `data`.
run_pglicko <- function(data,
                        init_r      = INIT_RATING,
                        init_sig    = best_sig_pglicko,
                        init_pw_sig = best_pwsig_pglicko,
                        rdmax       = RDMAX,
                        c_rd        = CVAL_PGLICKO) {
  run_pglicko_phase(data, NULL, init_sig, init_pw_sig, rdmax, c_rd)
}


# ==============================================================================
# 9. TRAIN PHASE — walk-forward over training data
# ==============================================================================
cat("\n── Train phase ─────────────────────────────────────────────────────────\n")

# ── 9a. Elo ──────────────────────────────────────────────────────────────────
elo_train_obj <- run_elo_playerratings(best_K_elo, train_df)

# ── 9b. P-Elo ────────────────────────────────────────────────────────────────
pelo_train_obj <- run_pelo(best_KK2, train_pelo,
                           init_ratings_pelo, init_pair_adv,
                           players_id_pelo,
                           new_player_initial_rating = INIT_RATING)

# ── 9c. Glicko ───────────────────────────────────────────────────────────────
glicko_train_obj <- run_glicko_playerratings(best_sig_glicko, train_df,
                                             init_r = INIT_RATING,
                                             cval   = CVAL_GLICKO,
                                             rdmax  = RDMAX)

# ── 9d. PGlicko ──────────────────────────────────────────────────────────────
pglicko_train_obj <- run_pglicko(train_df)


# ==============================================================================
# 10. TEST PHASE — continues from end-of-train state, optimized parameters
# ==============================================================================
cat("\n── Test phase ──────────────────────────────────────────────────────────\n")

# ── 10a. Elo ─────────────────────────────────────────────────────────────────
elo_test_obj <- run_elo_playerratings(best_K_elo, test_df,
                                      start_status = elo_train_obj$end_status)

# ── 10b. P-Elo ───────────────────────────────────────────────────────────────
pelo_test_obj <- run_pelo(best_KK2, test_pelo,
                          pelo_train_obj$ratings, pelo_train_obj$pair_adv,
                          players_id_pelo,
                          new_player_initial_rating = INIT_RATING)

pelo_test_auc <- calc_auc(test_df$Game_Result, pelo_test_obj$preds)
pelo_test_acc <- calc_acc(test_df$Game_Result, pelo_test_obj$preds)
pelo_test_ll  <- log_loss(test_df$Game_Result, pelo_test_obj$preds)

# ── 10c. Glicko ──────────────────────────────────────────────────────────────
glicko_test_obj <- run_glicko_playerratings(best_sig_glicko, test_df,
                                            start_status = glicko_train_obj$end_status,
                                            init_r = INIT_RATING,
                                            cval   = CVAL_GLICKO,
                                            rdmax  = RDMAX)

# ── 10d. PGlicko ─────────────────────────────────────────────────────────────
# Continues from the end-of-train state (ratings and pairwise table)
pglicko_test_obj <- run_pglicko_phase(test_df,
                                      pglicko_train_obj$end_state,
                                      best_sig_pglicko,
                                      best_pwsig_pglicko,
                                      RDMAX,
                                      CVAL_PGLICKO)


# ==============================================================================
# 11. SUMMARY TABLE  (test-set performance only)
# ==============================================================================
cat("\n── Summary ─────────────────────────────────────────────────────────────\n")

results <- data.frame(
  Model    = c("Elo", "PElo", "Glicko", "PGlicko"),
  Param1   = c(
    sprintf("K=%.4f",           best_K_elo),
    sprintf("K=%.4f",           best_KK2[1]),
    sprintf("init_sig=%.4f",    best_sig_glicko),
    sprintf("init_sig=%.4f",    best_sig_pglicko)
  ),
  Param2   = c(
    "",
    sprintf("K2=%.4f",          best_KK2[2]),
    "",
    sprintf("init_pw_sig=%.4f", best_pwsig_pglicko)
  ),
  Test_AUC = round(c(elo_test_obj$auc, pelo_test_auc, glicko_test_obj$auc, pglicko_test_obj$auc), 4),
  Test_Acc = round(c(elo_test_obj$acc, pelo_test_acc, glicko_test_obj$acc, pglicko_test_obj$acc), 4),
  Test_LL  = round(c(elo_test_obj$ll,  pelo_test_ll,  glicko_test_obj$ll,  pglicko_test_obj$ll),  6),
  stringsAsFactors = FALSE
)

print(results, row.names = FALSE)


# ==============================================================================
# 12. SAVE RESULTS
# ==============================================================================
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

write.csv(results,
          paste0(OUTPUT_DIR, "/optimization_results_", SCENARIO_NAME, ".csv"),
          row.names = FALSE)

# Train predictions
pred_train <- data.frame(
  Season       = train_df$Season,
  Player1      = train_df$Player1,
  Player2      = train_df$Player2,
  Game_Result  = train_df$Game_Result,
  Elo_pred     = elo_train_obj$preds,
  PElo_pred    = pelo_train_obj$preds,
  Glicko_pred  = glicko_train_obj$preds,
  PGlicko_pred = pglicko_train_obj$preds
)
write.csv(pred_train,
          paste0(OUTPUT_DIR, "/train_predictions_", SCENARIO_NAME, ".csv"),
          row.names = FALSE)

# Test predictions
pred_test <- data.frame(
  Season       = test_df$Season,
  Player1      = test_df$Player1,
  Player2      = test_df$Player2,
  Game_Result  = test_df$Game_Result,
  Elo_pred     = elo_test_obj$preds,
  PElo_pred    = pelo_test_obj$preds,
  Glicko_pred  = glicko_test_obj$preds,
  PGlicko_pred = pglicko_test_obj$preds
)
write.csv(pred_test,
          paste0(OUTPUT_DIR, "/test_predictions_", SCENARIO_NAME, ".csv"),
          row.names = FALSE)

cat(sprintf("\nResults saved to:\n  optimization_results_%s.csv\n  train_predictions_%s.csv\n  test_predictions_%s.csv\n",
            SCENARIO_NAME, SCENARIO_NAME, SCENARIO_NAME))


# ==============================================================================
# 13. HYPOTHETICAL VALUE BETTING  (test set only)
# ==============================================================================
cat("\n── Hypothetical Value Betting ──────────────────────────────────────────\n")

N_SIMS  <- 10000L
N_BETS  <- as.integer(nrow(test_df) * 0.5)

actual_outcomes <- test_df$Game_Result

elo_p    <- elo_test_obj$preds
pelo_p   <- pelo_test_obj$preds
glicko_p <- glicko_test_obj$preds
pg_p     <- pglicko_test_obj$preds

set.seed(4445)

cat("  Running PElo vs Elo       ...\n")
profits_PElo_vs_Elo     <- hypo_bet_scheme_value(elo_p,    pelo_p,   actual_outcomes, N_SIMS, N_BETS)
cat("  Running Elo vs PElo       ...\n")
profits_Elo_vs_PElo     <- hypo_bet_scheme_value(pelo_p,    elo_p,   actual_outcomes, N_SIMS, N_BETS)

cat("  Running Glicko vs Elo     ...\n")
profits_Glicko_vs_Elo   <- hypo_bet_scheme_value(elo_p,    glicko_p, actual_outcomes, N_SIMS, N_BETS)
cat("  Running Elo vs Glicko     ...\n")
profits_Elo_vs_Glicko   <- hypo_bet_scheme_value(glicko_p,    elo_p, actual_outcomes, N_SIMS, N_BETS)

cat("  Running PGlicko vs Elo    ...\n")
profits_PGlicko_vs_Elo  <- hypo_bet_scheme_value(elo_p,    pg_p,     actual_outcomes, N_SIMS, N_BETS)
cat("  Running Elo vs PGlicko    ...\n")
profits_Elo_vs_PGlicko  <- hypo_bet_scheme_value(pg_p,    elo_p,     actual_outcomes, N_SIMS, N_BETS)

cat("  Running Glicko vs PElo    ...\n")
profits_Glicko_vs_PElo  <- hypo_bet_scheme_value(pelo_p,   glicko_p, actual_outcomes, N_SIMS, N_BETS)
cat("  Running PElo vs Glicko    ...\n")
profits_PElo_vs_Glicko  <- hypo_bet_scheme_value(glicko_p,   pelo_p, actual_outcomes, N_SIMS, N_BETS)

cat("  Running PGlicko vs PElo   ...\n")
profits_PGlicko_vs_PElo <- hypo_bet_scheme_value(pelo_p,   pg_p,     actual_outcomes, N_SIMS, N_BETS)
cat("  Running PElo vs  PGlicko  ...\n")
profits_PElo_vs_PGlicko <- hypo_bet_scheme_value(pg_p,   pelo_p,     actual_outcomes, N_SIMS, N_BETS)

cat("  Running PGlicko vs Glicko   ...\n")
profits_PGlicko_vs_Glicko <- hypo_bet_scheme_value(glicko_p,   pg_p,     actual_outcomes, N_SIMS, N_BETS)
cat("  Running Glicko vs  PGlicko  ...\n")
profits_Glicko_vs_PGlicko <- hypo_bet_scheme_value(pg_p,   glicko_p,     actual_outcomes, N_SIMS, N_BETS)

profit_list <- list(
  "PElo vs Elo"     = profits_PElo_vs_Elo,
  "Elo vs PElo"     = profits_Elo_vs_PElo,
  "Glicko vs Elo"   = profits_Glicko_vs_Elo,
  "Elo vs Glicko"   = profits_Elo_vs_Glicko,
  "PGlicko vs Elo"  = profits_PGlicko_vs_Elo,
  "Elo vs PGlicko"  = profits_Elo_vs_PGlicko,
  "Glicko vs PElo"  = profits_Glicko_vs_PElo,
  "PElo vs Glicko"  = profits_PElo_vs_Glicko,
  "PGlicko vs PElo" = profits_PGlicko_vs_PElo,
  "PElo vs PGlicko" = profits_PElo_vs_PGlicko,
  "PGlicko vs Glicko" = profits_PGlicko_vs_Glicko,
  "Glicko vs PGlicko" = profits_Glicko_vs_PGlicko
)

rangeminmax <- range(unlist(profit_list), na.rm = TRUE)

betting_summary <- do.call(rbind, lapply(names(profit_list), function(nm) {
  x <- profit_list[[nm]]
  data.frame(
    Comparison    = nm,
    Mean_Profit   = round(mean(x),     4),
    SD_Profit     = round(sd(x),       4),
    Median_Profit = round(median(x),   4),
    P_Profit_Pos  = round(mean(x > 0), 4),
    stringsAsFactors = FALSE
  )
}))

cat("\n  Hypothetical Value Betting Summary\n")
print(betting_summary, row.names = FALSE)

pdf(paste0(OUTPUT_DIR, "/hypo_betting_", SCENARIO_NAME, ".pdf"),
    width = 9, height = 11, family = "serif")

op <- par(no.readonly = TRUE)
par(mfrow    = c(2, 1),
    mar      = c(5.5, 6.5, 7, 1.5),
    las      = 1,
    cex.main = 2.0,
    cex.lab  = 1.7,
    cex.axis = 1.5,
    lwd      = 1.6,
    mgp      = c(4, 1.2, 0),
    family   = "serif")

all_counts <- lapply(profit_list, function(x) hist(x, breaks = 25, plot = FALSE)$counts)
ymax <- max(unlist(all_counts), na.rm = TRUE)

plot_profit_hist <- function(x, label) {
  m  <- round(mean(x), 2)
  pp <- round(mean(x > 0) * 100, 1)
  hist(x,
       breaks = 25,
       col    = "lightgray",
       border = "black",
       main   = paste0(label, "\nmean = ", m, "  |  P(profit>0) = ", pp, "%"),
       xlab   = "Profit",
       ylab   = "Frequency",
       ylim   = c(0, ymax),
       xlim   = c(rangeminmax[1], rangeminmax[2]))
  abline(v = 0, col = "red", lty = 2, lwd = 2)
  box(lwd = 1.5)
}

for (nm in names(profit_list)) plot_profit_hist(profit_list[[nm]], nm)
plot.new()

par(op)
dev.off()
cat(sprintf("  Betting plot saved to: hypo_betting_%s.pdf\n", SCENARIO_NAME))

write.csv(betting_summary,
          paste0(OUTPUT_DIR, "/hypo_betting_summary_", SCENARIO_NAME, ".csv"),
          row.names = FALSE)

profits_df <- as.data.frame(profit_list)
colnames(profits_df) <- gsub(" ", "_", names(profit_list))
write.csv(profits_df,
          paste0(OUTPUT_DIR, "/hypo_betting_sims_", SCENARIO_NAME, ".csv"),
          row.names = FALSE)

cat(sprintf("  Betting results saved to:\n    hypo_betting_summary_%s.csv\n    hypo_betting_sims_%s.csv\n",
            SCENARIO_NAME, SCENARIO_NAME))

save.image(paste0(OUTPUT_DIR, "/optimize_rating_systems_", SCENARIO_NAME, ".RData"))