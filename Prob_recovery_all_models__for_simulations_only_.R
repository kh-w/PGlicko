library(dplyr)

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
setwd("C:/Users/woodr/OneDrive - Wayne State University/PhD Wayne State/Research/Mathematics/Shinki/Bayesian Pairwise Strengths/")

SCENARIO_NAME <- "5P"
OUTPUT_DIR    <- paste0("outputs_", SCENARIO_NAME, "_v2")   # outputs of the new Glicko / PGlicko runs

WP_FILE      <- paste0("datasets/simpledata_pairwise_glicko_compatible_", SCENARIO_NAME, "_wp.csv")
MODEL_COLS   <- c("Elo_pred", "PElo_pred", "Glicko_pred", "PGlicko_pred")
MODEL_LABELS <- c("Elo", "PElo", "Glicko", "PGlicko")

SERIES_COLORS <- c(True    = "black",
                   Elo     = "#E41A1C",
                   PElo    = "#377EB8",
                   Glicko  = "#4DAF4A",
                   PGlicko = "#984EA3")
SERIES_LTY <- c(True = 1, Elo = 2, PElo = 3, Glicko = 4, PGlicko = 5)
SERIES_LWD <- c(True = 2, Elo = 1.5, PElo = 1.5, Glicko = 1.5, PGlicko = 1.5)


# ==============================================================================
# 2. LOAD & COMBINE PREDICTIONS
# ==============================================================================
train_preds <- read.csv(paste0(OUTPUT_DIR, "/train_predictions_", SCENARIO_NAME, ".csv"),
                        stringsAsFactors = FALSE)
test_preds  <- read.csv(paste0(OUTPUT_DIR, "/test_predictions_",  SCENARIO_NAME, ".csv"),
                        stringsAsFactors = FALSE)

train_preds$Split <- "train"
test_preds$Split  <- "test"

full_preds <- rbind(train_preds, test_preds)

cat(sprintf("Train: %d rows  |  Test: %d rows  |  Total: %d rows\n",
            nrow(train_preds), nrow(test_preds), nrow(full_preds)))

# ==============================================================================
# 3. LOAD TRUE WIN PROBABILITIES  (optional)
# ==============================================================================
has_true_wp <- file.exists(WP_FILE)

if (has_true_wp) {
  wp_df <- read.csv(WP_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(wp_df)[1] <- "Rating_Period"
  
  if (nrow(wp_df) != nrow(full_preds)) {
    stop(sprintf(
      "Row mismatch: full_preds has %d rows but '%s' has %d rows.",
      nrow(full_preds), WP_FILE, nrow(wp_df)))
  }
  pair_labels <- colnames(wp_df)[-1]
  cat(sprintf("True WP loaded. Pairs: %d\n", length(pair_labels)))
} else {
  cat("No true WP file found — true WP series will be omitted from plots.\n")
  players_all <- sort(unique(c(full_preds$Player1, full_preds$Player2)))
  pair_labels <- unlist(lapply(seq_len(length(players_all) - 1), function(i)
    paste0(players_all[i], players_all[seq(i + 1, length(players_all))])))
}

all_seasons  <- sort(unique(full_preds$Season))
train_cutoff <- max(full_preds$Season[full_preds$Split == "train"])

cat(sprintf("Periods: %d  |  Pairs: %d  |  Train cutoff: period %d\n",
            length(all_seasons), length(pair_labels), train_cutoff))

if (has_true_wp) {
  full_preds$True_pred <- NA_real_
  for (i in 1:nrow(full_preds)) {
    sortXY <- sort(c(full_preds$Player1[i], full_preds$Player2[i]))
    col_name <- paste0(sortXY[1], sortXY[2])
    full_preds$True_pred[i] <- wp_df[i, col_name]
  }
}

full_pred_path <- paste0(OUTPUT_DIR, "/full_pred_", SCENARIO_NAME, ".csv")
write.csv(full_preds, full_pred_path)

# ===
# 3.5. Brier Score and breakdown
# ===

rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2, na.rm = TRUE))

library(verification)

result_elo <- brier(full_preds$Game_Result, full_preds$Elo_pred)
result_pelo <- brier(full_preds$Game_Result, full_preds$PElo_pred)
result_glick <- brier(full_preds$Game_Result, full_preds$Glicko_pred)
result_pglick <- brier(full_preds$Game_Result, full_preds$PGlicko_pred)
# result_oracle <- brier(full_preds$Game_Result, full_preds$True_pred)

brier_result <- data.frame(
  row.names = c("Brier Score", "Reliability (Lower is better)", "Resolution (Higher is better)", "Uncertainty (Higher = harder to predict)"),
  Elo = c(
    result_elo$bs,
    result_elo$bs.reliability,
    result_elo$bs.resol,
    result_elo$bs.uncert
  ),
  PElo = c(
    result_pelo$bs,
    result_pelo$bs.reliability,
    result_pelo$bs.resol,
    result_pelo$bs.uncert
  ),
  Glicko = c(
    result_glick$bs,
    result_glick$bs.reliability,
    result_glick$bs.resol,
    result_glick$bs.uncert
  ),
  PGlicko = c(
    result_pglick$bs,
    result_pglick$bs.reliability,
    result_pglick$bs.resol,
    result_pglick$bs.uncert
  )
)

if(has_true_wp){
  new_row <- data.frame(Elo = rmse(full_preds$Elo_pred, full_preds$True_pred), 
                        PElo = rmse(full_preds$PElo_pred, full_preds$True_pred), 
                        Glicko = rmse(full_preds$Glicko_pred, full_preds$True_pred), 
                        PGlicko = rmse(full_preds$PGlicko_pred, full_preds$True_pred), 
                        row.names = "allgamepred_RMSE")
  brier_result <- rbind(brier_result, new_row)
}
rm(new_row)

brier_result_path <- paste0(OUTPUT_DIR, "/brier_result_", SCENARIO_NAME, ".csv")
write.csv(brier_result, brier_result_path)

plot_pair <- function(A, B, full_pred_data){
  df <- full_pred_data[(full_pred_data$Player1 %in% c(A, B))&(full_pred_data$Player2 %in% c(A, B)),]
  for(i in 1:nrow(df)){
    if(all(sort(c(df$Player1[i], df$Player2[i]))!=c(df$Player1[i], df$Player2[i]))){
      p1 <- df$Player1[i]
      p2 <- df$Player2[i]
      df$Player1[i] <- p2
      df$Player2[i] <- p1
      df$Game_Result[i] <- 1-df$Game_Result[i]
      df$Elo_pred[i] <- 1-df$Elo_pred[i]
      df$PElo_pred[i] <- 1-df$PElo_pred[i]
      df$Glicko_pred[i] <- 1-df$Glicko_pred[i]
      df$PGlicko_pred[i] <- 1-df$PGlicko_pred[i]
      df$True_pred[i] <- 1-df$True_pred[i]
    }
  }
  y_range <- range(
    df$True_pred, df$Elo_pred, df$PElo_pred,
    df$Glicko_pred, df$PGlicko_pred,
    na.rm = TRUE
  )
  x <- 1:nrow(df)
  plot(x, df$True_pred,
       type  = "l", col = "black", lwd = 2,
       ylim  = c(y_range[1],y_range[2]+0.2),
       xlab  = "Observation", ylab = "Predicted Probability",
       main  = paste0("Player ", A, " vs ", B),
       cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.7)
  
  lines(x, df$Elo_pred,     col = "red",       lwd = 1, lty = 5)
  lines(x, df$PElo_pred,    col = "red",     lwd = 2, lty = 1)
  lines(x, df$Glicko_pred,  col = "royalblue",  lwd = 1, lty = 5)
  lines(x, df$PGlicko_pred, col = "royalblue", lwd = 2, lty = 1)
  
  legend("topright",
         legend = c("True", "Elo", "P-Elo", "Glicko", "P-Glicko"),
         col    = c("black", "red", "tomato", "royalblue", "dodgerblue"),
         lwd    = c(2, 1, 2, 1, 2),
         lty    = c(1, 5, 1, 5, 1),
         bty    = "n", cex = 0.6)
}

pdf(paste0(OUTPUT_DIR, "/probability_recovery_plots.pdf"), width = 16, height = 12)
par(mfrow = c(4, 4), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.5, 0))

for (pair in pair_labels) {
  plot_pair(substr(pair,1,1), substr(pair,2,2), full_preds)
}

dev.off()
