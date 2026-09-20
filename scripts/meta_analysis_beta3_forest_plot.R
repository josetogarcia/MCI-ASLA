# Block 3 (prognostic): resumen descriptivo, NO un meta-análisis.
# k=2, sin CI verificable en uno de los dos estudios -> no se calcula ningún
# modelo de efectos aleatorios ni se reportan I²/tau². Este script solo grafica
# los dos estudios y un marcador descriptivo no ponderado (dado, no derivado aquí).
#
# NOTA de columnas: este CSV usa el mismo esquema que Block 1/2
# (study, role, population, modality, n_total, n_pos, n_neg, auc, ci_lower,
# ci_upper, se_from_ci, source_note); acá "auc" es en realidad el F1.

dat <- read.csv("data/meta_analysis_beta3_data.csv", stringsAsFactors = FALSE)
stopifnot(nrow(dat) == 2)
print(dat[, c("study", "modality", "n_total", "auc", "ci_lower", "ci_upper")])

## --- Figure S3 ---
## forest() de metafor descarta filas donde no puede calcular varianza (CI=NA),
## así que Mirheidari (sin CI reportado) nunca se dibujaba con esa función.
## Acá se arma a mano con gráficos base, con control total por fila.
avg_est <- 0.7918
avg_lo  <- 0.7591
avg_hi  <- 0.8244

xlim <- range(c(dat$auc, dat$ci_lower, dat$ci_upper, avg_lo, avg_hi), na.rm = TRUE) + c(-0.05, 0.05)
study_rows <- nrow(dat):1  # de arriba hacia abajo, en el orden del CSV

plot(NA, xlim = xlim, ylim = c(-1.5, nrow(dat) + 1), yaxt = "n",
     xlab = "F1 (prognostic, single modality)", ylab = "", bty = "n")

for (i in seq_len(nrow(dat))) {
  y <- study_rows[i]
  has_ci <- !is.na(dat$ci_lower[i]) && !is.na(dat$ci_upper[i])
  if (has_ci) segments(dat$ci_lower[i], y, dat$ci_upper[i], y)
  points(dat$auc[i], y, pch = 15, cex = 1.3)
}

axis(2, at = study_rows, labels = paste0(dat$study, ifelse(is.na(dat$ci_lower), " *", "")),
     las = 1, tick = FALSE, cex.axis = 0.8, hadj = 1)
mtext("* 95% CI not reported", side = 1, line = 4, adj = 0, cex = 0.7)

## --- Marcador gris con trama: promedio descriptivo no ponderado (k=2) ---
## Tomado tal cual del reporte (no se recalcula acá) -- ver nota sobre la
## discrepancia con el promedio aritmético simple (0.772) antes de publicarlo.
row_avg <- 0
polygon(x = c(avg_lo, avg_est, avg_hi, avg_est),
        y = row_avg + c(0, 0.3, 0, -0.3),
        col = "grey85", border = "grey30", density = 20, angle = 45)

text(xlim[1], row_avg, pos = 4, cex = 0.8,
     "Descriptive average (unweighted, k=2) -- indicative only, not a formal meta-analysis")

## --- Tabla mínima ---
fmt_ci <- function(est, lo, hi) {
  ifelse(is.na(lo), sprintf("%.3f (n.r.)", est), sprintf("%.3f (%.3f-%.3f)", est, lo, hi))
}

summary_table <- data.frame(
  Study    = c(dat$study, "Descriptive average (unweighted model)"),
  Modality = c(dat$modality, NA),
  N_total  = c(dat$n_total, 2L),
  F1_CI    = c(fmt_ci(dat$auc, dat$ci_lower, dat$ci_upper),
               sprintf("%.4f (%.4f-%.4f)", avg_est, avg_lo, avg_hi)),
  Note     = c(dat$source_note, "Indicative only -- not a formal meta-analysis")
)
print(summary_table, row.names = FALSE)
