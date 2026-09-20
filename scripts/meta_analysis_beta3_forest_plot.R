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
## Cuadrados = F1 por estudio; Mirheidari no tiene CI reportado ("n.r."), por lo
## que su fila queda sin barra de intervalo (ci_lower/ci_upper = NA).
forest(dat$auc, ci.lb = dat$ci_lower, ci.ub = dat$ci_upper,
       slab = paste0(dat$study, ifelse(is.na(dat$ci_lower), " *", "")),
       rows = c(3, 2), ylim = c(-1.5, 5),
       xlab = "F1 (prognostic, single modality)", refline = NA,
       psize = 1.1, header = "Study", mlab = "")

mtext("* 95% CI not reported", side = 1, line = 4, adj = 0, cex = 0.7)

## --- Marcador gris con trama: promedio descriptivo no ponderado (k=2) ---
## Tomado tal cual del reporte (no se recalcula acá) -- ver nota sobre la
## discrepancia con el promedio aritmético simple (0.772) antes de publicarlo.
avg_est <- 0.7918
avg_lo  <- 0.7591
avg_hi  <- 0.8244
row_avg <- 0

polygon(x = c(avg_lo, avg_est, avg_hi, avg_est),
        y = row_avg + c(0, 0.3, 0, -0.3),
        col = "grey85", border = "grey30", density = 20, angle = 45)

text(par("usr")[1], row_avg, pos = 4, cex = 0.8,
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
