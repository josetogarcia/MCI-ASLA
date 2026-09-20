# Block 3 (prognostic, F1): resumen descriptivo, NO un meta-análisis.
# k=2, sin CI verificable en uno de los dos estudios -> no se calcula ningún
# modelo de efectos aleatorios ni se reportan I²/tau². Este script solo grafica
# los dos estudios y un marcador descriptivo no ponderado (dado, no derivado aquí).

dat <- read.csv("data/meta_analysis_beta3_data.csv", stringsAsFactors = FALSE)

## --- Figure S3 ---
## Cuadrados = F1 por estudio; Mirheidari no tiene CI reportado ("n.r."), por lo
## que su fila queda sin barra de intervalo (ci_lower/ci_upper = NA).
forest(dat$f1, ci.lb = dat$ci_lower, ci.ub = dat$ci_upper,
       slab = paste0(dat$study, ifelse(is.na(dat$ci_lower), " *", "")),
       rows = c(3, 2), ylim = c(-1.5, 5),
       xlab = "F1 (prognostic, single modality)", refline = NA,
       psize = 1.1, header = "Study", mlab = "")

mtext("* 95% CI not reported", side = 1, line = 4, adj = 0, cex = 0.7)

## --- Marcador gris con trama: promedio descriptivo no ponderado (k=2) ---
## Tomado tal cual del reporte (no se recalcula acá) -- ver nota en la respuesta
## sobre la discrepancia con el promedio aritmético simple (0.772).
avg_est <- 0.7918
avg_lo  <- 0.7591
avg_hi  <- 0.8244
row_avg <- 0

polygon(x = c(avg_lo, avg_est, avg_hi, avg_est),
        y = row_avg + c(0, 0.3, 0, -0.3),
        col = "grey85", border = "grey30", density = 20, angle = 45)

text(par("usr")[1], row_avg, pos = 4, cex = 0.8,
     "Descriptive average (unweighted, k=2) -- indicative only, not a formal meta-analysis")
