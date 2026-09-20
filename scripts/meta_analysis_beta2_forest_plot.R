# Meta-análisis ASLA (beta-amiloide): AUC pooled, efectos aleatorios (DerSimonian-Laird)
# Escala AUC directa -- ningún valor está cerca de 0/1, no se justifica transformación logit.

library(metafor)

dat <- read.csv("data/meta_analysis_beta2_data.csv", stringsAsFactors = FALSE)

# SE por método Hanley-McNeil (1982) para AUCs sin CI reportado pero con n y split conocidos
hanley_mcneil_se <- function(auc, n_pos, n_neg) {
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc^2) / (1 + auc)
  sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) + (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg))
}

## --- Análisis primario (k=3): solo estudios con varianza verificable desde un IC reportado ---
primary <- subset(dat, role == "primary")
m_primary <- rma(yi = auc, sei = se_from_ci, data = primary, method = "DL")
summary(m_primary)

## --- Sensitivity (k=4): + Hajjar acústico, varianza vía Hanley-McNeil ---
hajjar_acoustic <- subset(dat, study == "Hajjar_2023" & role == "sensitivity_acoustic")
hajjar_acoustic$se_from_ci <- hanley_mcneil_se(hajjar_acoustic$auc, hajjar_acoustic$n_pos, hajjar_acoustic$n_neg)

sens <- rbind(primary, hajjar_acoustic)
m_sens <- rma(yi = auc, sei = se_from_ci, data = sens, method = "DL")
summary(m_sens)

cat("I² primario (k=3):   ", round(m_primary$I2, 1), "%\n")
cat("I² sensitivity (k=4):", round(m_sens$I2, 1), "%\n")

## --- Figure S1 ---
## Cuadrados = las k=4 filas de "sens" con su propio IC95% (se_from_ci de cada una).
## Diamante azul = pooled primario (k=3); diamante naranja = pooled sensitivity (k=4).
## Se dejan 2 filas en blanco entre los estudios y los diamantes.
k <- nrow(sens)
forest(sens$auc, sei = sens$se_from_ci, slab = sens$study,
       rows = k:1 + 2, ylim = c(-2, k + 5),
       xlab = "AUC (amyloid positivity)", refline = NA,
       psize = 1.1, header = "Study", mlab = "")

addpoly(m_primary, row = 1,    mlab = "Primary (DL, k=3)",    col = "blue")
addpoly(m_sens,    row = -0.5, mlab = "Sensitivity (DL, k=4)", col = "orange")
