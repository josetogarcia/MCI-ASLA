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

## --- Análisis primario (k=3) ---
primary <- subset(dat, role == "primary")
m_primary <- rma(yi = auc, sei = se_from_ci, data = primary, method = "DL")
summary(m_primary)

## --- Sensitivity 1 (k=4): + Hajjar acústico, SE vía Hanley-McNeil ---
hajjar_acoustic <- subset(dat, study == "Hajjar_2023" & role == "sensitivity_acoustic")
hajjar_acoustic$se_from_ci <- hanley_mcneil_se(hajjar_acoustic$auc, hajjar_acoustic$n_pos, hajjar_acoustic$n_neg)

sens1 <- rbind(primary, hajjar_acoustic)
m_sens1 <- rma(yi = auc, sei = se_from_ci, data = sens1, method = "DL")
summary(m_sens1)

## --- Sensitivity 2: Hajjar NLP como punto narrativo, sin peso ---
# Split Aβ+/Aβ- del subgrupo no es verificable -> no se calcula SE, no entra al pool.
hajjar_nlp <- subset(dat, study == "Hajjar_2023" & role == "sensitivity_NLP")

## --- Comparación de los tres escenarios ---
comparison <- data.frame(
  scenario   = c("Primario (k=3)", "+ Hajjar acústico (k=4, ponderado)", "+ Hajjar NLP (no ponderado, referencia)"),
  auc_pooled = c(m_primary$b[1], m_sens1$b[1], NA),
  ci_lower   = c(m_primary$ci.lb, m_sens1$ci.lb, NA),
  ci_upper   = c(m_primary$ci.ub, m_sens1$ci.ub, NA),
  I2         = c(m_primary$I2, m_sens1$I2, NA),
  tau2       = c(m_primary$tau2, m_sens1$tau2, NA),
  auc_no_pool_note = c(NA, NA, hajjar_nlp$auc)
)
print(comparison, digits = 3)

## --- Forest plot: modelo de k=4 (primario + Hajjar acústico ponderado) ---
## + Hajjar NLP superpuesto como diamante abierto, fuera del pooling.
forest(m_sens1, slab = sens1$study, xlab = "AUC", refline = NA,
       ylim = c(-2, nrow(sens1) + 3),
       mlab = "Pooled DL: primario + Hajjar acústico (k=4)")

points(hajjar_nlp$auc, -1, pch = 5, cex = 1.4)
text(hajjar_nlp$auc, -1, "Hajjar NLP (no ponderado, split Aβ+/- no verificable)",
     pos = 4, cex = 0.75, offset = 1)
