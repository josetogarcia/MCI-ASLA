library(metafor)

# 1. Función Hanley & McNeil (1982)
hanley_mcneil_se <- function(auc, n_pos, n_neg) {
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc^2) / (1 + auc)
  sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) + (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg))
}

# 2. Carga directa desde tu archivo CSV
dat <- read.csv("data/meta_analysis_beta1_data.csv", stringsAsFactors = FALSE)

# 3. Derivar SE según la información disponible en las columnas
# Prioriza el IC 95% si existe; de lo contrario, usa Hanley & McNeil con n_pos y n_neg
dat$se_calc <- ifelse(
  !is.na(dat$ci_lower) & !is.na(dat$ci_upper),
  (dat$ci_upper - dat$ci_lower) / 3.92,
  hanley_mcneil_se(dat$auc, dat$n_pos, dat$n_neg)
)

# 4. Ajuste del modelo primario (DerSimonian-Laird) -- referencia I² = 96.78%
m_primary <- rma(yi = auc, sei = se_calc, data = dat, method = "DL")
summary(m_primary)

# 5. Forest plot
forest(m_primary,
       slab = dat$study,
       xlab = "AUC (MCI vs CU)",
       refline = 0.5,
       mlab = "Pooled AUC (DL)")

## --- Moderador: task_group, filtrado a los dos grupos comparables (k=9) ---
sub_task <- subset(dat, task_group %in% c("Free speech", "Incluye descripcion"))

# mods = ~ task_group - 1: sin intercepto -> cada coeficiente es directamente
# la media de AUC de ese grupo (no un contraste contra un grupo de referencia)
m_task <- rma(yi = auc, sei = se_calc, mods = ~ task_group - 1,
              method = "REML", data = sub_task)
summary(m_task)

# Medias por grupo + IC95%
coef(summary(m_task))[, c("estimate", "ci.lb", "ci.ub", "pval")]

# OJO con el QM de este modelo: al no tener intercepto, prueba H0: ambas medias = 0
# simultáneamente -- eso es trivial (ningún AUC es 0) y NO es un test de si los
# grupos difieren entre sí. El test correcto para "¿el moderador explica
# heterogeneidad significativa?" es comparar este modelo contra el nulo ajustado
# en el MISMO subset de k=9 (no contra el m_primary de arriba, que usa k distinto).
m_null_subset <- rma(yi = auc, sei = se_calc, method = "REML", data = sub_task)
anova(m_null_subset, m_task)

# I² residual del modelo con moderador vs. I² del modelo nulo global (96.78%)
cat("I² nulo (modelo completo):        ", round(m_primary$I2, 2), "%\n")
cat("I² nulo (mismo subset k=9):        ", round(m_null_subset$I2, 2), "%\n")
cat("I² residual (moderador task_group):", round(m_task$I2, 2), "%\n")
