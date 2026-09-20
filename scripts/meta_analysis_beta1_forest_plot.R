library(metafor)

## --- SE por Hanley & McNeil (1982), para AUCs sin CI reportado ---
hanley_mcneil_se <- function(auc, n_pos, n_neg) {
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc^2) / (1 + auc)
  sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) + (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg))
}

# Ajusta esta ruta al CSV correcto -- el código original apuntaba a dos archivos
# distintos ("beta1.description_data.csv" y "beta1_data.csv"); aquí se asume uno solo.
dat <- read.csv("data/meta_analysis_beta1_data.csv", stringsAsFactors = FALSE)

dat$se_calc <- ifelse(
  !is.na(dat$ci_lower) & !is.na(dat$ci_upper),
  (dat$ci_upper - dat$ci_lower) / 3.92,
  hanley_mcneil_se(dat$auc, dat$n_pos, dat$n_neg)
)

## --- Modelo global (diamante rojo) ---
m_overall <- rma(yi = auc, sei = se_calc, data = dat, method = "DL")

## --- Orden y filas para el forest plot con subgrupos ---
# Categorías con >=2 estudios se pooled (diamante azul); con 1 estudio
# (Metarugcheep_2022: phonemic fluency; Wang_2023: motor/DDK) se muestran
# solas, sin diamante -- exactamente lo que dice el pie de figura.
task_order <- c("Free speech", "Incluye descripcion", "phonemic fluency", "motor/DDK")
dat$task_group <- factor(dat$task_group, levels = task_order)
dat <- dat[order(dat$task_group), ]

n_by_group    <- table(dat$task_group)
pooled_groups <- names(n_by_group[n_by_group >= 2])

# Recorre las categorías de abajo hacia arriba asignando filas, dejando un
# hueco extra después de cada categoría que sí se pooled (para su diamante)
rows <- integer(nrow(dat))
cursor <- 1
subgroup_models <- list()
subgroup_row <- list()

for (g in rev(task_order)) {
  idx <- which(dat$task_group == g)
  if (length(idx) == 0) next
  rows[idx] <- cursor + rev(seq_along(idx)) - 1
  cursor <- cursor + length(idx)
  if (g %in% pooled_groups) {
    subgroup_models[[g]] <- rma(yi = auc, sei = se_calc, data = dat[idx, ], method = "DL")
    subgroup_row[[g]] <- cursor
    cursor <- cursor + 1
  }
  cursor <- cursor + 1
}

overall_row <- -1.5  # separado del resto, debajo de todas las categorías

## --- Forest plot ---
forest(dat$auc, sei = dat$se_calc, slab = dat$study,
       rows = rows, ylim = c(overall_row - 1.5, max(rows) + 3),
       xlab = "AUC (MCI vs. CU/HC)", refline = 0.5,
       psize = 1.1, header = "Study / Task category", mlab = "")

# Diamantes azules: subgrupos pooled
for (g in names(subgroup_models)) {
  addpoly(subgroup_models[[g]], row = subgroup_row[[g]],
          mlab = paste0("Subgroup pooled: ", g), col = "blue")
}

# Diamante rojo: estimado global (incluye las categorías de un solo estudio)
addpoly(m_overall, row = overall_row, mlab = "Overall (DL)", col = "red")
