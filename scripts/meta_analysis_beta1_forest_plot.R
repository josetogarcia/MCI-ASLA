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
# Categorías con >=2 estudios se pooled (diamante azul); con 1 estudio se
# muestran solas, sin diamante. El orden de las categorías (y sus nombres)
# se toma directo de los datos -- NO se hardcodea ningún string de
# task_group, para no depender de que coincida tildes/mayúsculas/espacios
# con lo que uno supone que dice el CSV.
cat("Categorías de task_group encontradas:\n")
print(table(dat$task_group, useNA = "ifany"))

n_by_group    <- sort(table(dat$task_group), decreasing = TRUE)  # grupos grandes primero
task_order    <- names(n_by_group)
pooled_groups <- names(n_by_group[n_by_group >= 2])
dat$task_group <- factor(dat$task_group, levels = task_order)
dat <- dat[order(dat$task_group), ]

# Recorre las categorías EN task_order (de arriba hacia abajo en la figura),
# asignando filas con un contador que solo baja -- así el diamante de cada
# subgrupo queda siempre inmediatamente bajo sus propios estudios, nunca
# encima ni pisando la fila de otra categoría.
counts          <- as.integer(n_by_group[task_order])
extra_per_group <- ifelse(task_order %in% pooled_groups, 1L, 0L)  # fila del diamante
total_rows      <- sum(counts) + sum(extra_per_group) + length(task_order)  # + huecos entre grupos

rows <- integer(nrow(dat))
subgroup_models <- list()
subgroup_row <- list()
row_ptr <- total_rows

for (g in task_order) {
  idx <- which(dat$task_group == g)
  n_g <- length(idx)
  if (n_g == 0) next
  rows[idx] <- row_ptr - seq_len(n_g) + 1
  row_ptr <- row_ptr - n_g
  if (g %in% pooled_groups) {
    subgroup_models[[g]] <- rma(yi = auc, sei = se_calc, data = dat[idx, ], method = "DL")
    subgroup_row[[g]] <- row_ptr       # el diamante va justo bajo el último estudio del grupo
    row_ptr <- row_ptr - 1
  }
  row_ptr <- row_ptr - 1               # hueco antes de la siguiente categoría
}

# Si dos estudios terminaron con la misma fila (p.ej. por un NA en task_group
# que no debería existir), mejor fallar acá con un mensaje claro que seguir
# y graficar dos labels superpuestos en silencio.
if (anyDuplicated(rows) > 0) {
  stop("Filas duplicadas en 'rows' -- revisa 'table(dat$task_group, useNA=\"ifany\")' arriba: ",
       "algún estudio puede tener NA o un valor de task_group inesperado.")
}

overall_row <- min(rows) - 3  # separado del resto, debajo de todas las categorías

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

## --- Tabla resumen mínima ---
fmt_est <- function(est, lb, ub) sprintf("%.2f [%.2f, %.2f]", est, lb, ub)

single_study_rows <- do.call(rbind, lapply(setdiff(task_order, pooled_groups), function(g) {
  r <- dat[dat$task_group == g, ]
  data.frame(
    Domain   = g,
    Model    = paste0("Single study (", r$study, ")"),
    k        = 1L,
    Estimate = fmt_est(r$auc, r$auc - 1.96 * r$se_calc, r$auc + 1.96 * r$se_calc),
    I2       = NA_character_
  )
}))

summary_table <- rbind(
  data.frame(Domain = "Block 1 (all task categories)", Model = "Overall (DL)",
             k = nrow(dat), Estimate = fmt_est(m_overall$b[1], m_overall$ci.lb, m_overall$ci.ub),
             I2 = sprintf("%.1f%%", m_overall$I2)),
  do.call(rbind, lapply(names(subgroup_models), function(g) {
    m <- subgroup_models[[g]]
    data.frame(Domain = g, Model = "Subgroup pooled (DL)", k = m$k,
               Estimate = fmt_est(m$b[1], m$ci.lb, m$ci.ub), I2 = sprintf("%.1f%%", m$I2))
  })),
  single_study_rows
)
print(summary_table, row.names = FALSE)
