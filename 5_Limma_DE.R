rm(list = ls())
library(limma)


setwd('/mnt/8tb/DARPA/Proteomics_R2')

tissues <- c('Plasma')
for (t in tissues)
{
  print(t)
  design <- read.csv(paste0('./limma/Cod_', t, '.csv'), row.names = 1, header = 1)
  exp_matrix <- read.csv(paste0('./limma/Exp_', t, '.csv'), row.names = 1, header = 1)

  cont_matrix <- makeContrasts(acute11p = acute11p - Normoxia,
                               acute8p = acute8p - Normoxia,
                               chronic11p = chronic11p - Normoxia,
                               chronic8p = chronic8p - Normoxia,
                               precondition11p = precondition11p - Normoxia,
                               precondition8p = precondition8p - Normoxia,
                               precondition8p11p = precondition8p - precondition11p,
                               precondition8pacute8p = precondition8p - acute8p,
                               levels = design)

  # Fit the expression matrix to a linear model
  fit <- lmFit(exp_matrix, design)
  # Compute contrast
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  # Bayes statistics of differential expression
  # *There are several options to tweak!*
  fit_contrast <- eBayes(fit_contrast)

  for (c in c("acute11p", "acute8p", "chronic11p", "chronic8p", "precondition11p", "precondition8p", "precondition8p11p", "precondition8pacute8p"))
  {
    result <- topTable(fit_contrast, coef = c, n = Inf, adjust = "BH")
    write.csv(result, paste0('./DE/', t, '/limma_', c, '.csv'))
  }
}
