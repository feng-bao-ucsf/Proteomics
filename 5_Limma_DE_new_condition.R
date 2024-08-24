rm(list = ls())
library(limma)


setwd('/mnt/8tb/DARPA/Proteomics_R2')

tissues <- c('Kidney', 'Liver', 'Muscle', 'Heart', 'Brain', 'Plasma')
for (t in tissues)
{
print(t)
  design <- read.csv(paste0('./limma/Cod_', t, '.csv'), row.names = 1, header = 1)
  exp_matrix <- read.csv(paste0('./limma/Exp_', t, '.csv'), row.names = 1, header = 1)

  cont_matrix <- makeContrasts(acute8pvsacute11p = acute8p - acute11p,
                               chronic8pvschronic11p = chronic8p - chronic11p,
                               acute8pvschronic8p = acute8p - chronic8p,
                               acute11pvschronic11p= acute11p - chronic11p,
                               levels = design)

  # Fit the expression matrix to a linear model
  fit <- lmFit(exp_matrix, design)
  # Compute contrast
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  # Bayes statistics of differential expression
  # *There are several options to tweak!*
  fit_contrast <- eBayes(fit_contrast)

  for (c in c("acute8pvsacute11p", "chronic8pvschronic11p", "acute8pvschronic8p", "acute11pvschronic11p"))
  {

  print(c)
    result <- topTable(fit_contrast, coef = c, n = Inf, adjust = "BH")
    write.csv(result, paste0('./DE/', t, '/limma_', c, '.csv'))
  }
}
