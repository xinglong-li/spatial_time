# Generate (approximate) posterior samples from models
# Warning - this takes up a LOT of memory. Recommended to use HPC

library(INLA)

m_samples = 1000
mod = readRDS('finalmod_naive.rds')

# Generate approximate posterior samples
samp <- inla.posterior.sample(n=m_samples, mod)

# reduce the size of samp by removing the reduntant information #
for (i in 1:m_samples)
{
  samp[[i]]$latent = samp[[i]]$latent[-c(grep('Predictor', rownames(samp[[i]]$latent), fixed = T)),]
}

saveRDS(samp, file = 'naive_samples.rds')
