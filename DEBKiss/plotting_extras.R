#HDPI figure code
# rethinking is not yet on CRAN but can be installed using devtools:

## install.packages(c("devtools","mvtnorm","loo","coda"), repos="https://cloud.r-project.org/",dependencies=TRUE)
## library(devtools)
## install_github("rmcelreath/rethinking")

## The code for the figure itself is 

library(coda)
library(rethinking)

pretty_posterior_plot <- function(samples, summary, ref.params, param, HDPI = 0.95, legend = TRUE, deBInfer_result = NULL){
  rethinking::dens(unlist(samples[,param]), show.HPDI = 0.95, main = param)
  abline(v = ref.params[param], lwd = 2, lty = 2)
  abline(v = summary$statistics[param, "Mean"])
  if(!is.null(deBInfer_result)){
    dprior <- paste("d", deBInfer_result$all.params[[param]]$prior, 
                    sep = "")
    plot.range <- seq(par("usr")[1], par("usr")[2], length.out = 100)
    prior.dens <- do.call(dprior, c(list(x = plot.range), 
                                    deBInfer_result$all.params[[param]]$hypers))
    lines(plot.range, prior.dens, col = 'red', lty = 3, lwd = 1.5)
  }
  if(legend) legend("topleft", legend = c("true parameter value", "posterior mean value", "95% HPDI"), lty = c(2,1, NA), pch = c(NA,NA,15), col = c("black", "black", rethinking::col.alpha("black",0.15)))
}

### where samples is a coda object, i.e. the samples slot of a debinfer_result object (i.e. debinfer_result$samples)

### summary is created using the summary method for the coda object (i.e. summary(debinfer_result$samples))

### ref.params is a named vector with the true values

### param is a string of the parameter to be plotted

#pdf("figs/earthworm_fixedEG_freeobs_estimate_vs_true.pdf")
#par(mfrow = c(3,3))
#for (pp in names(freeparams(fixedEG_freeobs1))[1:8]){
#    pretty_posterior_plot(fixedEG_freeobs_samples,
#                          fixedEG_freeobs_summary,
#                          ref.params = c(params_Lumter,
#                              c(sdlog.EWw = 0.05, sdlog.R = 0.1)),
#                          param = pp, legend = FALSE)
#}
#plot.new()
#legend("topleft",
#       legend = c("true parameter value", "posterior mean value", "95% HPDI"),
#       lty = c(2,1, NA), pch = c(NA,NA,15),
#       col = c("black", "black", rethinking::col.alpha("black",0.15)))
#dev.off()
