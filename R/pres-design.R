#' Check the simulation results for lower level designs in forward assembly. 
#' \code{pres.design} is only applicable when Mplus is used for simulations.
#' 
#' @param opt.results The object returned from forward assembly.
#' @param max.mk The maximum number of unique missing data
#'   patterns in the selected design. 
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export pres.design
#' @examples



pres.design <- function(
  opt.results,   #searching result from forward.opt
  max.mk         #maximum of missing slots within a pattern allowed
){
  sim.seq <- opt.results$design.order[max.mk]
  FNAME <- paste0("missing-", num.miss, "-sim.seq-", sim.seq, ".out")
  obj <- readModels(FNAME)
  return.obj <- list("Filename" = FNAME, "parameters" = obj$parameters)
  return(return.obj)
}

# Example
#pres.design(test,1)


### add wave missing simulation functions.
### to be done.


### user defined proportions
# seems difficult for users though if the possible number of missig is large...
# because for every level, they need to specify the prob vector

### potential issues with Mplus not being able to read lines longer than 90 characters




