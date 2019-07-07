#' To check the simulation results for lower level designs in forward selection

#' @return An object containing the information of the optimal PM design, with highest power for testing the
#' focal parameters, compared with other PM designs
#' @seealso \code{\link{simPM}} which wraps this function
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export pres.design
#' @examples



pres.design=function(
  opt.results,   #searching result from forward.opt
  max.mk         #maximum of missing slots within a pattern allowed
){
  sim.seq=opt.results$design.order[max.mk]
  FNAME=paste0("missing-",num.miss,"-sim.seq-",sim.seq,".out")
  obj=readModels(FNAME)
  return.obj=list("Filename"=FNAME,"parameters"=obj$parameters)
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




