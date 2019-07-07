#' To examine the lower level designs in forward selection

#' @return An object containing the information of the optimal PM design, with highest power for testing the
#' focal parameters, compared with other PM designs
#' @seealso \code{\link{simPM}} which wraps this function
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export pres.cost
#' @examples



pres.cost=function(opt.pattern,   # the obtained opt.pattern from the forward.opt command
                   costmx,        # should be exactly the same as
                   max.mk,
                   pc,
                   pd,
                   n,
                   k,
                   Time,
                   Time.complete){
  if (pd==0){
    if ((max.mk+1)==nrow(opt.pattern)){
      usematrix=opt.pattern
    }else{
      usematrix=opt.pattern[-(nrow(opt.pattern)-(max.mk+1)),]}
    ms.range=c((Time.complete*k+1):(Time*k))

    cost=sum(c(rep(n*(1-pc)/max.mk,max.mk),pc*n)*((1-usematrix[,ms.range])%*%costmx))
    opt.probs=c(rep(round((1-pc)/(nrow(usematrix)-1),6),max.mk),pc)
  }
  if(pd!=0){
    if ((max.mk+2)==nrow(opt.pattern)){
      usematrix=opt.pattern
    }else{
      usematrix=opt.pattern[-(nrow(opt.pattern)-(max.mk+2)),]}
    ms.range=c((Time.complete*k+1):(Time*k))
    cost=sum(c(rep(n*(1-pc-pd)/max.mk,max.mk),pc*n,pd*n)*((1-usematrix[,ms.range])%*%costmx))
    opt.probs=c(rep(round((1-pc-pd)/(nrow(usematrix)-2),6),max.mk),pc,pd)
  }
  obj=list("design.matrix"=usematrix,"probs"=opt.probs,"cost"=cost)

  return(obj)
}

## Example
# pres.cost(opt.pattern = test$opt.pattern,
#           costmx=rep(5,6),
#           max.mk=1,
#           pc=0.2,
#           pd=0.1,
#           n=200,
#           k=3,
#           Time=3,
#           Time.complete=1)


# get the simulation details of the design
# user may want to check on lower level options
# as all the output files have been saved in the working directory
# they can have a quick check by looking at the output file


### add wave missing simulation functions.
### to be done.


### user defined proportions
# seems difficult for users though if the possible number of missig is large...
# because for every level, they need to specify the prob vector

### potential issues with Mplus not being able to read lines longer than 90 characters




