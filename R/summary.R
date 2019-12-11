#' Extract the summary information of the simpm object.
#'
#' \code{summary.simpm} summarizes and extracts the important information about the optimal PM design.
#' 
#' @param object The output object returned by \code{\link{simPM}}.
#' @return The information about the optimal PM design.
#' 
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @method summary simpm
#' @export summary.simpm
#' @examples



summary.simpm <- function(object) {
  
  if (length(object) == 2) {
    focal.param <- object[[1]]$misc$focal.param
    opt.powers <- summaryParam(object[[1]]$opt.output)[focal.param, ]
    
    print("=================results summary================")
    print(object[[1]]$results)
    print("=================Optimal design=================")
    print(object[[1]]$opt.design)
    print("=================Optimal design for focal parameters=================")
    print(opt.powers)
    print("=================Optimal patterns===============")
    print(object[[1]]$opt.pattern)
    print("=================Optimal probs==================")
    print(object[[1]]$opt.probs)
    print("=================Optimal ns====================")
    print(object[[1]]$opt.ns)
    
  }else{
    focal.param <- object$misc$focal.param
    opt.powers <- summaryParam(object$opt.output)[focal.param, ]
    
    print("=================results summary================")
    print(object$results)
    print("=================Optimal design=================")
    print(object$opt.design)
    print("=================Optimal design for focal parameters=================")
    print(opt.powers)
    print("=================Optimal patterns===============")
    print(object$opt.pattern)
    print("=================Optimal probs==================")
    print(object$opt.probs)
    print("=================Optimal ns====================")
    print(object$opt.ns)
  }  
  
}

