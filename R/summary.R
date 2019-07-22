#' A summary function to extract the important information of the output object.


#' @return An object containing the information of the optimal PM design, with highest power for testing the
#' focal parameters, compared with other PM designs
#' 
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export summary.opt
#' @examples



summary.opt=function(object){
  
  focal.param=object$misc$focal.param
  opt.powers=summaryParam(object$opt.output)[focal.param,]
  
  
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

#summary.opt(re.ob)
