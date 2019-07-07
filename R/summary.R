#' A summary function to extract the important information of the output object.


#' @return An object containing the information of the optimal PM design, with highest power for testing the
#' focal parameters, compared with other PM designs
#' @seealso \code{\link{simPM}} which wraps this function
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export summary.opt
#' @examples



summary.opt=function(object){
  print("=================results summary================")
  print(object$results)
  print("=================Optimal design=================")
  print(object$opt.design)
  print("=================Optimal patterns===============")
  print(object$opt.pattern)
  print("=================Optimal probs==================")
  print(object$opt.probs)
  print("=================Optimal ns====================")
  print(object$opt.ns)
}

#summary.opt(re.ob)
