#' Plot the missing data patterns for the optimal PM design
#' @param object The simPM object
#' @param Time The total number of waves
#' @param k The number of observed variables at each wave
#' @param colbr colors for waves
#' @param col colors for complete vs. missing data
#' @import pheatmap
#' @import RColorBrewer
#' @export plotPM
#' @examples
#' \dontrun{
#' plotPM(wave.out,Time=5,k=3)
#' plotPM(indicator.out,Time=5,k=3)
#' plotPM(forward.out,Time=5,k=3)
#' }

plotPM=function(object,
                Time,
                k,
                colbr="PRGn",
                col=c("antiquewhite1","firebrick"),
                row.names=T,
                labels=T){
  data=object$opt.pattern
  if (row.names==T){
    row.names(data)=paste0(paste0("pat",1:nrow(data),":",sep=""),"n=",object$opt.ns)
  }
  annotation=data.frame(Wave=factor(rep(c(1:Time),each=k), labels = paste0("W",1:Time)))
  rownames(annotation) <- colnames(data)
  mat_colors <- list(Wave = brewer.pal(Time, colbr))
  names(mat_colors$Wave) <- paste0("W",1:Time)
  
  gaps=seq(k,(Time-1)*k,by=k)
  
  if (labels==T){
    pheatmap(data, scale = "none",col=col,
             cluster_rows = F, cluster_cols = F,legend=T,legend_breaks=c(0,1),legend_labels=c("complete","missing"),fontsize_row = 14, fontsize_col=20,fontsize=14,drop_levels=T,annotation = annotation,annotation_colors = mat_colors,angle_col=45,gaps_col=gaps)
  }
  
  if (labels==F){
    pheatmap(data, scale = "none",col=col,
             cluster_rows = F, cluster_cols = F,legend=T,legend_breaks=c(0,1),fontsize_row = 14,legend_labels=c("complete","missing"), fontsize_col=20,fontsize=14,drop_levels=T,angle_col=45,gaps_col=gaps)
  }
}
