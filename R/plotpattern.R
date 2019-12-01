#' Plot the missing data patterns for the optimal PM design.
#' 
#' \code{plotPM} plots the optimal PM design's missing patterns as a
#'   heatmap.
#' 
#' @param object The simPM object.
#' @param colbr Specify the colors for different waves. Default is "PRGn".
#' @param col Specify the colors for complete vs. missing data. 
#'    Default is c("antiquewhite1","firebrick").
#' @param labels logical scalar, indicating whether the label for waves is
#'    needed. Default is TRUE.
#' @param fontsize_col Specify the font size for the column labels. Default
#'    is 20.
#' @param fontsize_row Specify the font size for the row labels. Default is
#'    14.
#' @param fontsize Specify the font size for the legend. Default is 14.
#' @param angle_col Specify the angle of how the column labels are
#'    displayed. Default is 45. 
#' @param legend Logical scalar, indicating whether the legend is shown.
#'    Default is TRUE.
#' @param main Specify the plot title. 
#' @param ... Any additional arguments for \code{\link[pheatmap]{pheatmap}}.
#'
#' @import pheatmap
#' @import RColorBrewer
#' @export plotPM
#' 
#' @seealso \code{\link[pheatmap]{pheatmap}}
#' 
#' @examples
#' \dontrun{
#' plotPM(wave.out)
#' plotPM(indicator.out)
#' plotPM(forward.out,labels=F,col=c("gray96","gray35"),fontsize_row=26,fontsize=18,fontsize_col=26)
#' }

plotPM <- function(object,
                colbr="PRGn",
                col=c("antiquewhite1","firebrick"),
                row.names=T,
                labels=T,
                fontsize_col=20,
                fontsize_row=14,
                fontsize=14,
                angle_col=45,
                legend=T,
                main="",
                ...) {
  data <- object$opt.pattern
  
  Time <- object$misc$time
  k <- object$misc$k
  
  if (row.names==T) {
    row.names(data) <- paste0(paste0("pat", seq_len(nrow(data)),":",sep=""), "n=", object$opt.ns)
  }
  annotation <- data.frame(Wave = factor(rep(c(1:Time), each = k), labels = paste0("W",1:Time)))
  rownames(annotation) <- colnames(data)
  mat_colors <- list(Wave = brewer.pal(Time, colbr))
  names(mat_colors$Wave) <- paste0("W",1:Time)
  
  gaps <- seq(k, (Time-1)*k, by = k)
  
  if (labels==T) {
    pheatmap(data, scale = "none",col=col,
             cluster_rows = F, cluster_cols = F,legend = legend, legend_breaks = c(0,1), legend_labels = c("complete", "missing"), fontsize_row = fontsize_row, fontsize_col = fontsize_col, fontsize = fontsize, drop_levels = T, annotation = annotation, annotation_colors = mat_colors, angle_col = angle_col, gaps_col = gaps, main = main, ...)
  }
  
  if (labels==F) {
    pheatmap(data, scale = "none",col=col,
             cluster_rows = F, cluster_cols = F, legend = legend, legend_breaks = c(0, 1), fontsize_row = fontsize_row, legend_labels = c("complete", "missing"), fontsize_col = fontsize_col, fontsize = fontsize, drop_levels = T, angle_col = angle_col, gaps_col = gaps, main = main)
  }
}
