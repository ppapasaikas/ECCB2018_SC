require("plyr")

#' This is the function for selection of overdispersed genes adapted from:
#' https://github.com/10XGenomics/single-cell-3prime-paper/
#'
#' @param m  a (protentially sparse) gene x cells count matrix
#' @return a vector of normalized (robust Z-scores) dispersion values, one per gene.
select_variable_genes<-function(m) {
  df<-data.frame(mean=rowMeans(m+1/ncol(m)),cv=apply(m,1,sd)/rowMeans(m+1/ncol(m)),var=apply(m,1,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,unique(quantile(mean,seq(0.1,1,0.05),na.rm=TRUE) )  ,Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,(dispersion-bin_disp_median)/(bin_disp_mad+0.01) )
  return(df$dispersion_norm)
}