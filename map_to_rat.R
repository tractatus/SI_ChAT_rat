to.rat.stereotactic(dataset){
  dataset$ML<-dataset$ML*1.56
  dataset$DV<-dataset$DV*1.45
  return(dataset)
}

map.to.rat<-function(dataset, regi, col = NULL, ...){
  if(is.null(col)){
    col<-dataset$color
  }
  back.warp<-stereotactic.coordinates(dataset$ML, dataset$DV, regi, inverse =TRUE)
  
  scale.factor <- mean(dim(regi$transformationgrid$mx)/c(regi$transformationgrid$height, 
                                                         regi$transformationgrid$width))
  xMax <- max(c(regi$transformationgrid$mx, regi$transformationgrid$mxF), 
              na.rm = TRUE) * (1/scale.factor)
  xMin <- min(c(regi$transformationgrid$mx, regi$transformationgrid$mxF), 
              na.rm = TRUE) * (1/scale.factor)
  yMax <- max(c(regi$transformationgrid$my, regi$transformationgrid$myF), 
              na.rm = TRUE) * (1/scale.factor)
  yMin <- min(c(regi$transformationgrid$my, regi$transformationgrid$myF), 
              na.rm = TRUE) * (1/scale.factor)
  
  plot(c(xMin, xMax), c(yMin, yMax), ylim = c(yMax, yMin), 
       xlim = c(xMin, xMax), asp = 1, axes = F, xlab = "", ylab = "", col=0)
  
  line.type <- c(3, 1)[as.integer(regi$atlas$col == "#cccccc") + 1]
  line.type[1]<-1
  line.width<-rep(1, length(line.type))
  line.width[1]<-3
  fill.col <- c(gray(0.95), gray(0.85))[as.integer(regi$atlas$col == "#cccccc") + 1]
  fill.col[which(regi$atlas$col == "#aaaaaa")]<-"black"
  
  lapply(1:regi$atlas$numRegions, function(x) {
    polygon(regi$atlas$outlines[[x]]$xrT/scale.factor, 
            regi$atlas$outlines[[x]]$yrT/scale.factor, 
            border = 'black', lty = line.type[x], col = fill.col[x], lwd = line.width[x])
  })
  lapply(1:regi$atlas$numRegions, function(x) {
    polygon(regi$atlas$outlines[[x]]$xlT/scale.factor, 
            regi$atlas$outlines[[x]]$ylT/scale.factor, 
            border = 'black', lty = line.type[x], col = fill.col[x], lwd = line.width[x] )
  })
  
  # lapply(1:regi$atlas$numRegions, function(x) {
  #   polygon(regi$atlas$outlines[[x]]$xr/scale.factor, 
  #           regi$atlas$outlines[[x]]$yr/scale.factor, 
  #           border = 'orange', col='black', lty = line.type[x],  lwd = line.width[x])
  # })
  # lapply(1:regi$atlas$numRegions, function(x) {
  #   polygon(regi$atlas$outlines[[x]]$xl/scale.factor, 
  #           regi$atlas$outlines[[x]]$yl/scale.factor, 
  #           border = 'orange', lty = line.type[x], col='black', lwd = line.width[x] )
  # })
  
  
  index <- cbind(as.integer(round(back.warp$y)), as.integer(round(back.warp$x)))
  #index[,2]<- (index[,2]+ regi$centroidNorm[2])
  #index[,1]<- (index[,1]- regi$centroidNorm[1])
  index<-index*scale.factor
  
  remove<-unique(c(which(index[,1]<0 | index[,1]>dim(regi$transformationgrid$mx)[1] ), which(index[,2]<0 | index[,2]> dim(regi$transformationgrid$mx)[1] )))
  if(length(remove)>0){
    index<-index[-remove,]
    col<-col[-remove]
  }
  #index<-index*scale.factor
  
  xT <- (index[,2] + (regi$transformationgrid$mx[index] - index[,2]))
  yT <- (index[,1] + (regi$transformationgrid$my[index] - index[,1]))
  
  points(xT/scale.factor, yT/scale.factor, pch=16, col= col, ...)
  
}
