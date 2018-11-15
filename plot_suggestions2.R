
plot.suggestions2<-function (dataset, normalize.by = NULL, group = NULL, title = "Groups:", 
                             color = gray(0.4), exclude.below = 10, reduce.below = 100, 
                             exclude.regions = NULL, include.regions = NULL, xlab = "Cell count", 
                             log.scale = TRUE, bargraph = FALSE, fun = function(x) mean(x, 
                                                                                        na.rm = TRUE), ci.fun = function(x) c(fun(x) - 1.96 * 
                                                                                                                                se(x), fun(x) + 1.96 * se(x)), xlim = NULL, device = TRUE,  matching.string = NULL, ... ) 
{
  counts <- wholebrain::suggestions(dataset, exclude.below = exclude.below, 
                        reduce.below = reduce.below)
  if (!is.null(exclude.regions)) {
    remove <- which(row.names(counts) %in% exclude.regions)
    if (length(remove) > 0) {
      counts <- counts[-remove, ]
    }
  }
  if (!is.null(include.regions)) {
    dataset$acronym <- as.character(dataset$acronym)
    tableCount <- table(dataset$acronym, dataset$animal)
    to.be.included <- lapply(include.regions, function(x) {
      include<-tableCount[which(row.names(tableCount) %in% get.sub.structure(x)), ]
      ;
      if(length(dim(include))>1){
        include<-colSums(include)
      };
      return(include)
    })
    to.be.included <- do.call(rbind, to.be.included)
    row.names(to.be.included) <- include.regions
    counts <- rbind(counts, to.be.included)
  }
  
  if(is.null(matching.string)){
    matching.string = c("CTX", "CNU", "IB", "MB", "HB", "TH", "grey", "root", "VS", "fiber tracts")
  }
  
  dat<-data.frame(parent = unlist(lapply(row.names(counts), get.sup.structure, matching.string )), count = rowSums(counts))
  counts<-counts[order(dat$count, decreasing = TRUE),]
  dat<-dat[order(dat$count, decreasing = TRUE),]
  counts<-counts[order(dat$parent), ]
  
  if (!is.null(normalize.by)) {
    if (length(normalize.by) != ncol(counts)) {
      if(length(dim(normalize.by))>1){
        counts<-counts/normalization.table
      }else{
        normalize.by<-matrix(rep(normalize.by,  ncol(counts)), ncol= ncol(counts))
        for (j in 1:ncol(counts)) {
          counts[, j] <- counts[, j]/normalize.by[,j]
        }
      }
      #print(paste("Error: You have ", ncol(counts), " animals, but entered values ", 
      #           length(normalize.by), " for normalization."))
      #return()
    }else{
      for (j in 1:ncol(counts)) {
        counts[, j] <- counts[, j]/normalize.by[j]
      }
    }
  }
  if (log.scale) {
    counts <- log10(counts)
  }
  
  
  if (device) {
    quartz(width = 7.036585, height = 0.2099039 * nrow(counts))
  }
  layout(matrix(c(1, 1, 1, 2, 2, 2, 2), nrow = 1))
  par(mar = c(4, 0, 4, 0), ...)
  plot(rep(2.5, nrow(counts)), nrow(counts):1, col = 0, axes = F, 
       ylim = c(0.5, nrow(counts) + 0.5), ylab = "", xlab = "", 
       xlim = c(1, 5))
  mtext("Input region:", 3, cex = 0.9)
  for (i in 1:nrow(counts)) {
    regioncolor <- color.from.acronym(row.names(counts)[i])
    regioncolor <- adjustcolor(regioncolor, alpha.f = 0.2)
    y.lab <- (nrow(counts) + 1) - i
    polygon(c(1, 5, 5, 1), c(y.lab - 0.5, y.lab - 0.5, y.lab + 
                               0.5, y.lab + 0.5), col = regioncolor, border = FALSE)
    text(3, y.lab, name.from.acronym(row.names(counts)[i]), 
         cex = 0.9)
  }
  if (!is.null(group)) {
    par(xpd = TRUE)
    legend(3, -0.75, sort(unique(group)), pch = c(21), pt.bg = color, 
           title = title, bg = "white", horiz = TRUE, cex = 1.3, 
           xjust = 0.5)
    par(xpd = FALSE)
  }
  par(mar = c(4, 4, 4, 6), ...)
  zeros <- floor(range(counts[is.finite(counts)])[1])
  no.zero.values <- FALSE
  print(which(!is.finite(counts), arr.ind = TRUE))
  print(zeros)
  if (length(which(!is.finite(counts))) == 0) {
    no.zero.values <- TRUE
  }
  counts[!is.finite(counts)] <- zeros
  if (is.null(xlim)) {
    x.range <- c(floor(range(counts[is.finite(counts)])[1]), 
                 ceiling(range(counts[is.finite(counts)])[2]))
  }
  else {
    x.range <- xlim
  }
  
  plot(apply(counts, 1, max), nrow(counts):1 - 0.125, pch = 21, 
       bg = "white", ylim = c(0.5, nrow(counts) + 0.5), xlim = x.range, 
       xlab = "", axes = F, ylab = "", col = 0)
  for (i in 1:nrow(counts)) {
    regioncolor <- color.from.acronym(row.names(counts)[i])
    regioncolor <- adjustcolor(regioncolor, alpha.f = 0.15)
    y.lab <- (nrow(counts) + 1) - i
    polygon(c(x.range[1] - 1, x.range[2] + 1, x.range[2] + 
                1, x.range[1] - 1), c(y.lab - 0.5, y.lab - 0.5, y.lab + 
                                        0.5, y.lab + 0.5), col = regioncolor, border = FALSE)
  }
  if (log.scale) {
    log.range <- 10^seq(x.range[1], x.range[2])
    if (no.zero.values) {
      axis(1, at = seq(x.range[1], x.range[2]), las = 1, 
           labels = c(log.range))
      axis(3, at = seq(x.range[1], x.range[2]), las = 1, 
           labels = c(log.range))
    }
    else {
      axis(1, at = seq(x.range[1], x.range[2]), las = 1, 
           labels = c(0, log.range[-1]))
      axis(3, at = seq(x.range[1], x.range[2]), las = 1, 
           labels = c(0, log.range[-1]))
    }
    log.range <- unlist(lapply(1:(length(log.range) - 1), 
                               function(x) {
                                 seq(log.range[x], log.range[x + 1], by = log.range[x])
                               }))
    axis(1, at = log10(log.range), labels = FALSE)
    axis(3, at = log10(log.range), labels = FALSE)
    axis(2, at = nrow(counts):1, labels = row.names(counts), 
         las = 1)
    axis(4, at = nrow(counts):1, labels = row.names(counts), 
         las = 1)
    abline(v = log10(log.range), col = "lightblue")
    abline(h = 1:nrow(counts), lty = 2, col = "gray")
  }
  else {
    abline(h = 1:nrow(counts), lty = 2, col = "gray")
    axis(2, at = nrow(counts):1, labels = row.names(counts), 
         las = 1)
    axis(4, at = nrow(counts):1, labels = row.names(counts), 
         las = 1)
    abline(v = seq(x.range[1], x.range[2], length.out = 6), 
           col = "lightblue")
  }
  if (is.null(group)) {
    if (bargraph) {
      lapply(1:nrow(counts), function(x) {
        polygon(c(x.range[1], rep(fun(counts[x, ]), 2), 
                  x.range[1]), c(rep(nrow(counts) - x + 1, 4) + 
                                   c(-0.25, -0.25, 0.25, 0.25)), col = color)
        lines(ci.fun(counts[x, ]), rep(nrow(counts) - 
                                         x + 1, 2), lwd = 1.5)
        points(fun(counts[x, ]), nrow(counts) - x + 1, 
               pch = 21, bg = color, cex = 1.2)
      })
    }
    else {
      lapply(1:nrow(counts), function(x) {
        points(counts[x, ], rep(nrow(counts) - x + 1, 
                                ncol(counts)), pch = 21, bg = color, cex = 1.2)
      })
    }
  }
  else {
    k <- 1
    vertical <- seq(0.15, -0.15, length.out = length(unique(group)))
    for (j in sort(unique(group))) {
      if (bargraph) {
        lapply(1:nrow(counts), function(x) {
          polygon(c(x.range[1], rep(fun(counts[x, which(group == 
                                                          j)]), 2), x.range[1]), c(rep(nrow(counts) - 
                                                                                         x + 1 + vertical[k], 4) + c(-0.25, -0.25, 
                                                                                                                     0.25, 0.25)/length(unique(group))), col = color[k])
          lines(ci.fun(counts[x, which(group == j)]), 
                rep(nrow(counts) - x + 1 + vertical[k], 2), 
                lwd = 1.5)
          points(fun(counts[x, which(group == j)]), nrow(counts) - 
                   x + 1 + vertical[k], pch = 21, bg = color[k], 
                 cex = 1.2)
        })
      }
      else {
        lapply(1:nrow(counts), function(x) {
          points(counts[x, which(group == j)], rep(nrow(counts) - 
                                                     x + 1 + vertical[k], length(which(group == 
                                                                                         j))), pch = 21, bg = color[k], cex = 1.2)
        })
      }
      k <- k + 1
    }
  }
  box()
  if (log.scale) {
    if (!no.zero.values) {
      par(xpd = TRUE)
      polygon(c(mean(log10(log.range)[1:2]), log10(log.range)[2], 
                log10(log.range)[2], mean(log10(log.range)[1:2])), 
              c(-15, -15, nrow(counts) + 15, nrow(counts) + 
                  15), col = "white", border = "white")
      par(xpd = FALSE)
      abline(v = c(mean(log10(log.range)[1:2]), log10(log.range)[2]))
    }
  }
  mtext(xlab, 3, 2.2, cex = 0.8)
  mtext(xlab, 1, 2.2, cex = 0.8)
  if (log.scale) {
    counts <- 10^(counts)
  }
  return(counts)
}