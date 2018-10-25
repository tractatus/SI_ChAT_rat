library(wholebrain)
load('./data/all_rats.RData')

source('plot_suggestions2.R')

#remove section with dirty segmentation
datasets<-datasets[-which(datasets$animal == 'R2025' & datasets$AP == '-0.39'),]

#create group variable
datasets$group<-factor(datasets$animal)
levels(datasets$group)<-c('Medial', 'Lateral', 'Medial', 'Lateral')

#colors
colors<-list(
  lateral = rgb(175, 44, 120, maxColorValue =  255),
  medial = rgb(45, 135, 135, maxColorValue =  255)
)

#make initial suggestion for plot to then normalize on
counts<-plot.suggestions2(datasets, exclude.below = 10, reduce.below = 10, include.regions = c("ACA", "SNr", "LH"), matching.string = c("CTX", "CNU", "IB", "MB", "HB", "TH", "grey", "root", "VS", "fiber tracts") )

#make nornalization by brain region
normalization.table<-table(datasets$animal)
region.volume<-get.region.volume(row.names(counts))
normalization.table<-matrix(rep(normalization.table, length(region.volume)), ncol = length(normalization.table) )
normalization.table<-normalization.table*region.volume

#make plot
quartz(width = 6.497297, height = 6.508108)
counts<-plot.suggestions2(datasets, exclude.below = 10, reduce.below = 10, include.regions = c("ACA", "SNr", "LH"), xaxs = 'i', yaxs='i', mar= c(7.2,4,4,3.5), bargraph=TRUE, col = c(colors$lateral, colors$medial), device = FALSE, group = c('Medial', 'Lateral', 'Medial', 'Lateral'), normalize.by = normalization.table, log.scale = F, xlim=c(0,0.05), xlab = "Normalized cell count per mm^3" )
axis(1)
axis(3)