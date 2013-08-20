#############################################
###Function: summaryPlots.R
###  --------------
###momo.sander@googlemail.com
##############################################

infile <- "data/query_results.tab"


## @knitr load
# Loading.
library(ggplot2)
library(rcdk)
intable <- read.table(infile ,sep="\t", header=TRUE, na.strings = 'None')
mainFrame <- as.data.frame(intable)

## @knitr prepData
targetFrame <- mainFrame[mainFrame$assigned_target == assigned_target,]
targetFrame$lab <- sapply(targetFrame$pref_name, function(z) paste(strwrap(z, width = 20), collapse = "\n"))
targetFrame <- targetFrame[with(targetFrame, order(median)), ]
targetFrame$lab <- reorder(targetFrame$lab, targetFrame$median)
targetFrame$bool <- sapply(targetFrame$median, function(z) if(z == 4){'inactive'} else {'active'})

## @knitr plot
gen_frames <- function(data,z){
   plotFrame <- data[data$chembl_id ==  z,]
   plotFrame <- plotFrame[with(plotFrame, order(median)), ]
   plotFrame$lab <- reorder(plotFrame$lab, plotFrame$median)
   idx <- tail(unique(plotFrame$lab), n = 20)
   smi <- plotFrame$smiles[1]
   print(smi)
   m <- parse.smiles(as.character(smi))[[1]]
   img <- view.image.2d(m, 600,800)
   plot.new()
   title(unique(plotFrame$aggName))
   rasterImage(img, 0,0,1,1)
   print(unique(plotFrame$aggName))
   print(p %+% plotFrame[plotFrame$lab %in% idx,])
   print(bar %+% plotFrame[!duplicated(plotFrame[, 'pref_name']),])
}
plotFrame <- NULL

p <- ggplot(plotFrame)+
  geom_dotplot(aes(y=conc, x= lab,fill = bool ), binaxis = "y", stackdir = "center", binwidth = .1, dotsize = 2)+
  facet_wrap(~chembl_id)+
  xlab("")+
  ylab('pIC50')+
  coord_flip(ylim = c(3.5,11))+
  theme(legend.position = "none", axis.text.y = element_text(size = 11))

bar <- ggplot(plotFrame)+geom_bar(aes(x=bool, y = ..count.., fill = bool ))+theme(legend.position = "none", axis.text.y = element_text(size = 11))


## knitr plot_it
ids <- sort(unique(targetFrame$chembl_id))
plots <- sapply(ids, function(z) gen_frames(targetFrame, z))

## knitr data_sum
#plotFrame <- mainFrame[!duplicated(mainFrame[, c("chembl_id","assigned_target")]), ]
#ggplot(plotFrame)+geom_bar(aes(x = assigned_target, y = ..count..))+theme(axis.text.x = element_text(vjust = 0,angle = 45))




## @knitr lkp
library(plyr)
lkp <- targetFrame[!duplicated(targetFrame[,c('lab','target')]),c(11,2)]
lkp <- arrange(lkp, as.character(target))


