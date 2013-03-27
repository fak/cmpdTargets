#############################################
###Function: summaryPlots.R
###  --------------
###momo.sander@googlemail.com
##############################################

infile <- "data/query_results.tab"


## @knitr load
# Loading.
library(ggplot2)
intable <- read.table(infile ,sep="\t", header=TRUE, na.strings = 'None')
mainFrame <- as.data.frame(intable)

## @knitr prepData
mainFrame$lab <- abbreviate(mainFrame$target, minlength = 10)
mainFrame <- mainFrame[with(mainFrame, order(median)), ]
mainFrame$lab <- reorder(mainFrame$lab, mainFrame$median)

## @knitr plots
ids <- unique(mainFrame$chembl_id)
lapply(ids, function(z)
ggplot(mainFrame[mainFrame$chembl_id ==  z,])+geom_boxplot(aes(y=conc, x= lab))+geom_jitter(aes(y=conc, x=
 lab))+coord_flip(ylim = c(4.5,11))+facet_wrap(~name)
)


## @knitr oldschool
tt <- lapply(ids, function(z) 
ggplot(mainFrame[mainFrame$chembl_id ==  z,])+geom_boxplot(aes(y=conc, x= lab))+geom_jitter(aes(y=conc, x= lab))+coord_flip(ylim = c(4.5,11))+facet_wrap(~name)
)#+theme(axis.text.x=element_text(angle=-60, hjust = 0))

lapply(seq(length(tt)), function(z)
ggsave(tt[[z]], file = paste('visual/',z, '.pdf', sep = ''), useDingbats = F, width = 5, height = 3)
)

## @knitr legend
library(plyr)
lkp <- mainFrame[!duplicated(mainFrame[,c('lab','target')]),c(9,2)]
lkp <- arrange(lkp, as.character(lab))


