#############################################
###Function: compileKnitr.R
###  --------------
###momo.sander@googlemail.com
##############################################
library(knitr)

key <- gsub("-","", commandArgs()[8])
pubdir <- gsub("-","", commandArgs()[9])
path <- sprintf(key)
name <- strsplit(path, '\\.')[[1]][1]

knit(path)
system(sprintf("pandoc %s.md -s -f markdown+pipe_tables+table_captions -c buttondown.css -o %s.html", name, name))
system(sprintf("sed 's,<img src=\"\\([^\"]*\"\\) alt=\"plot of chunk plot\" />,<object data=\"\\1 type=\"image/svg+xml\"><\\/object>,g' %s.html > %s.html.sed", name, name))
system(sprintf("mv %s.html.sed %s.html",name, name))
system(sprintf("rsync -r  figure/%s_* %s/%s/wiki/figure", name, pubdir, name))
system(sprintf("rsync %s.html %s/%s/wiki/", name, pubdir, name))

