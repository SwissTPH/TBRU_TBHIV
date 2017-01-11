####R version: 3.2.2 (2015-08-14) -- "Fire Safety"####
####Copyright (C) 2015 The R Foundation for Statistical Computing###
####Platform: x86_64-pc-linux-gnu (64-bit)###
####Daniela Brites###
####date:15.12.2016
###Analysis of the relationship between number of singletons per strain and CD4 Cell counts###### 


#data  on number of singletons which are found in strain by cell counts in categories more or less than 250 Tcell/ul.
data <- read.table ("~/data/singletons_hiv_cell.counts_categorical.csv", header=T,sep="\t")

#test if there are more singletons in strains with lower than 250 TCell/ul 
wilcox.test (singletons~Cell_counts,data=data)


boxplot(singletons~Cell_counts,data=data,ylab="n° singletons per genome")

# number of singletons per genome and number of CDTCells per isolate.

data_sing_cells <- read.table ("~/data/singletons_cellCounts_HIV.csv",header=T,sep=",")
attach(data_sing_cells)
head(data_sing_cells)

plot(CD4.cells.ul,nøsingletons, ylab="n°singletons per genome", xlab="CD4 T-cells/ul")

#fit lowess 
fit_lowess <-lowess(CD4.cells.ul,nøsingletons)
lines(fit_lowess, col = 2)
