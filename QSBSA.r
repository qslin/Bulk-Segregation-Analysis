#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("ggplot2", lib="~/local/R_libs/") 
library("reshape2", lib="~/local/R_libs/")
library("labeling", lib="~/local/R_libs/")

d <- getwd()
A <- read.table(file = paste(d, "/", args[1], ".txt", sep = ""),header = TRUE,sep = "\t")

setEPS(paper="special",width = 10,height = 2)
postscript(file = paste(args[1], ".eps", sep = ""))
tu <- ggplot(data = A, aes(x=position, y=count)) + 
      geom_bar(stat="identity", position=position_dodge(), width = 1, color = "black", size = 0.001) +
      facet_grid(~contig, scales="free", space="free", switch = "x") +
      theme(legend.position = "none") +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(strip.background = element_rect(colour="white", fill="white")) +
      scale_y_continuous(expand = c(0,0), limits = c(0,as.numeric(args[2])), breaks = seq(0,as.numeric(args[2]),4)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) 
tu

