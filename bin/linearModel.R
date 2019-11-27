#!/usr/bin/env R

library(ggplot2)

args <- commandArgs(trailingOnly=T)
data <- read.table(args[1], header=F, na.strings = "NA", comment.char='')
colnames(data) <- c("reads","duplicates")

# calculate the linear regression model
linearModel <- lm(reads ~ duplicates, data=data)

# print statistics to file
modelSummary <- summary(linearModel)  
sink("Duplicates.txt")
print(modelSummary)
sink() 

# print scatter plot to file 
pdf("Duplicates.pdf")

p=ggplot(data, aes(reads, duplicates, color = reads)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  labs(x ="Total # of alignments", y = "Total # of duplicates") +
  theme( axis.line = element_line(colour = "darkblue", 
                                  size = 1, linetype = "solid")) +
  stat_smooth(method="lm", se=FALSE, color= "darkred") +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") 

p + theme(
  plot.title = element_text(color="darkred", size=14, face="bold"),
  axis.title.x = element_text(color="black", size=14),
  axis.title.y = element_text(color="black", size=14))

dev.off()