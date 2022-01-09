library(ggplot2)

df <- read.csv('../classification/classification_s=500.csv')
df[,1] <- NULL

pdf("allele_freq_corr_s=500.pdf")

ggplot(data=df, aes(x=cons_score , y=Allele_frequency, color=isPathogenic)) +
  xlab("Conservation Score")+ ylab("Allele Frequency")+ geom_point() + theme(legend.title = element_text(colour="Black", size=15)) 

dev.off()
