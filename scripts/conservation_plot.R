library(ggplot2)

df <- read.csv('../conservation/conservation_scores_s=500.tsv', sep='\t', h=T)

pdf("conservation_scores_s=500.pdf", width=25, height=8)

names(df) <- c("Position", "ConsensusAminoAcid" ,"ConservationScore")

ggplot(data=df, aes(x=Position , y=ConservationScore, color= ConservationScore)) +
  xlab("Position")+ ylab("Conservation Score")+ geom_point()+ geom_line(aes(color=ConservationScore)) +
  scale_color_gradient(low="Blue", high = "Red")+ theme(legend.title = element_text(colour="Black", size=15)) 

dev.off()
