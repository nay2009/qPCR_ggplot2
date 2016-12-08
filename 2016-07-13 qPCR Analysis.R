library(ggplot2)
library(gridExtra)
library(reshape2)

save.image(paste0(dir,"/",expNum, " W9 qPCR.RData"))

## Set environment variables---------
dir <- setwd("~/Google Drive/Weill Cornell/Research/Analysis/qPCR")
date <- paste0(Sys.Date())
expNum <- paste0("3-1A-3")

## Import raw qPCR data from Google Drive---------------------
qpcrdata <- read.delim("~/Google Drive/Weill Cornell/Research/Analysis/qPCR/3-1A-3 W9 qPCR.txt", stringsAsFactors=TRUE)

## Set variables as factors or numeric---------------------
qpcrdata$Replicate <- as.factor(qpcrdata$Replicate)
qpcrdata$Group <- as.factor(qpcrdata$Group)
qpcrdata$dCt <- qpcrdata$Ct-qpcrdata$Norm


## castqpcrdata <- dcast(qpcrdata[,-(8:9)], ... ~ Target, value.var = "dCt")
## normgene <- castqpcrdata[,"ACTB"]

## coldata <- castqpcrdata[,1:6]
## normct <- castqpcrdata[,-(1:6)]
## normct <- normct - normgene
## castqpcrdata <- cbind(coldata, normct)
## castqpcrdata <- castqpcrdata[,-grep("ACTB", colnames(castqpcrdata))]

## meltedqpcrdata <- melt(castqpcrdata, 
                       ## id.vars = colnames(coldata),
                       ## variable.name = "Target", 
                       ## value.name = "dCt")

meltedqpcrdata <- qpcrdata[,c(1:8,11)]


meltedqpcrdata <- meltedqpcrdata[!(is.na(meltedqpcrdata$dCt)),]
meltedqpcrdata <- meltedqpcrdata[-grep("IL9|ACTB", meltedqpcrdata$Target),]

meltedqpcrdata$PopulationDiet <- paste0(meltedqpcrdata$Population,"_", meltedqpcrdata$Diet)
meltedqpcrdata$PopulationDietGroup <- paste0(meltedqpcrdata$Population,"_", meltedqpcrdata$Diet,"_",meltedqpcrdata$Group)
meltedqpcrdata$TissuePopulationDiet <- paste0(meltedqpcrdata$Tissue,"_", meltedqpcrdata$Population,"_",meltedqpcrdata$Diet)
meltedqpcrdata$TissuePopulationDietGroup <- paste0(meltedqpcrdata$Tissue,"_", meltedqpcrdata$PopulationDietGroup)
meltedqpcrdata$TargetTissue <- paste0(meltedqpcrdata$Target," ", meltedqpcrdata$Tissue)

populationcolors <- c(Treg="#4f1310", 
                  ILC2="#212121", 
                  CD4T="#131e49")[meltedqpcrdata$Population]

populationdietcolors <- c(Treg_Control="#E39794", 
                      ILC2_Control="#d4d4d4", 
                      CD4T_Control="#b3c1fb",
                      Treg_Fat="#C7302A", 
                      ILC2_Fat="#707070", 
                      CD4T_Fat="#4266F6")[meltedqpcrdata$PopulationDiet]

mytheme <- theme(legend.position='top',
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = rel(1), angle = 90),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = rel(1.2)),
      panel.background = element_rect(color="#707070", linetype="solid", fill=NA))

cairo_pdf(paste0(expNum," W9 qPCR Analysis.pdf"), w=15, h=20)
ggplot(data=meltedqpcrdata, aes(x=PopulationDietGroup, y=(2^(-dCt)), fill=PopulationDiet, color=Population)) +
  facet_wrap(  ~ TargetTissue, scales="free_y", ncol=6) +
  stat_summary (fun.y = "mean", geom = "bar", position = position_dodge(0.9), width=0.8) +
  scale_fill_manual(values = populationdietcolors) +
  scale_color_manual(values = populationcolors ) +
  geom_jitter(position=position_jitterdodge(0.5), shape=21, fill="white", size=1.5) +
  labs(title=paste0(date," ",expNum," W9 qPCR Analysis")) +
  mytheme
dev.off()

prismdata <- dcast(meltedqpcrdata,  Tissue+PopulationDiet ~ Target + Group, mean , value.var = "dCt")
write.csv(prismdata, file = "3-1A-3 W9 Cast qPCR Analysis.csv")
