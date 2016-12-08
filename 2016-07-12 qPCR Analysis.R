library(ggplot2)
library(gridExtra)
library(reshape2)

save.image(paste0(dir,"/", date," ",expNum, " W8 qPCR.RData"))

## Set environment variables---------
dir <- setwd("~/Google Drive/Weill Cornell/Research/Analysis/qPCR")
date <- paste0(Sys.Date())
expNum <- paste0("3-B-1")

## Import raw qPCR data from Google Drive---------------------
qpcrdata <- read.delim("~/Google Drive/Weill Cornell/Research/Analysis/qPCR/3-1B-1 W8 qPCR.txt", stringsAsFactors=TRUE)

## Set variables as factors or numeric---------------------
qpcrdata$Rep <- as.factor(qpcrdata$Rep)
qpcrdata$Group <- as.factor(qpcrdata$Group)
qpcrdata$Ct <- as.numeric(qpcrdata$Ct)

castqpcrdata <- dcast(qpcrdata, ... ~ Target, value.var = "Ct")
normgene <- castqpcrdata[,"ACTB"]

coldata <- castqpcrdata[,1:6]
normct <- castqpcrdata[,-(1:6)]
normct <- normct - normgene
castqpcrdata <- cbind(coldata, normct)
castqpcrdata <- castqpcrdata[,-grep("ACTB", colnames(castqpcrdata))]

meltedqpcrdata <- melt(castqpcrdata, 
                       id.vars = colnames(coldata),
                       variable.name = "Target", 
                       value.name = "dCt")

meltedqpcrdata$PopulationDiet <- paste0(meltedqpcrdata$Population,"_", meltedqpcrdata$Diet)
meltedqpcrdata$PopulationDietGroup <- paste0(meltedqpcrdata$Population,"_", meltedqpcrdata$Diet,"_",meltedqpcrdata$Group)
meltedqpcrdata <- meltedqpcrdata[!(is.na(meltedqpcrdata$dCt)),]
meltedqpcrdata$TissueTarget <- paste0(meltedqpcrdata$Tissue," ", meltedqpcrdata$Target)

populationcolors <- c(Treg="#4f1310", 
                  ILC2="#212121", 
                  CD4T="#131e49")[meltedqpcrdata$Population]

populationdietcolors <- c(Treg_Control="#E39794", 
                      ILC2_Control="#d4d4d4", 
                      CD4T_Control="#b3c1fb",
                      Treg_Fat="#C7302A", 
                      ILC2_Fat="#707070", 
                      CD4T_Fat="#4266F6")[meltedqpcrdata$PopulationDiet]


cairo_pdf(paste0(expNum," W8 qPCR Analysis.pdf"), w=13, h=10)
ggplot(data=meltedqpcrdata[-grep("TGFB3", meltedqpcrdata$Target),], aes(x=PopulationDietGroup, y=(2^(-dCt)), fill=PopulationDiet, color=Population)) +
  facet_wrap(  ~ TissueTarget, scales="free_y", ncol=4) +
  stat_summary (fun.y = "mean", geom = "bar", position = position_dodge(0.9), width=0.8) +
  scale_fill_manual(values = populationdietcolors) +
  scale_color_manual(values = populationcolors ) +
  geom_jitter(position=position_jitterdodge(0.5), shape=21, fill="white", size=1.5) +
  labs(title=paste0(date," ",expNum," W8 qPCR Analysis")) +
  theme(legend.position='none',
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = rel(1.2)),
        panel.background = element_rect(color="#707070", linetype="solid", fill=NA))
dev.off()


prismdata <- dcast(meltedqpcrdata,  Tissue+PopulationDiet ~ Target + Group, mean , value.var = "dCt")
write.csv(prismdata, file = "3-1B-1 W8 Cast qPCR Analysis.csv")
