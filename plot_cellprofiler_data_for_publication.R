#import packages
library(tidyverse)
library(ggpol)
library(plyr)
#functions for ggplot plot layout adapted from: https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#ffffff"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.margin = unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#ffffff",fill="#ffffff"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#c00000","#1F497D","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#c00000","#6599D9","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}

#Load all data with the features per nucleus
read_csv("C:/Users/maart/Erasmus MC/source data/")

#convert the different conditions into factors for optimal plotting results
nuclei$Metadata_condition <- factor(nuclei$Metadata_condition,levels = c('noIR','2hIR','8hIR','24hIR'))
nuclei$Metadata_protein <- factor(nuclei$Metadata_protein,levels = c('WT','dDBD','dCTD','dDBDdCTD'))
nuclei$condition <- mapvalues(nuclei$Metadata_condition, from = c('noIR','2hIR','8hIR','24hIR'),to = c('no IR','2h','8h','24h'))

#Number of nuclei per protein
filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  count(vars = c("Metadata_protein"))

#Number of EdU positive cells
nuclei$Edu_pos <- nuclei$Intensity_IntegratedIntensity_EdU>500
table(nuclei$Edu_pos,nuclei$Metadata_protein,nuclei$Metadata_condition)

#Box-plot of foci numbers per nucleus plotted per timepoint
plot1 <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Children_RAD51_Count,x=condition,fill=Metadata_protein))+
  geom_boxplot(notch=TRUE)+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+ylab("Number of RAD51 foci/nucleus")+
  theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))+stat_summary(fun.y=mean, colour="white", geom="point", shape=1, size=2,show_guide = FALSE)
plot1
ggsave(plot1,filename = "RAD51 foci frequencies_boxplot.pdf",width = 10, units="in", height=6,useDingbats=F)


#Box-plot of foci numbers per nucleus plotted per genotype
plot1 <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Children_RAD51_Count,x=Metadata_protein,fill=Metadata_protein))+
  geom_boxplot(notch=TRUE)+
  facet_grid(.~condition)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+ylab("Number of RAD51 foci/nucleus")+
  theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size=14))+stat_summary(fun.y=mean, colour="white", geom="point",shape=1, size=2,show_guide = FALSE)
plot1

#bar plot with mean and SD per condition
plot2 <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  group_by(Metadata_protein,condition)%>%
  dplyr::summarise(mean_foci=mean(Children_RAD51_Count),sd_foci=sd(Children_RAD51_Count)) %>%
  mutate(prot_condition=paste(Metadata_protein,condition)) %>%
  ggplot(aes(y=mean_foci,x=prot_condition,fill=Metadata_protein))+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean_foci-sd_foci, ymax=mean_foci+sd_foci), width=.2,
                position=position_dodge(.9)) + scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("")+ylab("Number of RAD51 foci/nucleus")+theme(legend.position = "none")
plot2

#Plotting fraction of RAD51 positive cells with arbitrary cut-off
nuclei %>%
  filter(Intensity_IntegratedIntensity_EdU>500,) %>%
  dplyr::mutate(R51_positive=Children_RAD51_Count>25) %>%
  group_by(Metadata_condition,Metadata_protein,R51_positive) %>%
  dplyr::summarise(n=n())%>%
  dplyr::mutate(freq = n / sum(n)) %>%
  filter(R51_positive==TRUE) %>%
  ggplot(aes(y=freq,x=Metadata_condition,fill=Metadata_protein))+geom_bar(stat="identity")+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("")+ggtitle("cells > 25 foci")

nuclei %>%
  filter(Intensity_IntegratedIntensity_EdU>500) %>%
  dplyr::mutate(R51_positive=Children_RAD51_Count>30) %>%
  group_by(Metadata_condition,Metadata_protein,experiment,R51_positive) %>%
  dplyr::summarise(n=n())%>%
  dplyr::mutate(freq = n / sum(n)) %>%
  filter(R51_positive==TRUE) %>%
  group_by(Metadata_condition,Metadata_protein)%>%
  dplyr::summarise(mean=mean(freq),sd=sd(freq)/sqrt(n()))%>%
  ggplot(aes(y=mean,x=Metadata_condition,fill=Metadata_protein))+geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("")+ggtitle("cells > 25 foci")

#statistics comparing timepoints
statistical_data <- matrix(nrow = 4,ncol=4)
proteins <- c("no IR","2h","8h","24h")
for (i in 1:4){
  for (j in 1:4){
    group_a <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500,Metadata_protein=="WT",condition==proteins[i])
    group_b <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500,Metadata_protein=="WT",condition==proteins[j])
    statistical_data[i,j] <- t.test(group_a$Children_RAD51_Count,group_b$Children_RAD51_Count)$estimate[1]
  }
}
   
View(statistical_data)

#statistics compare variants
statistical_data <- matrix(nrow = 4,ncol=4)
proteins <- c("WT","dDBD","dCTD","dDBDdCTD")
for (i in 1:4){
  for (j in 1:4){
    group_a <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500,Metadata_protein==proteins[i],condition=="24h")
    group_b <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500,Metadata_protein==proteins[j],condition=="24h")
    statistical_data[i,j] <- t.test(group_a$Mean_RAD51_Intensity_IntegratedIntensity_MaskedGreen,group_b$Mean_RAD51_Intensity_IntegratedIntensity_MaskedGreen)$estimate[1]
  }
}

View(statistical_data)



#join the nuclei and RAD51 foci tables to plot the integrated intensity per focus
nuclei_R51 <- mutate(nuclei,Parent_Nuclei=ObjectNumber) %>%
    left_join(y=RAD51,by=c("Parent_Nuclei","ImageNumber","Metadata_FileLocation")) %>%
    filter(Intensity_IntegratedIntensity_EdU>500)
  
  
  plot3 <- ggplot(nuclei_R51,aes(y=Intensity_IntegratedIntensity_MaskedGreen,x=Metadata_condition.x,fill=Metadata_protein.x))+
    geom_boxplot(notch = T,outlier.shape = NA)+ facet_grid(.~Metadata_protein.x)+
    scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+ylab("Integrated intensity/focus")+ylim(0,20)+
    theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))
  plot3
  ggsave(plot3,filename = "RAD51 foci individual intensity_boxplot.pdf",width = 10, units="in", height=6,useDingbats=F)
  

#Plot mean RAD51 integrated density
plot4 <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Mean_RAD51_Intensity_IntegratedIntensity_MaskedGreen,x=Metadata_condition,fill=Metadata_protein))+geom_boxplot(notch=TRUE)+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+
  ylab("Mean RAD51 focus intensity/nucleus")+xlab("")+
  theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))
plot4
ggsave(plot4,filename = "RAD51 foci intensity_boxplot.pdf",width = 10, units="in", height=6,useDingbats=F)

#plot number of foci averaged per experimental replicate
mean_foci <- nuclei%>%
  filter(Intensity_IntegratedIntensity_EdU>500) %>%
  group_by(experiment,Metadata_protein,condition) %>%
  dplyr::summarize(RAD51_foci=mean(Children_RAD51_Count)) %>%
  group_by(Metadata_protein,condition) %>%
  dplyr::summarize(mean_foci=mean(RAD51_foci),sd_foci=sd(RAD51_foci),n=n())
#View(mean_foci)
mean_foci %>%
  ggplot(aes(y=mean_foci,x=condition)) + geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_foci-sd_foci, ymax=mean_foci+sd_foci), width=.2,
                position=position_dodge(.9)) +facet_grid(.~Metadata_protein)+ggtitle("Mean foci/experiment (n=3)")+
  scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+
  theme(legend.position = "none")

#fold plots
x <- nuclei%>%
  filter(Intensity_IntegratedIntensity_EdU>500) %>%
  group_by(experiment,Metadata_protein,condition) %>%
  dplyr::summarize(RAD51_foci=mean(Children_RAD51_Count)) %>%
  group_by(Metadata_protein,experiment) %>%
  dplyr::mutate(fold=RAD51_foci/dplyr::first(RAD51_foci)) %>%
  group_by(Metadata_protein,condition) %>%
  dplyr::summarize(mean_fold=mean(fold),sd_fold=sd(fold)/sqrt(n()),n=n())

library(clipr)
x %>%
ggplot(aes(x=condition,y=mean_fold,group=Metadata_protein,color=Metadata_protein))+
  geom_point()+geom_line()+geom_errorbar(aes(ymin=mean_fold-sd_fold, ymax=mean_fold+sd_fold))+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+
  theme(legend.position = "none")


x <- nuclei %>%
  group_by(experiment,Metadata_protein,condition) %>%
  dplyr::mutate(EdUpos = Intensity_IntegratedIntensity_EdU>500) %>%
  dplyr::summarize(EdU = (sum(EdUpos==TRUE)/n()))%>%
  group_by(Metadata_protein,condition)

#fraction of cells in S-phase
x %>%  
  dplyr::summarize(mean_EdU=mean(EdU),sd_EdU=sd(EdU),n=n()) %>%
  ggplot(aes(y=mean_EdU,x=condition,fill=Metadata_protein)) + geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_EdU-sd_EdU, ymax=mean_EdU+sd_EdU), width=.2,
                position=position_dodge(.9)) +facet_wrap(.~Metadata_protein,nrow = 1)+
  geom_point(data=x, aes(y=EdU, x=condition))+
  scale_colour_Publication()+ylab("% cells in S-phase")+ylim(0,1)+
  scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+
  theme(legend.position = "none")



