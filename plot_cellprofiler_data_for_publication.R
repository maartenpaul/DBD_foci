library(tidyverse)
library(ggpol)
library(plyr)
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



#nuclei <- read_csv("C:/Dropbox/Postdoc Genetics/R51 data/MyExpt_EdUNuclei_Nuclei.csv")
nuclei1 <- read_csv("D:/Dropbox/Postdoc Genetics/R51 data/DAPI_Foci_Nuclei.csv")
nuclei1$experiment <- "exp1"
nuclei2 <- read_csv("D:/OneDrive/Data2/MPexp1901_09 RAD51 foci/DAPI_Foci_Nuclei2.csv")
nuclei2$experiment <- "exp2"
nuclei3 <-  read_csv("D:/OneDrive/Data2/190527 ExpMP1905_05 RAD51 foci/DAPI_Foci_Nuclei.csv")
nuclei3$experiment <- "exp3"


nuclei <- rbind(nuclei1,nuclei2,nuclei3)
RAD51_1 <- read_csv("D:/Dropbox/Postdoc Genetics/R51 data/DAPI_Foci_RAD51.csv")
RAD51_2 <- read_csv("D:/OneDrive/Data2/MPexp1901_09 RAD51 foci/DAPI_Foci_RAD51.csv")
RAD51_3 <- read_csv("D:/OneDrive/Data2/190527 ExpMP1905_05 RAD51 foci/DAPI_Foci_RAD51.csv")


RAD51 <- rbind(RAD51_1,RAD51_2,RAD51_3)
#View(RAD51)

#nuclei <- mutate(nuclei,R51_positive=Children_RAD51_Count>20)

#nuclei$Metadata_protein[nuclei$Metadata_protein=='dCTDdDBD'] <- 'dDBDdCTD'

#plot(nuclei$Intensity_IntegratedIntensity_EdU,nuclei$Intensity_IntegratedIntensity_DAPI)

nuclei$Metadata_condition <- factor(nuclei$Metadata_condition,levels = c('noIR','2hIR','8hIR','24hIR'))
nuclei$Metadata_protein <- factor(nuclei$Metadata_protein,levels = c('WT','dDBD','dCTD','dDBDdCTD'))

nuclei$condition <- mapvalues(nuclei$Metadata_condition, from = c('noIR','2hIR','8hIR','24hIR'),to = c('no IR','2h','8h','24h'))
  

group_by(nuclei,Metadata_condition,Metadata_protein) %>%
  summarise(mean=mean(Children_RAD51_Count))

data <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  select(Protein=Metadata_protein,Condition=condition,RAD51_Count=Children_RAD51_Count,Mean_RAD51_Intensity=Mean_RAD51_Intensity_IntegratedIntensity_MaskedGreen)%>%
  mutate("Protein_Condition"=paste0(Protein,"_",Condition))%>%
  write_csv("data.csv")


#Foci frequencies
plot <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Children_RAD51_Count,x=condition,fill=Metadata_protein))+
  geom_boxjitter(errorbar.draw = TRUE,jitter.height = 0, jitter.width = 0.1,notch=TRUE,jitter.alpha=0.8,jitter.size=0.8,outlier.intersect=TRUE)+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+ylab("Number of RAD51 foci/nucleus")+
  theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))
plot


plotbox <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Children_RAD51_Count,x=condition,fill=Metadata_protein))+
  geom_boxplot(notch=TRUE)+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+ylab("Number of RAD51 foci/nucleus")+
  theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))+stat_summary(fun.y=mean, colour="white", geom="point", 
                                                                                                                                       shape=1, size=2,show_guide = FALSE)
plotbox

plotbox <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Children_RAD51_Count,x=Metadata_protein,fill=Metadata_protein))+
  geom_boxplot(notch=TRUE)+
  facet_grid(.~condition)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+ylab("Number of RAD51 foci/nucleus")+
  theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size=14))+stat_summary(fun.y=mean, colour="white", geom="point", 
                                                                                                                                       shape=1, size=2,show_guide = FALSE)
plotbox


filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  group_by(Metadata_protein,condition)%>%
  dplyr::summarise(mean_foci=mean(Children_RAD51_Count),sd_foci=sd(Children_RAD51_Count)) %>%
  mutate(prot_condition=paste(Metadata_protein,condition)) %>%
  ggplot(aes(y=mean_foci,x=prot_condition))+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean_foci-sd_foci, ymax=mean_foci+sd_foci), width=.2,
                position=position_dodge(.9)) 

#compare timepoints
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

#compare variants
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

group_a <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500,Metadata_protein=="WT",condition=="no IR")


group_b <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500,Metadata_protein=="WT",condition=="2h")

t.test(group_a$Children_RAD51_Count,group_b$Children_RAD51_Count)
  
  
  ggplot(aes(y=Children_RAD51_Count,x=condition,fill=Metadata_protein))+
  geom_boxplot(notch=TRUE)+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+ylab("Number of RAD51 foci/nucleus")+
  theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))+stat_summary(fun.y=mean, colour="white", geom="point", 
                                                                                                                                       shape=1, size=2,show_guide = FALSE)
plotbox

ggsave(plotbox,filename = "RAD51 foci frequencies_boxplot.pdf",width = 10, units="in", height=6,useDingbats=F)

#statistcs
filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  count(vars = c("Metadata_protein"))

#EdU positive
nuclei$Edu_pos <- nuclei$Intensity_IntegratedIntensity_EdU>500

table(nuclei$Edu_pos,nuclei$Metadata_protein,nuclei$Metadata_condition)


#Foci 
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

nuclei <- mutate(nuclei,Parent_Nuclei=ObjectNumber)

  nuclei_R51 <- left_join(x=nuclei,y=RAD51,by=c("Parent_Nuclei","ImageNumber","Metadata_FileLocation")) %>%
    filter(Intensity_IntegratedIntensity_EdU>500)
  
  
  plot_foci <- ggplot(nuclei_R51,aes(y=Intensity_IntegratedIntensity_MaskedGreen,x=Metadata_condition.x,fill=Metadata_protein.x))+
    geom_boxplot(notch = T,outlier.shape = NA)+ facet_grid(.~Metadata_protein.x)+
    scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+ylab("Integrated intensity/focus")+ylim(0,20)+
    theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))
  plot_foci
  ggsave(plot_foci,filename = "RAD51 foci individual intensity_boxplot.pdf",width = 10, units="in", height=6,useDingbats=F)
  

    left_join(x=RAD51,y=nuclei,by=c("Parent_Nuclei","ImageNumber")) %>%
    filter(x$Intensity_IntegratedIntensity_EdU>500)%>%
    ggplot(aes(y=Intensity_IntegratedIntensity_MaskedGreen,x=Metadata_condition.y,fill=Metadata_protein.y))+
    geom_boxplot()+ facet_grid(.~Metadata_protein.y)+ 
    scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+ylab("Integrated intensity/focus")+
    theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))
  
  
  #Distance from centroid
filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Mean_RAD51_Distance_Centroid_Nuclei,x=Metadata_condition))+geom_boxplot()+facet_grid(.~Metadata_protein)

#Mean RAD51 integrated density
plot2 <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Mean_RAD51_Intensity_MeanIntensity_MaskedGreen,x=Metadata_condition,fill=Metadata_protein))+geom_boxjitter(errorbar.draw = TRUE,jitter.height = 0, jitter.width = 0.1,notch=TRUE,jitter.alpha=0.8,jitter.size=0.8,outlier.intersect=TRUE)+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+
  ylab("Mean RAD51 focus intensity/nucleus")
  theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))
plot2
ggsave(plot2,filename = "RAD51 foci intensity.pdf",width = 10, units="in", height=6,useDingbats=F)

plot3 <- filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Mean_RAD51_Intensity_IntegratedIntensity_MaskedGreen,x=Metadata_condition,fill=Metadata_protein))+geom_boxplot(notch=TRUE)+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+
  ylab("Mean RAD51 focus intensity/nucleus")+xlab("")+
  theme(legend.position = "none")+ggtitle("BRCA2-HaloTag RAD51 foci IR (2Gy)")+ theme(plot.title = element_text(size=14))
plot3
ggsave(plot3,filename = "RAD51 foci intensity_boxplot.pdf",width = 10, units="in", height=6,useDingbats=F)

filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Mean_RAD51_AreaShape_Area,x=Metadata_condition))+geom_boxplot()+facet_grid(.~Metadata_protein)

filter(nuclei,Intensity_IntegratedIntensity_EdU>500) %>%
  ggplot(aes(y=Mean_RAD51_AreaShape_Eccentricity,x=Metadata_condition))+geom_boxplot()+facet_grid(.~Metadata_protein)

filter(nuclei,Intensity_IntegratedIntensity_EdU>500)


######
nuclei1 <- subset(nuclei,Metadata_protein=="dDBDdCTD"&Metadata_condition=="24hIR")

plot(log(nuclei1$Intensity_IntegratedIntensity_DAPI),log(nuclei1$Intensity_MeanIntensity_EdU),xlim=c(6,9))

nuclei$Intensity_IntegratedIntensity_DAPI)

#fold plot
mean_foci <- nuclei%>%
  filter(Intensity_IntegratedIntensity_EdU>500) %>%
  group_by(experiment,Metadata_protein,condition) %>%
  dplyr::summarize(RAD51_foci=mean(Children_RAD51_Count)) %>%
  group_by(Metadata_protein,condition) %>%
  dplyr::summarize(mean_foci=mean(RAD51_foci),sd_foci=sd(RAD51_foci),n=n())
write_clip(mean_foci)
#View(mean_foci)
mean_foci %>%
  ggplot(aes(y=mean_foci,x=condition)) + geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_foci-sd_foci, ymax=mean_foci+sd_foci), width=.2,
                position=position_dodge(.9)) +facet_grid(.~Metadata_protein)+ggtitle("Mean foci/experiment (n=3)")+
  scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+
  theme(legend.position = "none")


x <- nuclei%>%
  filter(Intensity_IntegratedIntensity_EdU>500) %>%
  group_by(experiment,Metadata_protein,condition) %>%
  dplyr::summarize(RAD51_foci=mean(Children_RAD51_Count)) %>%
  group_by(Metadata_protein,experiment) %>%
  dplyr::mutate(fold=RAD51_foci/dplyr::first(RAD51_foci)) %>%
  group_by(Metadata_protein,condition) %>%
  dplyr::summarize(mean_fold=mean(fold),sd_fold=sd(fold),n=n())

library(clipr)
write_clip(x)
%>%
ggplot(aes(x=condition,y=mean_fold,group=Metadata_protein,color=Metadata_protein))+
  geom_point()+geom_line()+geom_errorbar(aes(ymin=mean_fold-sd_fold, ymax=mean_fold+sd_fold))+
  facet_grid(.~Metadata_protein)+ scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+
  theme(legend.position = "none")


x <- nuclei %>%
  group_by(experiment,Metadata_protein,condition) %>%
  dplyr::mutate(EdUpos = Intensity_IntegratedIntensity_EdU>500) %>%
  dplyr::summarize(EdU = sum(EdUpos==TRUE)/n())%>%
  group_by(Metadata_protein,condition) %>%
  dplyr::summarize(mean_EdU=mean(EdU),sd_EdU=sd(EdU),n=n()) %>%
  ggplot(aes(y=mean_EdU,x=Metadata_protein)) + geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_EdU-sd_EdU, ymax=mean_EdU+sd_EdU), width=.2,
                position=position_dodge(.9)) +facet_wrap(.~condition)+
  scale_colour_Publication()+ylab("average fraction of EdU+ cells")+ylim(0,1)+
  scale_fill_Publication()+theme_Publication(base_size=16)+xlab("")+
  theme(legend.position = "none")
x
