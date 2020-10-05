# Function to visualize TIC x IT outliers per sigma group
# outliers identified based on the median absolute deviation (MAD) of TIC x IT values for each scan 

#Merve Öztoprak

TICxIT_outliers_visual <- function(light,heavy,analyte,replicate){
  ## First specify the packages of interest
  packages = c("plyr","dplyr","ggplot2","ggpubr","tidyverse","scales","cowplot",
               "reshape2")
  
  ## Now load or install&load all
  package.check <- lapply(packages,
                          FUN = function(x) {
                            if (!require(x, character.only = TRUE)) {
                              install.packages(x, dependencies = TRUE)
                              library(x, character.only = TRUE)
                            }
                          }
  )
  #read data
  data_light = read.csv(file=light,header = FALSE, na.strings = "NA", sep = ",", dec = ".", check.names=FALSE)
  data_heavy = read.csv(file=heavy,header=FALSE,na.strings = "NA",sep = ",",dec = ".",check.names=FALSE)
  data_light <- data_light[-(1:3),] %>% mutate_at(c(1:26, 28:32),as.numeric)#discard header with metadata and make values numeric
  data_heavy <- data_heavy[-(1:3),] %>% mutate_at(c(1:26, 28:32),as.numeric)#discard header with metadata and make values numeric
  
  #TICxIT correction Nsigma
  md_TICxIT_light<-median(data_light[,11]) #calculate median of TIC*IT
  md_TICxIT_heavy<-median(data_heavy[,11]) #calculate median of TIC*IT
  data_light$med_subt= abs(data_light[,11]-md_TICxIT_light) #inside function of MED=median(|xi-med(x)|) x=TIC*IT
  data_heavy$med_subt= abs(data_heavy[,11]-md_TICxIT_heavy) #inside function of MED=median(|xi-med(x)|) x=TIC*IT
  MADlight= median(data_light$med_subt)#outside funtion of MED=median(|xi-med(x)|) x=TIC*IT
  MADheavy= median(data_heavy$med_subt)#outside funtion of MED=median(|xi-med(x)|) x=TIC*IT
  k=1.4826 #constant defined for normally distributed data
  #outliers
  SDlight= k*MADlight
  SDheavy= k*MADheavy
  
  outlierslight_sub1sigma= subset(data_light,data_light$med_subt<1*SDlight)
  outliersheavy_sub1sigma= subset(data_heavy,data_heavy$med_subt<1*SDheavy)
  outlierslight_1sigma= subset(data_light,data_light$med_subt>1*SDlight&data_light$med_subt<2*SDlight)
  outliersheavy_1sigma= subset(data_heavy,data_heavy$med_subt>1*SDheavy&data_heavy$med_subt<2*SDheavy)
  outlierslight_2sigma= subset(data_light,data_light$med_subt>2*SDlight&data_light$med_subt<3*SDlight)
  outliersheavy_2sigma= subset(data_heavy,data_heavy$med_subt>2*SDheavy&data_heavy$med_subt<3*SDheavy)
  outlierslight_3sigma= subset(data_light,data_light$med_subt>3*SDlight)
  outliersheavy_3sigma= subset(data_heavy,data_heavy$med_subt>3*SDheavy)
  
  #TIC*IT visualizing median outliers
  #since TIC and IT are identical for mass light and heavy measured in same scan only outlier scans mass light are plotted
  #merges TIC*IT values for all sigma categories
  TICxIT_plot_light_up3sigma= data.frame(cbind(outlierslight_3sigma$V3,outlierslight_3sigma$V11))
  TICxIT_plot_light_sub3sigma= data.frame(cbind(outlierslight_2sigma$V3,outlierslight_2sigma$V11))
  TICxIT_plot_light_sub2sigma= data.frame(cbind(outlierslight_1sigma$V3,outlierslight_1sigma$V11))
  TICxIT_plot_light_sub1sigma= data.frame(cbind(outlierslight_sub1sigma$V3,outlierslight_sub1sigma$V11))
  colnames(TICxIT_plot_light_up3sigma)<- c("Retention Time (min)","up3_TIC")
  colnames(TICxIT_plot_light_sub3sigma)<- c("Retention Time (min)","sub3_TIC")
  colnames(TICxIT_plot_light_sub2sigma)<- c("Retention Time (min)","sub2_TIC")
  colnames(TICxIT_plot_light_sub1sigma)<- c("Retention Time (min)","sub1_TIC")
  merged_TICxIT_light= merge(TICxIT_plot_light_up3sigma,TICxIT_plot_light_sub3sigma,by=c("Retention Time (min)"),all=TRUE) #merges have to be nested since only 2 data frames permissable per merge
  merged_TICxIT_light= merge(merged_TICxIT_light,TICxIT_plot_light_sub2sigma,by=c("Retention Time (min)"),all=TRUE)
  merged_TICxIT_light= merge(merged_TICxIT_light,TICxIT_plot_light_sub1sigma,by=c("Retention Time (min)"),all=TRUE)
  melted_TICxIT_light= melt(merged_TICxIT_light,id.vars="Retention Time (min)")#melt for plotting
  melted_TICxIT_light= na.omit(melted_TICxIT_light)#remove NA rows
  melted_TICxIT_light= melted_TICxIT_light[order(melted_TICxIT_light$"Retention Time (min)"),]#sorts in ascending scan number order
  colnames(melted_TICxIT_light)<-c("Retention Time (min)","sigma_group","TIC*IT") 
  
  #scatter plot of TIC*IT values with different sigma groups in different colors
  TICxIT_outliers_plot <- ggplot(data=melted_TICxIT_light,aes(x=melted_TICxIT_light$"Retention Time (min)",y=melted_TICxIT_light$"TIC*IT",col=as.factor(melted_TICxIT_light$"sigma_group")))+
    geom_point() + 
    scale_x_continuous(breaks=seq(0,70,5))+
    theme_classic()+
    labs(col="TIC*IT by Sigma factor", x="Retention Time (min)", y="TICxIT")+
    scale_color_manual(labels = c(">3sigma", "<3sigma","<2sigma","<1sigma"), values = c("#CC79A7","#009E73","#E69F00","#000000"))+
    theme(legend.background = element_rect(fill="white"))+
    theme(legend.position = "top")
  TICxIT_outliers_plot_zoom <- TICxIT_outliers_plot +
    coord_cartesian(xlim = c(5, 15))+
    theme(legend.position = "none")
  TICxIT_outliers_plot_zoom2 <- TICxIT_outliers_plot +
    coord_cartesian(xlim = c(15, 70))+
    theme(legend.position = "none")+
    scale_y_continuous(labels = label_scientific(),breaks =seq(100000,300000,50000),limits=c(100000,300000))
  combined_TICxIT_plots= ggarrange(TICxIT_outliers_plot,# First row with first plot
                                   ggarrange(TICxIT_outliers_plot_zoom, TICxIT_outliers_plot_zoom2, ncol = 2, labels = c("B", "C")), # Second row plots
                                   nrow = 2, labels = "A"       # Label of the first plot
  ) 
  combined_TICxIT_plots  
  
  #outlier histogram plots
  outlierslight_all_sigma <-data.frame(sigma=c("<1sigma","1sigma","2sigma","3sigma and >"),Nscans=c(nrow(outlierslight_sub1sigma),nrow(outlierslight_1sigma),nrow(outlierslight_2sigma),nrow(outlierslight_3sigma)))
  #histogram of outliers per sigma bins
  sigma_plot <- ggplot(data=outlierslight_all_sigma,aes(x=sigma,y=Nscans))+
    geom_bar(stat="identity")+
    geom_text(aes(label=Nscans),vjust=-0.3, size=3)+
    scale_y_continuous(breaks=seq(0,10000,1000))+
    expand_limits(y=10000)+
    theme_classic()+
    theme(axis.title.x=element_blank())+
    theme(plot.title = element_text(hjust=0.5, vjust=3))+
    ggtitle("TIC x IT filtered for MAD")+
    theme(text = element_text(size=10),rect = element_rect(fill = "transparent"),axis.text.x = element_text(angle=45, hjust=1)) # all rectangles
  
  #TIC with sub1 sigma in grey, sub2sigma in green, sub3 sigma in blue and >3sigma in red
  #merges TIC values for all sigma categories
  TIC_plot_light_up3sigma= data.frame(cbind(outlierslight_3sigma$V3,outlierslight_3sigma$V9))
  TIC_plot_light_sub3sigma= data.frame(cbind(outlierslight_2sigma$V3,outlierslight_2sigma$V9))
  TIC_plot_light_sub2sigma= data.frame(cbind(outlierslight_1sigma$V3,outlierslight_1sigma$V9))
  TIC_plot_light_sub1sigma= data.frame(cbind(outlierslight_sub1sigma$V3,outlierslight_sub1sigma$V9))
  colnames(TIC_plot_light_up3sigma)<- c("Retention Time (min)","up3_TIC")
  colnames(TIC_plot_light_sub3sigma)<- c("Retention Time (min)","sub3_TIC")
  colnames(TIC_plot_light_sub2sigma)<- c("Retention Time (min)","sub2_TIC")
  colnames(TIC_plot_light_sub1sigma)<- c("Retention Time (min)","sub1_TIC")
  merged_TIC_light= merge(TIC_plot_light_up3sigma,TIC_plot_light_sub3sigma,by=c("Retention Time (min)"),all=TRUE) #merges have to be nested since only 2 data frames permissable per merge
  merged_TIC_light= merge(merged_TIC_light,TIC_plot_light_sub2sigma,by=c("Retention Time (min)"),all=TRUE)
  merged_TIC_light= merge(merged_TIC_light,TIC_plot_light_sub1sigma,by=c("Retention Time (min)"),all=TRUE)
  melted_TIC_light= melt(merged_TIC_light,id.vars="Retention Time (min)")#melt for plotting
  melted_TIC_light= na.omit(melted_TIC_light)#remove NA rows
  melted_TIC_light= melted_TIC_light[order(melted_TIC_light$"Retention Time (min)"),]#sorts in ascending scan number order
  colnames(melted_TIC_light)<-c("Retention Time (min)","sigma_group","TIC")
  
  #scatter plot of TIC values with different sigma groups in different colors
  ma= max(melted_TIC_light$"TIC")
  mi= min(melted_TIC_light$"TIC")
  TIC_plot <- ggplot(data=melted_TIC_light,aes(x=melted_TIC_light$"Retention Time (min)",y=melted_TIC_light$"TIC",col=as.factor(melted_TIC_light$"sigma_group")))+
    geom_point() + 
    scale_x_continuous(breaks=seq(0,70,5))+
    theme_classic()+
    scale_color_manual(labels = c(">3sigma", "<3sigma","<2sigma","<1sigma"), values = c("#CC79A7","#009E73","#E69F00","#000000"))+
    labs(col="TIC x IT filtered for MAD", x="Retention Time (min)", y="Total Ion Current (TIC)")+
    ylim(mi,ma)+ 
    annotation_custom(ggplotGrob(sigma_plot), xmin = 40, xmax = 70,ymin = 0.3*ma,ymax = ma)+
    theme(legend.position = "none")
  TIC_plot_wrapped <- ggplot(data=melted_TIC_light,aes(x=melted_TIC_light$"Retention Time (min)",y=melted_TIC_light$"TIC",col=as.factor(melted_TIC_light$"sigma_group")))+
    geom_point() + 
    scale_x_continuous(breaks=seq(0,70,5))+
    theme_classic()+
    scale_color_manual(labels = c(">3sigma", "<3sigma","<2sigma","<1sigma"), values = c("#CC79A7","#009E73","#E69F00","#000000"))+
    labs(col="TIC x IT filtered for MAD", x="Retention Time (min)", y="Total Ion Current (TIC)")+
    theme(legend.position = "none")+
    facet_wrap(~melted_TIC_light$"sigma_group",nrow=2)
  combined_TIC_plots <-ggarrange(TIC_plot,TIC_plot_wrapped,
                                 nrow = 2, labels = c("D", "E"))
  combined_TIC_plots
  
  #combined plots
  Full_combined_TICxIT_plots<-ggarrange(combined_TICxIT_plots,combined_TIC_plots,ncol = 2)
  Full_combined_TICxIT_plots
  
  #save
  dir.create(paste0("TICxIT_MAD_figs"), showWarnings = FALSE) #stops warnings if folder already exists
  ggsave(plot = Full_combined_TICxIT_plots,
         filename=file.path("TICxIT_MAD_figs", file=paste(analyte,replicate,"TICxIT_MAD_cleaned_plots.tiff")),
         width=300, height=200, units="mm", dpi=300, compression = "lzw")
}