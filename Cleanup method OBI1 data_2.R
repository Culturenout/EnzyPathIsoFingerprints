#General function to take OBIone output and produce R weighted values 

data_light = read.csv(file="1_light_MEP_mixed_raw.csv",header = FALSE, na.strings = "NA", sep = ",", dec = ".", check.names=FALSE)
data_heavy = read.csv(file="1_heavy_MEP_mixed_raw.csv",header=FALSE,na.strings = "NA",sep = ",",dec = ".",check.names=FALSE)

#sub in output and col names light and heavy for mass in function

###generalized function###
Corrected_Orbi_output<- function(filename1,filename2,analyte,replicate,turntime2){
  
  ## First specify the packages of interest
  packages = c("plyr","dplyr","ggplot2","ggpubr","tidyverse","scales","cowplot",
               "reshape2", "viridis")
  
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
  data_light = read.csv(file=filename1,header = FALSE, na.strings = "NA", sep = ",", dec = ".", check.names=FALSE)
  data_heavy = read.csv(file=filename2,header=FALSE,na.strings = "NA",sep = ",",dec = ".",check.names=FALSE)
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
  
  #output csv corrected for <1, <2 sigma TIC*IT
  clean_data_light_sub1sigma<-subset(data_light,data_light$med_subt<1*SDlight)#goal
  clean_data_heavy_sub1sigma<-subset(data_heavy,data_heavy$med_subt<1*SDheavy)#goal
  clean_data_light_sub2sigma<-subset(data_light,data_light$med_subt<2*SDlight)#goal
  clean_data_heavy_sub2sigma<-subset(data_heavy,data_heavy$med_subt<2*SDheavy)#goal
  clean_data_light_sub1sigma<-clean_data_light_sub1sigma[,c(-1,-34)]#remove med_subt and unnecessary column
  clean_data_heavy_sub1sigma<-clean_data_heavy_sub1sigma[,c(-1,-34)]#remove med_subt and unnecessary column
  clean_data_light_sub2sigma<-clean_data_light_sub2sigma[,c(-1,-34)]#remove med_subt and unnecessary column
  clean_data_heavy_sub2sigma<-clean_data_heavy_sub2sigma[,c(-1,-34)]#remove med_subt and unnecessary column
  colnames(clean_data_light_sub1sigma)<-c("Measured Mass:","Ret. Time:","Scan Number:","Deviation [ppm]:","Deviation [mmu]:","Abs. Intensity:","Rel. Intensity:",
                                      "TIC:","IT [ms]:","TIC*IT:","Elapsed scan time [ms]:","Target:","FT Resolution:","RF [V]:","Lockmass correction [ppm]:","Lockmass #1 [m/z]:",
                                      "Lockmass #2 [m/z]:","Lockmass #3 [m/z]:","Number of Packets:","A","B","C","Vt","Mode","Peak Noise","Peak Flags","Peak Resolution","Peak Baseline",
                                      "Ion Injection Time (ms):","Max. Ion Time (ms):","Repeller Voltage (V)","CI reagent gas port:")
  colnames(clean_data_heavy_sub1sigma)<-c("Measured Mass:","Ret. Time:","Scan Number:","Deviation [ppm]:","Deviation [mmu]:","Abs. Intensity:","Rel. Intensity:",
                                      "TIC:","IT [ms]:","TIC*IT:","Elapsed scan time [ms]:","Target:","FT Resolution:","RF [V]:","Lockmass correction [ppm]:","Lockmass #1 [m/z]:",
                                      "Lockmass #2 [m/z]:","Lockmass #3 [m/z]:","Number of Packets:","A","B","C","Vt","Mode","Peak Noise","Peak Flags","Peak Resolution","Peak Baseline",
                                      "Ion Injection Time (ms):","Max. Ion Time (ms):","Repeller Voltage (V)","CI reagent gas port:")
  colnames(clean_data_light_sub2sigma)<-c("Measured Mass:","Ret. Time:","Scan Number:","Deviation [ppm]:","Deviation [mmu]:","Abs. Intensity:","Rel. Intensity:",
                                      "TIC:","IT [ms]:","TIC*IT:","Elapsed scan time [ms]:","Target:","FT Resolution:","RF [V]:","Lockmass correction [ppm]:","Lockmass #1 [m/z]:",
                                      "Lockmass #2 [m/z]:","Lockmass #3 [m/z]:","Number of Packets:","A","B","C","Vt","Mode","Peak Noise","Peak Flags","Peak Resolution","Peak Baseline",
                                      "Ion Injection Time (ms):","Max. Ion Time (ms):","Repeller Voltage (V)","CI reagent gas port:")
  colnames(clean_data_heavy_sub2sigma)<-c("Measured Mass:","Ret. Time:","Scan Number:","Deviation [ppm]:","Deviation [mmu]:","Abs. Intensity:","Rel. Intensity:",
                                      "TIC:","IT [ms]:","TIC*IT:","Elapsed scan time [ms]:","Target:","FT Resolution:","RF [V]:","Lockmass correction [ppm]:","Lockmass #1 [m/z]:",
                                      "Lockmass #2 [m/z]:","Lockmass #3 [m/z]:","Number of Packets:","A","B","C","Vt","Mode","Peak Noise","Peak Flags","Peak Resolution","Peak Baseline",
                                      "Ion Injection Time (ms):","Max. Ion Time (ms):","Repeller Voltage (V)","CI reagent gas port:")
  print(paste("Datapoints Raw data 1)light 2)heavy:",NROW(data_light),NROW(data_heavy)))
  print(paste("Datapoints TICxIT_cleaned Cleaned up data [(|xi-med(x)|)< 1sigma] 1)light 2)heavy:",NROW(clean_data_light_sub1sigma),NROW(clean_data_heavy_sub1sigma)))
  print(paste("Datapoints TICxIT_cleaned Cleaned up data [(|xi-med(x)|)< 2sigma] 1)light 2)heavy:",NROW(clean_data_light_sub2sigma),NROW(clean_data_heavy_sub2sigma)))
  print(paste("Datapoints TICxIT_cleaned Outliers [(|xi-med(x)|)< 1sigma 1)light 2)heavy:",NROW(outlierslight_sub1sigma),NROW(outliersheavy_sub1sigma)))
  print(paste("Datapoints TICxIT_cleaned Outliers [(|xi-med(x)|)< 2sigma 1)light 2)heavy:",NROW(outlierslight_1sigma),NROW(outliersheavy_1sigma)))
  
  #merge data tables
  data_light_sub1sigma= cbind(clean_data_light_sub1sigma$"Ret. Time:",clean_data_light_sub1sigma$"Scan Number:",clean_data_light_sub1sigma$"Abs. Intensity:",clean_data_light_sub1sigma$"Peak Noise")
  data_heavy_sub1sigma= cbind(clean_data_heavy_sub1sigma$"Ret. Time:",clean_data_heavy_sub1sigma$"Scan Number:",clean_data_heavy_sub1sigma$"Abs. Intensity:",clean_data_heavy_sub1sigma$"Peak Noise")
  data_light_sub2sigma= cbind(clean_data_light_sub2sigma$"Ret. Time:",clean_data_light_sub2sigma$"Scan Number:",clean_data_light_sub2sigma$"Abs. Intensity:",clean_data_light_sub2sigma$"Peak Noise")
  data_heavy_sub2sigma= cbind(clean_data_heavy_sub2sigma$"Ret. Time:",clean_data_heavy_sub2sigma$"Scan Number:",clean_data_heavy_sub2sigma$"Abs. Intensity:",clean_data_heavy_sub2sigma$"Peak Noise")
  colnames(data_light_sub1sigma)<- c("Ret. Time:","Scan Number:","light Abs. Intensity:","light Peak Noise")
  colnames(data_heavy_sub1sigma)<- c("Ret. Time:","Scan Number:","heavy Abs. Intensity:","heavy Peak Noise")
  colnames(data_light_sub2sigma)<- c("Ret. Time:","Scan Number:","light Abs. Intensity:","light Peak Noise")
  colnames(data_heavy_sub2sigma)<- c("Ret. Time:","Scan Number:","heavy Abs. Intensity:","heavy Peak Noise")
  
  merged_data_sub1sigma= merge(data_light_sub1sigma,data_heavy_sub1sigma,By="Scan Number:")
  merged_data_sub2sigma= merge(data_light_sub2sigma,data_heavy_sub2sigma,By="Scan Number:")
  merged_data_sub1sigma$count_light= (merged_data_sub1sigma$"light Abs. Intensity:"/merged_data_sub1sigma$"light Peak Noise")*4.24
  merged_data_sub1sigma$count_heavy= (merged_data_sub1sigma$"heavy Abs. Intensity:"/merged_data_sub1sigma$"heavy Peak Noise")*4.24
  merged_data_sub2sigma$count_light= (merged_data_sub2sigma$"light Abs. Intensity:"/merged_data_sub2sigma$"light Peak Noise")*4.24
  merged_data_sub2sigma$count_heavy= (merged_data_sub2sigma$"heavy Abs. Intensity:"/merged_data_sub2sigma$"heavy Peak Noise")*4.24
  merged_data_sub1sigma$NL_ratios= (merged_data_sub1sigma$"heavy Abs. Intensity:"/merged_data_sub1sigma$"light Abs. Intensity:")
  merged_data_sub2sigma$NL_ratios= (merged_data_sub2sigma$"heavy Abs. Intensity:"/merged_data_sub2sigma$"light Abs. Intensity:")
  merged_data_sub1sigma$count_ratios= (merged_data_sub1sigma$count_heavy/merged_data_sub1sigma$count_light)
  merged_data_sub2sigma$count_ratios= (merged_data_sub2sigma$count_heavy/merged_data_sub2sigma$count_light)
  merged_data_sub1sigma$weighted_count_light= (merged_data_sub1sigma$count_light*merged_data_sub1sigma$"light Abs. Intensity:")/sum(merged_data_sub1sigma$"light Abs. Intensity:")
  merged_data_sub1sigma$weighted_count_heavy= (merged_data_sub1sigma$count_heavy*merged_data_sub1sigma$"heavy Abs. Intensity:")/sum(merged_data_sub1sigma$"heavy Abs. Intensity:")
  merged_data_sub2sigma$weighted_count_light= (merged_data_sub2sigma$count_light*merged_data_sub2sigma$"light Abs. Intensity:")/sum(merged_data_sub2sigma$"light Abs. Intensity:")
  merged_data_sub2sigma$weighted_count_heavy= (merged_data_sub2sigma$count_heavy*merged_data_sub2sigma$"heavy Abs. Intensity:")/sum(merged_data_sub2sigma$"heavy Abs. Intensity:")
  merged_data_sub1sigma$R= merged_data_sub1sigma$weighted_count_heavy/merged_data_sub1sigma$weighted_count_light
  merged_data_sub2sigma$R= merged_data_sub2sigma$weighted_count_heavy/merged_data_sub2sigma$weighted_count_light
  
  #plot TIC*IT corrected NL 
  data.frame(merged_data_sub1sigma)
  data.frame(merged_data_sub2sigma)
  naive_NL_1_sigma_plot <- ggplot(data=merged_data_sub1sigma,aes(y=merged_data_sub1sigma[,9],x=merged_data_sub1sigma[,2]))+
    geom_point() + 
    scale_x_continuous(breaks=seq(0,15000,2000))+
    theme_classic()+
    ggtitle("TIC*IT MAD <1sigma filtered NL ratio")+
    xlab("Scan Number:")+
    ylab("NL_ratios")+
    theme(legend.position = "none")
  naive_NL_2_sigma_plot <- ggplot(data=merged_data_sub2sigma,aes(y=merged_data_sub2sigma[,9],x=merged_data_sub2sigma[,2]))+
    geom_point() + 
    scale_x_continuous(breaks=seq(0,15000,2000))+
    theme_classic()+
    ggtitle("TIC*IT MAD <2sigma filtered NL ratio")+
    xlab("Scan Number:")+
    ylab("NL_ratios")+
    theme(legend.position = "none")
  combined_naive_NL_plot <-ggarrange(naive_NL_1_sigma_plot,naive_NL_2_sigma_plot,
                                     ncol = 2, labels = c("A", "B"))
  #only include scans 1 min after valve turn time
  Analysis_window_1sigma= subset(merged_data_sub1sigma,merged_data_sub1sigma$"Ret. Time:">turntime2+1)
  Analysis_window_2sigma= subset(merged_data_sub2sigma,merged_data_sub2sigma$"Ret. Time:">turntime2+1)
  NL_1_sigma_plot <- ggplot(data=Analysis_window_1sigma,aes(y=Analysis_window_1sigma[,9],x=Analysis_window_1sigma[,2]))+
    geom_point() + 
    scale_x_continuous(breaks=seq(0,15000,2000))+
    theme_classic()+
    ggtitle("TIC*IT MAD <1sigma filtered NL ratio only valve turn time + 1 min")+
    xlab("Scan Number:")+
    ylab("NL_ratios")+
    theme(legend.position = "none")
  NL_2_sigma_plot <- ggplot(data=Analysis_window_2sigma,aes(y=Analysis_window_2sigma[,9],x=Analysis_window_2sigma[,2]))+
    geom_point() + 
    scale_x_continuous(breaks=seq(0,15000,2000))+
    theme_classic()+
    ggtitle("TIC*IT MAD <2sigma filtered NL ratio only valve turn time + 1 min")+
    xlab("Scan Number:")+
    ylab("NL_ratios")+
    theme(legend.position = "none")
  combined_NL_plot <-ggarrange(NL_1_sigma_plot,NL_2_sigma_plot,
                               ncol = 2, labels = c("D", "E"))
  #combined plots
  Full_combined_NL_plots<-ggarrange(combined_naive_NL_plot,combined_NL_plot,nrow = 2)
  Full_combined_NL_plots
  
  #save
  dir.create(paste0("NL_figs"), showWarnings = FALSE) #stops warnings if folder already exists
  ggsave(plot = Full_combined_NL_plots,
         filename=file.path("NL_figs", file=paste(analyte,replicate,"NL_plots.tiff")),
         width=300, height=200, units="mm", dpi=300, compression = "lzw")
  
  #save data tables filtered for TICxIT MAD outliers and scans 1 min after valve turn time
  dir.create(paste0("1sigma_TICxIT_and_time_filtered_data"), showWarnings = FALSE) #stops warnings if folder already exists
  dir.create(paste0("2sigma_TICxIT_and_time_filtered_data"), showWarnings = FALSE) #stops warnings if folder already exists
  filtered_light_1sigma=subset(clean_data_light_sub1sigma,clean_data_light_sub1sigma[,2]>turntime2+1)
  filtered_light_2sigma=subset(clean_data_light_sub2sigma,clean_data_light_sub2sigma[,2]>turntime2+1)
  filtered_heavy_1sigma=subset(clean_data_heavy_sub1sigma,clean_data_heavy_sub1sigma[,2]>turntime2+1)
  filtered_heavy_2sigma=subset(clean_data_heavy_sub2sigma,clean_data_heavy_sub2sigma[,2]>turntime2+1)
  write.csv(filtered_light_1sigma,file.path(paste0("1sigma_TICxIT_and_time_filtered_data"),file=paste(analyte,replicate,"light_1sigma_filtered_raw_data.csv")),row.names = FALSE)
  write.csv(filtered_light_2sigma,file.path(paste0("2sigma_TICxIT_and_time_filtered_data"),file=paste(analyte,replicate,"light_2sigma_filtered_raw_data.csv")),row.names = FALSE)
  write.csv(filtered_heavy_1sigma,file.path(paste0("1sigma_TICxIT_and_time_filtered_data"),file=paste(analyte,replicate,"heavy_1sigma_filtered_raw_data.csv")),row.names = FALSE)
  write.csv(filtered_heavy_2sigma,file.path(paste0("2sigma_TICxIT_and_time_filtered_data"),file=paste(analyte,replicate,"heavy_2sigma_filtered_raw_data.csv")),row.names = FALSE)
  
  #calculate summary values and weighted average R
  sd.p=function(x){sd(x)*sqrt((length(x)-1)/length(x))} # dtandard deviation of the population function
  
  Average_R_1sigma=mean(Analysis_window_1sigma$count_ratios)#Average R
  StdDevR_1sigma=sd.p(Analysis_window_1sigma$NL_ratios)
  StdErr_1sigma=StdDevR_1sigma/sqrt(NROW(Analysis_window_1sigma$NL_ratios))
  RelStdErr_1sigma=StdErr_1sigma/Average_R_1sigma
  Avrg_Weigh_R_1sigma=sum(Analysis_window_1sigma$weighted_count_heavy)/sum(Analysis_window_1sigma$weighted_count_light)
  Shot_Noise_Limit_1sigma=sqrt(1/sum(Analysis_window_1sigma$count_heavy))
  RelStdErr_NL_div_ShotNoise_1sigma=RelStdErr_1sigma/Shot_Noise_Limit_1sigma
  perC_1sigma=Avrg_Weigh_R_1sigma/6 #calculate div6 column(=per carbon Weighted R value)
  
  Average_R_2sigma=mean(Analysis_window_2sigma$count_ratios)#Average R
  StdDevR_2sigma=sd.p(Analysis_window_2sigma$NL_ratios)
  StdErr_2sigma=StdDevR_2sigma/sqrt(NROW(Analysis_window_2sigma$NL_ratios))
  RelStdErr_2sigma=StdErr_2sigma/Average_R_2sigma
  Avrg_Weigh_R_2sigma=sum(Analysis_window_2sigma$weighted_count_heavy)/sum(Analysis_window_2sigma$weighted_count_light)
  Shot_Noise_Limit_2sigma=sqrt(1/sum(Analysis_window_2sigma$count_heavy))
  RelStdErr_NL_div_ShotNoise_2sigma=RelStdErr_2sigma/Shot_Noise_Limit_2sigma
  perC_2sigma=Avrg_Weigh_R_2sigma/6 #calculate div6 column(=per carbon Weighted R value)
  
  #Bulkcorrection
  #BulkC for Std
  x=0.0108998945729143
  #BulkC for MEP
  y=0.0107463731505
  #BulkC for MVA
  z=0.0109599248993
  Bulkcorr_1sigma=ifelse(analyte=="Std",(6*perC_1sigma)-(5*x),
                         ifelse(analyte=="MEP",(6*perC_1sigma)-(5*y),
                                ifelse(analyte=="MVA",(6*perC_1sigma)-(5*z),0
                                )
                         )
  )
  Bulkcorr_2sigma=ifelse(analyte=="Std",(6*perC_2sigma)-(5*x),
                         ifelse(analyte=="MEP",(6*perC_2sigma)-(5*y),
                                ifelse(analyte=="MVA",(6*perC_2sigma)-(5*z),0
                                )
                         )
  )
  
  Weighted_R_values_1sigma<-data.frame(ID=analyte,Replicate=replicate,Average_R_1sigma,StdDevR_1sigma,StdErr_1sigma,RelStdErr_1sigma,Avrg_Weigh_R_1sigma,Shot_Noise_Limit_1sigma,RelStdErr_NL_div_ShotNoise_1sigma,perC_1sigma,Bulkcorr_1sigma)
  Weighted_R_values_2sigma<-data.frame(ID=analyte,Replicate=replicate,Average_R_2sigma,StdDevR_2sigma,StdErr_2sigma,RelStdErr_2sigma,Avrg_Weigh_R_2sigma,Shot_Noise_Limit_2sigma,RelStdErr_NL_div_ShotNoise_2sigma,perC_2sigma,Bulkcorr_2sigma)
  
  dir.create(paste0("Weighted R values calc 1sigma"), showWarnings = FALSE) #stops warnings if folder already exists
  dir.create(paste0("Weighted R values calc 2sigma"), showWarnings = FALSE) #stops warnings if folder already exists
  write.csv(Weighted_R_values_1sigma, file.path(paste0("Weighted R values calc 1sigma"), file=paste(analyte,replicate,"weighted_R_calc_1sigma.csv")), row.names = FALSE)
  print(head(Weighted_R_values_1sigma))
  write.csv(Weighted_R_values_2sigma, file.path(paste0("Weighted R values calc 2sigma"), file=paste(analyte,replicate,"weighted_R_calc_2sigma.csv")), row.names = FALSE)
  
  ###make summarized data table###
  setwd("C:/Users/moztoprak/OneDrive - NIOZ/PhD/2020/Phytane project/Statistics/Mixed day measurements/Weighted R values calc 1sigma")
  myfiles<-list.files(path = ".", pattern = "csv")
  All_measurements<-ldply(myfiles,read.csv)
  All_measurements %>% arrange(ID)
  dir.create(paste0("summarized values 1sigma"), showWarnings = FALSE) #stops warnings if folder already exists
  write.csv(All_measurements,file.path(paste0("summarized values 1sigma"), file=paste("Weighted_R_summerized_1sigma.csv")), row.names = FALSE)
  
  setwd("C:/Users/moztoprak/OneDrive - NIOZ/PhD/2020/Phytane project/Statistics/Mixed day measurements/Weighted R values calc 2sigma")
  myfiles<-list.files(path = ".", pattern = "csv")
  All_measurements<-ldply(myfiles,read.csv)
  All_measurements %>% arrange(ID)
  dir.create(paste0("summarized values 2sigma"), showWarnings = FALSE) #stops warnings if folder already exists
  write.csv(All_measurements,file.path(paste0("summarized values 2sigma"), file=paste("Weighted_R_summerized_2sigma.csv")), row.names = FALSE)
  
  setwd("C:/Users/moztoprak/OneDrive - NIOZ/PhD/2020/Phytane project/Statistics/Mixed day measurements")
}


# look up and set directory
setwd("C:/Users/moztoprak/OneDrive - NIOZ/PhD/2020/Phytane project/Statistics/Mixed day measurements")
list.files(path = ".", pattern = "raw")

Corrected_Orbi_output("1_light_Phy_std_mixed_raw.csv","1_heavy_Phy_std_mixed_raw.csv","Std","1",10)
Corrected_Orbi_output("1_light_MEP_mixed_raw.csv","1_heavy_MEP_mixed_raw.csv","MEP","1",10)
Corrected_Orbi_output("1_light_MVA_mixed_raw.csv","1_heavy_MVA_mixed_raw.csv","MVA","1",10)

Corrected_Orbi_output("2_light_Phy_std_mixed_raw.csv","2_heavy_Phy_std_mixed_raw.csv","Std","2",10)
Corrected_Orbi_output("2_light_MEP_mixed_raw.csv","2_heavy_MEP_mixed_raw.csv","MEP","2",10)
Corrected_Orbi_output("2_light_MVA_mixed_raw.csv","2_heavy_MVA_mixed_raw.csv","MVA","2",10)

Corrected_Orbi_output("3_light_Phy_std_mixed_raw.csv","3_heavy_Phy_std_mixed_raw.csv","Std","3",10)
Corrected_Orbi_output("3_light_MEP_mixed_raw.csv","3_heavy_MEP_mixed_raw.csv","MEP","3",10)
Corrected_Orbi_output("3_light_MVA_mixed_raw.csv","3_heavy_MVA_mixed_raw.csv","MVA","3",10)

##### add average for every column row for every analyte to summary output table maybe just do it in excel?
### standardize code somehow take the working directory out at the end