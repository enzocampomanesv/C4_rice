library(ggplot2)
# devtools::install_github('Mikata-Project/ggthemr')
library(ggthemr)
library(RColorBrewer)
library(cowplot)
library(zoo)
library(tidyverse)
library(zeallot)
ggthemr("fresh")

setwd("/home/jovyan/private/C4_Rice/")
rm(list=ls(all=TRUE))
all_rcp <- c('baseline','SSP1-2.6','SSP2-4.5','SSP4-6.0','SSP5-8.5')

data_prep <- function(rcp, year){
  if(rcp=="baseline"){
    c3 <- read.csv(paste0("Data/Model results/Crop_yield_final/",rcp,"/latitudinal/C3_",rcp,"_riceOnly.csv"))
    c4 <- read.csv(paste0("Data/Model results/Crop_yield_final/",rcp,"/latitudinal/C4_",rcp,"_riceOnly.csv"))
  }else{
    fname_rcp = gsub("[.-]", "", rcp)
    print(substr(fname_rcp,5,6))
    c3 <- read.csv(paste0("Data/Model results/Crop_yield_final/",fname_rcp,"/latitudinal/C3_",year,"_",substr(fname_rcp,5,6),"_riceOnly.csv"))
    c4 <- read.csv(paste0("Data/Model results/Crop_yield_final/",fname_rcp,"/latitudinal/C4_",year,"_",substr(fname_rcp,5,6),"_riceOnly.csv"))
  }
  c3 <- c3[c3$avg_yld!="#N/A",]
  c4 <- c4[c4$avg_yld!="#N/A",]
  names(c3)[3] <- "yld"
  names(c4)[3] <- "yld"
  c3$type <- 'C3'
  c4$type <- 'C4'
  
  c3$yld <- as.double(c3$yld)
  c4$yld <- as.double(c4$yld)
  c3$sd_yld <- as.double(c3$sd_yld)
  c4$sd_yld <- as.double(c4$sd_yld)
  
  combined <- rbind(c3,c4)
  combined$sd_min <- combined$yld-combined$sd_yld
  combined$sd_max <- combined$yld+combined$sd_yld
  return(list(combined,c3,c4))
}

plot_all_yield <- function(rcp_list, year){
  plot_yield <- function(rcp, year){
    data_list <- data_prep(rcp,year)
    combined <- data_list[[1]]
    
    latxYield_combined <- combined %>%
      ggplot(aes(x=Latitude,y=rollmean(yld, 100,na.pad=TRUE)))+
      geom_line(aes(color=type),size=0.3, alpha=1)+
      geom_ribbon(aes(ymin=sd_min,ymax=sd_max,fill=type),alpha=0.2)+
      geom_vline(xintercept=0, color="black", size=1, alpha=0.3)+
      scale_y_continuous(limits=c(0,30), breaks=(seq(from=0,to=30,by=10)))+
      scale_x_continuous(limits=c(-15,60),breaks=(seq(from=-15,to=60,by=15)))+
      ylab("Yield (t/ha)")+
      xlab("Latitude")+
      coord_flip()+
      theme(legend.title= element_blank(),
            legend.position = c(0.8, 0.7))
    
    if(rcp=="baseline"){
      latxYield_combined <- latxYield_combined+
        ggtitle(paste0("Baseline"))
    }else{
      # rcp <- paste0('SSP ',substr(rcp,5,5),'.',substr(rcp,6,6))
      latxYield_combined <- latxYield_combined+
        ggtitle(paste0(rcp," (",year, ")"))
    }
    latxYield_combined <- latxYield_combined+
      theme(plot.title = element_text(size=12))
    return(latxYield_combined)
  }
  
  data_plots <- lapply(rcp_list, plot_yield, year=year)
  
  #adjust axis titles in each plot
  data_plots[[2]] <- data_plots[[2]]+
    theme(
      axis.title.y=element_blank(),
      axis.text.y=element_blank())+
    guides(fill = "none", color="none")
  
  data_plots[[3]] <- data_plots[[3]]+
    theme(
      axis.title.y=element_blank(),
      axis.text.y=element_blank())+
    guides(fill = "none", color="none")
  data_plots[[4]] <- data_plots[[4]]+
    guides(fill = "none", color="none")
  data_plots[[5]] <- data_plots[[5]]+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank())+
    guides(fill = "none", color="none")
  
  data_list <- data_prep('baseline',year)
  combined <- data_list[[1]]
  num_pix_plot <- combined %>%
    ggplot(aes(x=Latitude,y=rollmean(num_pix, 100,na.pad=TRUE)))+
    geom_bar(stat='identity',width=0.2,show.legend = FALSE)+
    ylab("Number of pixels")+
    coord_flip()+
    ggtitle(paste0("Pixel count"))+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          plot.title = element_text(size=12))
  
  return(plot_grid(data_plots[[1]],data_plots[[2]],
            data_plots[[3]],data_plots[[4]],
            data_plots[[5]],num_pix_plot,ncol=3))
}

plot_all_yield(all_rcp, 2050)
ggsave('./Data/Model results/Crop_yield_final/plots/latPlot_2050_Topt5_v2.png',
       dpi=300)
plot_all_yield(all_rcp, 2099)
ggsave('./Data/Model results/Crop_yield_final/plots/latPlot_2099_Topt5_v2.png',
       dpi=300)


# d_list <- data_prep('SSP585', 2099)
# 
# combined <- d_list[[1]]
# c3 <- d_list[[2]]
# c4 <- d_list[[3]]

calculate_diff <- function(rcp, year){
  df_list <- data_prep(rcp, year)
  combined <- df_list[[1]]
  c3 <- df_list[[2]]
  c4 <- df_list[[3]]
  names(c3)[3] <- "c3_yld"
  names(c4)[3] <- "c4_yld"
  names(c3)[4] <- "c3_sd_yld"
  names(c4)[4] <- "c4_sd_yld"
  diff <- left_join(c3[,c(2:4)],c4[,c(2:4)], by="Latitude")
  diff$diff <- diff$c4_yld-diff$c3_yld
  diff$sd_min <- ifelse(diff$diff < 0, diff$diff-diff$c3_sd_yld, diff$diff-diff$c4_sd_yld)
  diff$sd_max <- ifelse(diff$diff < 0, diff$diff+diff$c4_sd_yld, diff$diff+diff$c3_sd_yld)
  
  return(diff)
}
define_lat_change <- function(df){
  for(i in seq((nrow(df)-1))){
    if((df[i,c("diff")] < 0) & (df[i+1,c("diff")] >= 0)){
      return(df[i,1])
    }
  }
  return(0)
}

# diff1 <- calculate_diff('baseline', 2050)
# lat_ch <- define_lat_change(diff)

lat_change <- data.frame(matrix(nrow=(length(all_rcp)*2),ncol=3))
names(lat_change) <- c('yr', 'rcp', 'latitude_of_shift')
count <- 1
for(yr in c(2050,2099)){
  for(rcp in all_rcp){
    temp_df <- calculate_diff(rcp, yr)
    lat_change[count,] <- c(as.numeric(yr),rcp, as.numeric(define_lat_change(temp_df)))
    count <- count+1
  }
}
write.csv(lat_change, './Data/Model results/Crop_yield_final/plots/lat_change_topt5.csv', row.names = FALSE)

############ DIFFERENCE PLOT #############
# plot_diff <- function(rcp, year){
#   diff <- calculate_diff(rcp, year)
#   # df_list <- data_prep(rcp, year)
#   # combined <- df_list[[1]]
#   # c3 <- df_list[[2]]
#   # c4 <- df_list[[3]]
#   # names(c3)[3] <- "c3_yld"
#   # names(c4)[3] <- "c4_yld"
#   # 
#   # diff <- left_join(c3[,c(2:3)],c4[,c(2,3)], by="Latitude")
#   # diff$diff <- diff$c4_yld-diff$c3_yld
#   lat_change <- define_lat_change(diff)
#   diff_plot <- diff %>%
#     ggplot(aes(x=Latitude,y=rollmean(diff, 100,na.pad=TRUE)))+
#     geom_ribbon(aes(ymin=sd_min,ymax=sd_max), fill='dark green',alpha=0.2, show.legend=FALSE)+
#     geom_line(size=0.5,color = 'dark green', n=2)+
#     scale_y_continuous(limits=c(-12,12), breaks=(seq(from=-12,to=12,by=3)))+
#     scale_x_continuous(limits=c(-15,60),breaks=(seq(from=-15,to=60,by=15)))+
#     geom_hline(yintercept=0, color="black", linetype="dashed", size=1, alpha=0.8)+
#     # geom_vline(xintercept=lat_change, color="red", size=0.4, alpha=0.7)+
#     ylab("Yield difference (t/ha)")+
#     xlab("Latitude")+
#     coord_flip()+
#     ggtitle(paste0(rcp,"-",year))
#   
#   if(rcp=="baseline"){
#     diff_plot <- diff_plot+
#       ggtitle(paste0("Baseline"))
#   }else{
#     rcp <- paste0('RCP ',substr(rcp,5,5),'.',substr(rcp,6,6))
#     diff_plot <- diff_plot+
#       ggtitle(paste0(rcp,"-",year))
#   }
#   return(list(diff_plot, diff))
# }
# c(plt, diff) %<-% plot_diff('SSP585', 2050)
plot_all_diff <- function(rcp_list, year){
  plot_diff <- function(rcp, year){
    diff <- calculate_diff(rcp, year)
    lat_change <- define_lat_change(diff)
    diff_plot <- diff %>%
      ggplot(aes(x=Latitude,y=rollmean(diff, 100,na.pad=TRUE)))+
      geom_ribbon(aes(ymin=sd_min,ymax=sd_max), fill='dark green',alpha=0.2, show.legend=FALSE)+
      geom_line(size=0.5,color = 'dark green', n=2)+
      scale_y_continuous(limits=c(-6,6), breaks=(seq(from=-6,to=6,by=3)))+
      scale_x_continuous(limits=c(-15,60),breaks=(seq(from=-15,to=60,by=15)))+
      geom_hline(yintercept=0, color="black", linetype="dashed", size=1, alpha=0.8)+
      # geom_vline(xintercept=lat_change, color="red", size=0.4, alpha=0.7)+
      ylab("Yield difference (t/ha)")+
      xlab("Latitude")+
      coord_flip()+
      ggtitle(paste0(rcp,"-",year))

    if(rcp=="baseline"){
      diff_plot <- diff_plot+
        ggtitle(paste0("Baseline"))
    }else{
      # rcp <- paste0('RCP ',substr(rcp,5,5),'.',substr(rcp,6,6))
      diff_plot <- diff_plot+
        ggtitle(paste0(rcp," (",year,")"))
    }
    write.csv(diff, paste0('./Data/Model results/Crop_yield_final/plots/diff_plot_',rcp,'_',year,'_topt5.csv'), row.names = FALSE)
    return(diff_plot)
  }
  data_plots <- lapply(rcp_list, plot_diff, year=year)
  data_plots[[2]] <- data_plots[[2]]+
    theme(
      axis.title.y=element_blank(),
      axis.text.y=element_blank())+
    guides(fill = "none", color="none")
  
  data_plots[[3]] <- data_plots[[3]]+
    theme(
      axis.title.y=element_blank(),
      axis.text.y=element_blank())+
    guides(fill = "none", color="none")
  data_plots[[4]] <- data_plots[[4]]+
    guides(fill = "none", color="none")
  data_plots[[5]] <- data_plots[[5]]+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank())+
    guides(fill = "none", color="none")
  
  data_list <- data_prep('baseline',year)
  combined <- data_list[[1]]
  num_pix_plot <- combined %>%
    ggplot(aes(x=Latitude,y=rollmean(num_pix, 100,na.pad=TRUE)))+
    geom_bar(stat='identity',width=0.2,show.legend = FALSE)+
    ylab("Number of pixels")+
    coord_flip()+
    ggtitle(paste0("Pixel count"))+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          plot.title = element_text(size=12))
  
  return(plot_grid(data_plots[[1]],data_plots[[2]],
                   data_plots[[3]],data_plots[[4]],
                   data_plots[[5]],num_pix_plot,ncol=3))
}

plot_all_diff(all_rcp, 2050)
ggsave('./Data/Model results/Crop_yield_final/plots/diffPlot_2050_latLineChange_Topt5_v2.png',
       width=9,height=8,units="in", dpi=300)
plot_all_diff(all_rcp, 2099)
ggsave('./Data/Model results/Crop_yield_final/plots/diffPlot_2099_latLineChange_Topt5_v2.png',
       width=9,height=8,units="in", dpi=300)

