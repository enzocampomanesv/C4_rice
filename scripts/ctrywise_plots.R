library(ggplot2)
library(ggthemr)
library(RColorBrewer)
library(cowplot)
library(sf)
library(raster)
library(tidyverse)
library(rgdal)
ggthemr("fresh")

setwd("/home/jovyan/private/C4_Rice/")
rm(list=ls(all=TRUE))
rcp_list <- c('baseline','SSP1-2.6','SSP2-4.5','SSP4-6.0','SSP5-8.5')
tp_list <- c('C3','C4')
yr_list <- c('2050','2099')
ctry_list <- c('BGD','CHN','IND','IDN','MMR','PHL','THA','VNM')
# crop_yld_shp <- function(yld, shp, out_name){
#   print(paste("Writing", out_name))
#   yld_ras <- raster(yld_name)
#   print("Cropping")
#   yld_crop <- raster::crop(yld_ras, shp)
#   print("Masking")
#   yld_crop <- raster::mask(yld_crop,shp)
#   yld_crop[yld_crop==0] <- NA
#   writeRaster(yld_crop, out_name, overwrite=TRUE)
# }
# 
# for(ctry in ctry_list){
#   ctry_bound <- c(list.files(paste0("Data/Country_mask"),
#                              recursive = F,
#                              full.names = T,
#                              pattern = paste0("(",as.character(ctry),").*\\.shp$")))
#   # print(ctry_bound)
#   ctry_shp <- readOGR(ctry_bound)
#   for(rcp in rcp_list){
#     for(tp in tp_list){
#       if(rcp == "baseline"){
#         yld_name <- paste0('Data/Model results/Crop_yield/',
#                            rcp,'/rice_only/Avg_yld_',tp,'_riceOnly_rs.tif')
#         out_name <- paste0('Data/Model results/Crop_yield/',
#                            rcp,'/countrywise/Avg_yld_',tp,'_',ctry,'_riceOnly_rs.tif')
#         crop_yld_shp(yld_name,ctry_shp,out_name)
#       }else{
#         for(yr in yr_list){
#           rcp_name = paste0(substr(rcp,5,5),substr(rcp,6,6))
#           yld_name <- paste0('Data/Model results/Crop_yield/',
#                              rcp,'/rice_only/Avg_yld_',tp,'_',
#                              yr,'_',rcp_name,'_riceOnly_rs.tif')
#           out_name <- paste0('Data/Model results/Crop_yield/',
#                              rcp,'/countrywise/Avg_yld_',tp,'_',
#                              yr,'_',rcp_name,'_',ctry,'_riceOnly_rs.tif')
#           crop_yld_shp(yld_name,ctry_shp,out_name)
#         }
#       }
#     }
#   }
# }

# ctry_yld_df <- data.frame(matrix(nrow=0,
#                                  ncol=5))
# names(ctry_yld_df) <- c('ctry','rcp','tp','yr','sum_yld')
# count <- 1
# for(ctry in ctry_list){
#   for(tp in tp_list){
#     for(rcp in rcp_list){
#       rcp_name = paste0(substr(rcp,5,5),substr(rcp,6,6))
#       if(rcp == "baseline"){
#         yld_name <- paste0('Data/Model results/Crop_yield/',
#                            rcp,'/countrywise/Avg_yld_',tp,'_',ctry,'_riceOnly_rs.tif')
#         print(paste(ctry,rcp,tp))
#         yld_ras <- raster(yld_name)
#         ctry_yld_df[count,c(1:4)] <- c(ctry,rcp,tp,'base')
#         ctry_yld_df[count,5] <- cellStats(yld_ras, sum)
#         count <- count+1
#       }else{
#         for(yr in yr_list){
#           yld_name <- paste0('Data/Model results/Crop_yield/',
#                              rcp,'/countrywise/Avg_yld_',tp,'_',
#                              yr,'_',rcp_name,'_',ctry,'_riceOnly_rs.tif')
#           print(paste(ctry,rcp,tp,yr))
#           yld_ras <- raster(yld_name)
#           ctry_yld_df[count,c(1:4)] <- c(ctry,rcp,tp,yr)
#           ctry_yld_df[count,5] <- cellStats(yld_ras, sum)
#           count <- count+1
#         }
#       }
#     }
#   }
# }
# write.csv(ctry_yld_df,'Data/Model results/Crop_yield/countrywise_yld.csv')

#calculate differences in yield
ctry_yld_df <- read.csv('Data/Model results/Crop_yield/countrywise_yld_v2.csv')
ctry_yld_df[,1] <- NULL
ctry_yld_df$yld <- ctry_yld_df$sum_yld/ctry_yld_df$num_pix
yld_diff_df <- data.frame(matrix(nrow=0,
                                 ncol=4))
names(yld_diff_df) <- c('ctry','rcp','yr','diff_yld')
count <- 1
for(ctry in ctry_list){
  for(rcp in rcp_list){
    if(rcp == "baseline"){
      diff <- ctry_yld_df[((ctry_yld_df$ctry==ctry)&
                             (ctry_yld_df$rcp==rcp)&
                             (ctry_yld_df$yr=="base")&
                             (ctry_yld_df$tp=="C4")),c('yld')]-
              ctry_yld_df[((ctry_yld_df$ctry==ctry)&
                             (ctry_yld_df$rcp==rcp)&
                             (ctry_yld_df$yr=="base")&
                             (ctry_yld_df$tp=="C3")),c('yld')]
      print(paste(ctry,rcp,diff))
      
      yld_diff_df[count,] <- c(ctry,rcp,'base',diff)
      count <- count+1
    }else{
      for(yr in yr_list){
        diff <- ctry_yld_df[((ctry_yld_df$ctry==ctry)&
                               (ctry_yld_df$rcp==rcp)&
                               (ctry_yld_df$yr==yr)&
                               (ctry_yld_df$tp=="C4")),c('sum_yld')]-
          ctry_yld_df[((ctry_yld_df$ctry==ctry)&
                         (ctry_yld_df$rcp==rcp)&
                         (ctry_yld_df$yr==yr)&
                         (ctry_yld_df$tp=="C3")),c('sum_yld')]
        print(paste(ctry,rcp,yr,diff))
        
        yld_diff_df[count,] <- c(ctry,rcp,yr,diff)
        count <- count+1
      }
    }
  }
}

# diff_plot_ctry <- diff[diff$ctry %in% c('CHN','IND','BGD'),] %>%
#   ggplot(aes(x=Latitude,y=rollmean(diff, 50,na.pad=TRUE),color=ctry))+
#   geom_line(size=1)+
#   scale_y_continuous(limits=c(-1,3), breaks=(seq(from=-1,to=3,by=1)))+
#   scale_x_continuous(limits=c(-20,60),breaks=(seq(from=-20,to=60,by=10)))+
#   geom_hline(yintercept=0, color="black", linetype="dashed", size=1, alpha=0.8)+
#   geom_vline(xintercept=30, color="red", size=0.4, alpha=0.7)+
#   ylab("Yield")+
#   xlab("Latitude")+
#   # ylim(-3,3)+
#   coord_flip()+
#   ggtitle("Difference in yield across latitudes")
# 
# diff_plot_ctry


# ctry_yld_agg <- aggregate(diff[,c(2,3,5,6)], by=list(diff$ctry), sum)
# ctry_yld_diff_agg <- aggregate(diff[,c(7)], by=list(diff$ctry), mean)
# names(ctry_yld_diff_agg)[2] <- 'avg_diff'
# 
# ctry_yld_agg_join <- left_join(ctry_yld_agg,ctry_yld_diff_agg, by="Group.1")
# names(ctry_yld_agg_join)[1] <- 'ctry'

yld_diff_df$ctry_f <- ifelse(yld_diff_df$ctry=="BGD","Bangladesh",
                             ifelse(yld_diff_df$ctry=="CHN","China",
                              ifelse(yld_diff_df$ctry=="IND","India",
                               ifelse(yld_diff_df$ctry=="IDN","Indonesia",
                                ifelse(yld_diff_df$ctry=="MMR","Myanmar",
                                 ifelse(yld_diff_df$ctry=="PHL","Philippines",
                                  ifelse(yld_diff_df$ctry=="THA","Thailand",
                                   ifelse(yld_diff_df$ctry=="VNM","Viet Nam","N/A"))))))))
yld_diff_df[yld_diff_df$rcp=="baseline",c('rcp')] <- "Baseline"
yld_diff_df$diff_yld <- as.numeric(yld_diff_df$diff_yld)
yld_diff_df$diff_yld_M <- as.numeric(yld_diff_df$diff_yld)/1000000

yld_diff_df <- read.csv('Data/Model results/Crop_yield/ctry_yield_diff.csv')

fao_prod <- read.csv("Data/FAOSTAT_data_rice.csv")
fao_prod <- fao_prod[,c("Area","Value")]
names(fao_prod)[1] <- "ctry_f"
fao_prod$valueM <- fao_prod$Value/1000000

yld_fao_join <- left_join(yld_diff_df,fao_prod,by="ctry_f")

# yld_fao_join <- yld_fao_join[order(yld_fao_join$valueM, decreasing=TRUE),]
# x=reorder(feat, imp),y=imp
# ctry_yld_plot <- yld_fao_join %>%
#   ggplot(aes(x=reorder(ctry_f,valueM), y=avg_diff,fill=ctry_f))+ #fill should be by RCP scenario
#   geom_bar(stat='identity',width=0.7,show.legend = FALSE)+
#   xlab("Country")+
#   ylab("Average of yield difference")+
#   geom_hline(yintercept=0, color="black", size=1, alpha=0.8)+
#   coord_flip()+
#   ggtitle("Yield difference (C4-C3)")

create_yld_diff_plot <- function(yld_df, fao_df, year, mn, mx, by){
  ctry_yld_plot <- yld_df[yld_df$yr %in% c('base',year),] %>%
    ggplot(aes(x=reorder(ctry_f,valueM), y=diff_yld,fill=rcp))+ #fill should be by RCP scenario
    geom_bar(position="dodge",stat='identity',width=0.7)+
    scale_fill_brewer(palette = "RdYlGn", direction = -1)+
    scale_y_continuous(breaks=seq(from=mn, to = mx, by=by), limits=c(mn,mx))+
    # ylim(-8.5,1)+
    xlab("Country")+
    ylab("Yield difference (in t/ha)")+
    geom_hline(yintercept=0, color="black", size=1, alpha=0.8)+
    coord_flip()+
    labs(fill='SSP')+
    ggtitle(paste0("Yield difference (C4-C3) in ",year))+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank())
  ctry_prod_plot <- fao_df %>%
    ggplot(aes(x=reorder(ctry_f,valueM), y=valueM))+ #fill should be by RCP scenario
    geom_bar(stat='identity',width=0.7, fill="orange2")+
    xlab("Country")+
    ylab("Production (in million T)")+
    coord_flip()+
    ggtitle("Rice production by country (FAO)")
  return(plot_grid(ctry_prod_plot,ctry_yld_plot,align="h"))
}

range(yld_diff_df$diff_yld)
diff_2050 <- create_yld_diff_plot(yld_fao_join, fao_prod, '2050', -2.5, 1.5, 1)
diff_2050
ggsave(paste0('./Data/Model results/Crop_yield/plots/ctry_yldDiffxFAOprod_baseline_','2050','_v4.png'),
       width=10,height=6,units="in")
diff_2099 <- create_yld_diff_plot(yld_fao_join, fao_prod, '2099', -2.5, 1.5, 1)
diff_2099
ggsave(paste0('./Data/Model results/Crop_yield/plots/ctry_yldDiffxFAOprod_baseline_','2099','_v4.png'),
       width=10,height=6,units="in")

write.csv(ctry_yld_df, './Data/Model results/Crop_yield/ctry_yield.csv')
write.csv(yld_diff_df, './Data/Model results/Crop_yield/ctry_yield_diff.csv')

# ctry_yld_plot



# ctry_yld_plot <- ctry_yld_agg_join %>%
#   ggplot(aes(x=ctry_f, y=avg_diff,fill=ctry_f))+ #fill should be by RCP scenario
#   geom_bar(stat='identity',width=0.7,show.legend = FALSE)+
#   xlab("Country")+
#   ylab("Average of yield difference")+
#   coord_flip()+
#   ggtitle("Yield difference (C4-C3)")
# 
# ctry_prod_plot <- fao_prod %>%
#   ggplot(aes(x=Area, y=valueM))+ #fill should be by RCP scenario
#   geom_bar(stat='identity',width=0.7, fill="orange2")+
#   xlab("Country")+
#   ylab("Production (in million T)")+
#   coord_flip()+
#   ggtitle("Rice production by country (FAO)")+
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank())

# ctry_yld_plot
# ctry_prod_plot
