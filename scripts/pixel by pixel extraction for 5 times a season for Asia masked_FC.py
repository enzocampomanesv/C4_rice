"""
Created on Sun Mar 21 00:41:50 2021

@author: Deepak
"""

import os,sys
from osgeo import gdal
from datetime import date
from calendar import monthrange
import numpy as np
from statistics import mean 



'******************************************************************************'

# This script is used to extract the parameters for the model PEMOC for  year 1982 to 2015

wd = "/data/private/C4_Rice"

##############################################
#Function to find SOS and EOS dates from image file from timesat (The files have been created using the NDVI images)
#there are 24 NDVI image files in a year (2 images per month)

def date_find(d,no,yy):
    img_no=int(d-no*24)
    
    if img_no ==0 and d-no*24>0:
        mo =1
    elif img_no ==1 or img_no ==2 :
        mo =1
    
    elif img_no ==3 or img_no ==4 :
        mo =2
    elif img_no ==5 or img_no ==6 :
        mo =3
    elif img_no ==7 or img_no ==8 :
        mo =4
    elif img_no ==9 or img_no ==10 :
        mo =5
    elif img_no ==11 or img_no ==12 :
        mo =6
    elif img_no ==13 or img_no ==14 :
        mo =7
    elif img_no ==15 or img_no ==16 :
        mo = 8
    elif img_no ==17 or img_no ==18 :
        mo =9
    elif img_no ==19 or img_no ==20 :
        mo =10
    elif img_no ==21 or img_no ==22 :
        mo =11
    else:
        mo =12
    
    
    remain = d-(no*24+img_no)


    no_days= monthrange(yy,mo)
    
    if img_no%2==0:
        plant_day=int(15+(no_days[1]-15)*remain)
    else:
        plant_day=int(15*remain)
    if plant_day==0:
        plant_day +=1    
        
    return(plant_day,mo,yy)   

'****************************************************************************************'
#function to extract CO2 values
def co2_extr(dmy_list,r,c):
    
    #Co2 file path stored in the server
    dir_name=wd+'/Data/Carbon/monthly/resampled/'
    co2_val=[]
    
    for dmy in dmy_list:
        try:
            month = dmy[1]
            year = dmy[2]
            fname='CO2_'+str(year)+'_'+str(month)+'_rs.tif' 
            dataset = gdal.Open(dir_name+fname)   
            band = dataset.GetRasterBand(1)
            data = band.ReadAsArray(0, 0, cols, rows)
            val=data[row][col]
            if val>0 :
                co2_val.append(val)
            dataset=band=data=None
            #print(mean(co2_val))
            #print(str(month)+' '+str(year))
        except:
            continue
    try:            
        return [mean(co2_val), max(co2_val)]
    except:
        return [0,0]

'********************************************************************************'
#The code will run for 34 years starting from 1982 to 2015.

yr_n=1 #yr_n should be given accoring to year if the code wants to run separate for year; 1 for 1982, 2 for 1983  and go on 

# register all of the GDAL drivers
gdal.AllRegister()

# In the following line gdal opens the image for masking. The mask image containg pixels only in South and East Asia
ref_ds = gdal.Open(wd+"/Data/Asia mask/asia_mask.tif")
band_ref= ref_ds.GetRasterBand(1)


rows = ref_ds.RasterYSize
cols = ref_ds.RasterXSize
ref_projection=ref_ds.GetProjection()
ref_geotrans=ref_ds.GetGeoTransform()
xOrigin = ref_geotrans[0]
yOrigin = ref_geotrans[3]
pixelWidth = ref_geotrans[1]
pixelHeight = -ref_geotrans[5]

mask = band_ref.ReadAsArray(0, 0, cols, rows,buf_type=gdal.GDT_Int16)

#Specify the start and end of period under study


for yr in  range(1988,2013):
    #specify the path of
    path =wd+'/Data/Season Extraction results/' + str(yr)
    os.chdir(path)
    sos_name = 'SOS_'+str(yr)+'_season1.img' 
    eos_name = 'EOS_'+str(yr)+'_season1.img'
    
    dataset_sos = gdal.Open(sos_name)   
    band_sos = dataset_sos.GetRasterBand(1)
    data_sos = band_sos.ReadAsArray(0, 0, cols, rows,buf_type=gdal.GDT_Float32)

    dataset_eos = gdal.Open(eos_name)   
    band_eos = dataset_eos.GetRasterBand(1)
    data_eos = band_eos.ReadAsArray(0, 0, cols, rows,buf_type=gdal.GDT_Float32)    
              
#Masking sos values with the mask image
    data_sos=mask*data_sos
#Creating empty arrays for saving parameters

    co2_array=np.empty([rows,cols])

#Loop to start processing

    for row in range(rows):
        for col in range(cols):
            #checking for season 
            try:
            
                if data_sos[row][col]<=yr_n*24 or data_eos[row][col]==0:
                    sos_date= 0
                    eos_date= 0
                elif data_sos[row][col]>= data_eos[row][col]  :
                    sos_date= 0
                    eos_date= 0 
                elif data_eos[row][col]-data_sos[row][col] < 4 : 
                    sos_date= 0
                    eos_date= 0
                else:
                    if data_eos[row][col]>(yr_n+1)*24 and yr!=2015 : # data is only available up to 2015, therefore EOS of 2015 can not be used
                        new_n=yr_n+1
                        new_yr=yr+1
                        
                        # finding 5 dates to extract parameters from a season (start date,quartile date, mid date, three quarter date and end date)
                        sos_date=date_find(data_sos[row][col],yr_n,yr)
                        eos_date=date_find(data_eos[row][col],new_n,new_yr)
                        
                        mid_val=(data_sos[row][col]+data_eos[row][col])/2
                        q_val = mid_val-((data_eos[row][col]-data_sos[row][col])/4)
                        three_q_val= mid_val+((data_eos[row][col]-data_sos[row][col])/4)
                        
                        if mid_val>(yr_n+1)*24:
                            mid_date= date_find(mid_val,new_n,new_yr)
                        else:
                            mid_date= date_find(mid_val,yr_n,yr)
                        if q_val>(yr_n+1)*24:
                            q_date= date_find(q_val,new_n,new_yr)
                        else:
                            q_date= date_find(q_val,yr_n,yr)
                        if three_q_val>(yr_n+1)*24:
                            three_q_date= date_find(three_q_val,new_n,new_yr)
                        else:
                            three_q_date= date_find(three_q_val,yr_n,yr)    
                                             
                    elif data_eos[row][col]>(yr_n+1)*24 and yr==2015 : # data is only available up to 2015, therefore EOS of 2015 can not be used
                        sos_date= 0
                        eos_date= 0
                    else:
                        sos_date=date_find(data_sos[row][col],yr_n,yr)
                        eos_date=date_find(data_eos[row][col],yr_n,yr)
                        
                        mid_val=(data_sos[row][col]+data_eos[row][col])/2
                        mid_date= date_find(mid_val,yr_n,yr)
                        
                        q_val = mid_val-((data_eos[row][col]-data_sos[row][col])/4)
                        q_date= date_find(q_val,yr_n,yr)
                        
                        three_q_val= mid_val+((data_eos[row][col]-data_sos[row][col])/4)
                        three_q_date= date_find(three_q_val,yr_n,yr)
                        
            except KeyboardInterrupt:
                sys.exit()

#from here >> extraction of parameters
# if there is no sos and eos then no extraction
            if sos_date==0 and eos_date==0:
                co2_array[row][col]=0
                #print("No extraction")
                continue
            else:
                if row==0 and col==0 : # for the first pixel
                    dates_extraction=[sos_date,q_date,mid_date,three_q_date,eos_date]
                    avg_co2 = co2_extr(dates_extraction, row, col)
                    co2_array[row][col] = avg_co2[0]
                    print(str(yr)+ ' ' + str(row) +' '+ str(col) + ' ' + str(avg_co2[0]))
                    #print('Parameters extracted for row = '+ str(row)+' column = '+ str(col)+' for year '+str(yr))
                    continue

                elif row==0 and col!=0:
                    if data_sos[row][col]==data_sos[row][col-1]:
                        co2_array[row][col] = avg_co2[row][col-1]
                        continue
                    else:
                        dates_extraction=[sos_date,q_date,mid_date,three_q_date,eos_date]             
                        avg_co2 = co2_extr(dates_extraction, row, col)
                        co2_array[row][col] = avg_co2[0]
                        print(str(yr)+ ' ' + str(row) +' '+ str(col) + ' ' + str(avg_co2[0]))               
                        #print('Parameters extracted for row = '+ str(row)+' column = '+ str(col)+' for year '+str(yr))
                        continue
                elif row!=0 and col==0:
                    if data_sos[row][col]==data_sos[row-1][col]:
                        co2_array[row][col] = avg_co2[row-1][col]
                        continue
                    else:
                        dates_extraction=[sos_date,q_date,mid_date,three_q_date,eos_date]
                        avg_co2 = co2_extr(dates_extraction, row, col)
                        co2_array[row][col] = avg_co2[0]
                        print(str(yr)+ ' ' + str(row) +' '+ str(col) + ' ' + str(avg_co2[0]))
                        #print('Parameters extracted for row = '+ str(row)+' column = '+ str(col)+' for year '+str(yr))
                        continue
                elif row!=0 and col!=0:
                    if data_sos[row][col]==data_sos[row][col-1]:
                        co2_array[row][col] = co2_array[row][col-1]
                        continue
                    elif data_sos[row][col]==data_sos[row-1][col]:
                        co2_array[row][col] = co2_array[row-1][col]
                        continue
                    else:
                        dates_extraction=[sos_date,q_date,mid_date,three_q_date,eos_date]

                        avg_co2 = co2_extr(dates_extraction, row, col)
                        co2_array[row][col] = avg_co2[0]
                        print(str(yr)+ ' ' + str(row) +' '+ str(col) + ' ' + str(avg_co2[0]))      
                        #print('Parameters extracted for row = '+ str(row)+' column = '+ str(col)+' for year '+str(yr))
                        continue
           
#################################################################################   
    #Tiff file creation
    os.chdir(path)
    co2_file_name='Average_CO2_'+str(yr)+'.tiff'
    
    # register all of the GDAL drivers
    gdal.AllRegister()
    driver = gdal.GetDriverByName('GTiff')
    
    co2_ds = driver.Create(co2_file_name, cols, rows, 1, gdal.GDT_Float32)
    co2_ds.SetProjection(ref_projection)
    co2_ds.SetGeoTransform(ref_geotrans)
    co2_band=co2_ds.GetRasterBand(1)
    co2_band.WriteArray(co2_array)
    co2_band.FlushCache()
    co2_band.ComputeStatistics(False)
    
#############################################    
    print('Images created for year '+ str(yr))
    
    co2_ds=None
    yr_n+=1
    
           
print('Process finished')  
     