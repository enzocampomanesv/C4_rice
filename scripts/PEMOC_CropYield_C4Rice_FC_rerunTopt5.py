# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
from osgeo import gdal
import sys

c4 = int(sys.argv[1])
projected = int(sys.argv[2])
year_list=['2050','2099']
rcp_list=['SSP126','SSP245','SSP460','SSP585']

# Function to calculate crop yield using PEMOC
def calculate_GPP(avg_ndvi,max_ndvi,temp,vpd,swn,co2,sos, eos, c4):
        ax=1.45324
        bx=0.18079
        cx=0.22933
        dx=1.00
        ex=0.23264
        hi=0.55
        mc=0.14
        r=0.9
        f=1/0.45
        
        gpp_fact=0.63 # convert GPPmax to GPP
        resp_fact=0.5
        Topt=25+5 #Topt +5 = 30
        #avg_crop_grow_period=120 #use abs(avg_eos-avg_sos)*15 days
        avg_crop_grow_period = (abs(eos-sos))*15
        factor=(f*r*hi)/(1-mc)

        fpar=ax*avg_ndvi-bx
        fpar_max=ax*max_ndvi-bx
        fa=fpar/fpar_max
        ft=(1.1814/((1+np.exp(cx*(Topt-10-temp)))*(1+np.exp(dx*(temp-10-Topt)))))
        fm=(1-ex*np.log(vpd))
        par=2.3*12.011*swn

        if(c4==1):
            efficiency = 0.06 #c4_efficiency = 0.06
        else: #c3 GPP calc
            mxquant = 0.08
            CO2_conc = 0.9*co2 #380 should be CO2 raster
            tau=2600*(0.57**((temp-25)/10))
            gamma = 209000/(2*tau)
            efficiency = mxquant*((CO2_conc-gamma)/(CO2_conc+2*gamma))
        
        output= (factor*avg_crop_grow_period*gpp_fact*resp_fact*efficiency*fpar*fa*ft*fm*par)/100
        return output

wd_path="Data/"
avg_NDVI_filename=wd_path+"Average and SD parameter values/Avg_NDVI.tiff"
max_NDVI_filename=wd_path+"Average and SD parameter values/Average_max_NDVI.tiff"

for year in year_list:
    for rcp in rcp_list:
        rcp_name = rcp[-2]+rcp[-1]
        rcp_num = rcp[-2]+"."+rcp[-1]
        print(year+" "+rcp_name+" "+rcp_num)
        if projected == 1:
            avg_temp_filename=wd_path+"Temp/"+year+"/"+rcp_num+"/avg_temp_" + year+"_"+rcp_name+".tiff"
            max_temp_filename=wd_path+"Temp/"+year+"/"+rcp_num+"/max_temp_" + year+"_"+rcp_name+".tiff"
            min_temp_filename=wd_path+"Temp/"+year+"/"+rcp_num+"/min_temp_" + year+"_"+rcp_name+".tiff"
            avg_co2_filename=wd_path+"Carbon/projections/"+rcp+"/annual/"+year+"/"+"CO2_proj_"+rcp+"_"+year+".tif"
        else:
            avg_temp_filename=wd_path+"Average and SD parameter values/Avg_temp.tiff"
            max_temp_filename=wd_path+"Average and SD parameter values/Average_max_temp.tiff"
            min_temp_filename=wd_path+"Average and SD parameter values/Average_min_temp.tiff"
            avg_co2_filename=wd_path+"Average and SD parameter values/Avg_CO2.tiff"

        avg_swn_filename=wd_path+"Average and SD parameter values/Avg_swn.tiff"

        avg_vap_filename=wd_path+"Average and SD parameter values/Avg_vap.tiff"

        avg_eos_filename=wd_path+"Average and SD parameter values/Average_eos.tiff"
        avg_sos_filename=wd_path+"Average and SD parameter values/Average_sos.tiff"

        ref_ds = gdal.Open(wd_path+"Asia mask/asia_mask.tif")
        band_ref= ref_ds.GetRasterBand(1)
        rows = ref_ds.RasterYSize
        cols = ref_ds.RasterXSize
        ref_projection=ref_ds.GetProjection()
        ref_geotrans=ref_ds.GetGeoTransform()
        xOrigin = ref_geotrans[0]
        yOrigin = ref_geotrans[3]
        pixelWidth = ref_geotrans[1]

        avg_NDVI_ds=gdal.Open(avg_NDVI_filename)
        avg_NDVI= avg_NDVI_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

        max_NDVI_ds=gdal.Open(max_NDVI_filename)
        max_NDVI=max_NDVI_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

        avg_temp_ds=gdal.Open(avg_temp_filename)
        avg_T=avg_temp_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

        max_temp_ds=gdal.Open(max_temp_filename)
        max_T=max_temp_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

        min_temp_ds=gdal.Open(min_temp_filename)
        min_T=min_temp_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

        avg_vap_ds=gdal.Open(avg_vap_filename)
        avg_vap=avg_vap_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

        avg_swn_ds=gdal.Open(avg_swn_filename)
        avg_swn=avg_swn_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

        avg_co2_ds=gdal.Open(avg_co2_filename)
        avg_co2=avg_co2_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

        avg_eos_ds=gdal.Open(avg_eos_filename)
        avg_eos=avg_eos_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)
        avg_sos_ds=gdal.Open(avg_sos_filename)
        avg_sos=avg_sos_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

        svp_max=0.6108*np.exp((17.27*max_T)/(max_T+237.3)) #Buck's equation
        svp_min=0.6108*np.exp((17.27*min_T)/(min_T+237.3))
        svp=(svp_max+svp_min)/2

        vpd=svp-avg_vap

        avg_yield_array= np.empty([rows,cols])

        # values for different parameters: Temperature =1, Vapour pressure= 2, radiation= 3, NDVI =4

        for row in range(rows):
            for col in range(cols):
                if np.isnan(avg_NDVI[row][col])==True:
                    avg_yield_array[row][col]=np.nan
                    continue
                elif np.isnan(max_NDVI[row][col])==True:
                    avg_yield_array[row][col]=np.nan
                    continue        
                elif np.isnan(avg_T[row][col])==True:
                    avg_yield_array[row][col]=np.nan
                    continue  
                
                elif np.isnan(avg_vap[row][col])==True:
                    avg_yield_array[row][col]=np.nan
                    continue 
                
                elif np.isnan(avg_swn[row][col])==True:
                    avg_yield_array[row][col]=np.nan
                    continue 
                
                elif np.isnan(vpd[row][col])==True:
                    avg_yield_array[row][col]=np.nan
                    continue         
        
                elif row!=0 and avg_NDVI[row][col]== avg_NDVI[row-1][col] and avg_T[row][col]==avg_T[row-1][col] and avg_vap[row][col]==avg_vap[row-1][col] and avg_swn[row][col]==avg_swn[row-1][col] and avg_co2[row][col]==avg_co2[row-1][col]:
                    avg_yield_array[row][col]=avg_yield_array[row-1][col]

                elif col!=0 and avg_NDVI[row][col]== avg_NDVI[row][col-1] and avg_T[row][col]==avg_T[row][col-1] and avg_vap[row][col]==avg_vap[row][col-1] and avg_swn[row][col]==avg_swn[row][col-1] and avg_co2[row][col]==avg_co2[row][col-1]:
                    avg_yield_array[row][col]=avg_yield_array[row][col-1]

                else:
                    crop_yield = calculate_GPP(avg_NDVI[row][col],max_NDVI[row][col],avg_T[row][col],vpd[row][col],avg_swn[row][col],avg_co2[row][col],avg_sos[row][col],avg_eos[row][col], c4)
                    avg_yield_array[row][col]=crop_yield

                    print(f'yield calculated for row {row} and column {col} with {crop_yield}', end='\r' )
        
        if c4==1:
            tp = "C4"
        else:
            tp = "C3"

        if projected == 1:
            img_name=wd_path+'Model results/Crop_yield/'+rcp+'/Avg_yld_'+tp+'_'+year+'_'+rcp_name+'_rice_Topt5.tif'
        else:
            img_name=wd_path+'Model results/Crop_yield/baseline/Avg_yld_'+tp+'_rice_Topt5.tif'
        
        gdal.AllRegister()
        driver = gdal.GetDriverByName('GTiff')

        #creating a tiff for the yield
        yield_ds = driver.Create(img_name, cols, rows,1 , gdal.GDT_Float32)
        yield_band=yield_ds.GetRasterBand(1)
        yield_band.WriteArray(avg_yield_array)
        yield_band.SetDescription('Yield')

        yield_ds.SetProjection(ref_projection)
        yield_ds.SetGeoTransform(ref_geotrans)
        yield_ds=None

print('Process finished')