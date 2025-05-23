"""
@author: Deepak
@modified by: Enzo
"""
import sys
import numpy as np
from osgeo import gdal

c4_c3 = int(sys.argv[1])

# Function to perform sensitivity analysis. Sample array contains 10000 random values will be created for each parameter
def find_sensitivity(avg_ndvi,max_ndvi,avg_temp,sd_temp,avg_vap,sd_vap,avg_swn,sd_swn, avg_co2, sd_co2, avg_eos, avg_sos, c4):
    def calculate_GPP(avg_ndvi, max_ndvi, temp, vpd, swn, co2, eos, sos, c4):
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
        Topt=25+5
        #avg_crop_grow_period=120 #use abs(avg_eos-avg_sos)*15 days
        avg_crop_grow_period=(abs(eos-sos))*15
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
        
        # print(np.mean(efficiency))
        output= (factor*avg_crop_grow_period*gpp_fact*resp_fact*efficiency*fpar*fa*ft*fm*par)/100
        return output
    sample_temp_A=np.random.uniform(low=avg_temp-sd_temp,high=avg_temp+sd_temp,size=(10000,))
    sample_temp_B=np.random.uniform(low=avg_temp-sd_temp,high=avg_temp+sd_temp,size=(10000,))
    
    sample_vap_A=np.random.uniform(low=avg_vap-sd_vap,high=avg_vap+sd_vap,size=(10000,))
    sample_vap_B=np.random.uniform(low=avg_vap-sd_vap,high=avg_vap+sd_vap,size=(10000,))
    
    sample_swn_A=np.random.uniform(low=avg_swn-sd_swn,high=avg_swn+sd_swn,size=(10000,))
    sample_swn_B=np.random.uniform(low=avg_swn-sd_swn,high=avg_swn+sd_swn,size=(10000,))

    sample_co2_A=np.random.uniform(low=avg_co2-sd_co2,high=avg_co2+sd_co2,size=(10000,))
    sample_co2_B=np.random.uniform(low=avg_co2-sd_co2,high=avg_co2+sd_co2,size=(10000,))
    
    #Buck's equation-  calculation of vpd
    def bucks(Ta, Vap):
        max_T=np.amax(Ta)
        min_T=np.amin(Ta)
        svp_max=0.6108*np.exp((17.27*max_T)/(max_T+237.3)) 
        svp_min=0.6108*np.exp((17.27*min_T)/(min_T+237.3))
        svp=(svp_max+svp_min)/2
        return(svp-Vap)

    sample_vpd_A = bucks(sample_temp_A,sample_vap_A)
    sample_vpd_B = bucks(sample_temp_B,sample_vap_B)
    sample_vpd_A[sample_vpd_A<=0] = 0.01
    sample_vpd_B[sample_vpd_B<=0] = 0.01
    
    # Set A
    output_A = calculate_GPP(avg_ndvi, max_ndvi, sample_temp_A, sample_vpd_A, sample_swn_A, sample_co2_A, avg_eos, avg_sos, c4)
    # A (temp and VPD bec VPD affects temp), B(SWN, CO2)
    output_B_except_temp = calculate_GPP(avg_ndvi, max_ndvi, sample_temp_A, sample_vpd_A, sample_swn_B,sample_co2_B, avg_eos, avg_sos, c4)
    # A (VPD only), B (Temp, CO2, SWN)
    output_B_except_vpd= calculate_GPP(avg_ndvi, max_ndvi, sample_temp_B, sample_vpd_A, sample_swn_B,sample_co2_B, avg_eos, avg_sos, c4)
    # A (SWN only), B (Temp, VPD, CO2)
    output_B_except_swn= calculate_GPP(avg_ndvi, max_ndvi, sample_temp_B, sample_vpd_B, sample_swn_A,sample_co2_B, avg_eos, avg_sos, c4)
    # A (CO2), B (Temp, VPD, SWN)
    output_B_except_co2= calculate_GPP(avg_ndvi, max_ndvi, sample_temp_B, sample_vpd_B, sample_swn_B, sample_co2_A, avg_eos, avg_sos, c4)
    
    #variance of output from parameter set A
    output_mean_A=np.mean(output_A)
    D_A=np.var(output_A)
    if D_A == 0:
        s_temp = 63.75
        s_vpd = 63.75
        s_swn = 63.75
        s_co2 = 63.75
        max_ind = -1
        co2_v_temp = 0
        total = s_temp+s_vpd+s_swn+s_co2
    else:
        #variance with Temp
        temp_abs=np.absolute(output_A*output_B_except_temp-output_mean_A**2)
        D_temp= np.mean(temp_abs)
        #variance with VPD
        vpd_abs=np.absolute(output_A*output_B_except_vpd-output_mean_A**2)
        D_vpd= np.mean(vpd_abs)
        #variance with SWN
        swn_abs=np.abs(output_A*output_B_except_swn-output_mean_A**2)
        D_swn= np.mean(swn_abs)
        #variance with CO2
        co2_abs=np.abs(output_A*output_B_except_co2-output_mean_A**2)
        D_co2= np.mean(co2_abs)

        #Sensitivity index of Temp
        s_temp= D_temp/D_A
        #Sensitivity index of VPD
        s_vpd= D_vpd/D_A 
        #Sensitivity index of SWN
        s_swn= D_swn/D_A
        #Sensitivity index of CO2
        s_co2= D_co2/D_A

        total=s_temp+s_vpd+s_swn+s_co2
        # print(s_temp,s_vpd,s_swn,s_co2)
        s_temp = s_temp*255
        s_vpd = s_vpd*255
        s_swn = s_swn*255
        s_co2 = s_co2*255
        
        #Get variable with highest sensitivity
        res = [(s_temp)/total,(s_vpd)/total,(s_swn)/total,(s_co2)/total]
        max_res = [(s_temp)/total,(s_vpd)/total,(s_swn)/total,(s_co2)/total]
        max_sens = max(res)
        max_ind = res.index(max_sens)+1
        co2_v_temp = 0
        if s_co2 > s_temp:
            co2_v_temp = 1
        else:
            co2_v_temp = -1
        # max_res.remove(max_sens)
        # max_sens2 = max(max_res)
        # max2_ind = res.index(max_sens2)+1

    return [s_temp/total,s_vpd/total,s_swn/total,s_co2/total,max_ind,co2_v_temp]


wd_path="/home/jovyan/private/C4_Rice/Data/"
avg_NDVI_filename=wd_path+"Average and SD parameter values/Avg_NDVI.tiff"
max_NDVI_filename=wd_path+"Average and SD parameter values/Average_max_NDVI.tiff"
sd_NDVI_filename=wd_path+"Average and SD parameter values/SD_NDVI.tiff"


avg_temp_filename=wd_path+"Average and SD parameter values/Avg_temp.tiff"
max_temp_filename=wd_path+"Average and SD parameter values/Average_max_temp.tiff"
min_temp_filename=wd_path+"Average and SD parameter values/Average_min_temp.tiff"
sd_temp_filename=wd_path+"Average and SD parameter values/SD_temp.tiff"

avg_swn_filename=wd_path+"Average and SD parameter values/Avg_swn.tiff"
sd_swn_filename=wd_path+"Average and SD parameter values/SD_swn.tiff"

avg_vap_filename=wd_path+"Average and SD parameter values/Avg_vap.tiff"
sd_vap_filename=wd_path+"Average and SD parameter values/SD_vap.tiff"

avg_co2_filename=wd_path+"Average and SD parameter values/Avg_CO2.tiff"
sd_co2_filename=wd_path+"Average and SD parameter values/SD_CO2.tiff"

avg_eos_filename=wd_path+"Average and SD parameter values/Average_eos.tiff"
avg_sos_filename=wd_path+"Average and SD parameter values/Average_sos.tiff"

mask_tif = wd_path+"Asia_mask/asia_mask.tif"
print(mask_tif)
ref_ds = gdal.Open(mask_tif) # The mask image containg pixels only in South and East Asia
band_ref= ref_ds.GetRasterBand(1)
rows = ref_ds.RasterYSize
cols = ref_ds.RasterXSize
ref_projection=ref_ds.GetProjection()
ref_geotrans=ref_ds.GetGeoTransform()
xOrigin = ref_geotrans[0]
yOrigin = ref_geotrans[3]
xmax = xOrigin + ref_geotrans[1] * ref_ds.RasterXSize
ymin = yOrigin + ref_geotrans[5] * ref_ds.RasterYSize
asia_mask = band_ref.ReadAsArray(0,0, cols, rows)
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

sd_temp_ds=gdal.Open(sd_temp_filename)
sd_T=sd_temp_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

avg_vap_ds=gdal.Open(avg_vap_filename)
avg_vap=avg_vap_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

sd_vap_ds=gdal.Open(sd_vap_filename)
sd_vap=sd_vap_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

avg_swn_ds=gdal.Open(avg_swn_filename)
avg_swn=avg_swn_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

sd_swn_ds=gdal.Open(sd_swn_filename)
sd_swn=sd_swn_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

avg_co2_ds=gdal.Open(avg_co2_filename)
avg_co2=avg_co2_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

sd_co2_ds=gdal.Open(sd_co2_filename)
sd_co2=sd_co2_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

def normalize(data):
    data_min = np.nanmin(data)
    data_max = np.nanmax(data)
    normalized_data = (data - data_min) / (data_max - data_min)
    return normalized_data
# Normalize each sample
avg_T = normalize(avg_T.copy())
avg_vap = normalize(avg_vap.copy())
avg_swn = normalize(avg_swn.copy())
avg_co2 = normalize(avg_co2.copy())
sd_T = normalize(sd_T.copy())
sd_vap = normalize(sd_vap.copy())
sd_swn = normalize(sd_swn.copy())
sd_co2 = normalize(sd_co2.copy())

avg_eos_ds=gdal.Open(avg_eos_filename)
avg_eos=avg_eos_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)
avg_sos_ds=gdal.Open(avg_sos_filename)
avg_sos=avg_sos_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

mask_sense_array = np.empty([rows,cols])
temp_sense_array= np.empty([rows,cols])
vap_sense_array= np.empty([rows,cols])
rad_sense_array= np.empty([rows,cols])
co2_sense_array= np.empty([rows,cols])
max_sense_array= np.empty([rows,cols])
max2_sense_array= np.empty([rows,cols])

for row in range(rows):
    for col in range(cols):
        if np.isnan(avg_NDVI[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue
        elif np.isnan(max_NDVI[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue        
        elif np.isnan(avg_T[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue  
        elif np.isnan(sd_T[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue         
        elif np.isnan(avg_vap[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue 
        elif np.isnan(sd_vap[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue        
        elif np.isnan(avg_swn[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue 
        elif np.isnan(sd_swn[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue 
        elif np.isnan(avg_co2[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue
        elif np.isnan(sd_co2[row][col])==True:
            temp_sense_array[row][col]=np.nan
            vap_sense_array[row][col]=np.nan
            rad_sense_array[row][col]=np.nan
            co2_sense_array[row][col]=np.nan
            max_sense_array[row][col]=np.nan
            max2_sense_array[row][col]=np.nan
            continue        
        elif row!=0 and avg_NDVI[row][col]== avg_NDVI[row-1][col] and avg_T[row][col]==avg_T[row-1][col] and avg_vap[row][col]==avg_vap[row-1][col] and avg_swn[row][col]==avg_swn[row-1][col] and avg_co2[row][col]==avg_co2[row-1][col]:
            temp_sense_array[row][col]=temp_sense_array[row-1][col]
            vap_sense_array[row][col]=vap_sense_array[row-1][col]
            rad_sense_array[row][col]=rad_sense_array[row-1][col]
            co2_sense_array[row][col]=co2_sense_array[row-1][col]
            max_sense_array[row][col]=max_sense_array[row-1][col]
            max2_sense_array[row][col]=max2_sense_array[row-1][col]
        elif col!=0 and avg_NDVI[row][col]== avg_NDVI[row][col-1] and avg_T[row][col]==avg_T[row][col-1] and avg_vap[row][col]==avg_vap[row][col-1] and avg_swn[row][col]==avg_swn[row][col-1] and avg_co2[row][col]==avg_co2[row][col-1]:
            temp_sense_array[row][col]=temp_sense_array[row][col-1]
            vap_sense_array[row][col]=vap_sense_array[row][col-1]
            rad_sense_array[row][col]=rad_sense_array[row][col-1]
            co2_sense_array[row][col]=co2_sense_array[row-1][col]
            max_sense_array[row][col]=max_sense_array[row-1][col]
            max2_sense_array[row][col]=max2_sense_array[row-1][col]
        else: # if all input rasters have non-NA values
            # print(sd_T[row][col],sd_vap[row][col],sd_swn[row][col])
            sensitive_parameter = find_sensitivity(avg_NDVI[row][col],max_NDVI[row][col],avg_T[row][col],sd_T[row][col],avg_vap[row][col],sd_vap[row][col],avg_swn[row][col],sd_swn[row][col],avg_co2[row][col], sd_co2[row][col],avg_eos[row][col], avg_sos[row][col],c4_c3)
            temp_sense_array[row][col]=sensitive_parameter[0]*asia_mask[row][col]
            vap_sense_array[row][col]=sensitive_parameter[1]*asia_mask[row][col]
            rad_sense_array[row][col]=sensitive_parameter[2]*asia_mask[row][col]
            co2_sense_array[row][col]=sensitive_parameter[3]*asia_mask[row][col]
            max_sense_array[row][col]=sensitive_parameter[4]*asia_mask[row][col]
            max2_sense_array[row][col]=sensitive_parameter[5]*asia_mask[row][col]
            print(f'Index calculated for row {row} and column {col}: {sensitive_parameter}', end = '\r' )

if c4_c3==1:
    img_name=wd_path+'Model results/Sensitivity/C4_sensitivity_withCO2_Topt5.tiff'
    out_name = wd_path+'Model results/Sensitivity/C4_sensitivity_withCO2_Topt5_clip.tiff'
else:
    img_name=wd_path+'Model results/Sensitivity/C3_sensitivity_withCO2_Topt5.tiff'
    out_name = wd_path+'Model results/Sensitivity/C3_sensitivity_withCO2_Topt5_clip.tiff'
gdal.AllRegister()
driver = gdal.GetDriverByName('GTiff')

#creating a tiff with three layers
sensitivity_ds = driver.Create(img_name, cols, rows,6 , gdal.GDT_Int16)
temp_band=sensitivity_ds.GetRasterBand(1)
temp_band.WriteArray(temp_sense_array)
temp_band.SetDescription('Temperature')

vap_band=sensitivity_ds.GetRasterBand(2)
vap_band.WriteArray(vap_sense_array)
vap_band.SetDescription('VPD')

rad_band=sensitivity_ds.GetRasterBand(3)
rad_band.WriteArray(rad_sense_array)
rad_band.SetDescription('Radiation')

co2_band=sensitivity_ds.GetRasterBand(4)
co2_band.WriteArray(co2_sense_array)
co2_band.SetDescription('CO2')

max_band=sensitivity_ds.GetRasterBand(5)
max_band.WriteArray(max_sense_array)
max_band.SetDescription('MaxSens')

max2_band=sensitivity_ds.GetRasterBand(6)
max2_band.WriteArray(max2_sense_array)
max2_band.SetDescription('CO2_v_T')

sensitivity_ds.SetProjection(ref_projection)
sensitivity_ds.SetGeoTransform(ref_geotrans)
sensitivity_ds=None

print(img_name)
print('Process finished')