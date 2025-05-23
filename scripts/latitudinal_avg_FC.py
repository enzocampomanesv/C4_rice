from osgeo import gdal
import numpy as np
import csv
import sys
fields =['Sl No','Latitude','avg_yld','sd_yld','num_pix']

c4 = int(sys.argv[1])
if c4==1:
    tp = "C4"
else:
    tp = "C3"
base = int(sys.argv[2])
if(base==1):
    year_list=['2099']
    rcp_list=['baseline']
else:
    year_list=['2050','2099']
    rcp_list=['SSP126','SSP245','SSP460','SSP585']
#change output filename for c3 and c4
#ctry=['BGD','CHN','IND','IDN','MMR','PHL','THA','VNM']
wd_path="Data/"
for year in year_list:
    for rcp in rcp_list:
        rcp_name = rcp[-2]+rcp[-1]
        rcp_num = rcp[-2]+"."+rcp[-1]
        if(base==1):
            out_filename= wd_path+'Model results/Crop_yield_final/'+rcp+'/latitudinal/'+tp+'_'+rcp+'_riceOnly.csv'
        else:
            out_filename= wd_path+'Model results/Crop_yield_final/'+rcp+'/latitudinal/'+tp+'_'+year+'_'+rcp_name+'_riceOnly.csv'
        with open(out_filename,'w',newline='') as csvfile:
            csvwriter = csv.writer(csvfile) 
            csvwriter.writerow(fields)
        csvfile.close     

        gdal.AllRegister()

        # Change file name for input => c3 or c4
        if(base==1):
            img_ds = gdal.Open(wd_path+'Model results/Crop_yield/'+rcp+'/Avg_yld_'+tp+'_rice_rs_Topt5.tif')
        else:
            img_ds = gdal.Open(wd_path+'Model results/Crop_yield/'+rcp+'/Avg_yld_'+tp+'_'+year+'_'+rcp_name+'_rice_rs_Topt5.tif')
        asia_rice = gdal.Open(wd_path+"Asia mask/asia-rice-extent.tif")
        rows = img_ds.RasterYSize
        cols = img_ds.RasterXSize
        ref_projection=img_ds.GetProjection()
        ref_geotrans=img_ds.GetGeoTransform()
        xOrigin = ref_geotrans[0]
        yOrigin = ref_geotrans[3]
        pixelWidth = ref_geotrans[1]
        pixelHeight = -ref_geotrans[5]

        img_data=img_ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows,buf_type=gdal.GDT_Float32)
        rice_mask = asia_rice.GetRasterBand(1).ReadAsArray(0, 0, cols, rows,buf_type=gdal.GDT_Float32)
        img_data = img_data*rice_mask

        gdal.AllRegister()
        driver = gdal.GetDriverByName('GTiff')

        print("Saving "+rcp+" "+tp)
        if(base==1):
            img_name = wd_path+'Model results/Crop_yield_final/'+rcp+'/rice_only/Avg_yld_'+tp+'_riceOnly_rs.tif'
        else:
            img_name = wd_path+'Model results/Crop_yield_final/'+rcp+'/rice_only/Avg_yld_'+tp+'_'+year+'_'+rcp_name+'_riceOnly_rs.tif'        
        yield_ds = driver.Create(img_name, cols, rows,1 , gdal.GDT_Float32)
        yield_ds.SetProjection(ref_projection)
        yield_ds.SetGeoTransform(ref_geotrans)
        yield_band=yield_ds.GetRasterBand(1)
        yield_band.WriteArray(img_data)
        yield_band.SetDescription('Yield')

        yield_ds=None
        no=1

        for row in range(rows):

            lat= yOrigin-(row*pixelHeight)
            lat_array=img_data[row]
            lat_array_nan=np.where(lat_array<=0, np.nan,lat_array)

            pixel_no=lat_array_nan.size - np.count_nonzero(np.isnan(lat_array_nan))

            if pixel_no==0:
                avg_yield='#N/A'
                sd_yield='#N/A'
            else:
                avg_yield=np.nanmean(lat_array_nan) 
                sd_yield=np.nanstd(lat_array_nan)

            row_list=[no,round(lat,4),avg_yield,sd_yield,pixel_no]
            
            with open(out_filename, 'a',newline='') as csvfile:  
                # creating a csv writer object  
                csvwriter = csv.writer(csvfile) 
                csvwriter.writerow(row_list)
            csvfile.close()
            print(f'{row_list}', end='\r')
            no+=1
        img_ds=None

print('Process finished')
