#!/usr/bin/env python
"""
Script to calculate the mean value for each band from a generic input image, for each polygon in an input shapefile

Does the same kind of thing as Starspan - but only includes pixels whose centre co-ordinate is within the shapefile

Will later modify to include other stats if necessary i.e. standard deviation, median etc.


Written by: Jasmine Muir at University of New England
Date: /7/2016

Usage: mean_poly_generic.py <input_image> <input_shapefile> <out_csv>
Example: ./mean_poly.py S15_8March2015_MidNSW_ortho_TOA_LS boundaries_subset_gda_remerge.shp outstats.csv
"""
from __future__ import print_function, division
import sys
import os
import math
import random
from osgeo import gdal
from osgeo import ogr
import numpy as np
import numpy.ma as ma


def mksingleshape(shapefile):
    shapebasename = shapefile.split('.')[0]
    #Takes an input shapefile and creates a new shapefile for each feature
    outshapeList = []
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapefile, 0)
    layer = dataSource.GetLayer()
    proj = layer.GetSpatialRef()
    count = 1
    for feature in layer:
        geom = feature.GetGeometryRef()
        name = str(feature.GetField("BlocID")).strip()
        crop = "_".join(str(feature.GetField("Crop")).strip().split(" "))		
        outshape = "%s_out%s_%s.shp" %(shapebasename,name,crop)
        count += 1
        outshapeList.append(outshape)
        if not os.path.exists(outshape):
            # Create the output shapefile
            outDataSource = driver.CreateDataSource(outshape)
            outLayer = outDataSource.CreateLayer(outshape.split('.')[0], proj, geom_type=ogr.wkbPolygon)      
            outLayer.CreateFeature(feature)
            outDataSource.Destroy()
    dataSource.Destroy()
          
    return outshapeList



def writeOutImg(inputArray, outfile, n, m, c, TLX, TLY, nulVal, proj, dType):
    # Write the output DEM into an image file with GDAL
    nBands = 1  
    drvr = gdal.GetDriverByName('HFA')
    ds = drvr.Create(outfile, n, m, nBands, dType, ['COMPRESS=YES'])
    band = ds.GetRasterBand(1)

    band.WriteArray(inputArray)
    ds.SetGeoTransform((TLX, c, 0, TLY, 0, -c))
    ds.SetProjection(proj)
    progress=cuiprogress.CUIProgressBar() 
    calcstats.calcStats(ds, progress, ignore=nulVal)
    del ds

def getStats(infile):
    ds = gdal.Open(infile)
    band = np.array(ds.GetRasterBand(1).ReadAsArray())
    bandMasked = np.ma.array(band, mask = (band==-999))
    bandMaskedFill = np.ma.filled(bandMasked,np.nan)
    mean = np.nanmean(bandMaskedFill.flatten())
    sd = np.nanstd(bandMaskedFill.flatten())

    per5 = np.nanpercentile(bandMaskedFill.flatten(),5)
    per25 = np.nanpercentile(bandMaskedFill.flatten(),25)
    per50 = np.nanpercentile(bandMaskedFill.flatten(),50)
    per75 = np.nanpercentile(bandMaskedFill.flatten(),75)
    per95 = np.nanpercentile(bandMaskedFill.flatten(),95)
    #print(per5)
    stats ="%s,%s,%s,%s,%s,%s,%s" %(per5,per25,per50,per75,per95,mean,sd)
    #print(stats)
    del ds
    return stats    
      
    
def main():

    shapefile = sys.argv[1]    
    outfile = sys.argv[2]
    inImageList = sys.argv[3:]
    #sys.exit()
    
    if os.path.exists(outfile):
        os.remove(outfile)
    outstats = open(outfile,'w')
	#Write the output header    
    outstats.write("Outshape,Image,Date,CropType,per5,per25,per50,per75,per95,mean,sd\n")
     
    
    #create single feature shapefiles from the input shapefile    
    (outshapeList) = mksingleshape(shapefile)
    print(outshapeList)
    
    for image in inImageList:
	    #Check if shapefiles in image extent
        ds = gdal.Open(image)
        geotransform = ds.GetGeoTransform()
        m = ds.RasterXSize
        n = ds.RasterYSize
        TLX = geotransform[0]
        TLY = geotransform[3]
        c = geotransform[1]
        BRX = float(TLX) + (float(c) * float(m))
        BRY = float(TLY) - (float(c) * float(n))
        nBands = ds.RasterCount
        if len(image.split('_'))>3:
            date = image.split('_')[-4][:8]
        else:
            date = image.split('_')[-2][:8]		
        print(date)
        #continue		
    
       
        #use gdal_rasterize to create output image for each single feature shapefile
        outList = []
        for outshape in outshapeList:
            print(image,outshape)
            driver = ogr.GetDriverByName("ESRI Shapefile")
            dataSource = driver.Open(outshape, 0)
            layer = dataSource.GetLayer()
            extent = layer.GetExtent()			
            croptype = outshape.split(".")[0].split("_")[-1]
            print(croptype)			
            xmin = np.floor(float(extent[0]))
            ymin = np.floor(float(extent[2])) - 2*c 
            xmax = np.ceil(float(extent[1])) + 2*c
            ymax = np.ceil(float(extent[3]))
            if xmin > TLX and xmax < BRX and ymin > BRY and ymax < TLY:
                #need to make a fresh copy of the raster to burn into - now subseting the main raster to each polygon extent
                outraster = image.split('.')[0] + '_' + outshape.split('/')[-1].split('.')[0] + '.tif'
                print(outraster)
        
                if os.path.exists(outraster):
                    os.remove(outraster)
                if not os.path.exists(outraster):
                    cmd = "gdal_translate -projwin %s %s %s %s -a_nodata -999 %s %s" %(xmin,ymax,xmax,ymin,image,outraster)
                    os.system(cmd)
        
                    cmd = "gdal_rasterize -b 1 -i -burn -999 -l %s %s %s" %(outshape.split('.')[0], outshape, outraster)
                    os.system(cmd)                    
                    stats = getStats(outraster)
                    outline = "%s,%s,%s,%s,%s\n"%(outshape.split('.')[0], image, date, croptype, stats)
                    outstats.write(outline)
                    print(outline)                      
            else:
                print("Shapefile %s not within imagery extent" %outshape)
        
    outstats.close()
        
if __name__ == "__main__":
    main()
