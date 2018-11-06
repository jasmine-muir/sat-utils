#!/usr/bin/env python
"""
Unzip the Sentinel 2 zipfiles downloaded from Google Cloud Public bucket.
Uses a search pattern in the current directory to fild all zip files i.e. *.zip or they can be listed on the cmd line
Create indicies: MSAVI, redEdgeNDVI, NDWI - also resample red edge and swir bands from 20m to 10m to match nir and red band resolution
Note: It's probably not the best practise to downsample - but will see how it looks after testing

Written by: Jasmine Muir at Astron
Date: 6/11/2018

Usage: create_indicies.py <input_zip_file_list>
Example: ./create_indicies.py *.zip
"""
from __future__ import print_function
import sys
import os
import argparse
import glob

from osgeo import gdal
from osgeo import ogr
import numpy as np
import numpy.ma as ma


def writeOutImg(inputArray, outfile, n, m, c, TLX, TLY, nulVal, proj, dType):
    # Write the output DEM into an image file with GDAL
    nBands = 1  
    drvr = gdal.GetDriverByName('GTiff')
    ds = drvr.Create(outfile, n, m, nBands, dType, ['COMPRESS=LZW'])
    band = ds.GetRasterBand(1)    
    band.WriteArray(inputArray)
    band.SetNoDataValue(nulVal)
    ds.SetGeoTransform((TLX, c, 0, TLY, 0, -c))
    ds.SetProjection(proj)

    del ds
    

def doMSAVI(a, b, out):
    dsa = gdal.Open(a)
    dsb = gdal.Open(b)
    #Find row and column number
    m = dsa.RasterXSize
    n = dsa.RasterYSize
    #Find input image origin (TLX,TLY) and pixel size
    geotransform = dsa.GetGeoTransform()
    TLX = geotransform[0]
    TLY = geotransform[3]
    c = geotransform[1] 
    proj = dsa.GetProjection()  
    nulVal = dsa.GetRasterBand(1).GetNoDataValue()
    nulOutVal = -999
    
    band1 = np.array(dsa.GetRasterBand(1).ReadAsArray()).astype(np.float32)
    band2 = np.array(dsb.GetRasterBand(1).ReadAsArray()).astype(np.float32)
    print(np.min(band1),np.max(band1),np.min(band2),np.max(band2))    
    band1_norm = band1/np.max(band1)
    band2_norm = band2/np.max(band2)
    print(np.max(band1_norm))
    print(np.max(band2_norm))
    
    ratio = (2*band2_norm+1-np.sqrt((2*band2_norm+1)**2 - 8*(band2_norm-band1_norm)))/2
    ratio = np.where(np.logical_or(np.isnan(ratio),ratio==nulVal),nulOutVal, ratio)
    print(np.mean(ratio))
    writeOutImg(ratio, out, m, n, c, TLX, TLY, nulOutVal, proj, gdal.GDT_Float32)    
    
def doRatioDiff(a, b, out):
    dsa = gdal.Open(a)
    dsb = gdal.Open(b)
    #Find row and column number
    m = dsa.RasterXSize
    n = dsa.RasterYSize
    #Find input image origin (TLX,TLY) and pixel size
    geotransform = dsa.GetGeoTransform()
    TLX = geotransform[0]
    TLY = geotransform[3]
    c = geotransform[1] 
    proj = dsa.GetProjection()  
    nulVal = dsa.GetRasterBand(1).GetNoDataValue()
    nulOutVal = -999

    band1 = np.array(dsa.GetRasterBand(1).ReadAsArray()).astype(np.float32)
    band2 = np.array(dsb.GetRasterBand(1).ReadAsArray()).astype(np.float32)

    ratio = (band1-band2)/(band1+band2)
    ratio = np.where(np.logical_or(np.isnan(ratio),ratio==nulVal),nulOutVal, ratio)
    print(np.mean(ratio))
    writeOutImg(ratio, out, m, n, c, TLX, TLY, nulOutVal, proj, gdal.GDT_Float32)
                    
    
def main():

    #Set up inputs
	zipfileList = sys.argv[1:]
	#10m pixels
	redBand = 4
	NIRBand = 8
	#20m Pixels
	redEdgeBand = 5
	SWIRBand = 11
	
	for zipfile in zipfileList:
	    unzipfolder = zipfile.split('.')[0]
	    if not os.path.exists(unzipfolder):
                cmd = "unzip %s -d %s" %(zipfile,unzipfolder)
                os.system(cmd)

	    for filename in glob.iglob("%s/*/GRANULE/*/IMG_DATA/*.jp2" %unzipfolder, recursive=True):
                try:    
                    band = int(filename.split('/')[-1].split('.')[0][-2:])
                    print(band)
                    if band == redBand:
                        redBandFile = filename
                    elif band == NIRBand:
                        nirBandFile = filename
                    elif band == redEdgeBand:
                        redEdgeBandFile = filename
                    elif band == SWIRBand:
                        SWIRBandFile = filename
                except:
                    print("Band not an integer")	
            #Resample 20m images to
	    redEdgeBandFileResample = redEdgeBandFile.split('/')[-1].split('.')[0] + '_resample.tif'
	    if not os.path.exists(redEdgeBandFileResample):
	        cmd = "gdalwarp -tr 10 10 -r cubic %s %s" %(redEdgeBandFile,redEdgeBandFileResample)
	        os.system(cmd)
	    print(SWIRBandFile)
	    SWIRBandFileResample = SWIRBandFile.split('/')[-1].split('.')[0] + '_resample.tif'
	    print(SWIRBandFileResample)	    
	    if not os.path.exists(SWIRBandFileResample):	   
	        cmd = "gdalwarp -tr 10 10 -r cubic %s %s" %(SWIRBandFile,SWIRBandFileResample)
	        os.system(cmd)
	    msaviFile = "_".join(redEdgeBandFile.split('/')[-1].split('.')[0].split('_')[:-1]) + '_msavi.tif'
	    ndwiFile = "_".join(redEdgeBandFile.split('/')[-1].split('.')[0].split('_')[:-1]) + '_ndwi.tif'
	    rendviFile = "_".join(redEdgeBandFile.split('/')[-1].split('.')[0].split('_')[:-1]) + '_rendvi.tif'
	    
	    if not os.path.exists(msaviFile):
	        print("Creating %s" %(msaviFile))
	        doMSAVI(redBandFile,nirBandFile,msaviFile)
	    if not os.path.exists(ndwiFile):
	        print("Creating %s" %(ndwiFile))
	        doRatioDiff(nirBandFile,SWIRBandFileResample,ndwiFile)
	    if not os.path.exists(rendviFile):		
	        print("Creating %s" %(rendviFile))
	        doRatioDiff(nirBandFile,redEdgeBandFileResample,rendviFile)		

if __name__ == "__main__":
    main()
