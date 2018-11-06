#!/usr/bin/env python

from __future__ import print_function
import sys
import os

from osgeo import gdal

import numpy as np
from scipy.ndimage import median_filter

from rios import cuiprogress
from rios import calcstats
from rios import rat


def writeOutImg(inputArray, outfile, n, m, c, TLX, TLY, nulVal, proj, dType):
    # Write the output DEM into an image file with GDAL
    nBands = 1  
    drvr = gdal.GetDriverByName('GTiff')
    ds = drvr.Create(outfile, n, m, nBands, dType, ['COMPRESS=LZW', 'TFW=YES'])
    band = ds.GetRasterBand(1)
    
    band.WriteArray(inputArray)
    ds.SetGeoTransform((TLX, c, 0, TLY, 0, -c))
    ds.SetProjection(proj)
    progress=cuiprogress.CUIProgressBar() 
    calcstats.calcStats(ds, progress, ignore=nulVal)
    del ds


def creatColorTable():
    
    colorArray = np.array([
        [0,240,240,240,0],
        [1,255,0,0,0],
        [2,255,165,0,1],
        [3,255,255,0,1],
        [4,0,255,0,1],
        [5,0,255,255,1],
        [6,0,0,255,1],
        [7,176,48,96,1],
        [8,255,0,255,1]]
    )            

    return colorArray
    
def medianFilter(infile, nullVal, winsize, outfile):  
    
    
    ds = gdal.Open(infile)
    
    imageArray = np.array(ds.GetRasterBand(int(1)).ReadAsArray())
    
    m = ds.RasterXSize
    n = ds.RasterYSize

    #Find input image origin (TLX,TLY) and pixel size
    geotransform = ds.GetGeoTransform()
    TLX = geotransform[0]
    TLY = geotransform[3]
    c = geotransform[1]
    nulVal = 0
    proj = ds.GetProjection() 
    
    medianFilterArray = median_filter(imageArray, size=winsize)
    writeOutImg(medianFilterArray, outfile, m, n, c, TLX, TLY, nulVal, proj, gdal.GDT_Byte)
    #Create color table
    colorArray = creatColorTable()
    rat.setColorTable(outfile,colorArray)
    

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    winsize = sys.argv[3]
    
    medianFilter(infile, 0, int(winsize), outfile)
    
    



        
if __name__ == "__main__":
    main()


