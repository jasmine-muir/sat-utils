#!/usr/bin/env python
"""
Retrieve all files for a single tile from Google Cloud Storage (GCS), and save as a single zip file 
with a sensible name. 

Does not change the names of the files from those present on the bucket storage. 

"""
from __future__ import print_function, division

import sys
import os
import argparse
import zipfile
import tempfile
import subprocess


def getCmdargs():
    """
    Get command line arguments
    """
    p = argparse.ArgumentParser(description="""
        Retrieve all the files for a single tile/date from Google Cloud Storage, and store in a 
        single zip with a sensible name. 
    """)
    p.add_argument("--tile", help="Name of Sentinel-2 tile e.g. 56JMQ")
    p.add_argument("--date", help=("Acquisition date of desired imagery. "+
        "Specify as yyyy-mm-dd, in UTC time zone, as given by output of s2gcs_catalogtile.py "+
        "e.g. 2017-01-05"))
    p.add_argument("--verbose", default=False, action="store_true",
        help="Print informational messages as it goes")
    cmdargs = p.parse_args()
    return cmdargs


def mainRoutine():
    """
    Main routine
    """
    cmdargs = getCmdargs()

    # The prefix of the Google bucket
    bucketPrefix = "gs://gcp-public-data-sentinel-2/tiles"
    
    tilePrefixStr = makeTilePrefix(bucketPrefix, cmdargs)
    cmdWords = ["gsutil", "ls", "-R", tilePrefixStr]
    if cmdargs.verbose:
        print("Querying tile", tilePrefixStr)
    proc = subprocess.Popen(cmdWords, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = proc.communicate()
    if len(stderr) > 0 or len(stdout) == 0:
        print("Unable to query tile at '{}'".format(tilePrefixStr), file=sys.stderr)
        print("stderr:", stderr)
        sys.exit(1)
        
    filelist = [line.strip() for line in stdout.split('\n') if len(line) > 0]
    filelist = [fn for fn in filelist if not fn.endswith('/:')]
    
    zipfileName = makeZipfileName(cmdargs)
    buildArchive(zipfileName, filelist, tilePrefixStr, cmdargs.verbose)


def makeTilePrefix(bucketPrefix, cmdargs):
    """
    Make the prefix string for the selected tile/date. This is the prefix of the bucket key 
    string which identifies the "subdirectory" containing all files for this tile/date. 
    
    """
    zone = int(cmdargs.tile[:-3])
    latBand = cmdargs.tile[-3].upper()
    tileColRow = cmdargs.tile[-2:].upper()
    
    tileTopdir = "{}/{:02}/{}/{}".format(bucketPrefix, zone, latBand, tileColRow)
    
    cmdWords = ["gsutil", "ls", tileTopdir]
    proc = subprocess.Popen(cmdWords, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = proc.communicate()
    
    if len(stderr) > 0:
        print("Unable to list {}".format(tileTopdir), file=sys.stderr)
        print("stderr:", stderr)
        sys.exit(1)
    
    subdirList = [line.strip() for line in stdout.split('\n') if len(line) > 0]
    subdirList = [subdir[:-1] if subdir.endswith('/') else subdir for subdir in subdirList]
    
    dateStr = cmdargs.date.replace('-', '')
    matchingSubdirList = [subdir for subdir in subdirList if os.path.basename(subdir).split('_')[2][:8] == dateStr]
    if len(matchingSubdirList) == 0:
        print("Unable to find date '{}' for tile '{}'".format(cmdargs.date, cmdargs.tile), file=sys.stderr)
        sys.exit(1)
    elif len(matchingSubdirList) > 1:
        print("Found {} matches for date '{}' and tile '{}'".format(cmdargs.date, cmdargs.tile), file=sys.stderr)
        sys.exit(1)
    
    prefixStr = matchingSubdirList[0]
    return prefixStr


def makeZipfileName(cmdargs):
    """
    Return the name of the zip file to be created. 
    """
    tile = cmdargs.tile.lower()
    date = cmdargs.date.replace('-', '')
    zipfileName = "s2_{}_{}.zip".format(tile, date)
    return zipfileName


def buildArchive(zipfileName, filelist, tilePrefixStr, verbose):
    """
    Build a zipfile archive of the given list of files from GCS. 
    
    """
    zf = zipfile.ZipFile(zipfileName, 'w')
    subdirName = os.path.basename(tilePrefixStr)
    prefixLen = len(tilePrefixStr)
    if verbose:
        print("Downloading", len(filelist), "files")
    for filename in filelist:
        basename = os.path.basename(filename)
        (fd, tmpName) = tempfile.mkstemp(prefix=basename, dir='.')
        os.close(fd)
        
        if verbose:
            print("File", basename)
        cmdWords = ["gsutil", "cp", filename, tmpName]
        proc = subprocess.Popen(cmdWords, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = proc.communicate()
        
        if os.path.exists(tmpName):
            archiveName = "{}{}".format(subdirName, filename[prefixLen:])
            zf.write(tmpName, archiveName)
            os.remove(tmpName)
    zf.close()
    

if __name__ == "__main__":
    mainRoutine()
