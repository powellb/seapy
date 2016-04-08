import numpy as np
import re
import os
import seapy

def get_modis(url,outdir,grid):
    #Grid limits
    g = seapy.model.asgrid(grid)
    g_s, g_n = np.min(g.lat_rho), np.max(g.lat_rho)
    g_w, g_e = np.min(g.lon_rho), np.max(g.lon_rho)
    
    #Get file list
    print('Accessing '+url)
    f = seapy.download.Download(url).open_file()
    lines = f.split('\r\n')
    
    #Split file list into .xml and .bz2, ignoring all other files
    xmlfiles = []
    bzfiles = []
    for l in lines:
        if l.endswith('.xml'):
            xmlfiles.append(l.split()[-1])
        elif l.endswith('.bz2'):
            bzfiles.append(l.split()[-1])

    #Search xml files for data in grid region, download data if it exists
    print('Fetching list of files containing data in region...')
    flist=[]
    Nfiles=0
    Nfailedconnect=0
    for x,b in zip(xmlfiles,bzfiles):
        xmlstr = Download(url+x).open_file()
        if xmlstr:
            m_w = float(re.compile('<Westernmost_Longitude>(.*?)</Westernmost_Longitude>').findall(xmlstr)[0])
            m_e = float(re.compile('<Easternmost_Longitude>(.*?)</Easternmost_Longitude>').findall(xmlstr)[0])
            m_n = float(re.compile('<Northernmost_Latitude>(.*?)</Northernmost_Latitude>').findall(xmlstr)[0])
            m_s = float(re.compile('<Southernmost_Latitude>(.*?)</Southernmost_Latitude>').findall(xmlstr)[0])
            if m_w < g_e and m_e > g_w and m_n > g_s and m_s < g_n:
                Nfiles += 1
                flist.append([url+b,outdir+b])
        else:
            Nfailedconnect += 1
    print('%d files found with data in region and %d/%d unsuccessful connections' %(Nfiles,Nfailedconnect,len(xmlfiles))
    print('Downloading files...')
    Nday=0
    Nnight=0
    for f in seapy.progressbar.progress(flist):
      seapy.download.Download(f[0],f[1]).get_file()
      
    print('Extracting files...')
    for f in seapy.progressbar.progress(flist):
      if os.path.isfile(f[1]):
          os.system('bzip2 -d '+f[1])
          ofile=f[1].rstrip('.bz2')
          if '_N-' in ofile:
            os.system('ncks -O -v lat,lon,sea_surface_temperature4,SSES_bias_error4,SSES_standard_deviation_error4,proximity_confidence4 '+\
                ofile +' '+ofile + ' > /dev/null')
            Nnight += 1
          else:
            os.system('ncks -O -v lat,lon,sea_surface_temperature,SSES_bias_error,SSES_standard_deviation_error,proximity_confidence'+\
                ofile+' '+ofile +' > /dev/null')
            Nday += 1
    print('%d day files and %d night files obtained' %(Nday,Nnight))