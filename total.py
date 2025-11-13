from vos import Client
import argparse
import time
import numpy as np
import os
from astropy.io import fits
import re
from astropy.coordinates import SkyCoord
import pandas as pd
from astropy.table import hstack,Table
from astropy import units as u
from astropy.wcs import WCS
from scipy.ndimage.filters import median_filter
import shutil
import logging
import subprocess
import matplotlib.pyplot as plt

# Set up logging configuration
log_filename = 'download_log.txt'
logging.basicConfig(filename=log_filename, level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')


# this function checks if the galaxy is in an existing tile and returns the name of the tile:

def in_UNIONS(RA_input, DEC_input, listing, v=True):
    # units of input in degrees
    # listing is the name of the file containing the list of tiles
    
    # i dont remember what this command was for CHECK IF DELETE
    #os.system('export PATH="${HOME}/.local/bin:${PATH}"')

    tiles_file_location = ''
    tiles_file_name = listing
    tiles = np.loadtxt(tiles_file_location + tiles_file_name, dtype=str)
    
    xxx_cen = 0
    yyy_cen = 0

    # separate name in set of xxx and yyy int:
    
    xxx = np.empty(len(tiles), dtype=int)
    yyy = np.empty(len(tiles), dtype=int)
    distance = np.empty(len(tiles), dtype=float)
    DEC_cen = np.empty(len(tiles),dtype=float)
    RA_cen = np.empty(len(tiles),dtype=float)
    for ii in range(len(tiles)):
        xxx[ii], yyy[ii] = np.array(tiles[ii].split('.'), dtype=int)
        
        DEC_cen[ii] = (yyy[ii] / 2) - 90
        RA_cen[ii] = xxx[ii] / (2 * np.cos(DEC_cen[ii] * np.pi / 180))

        distance[ii] = np.sqrt(np.power((DEC_input-DEC_cen[ii]),2)+np.power((RA_input-RA_cen[ii]),2))
      #  print(RA_cen,DEC_cen,RA_input,DEC_input,distance)

    ch = np.where(distance==min(distance))
    
    xxx_cen = xxx[ch]
    yyy_cen = yyy[ch]
    d = distance[ch]
    check_if_tile = False    
    
    if(min(distance)<0.5):
       check_if_tile = True       
       
    return check_if_tile,xxx_cen,yyy_cen

# this function makes use of the previous one to create folder named CFIS.XXX.YYY and download the tiles in whichever filter exists
def download(coords_ra,coords_dec):
    # ra dec units in deg
    coords_ra = float(coords_ra)
    coords_dec = float(coords_dec)

    # first check is in DR5
    # commented lines for g and i, should rethink this
    
    check = in_UNIONS(coords_ra, coords_dec,listing='UNIONS_CFIS_LSB_r_DR5.txt')
    if(check[0]==False):
        logging.info(f'{coords_ra,coords_dec} not in UNIONS')

    if(check[0]):
        xxx_input = int(check[1])
        yyy_input = int(check[2])

        xxx_input = '{:03d}'.format(xxx_input)
        yyy_input = '{:03d}'.format(yyy_input)


        
        if yyy_input==0:
            yyy_input = '000'
        folder = '../CFIS.{}.{}'.format(xxx_input,yyy_input)
        os.makedirs(folder)
        
        file_u = 'CFIS.{}.{}.u.fits'.format(xxx_input,yyy_input)
        file_r = 'CFIS.{}.{}.r.fits'.format(xxx_input,yyy_input)
        file_g = 'calexp-CFIS_{}_{}.fits'.format(xxx_input,yyy_input)

        try:
            # Downloading files
            dw_u = 'vcp vos:cfis/tiles_DR5/{} {}'.format(file_u, folder)
            dw_r = 'vcp vos:cfis/tiles_DR5/{} {}'.format(file_r, folder)            
            dw_g = 'vcp vos:cfis/whigs/stack_images_CFIS_scheme/{} {}'.format(file_g,folder)
            
            os.system(dw_u)
            os.system(dw_r)
            os.system(dw_g)
            
            new_g = '{}/CFIS.{}.{}.g.fits'.format(folder,xxx_input,yyy_input)
            old_g = '{}/{}'.format(folder,file_g)
            r_img = '{}/{}'.format(folder,file_r)
            
            hdg = fits.open(old_g)
            hdr = fits.open(r_img)
            
            hdr[0].data = hdg[1].data
            hdr[0].header['FILTER'] = 'G'
            hdr[0].header['GAIN'] = 1

            hdr.writeto(new_g,overwrite=True)
            
            hdr.close()
            hdg.close()
            
            # Check if both files exist
            if os.path.exists(os.path.join(folder, file_u)) and os.path.exists(os.path.join(folder, file_r)):
                print('Downloaded')
                logging.info(f'Downloaded {file_u,file_r}')
            else:
                raise FileNotFoundError("One or both files not found.")

        except Exception as e:
            print(f'Failed download: {e}')

            # If files don't exist, remove the entire folder
            if not (os.path.exists(os.path.join(folder, file_u)) and os.path.exists(os.path.join(folder, file_r))):
                try:
                    shutil.rmtree(folder)
                    print(f'The folder {folder} and its contents have been deleted.')
                    logging.info(f'The folder {folder} and its contents have been deleted.')
                except Exception as e:
                    print(f'Error deleting folder: {e}')



def detect(folder,ra_sky,dec_sky,gal,dist):
    gal = gal
    d = dist
    rad = np.arcsin((100*np.power(10,3))/(d*np.power(10,6)))*206265
    rad = rad/3600/2

    # detects in the r' filter
    pattern = r"CFIS\.(\d+)\.(\d+)\.([urgiz]).fits"
    files = os.listdir(folder)
    fits_files = [file for file in files if file.endswith('.fits') and file.startswith('CFIS')]
    os.chdir(folder)
    for fits_file in fits_files:
        file_path = os.path.join(fits_file)

        match = re.match(pattern,fits_file)
        xxx = int(match.group(1))
        yyy = int(match.group(2))
        xxx = '{:03d}'.format(xxx)
        yyy = '{:03d}'.format(yyy)
        if 'r' in fits_file:
            rfile = fits_file
            with fits.open(file_path) as hdul:
                flt = str(hdul[0].header['FILTER'])
                gain = str(hdul[0].header['GAIN'])
                seeing = str(hdul[0].header['IQFINAL'])
                zp = str(30)
                pxsc = str(0.187)
                
                # first we run it over the entire tile
                
                w = WCS(hdul[0].header)
                fname = 'gauss_4.0_7x7.conv'
                cat1full = '{}_{}.{}.{}.cat'.format(flt[0],gal,xxx,yyy)

                os.system('cp ../unions_gcs/default* .')
                os.system('cp ../unions_gcs/gauss* .')
                os.system('cp ../unions_gcs/unions* .')
                
                
                

                hlf1 = hdul[0].data[0:5000]
                hlf2 = hdul[0].data[5000:10000]
                
                final_image = '{}_fin_{}.{}.fits'.format(flt[0],xxx,yyy)

                if final_image not in files:
                    filtereda = median_filter(hlf1,size=(50))
                    restaa = hlf1 - filtereda
                    filtered2a = median_filter(restaa,size=(8))
                    rest2a = restaa - filtered2a
                    filteredb = median_filter(hlf2,size=(50))
                    restab = hlf2 - filteredb
                    filtered2b = median_filter(restab,size=(8))
                    rest2b = restab - filtered2b
                    hdul[0].header['PGC'] = gal
                    hdul[0].data = np.concatenate((rest2a,rest2b),axis=0)
                    hdul.writeto(final_image)
                
                pt1 = 'sex '+final_image+' -CATALOG_NAME '+cat1full+' -SEEING_FWHM '+seeing+' -DETECT_THRESH 1.2 -FILTER_NAME '+fname
                pt2 = ' -MAG_ZEROPOINT '+zp+' -GAIN '+gain+' -PIXEL_SCALE '+pxsc+' -PHOT_APERTURES 9 '
                sx = pt1+pt2
                logging.info(sx)
                os.system(sx)


                
    for fits_file in fits_files:
        ufile = fits_file
        file_path = os.path.join(fits_file)
        try:
            with fits.open(file_path) as hdul:
                try:
                    flt = str(hdul[0].header['FILTER'])
                    gain = str(hdul[0].header['GAIN'])
                    print(flt)
                    zp = str(30)
                    pxsc = str(0.187)
                except:
                    flt = str(hdul[1].header['HIERARCH FPA.FILTERID'])
                    gain = str(hdul[1].header['HIERARCH CELL.GAIN'])
                    seeing = str(2)
                    zp = str(30)
                    pxsc = str(0.258)
                    
                # first we run it over the entire tile
                
                w = WCS(hdul[0].header)
                fname = 'gauss_4.0_7x7.conv'
                cat2full = '{}_{}.{}.{}.cat'.format(flt[0],gal,xxx,yyy)


                hlf1 = hdul[0].data[0:5000]
                hlf2 = hdul[0].data[5000:10000]
                
                final_image2 = '{}_fin_{}.{}.fits'.format(flt[0],xxx,yyy)

                if final_image not in files:
                    filtereda = median_filter(hlf1,size=(50))
                    restaa = hlf1 - filtereda
                    filtered2a = median_filter(restaa,size=(8))
                    rest2a = restaa - filtered2a
                    filteredb = median_filter(hlf2,size=(50))
                    restab = hlf2 - filteredb
                    filtered2b = median_filter(restab,size=(8))
                    rest2b = restab - filtered2b
                    hdul[0].header['PGC'] = gal
                    hdul[0].data = np.concatenate((rest2a,rest2b),axis=0)
                    hdul.writeto(final_image2)
                
#                 pt1 = 'sex '+final_image+' -CATALOG_NAME '+cat2full+' -SEEING_FWHM '+seeing+' -DETECT_THRESH 1.2 -FILTER_NAME '+fname
#                 pt2 = ' -MAG_ZEROPOINT '+zp+' -GAIN '+gain+' -PIXEL_SCALE '+pxsc+' -PHOT_APERTURES 9 '
#                 sx = pt1+pt2
#                 logging.info(sx)
#                 os.system(sx)
            
        except Exception as e:
            print('These are not the files you are looking for')


        match = re.match(pattern,fits_file)
        xxx = int(match.group(1))
        yyy = int(match.group(2))
        yyy = '{:03d}'.format(yyy)
        final_image = 'r_fin_{}.{}.fits'.format(flt[0],xxx,yyy)
        r_image = 'R_fin_{}.{}.{}-{}.fits'.format(xxx,yyy,ra_sky,dec_sky)
        

        if os.path.exists(cat2full):
            print('Already run')
        else:
            pt1 = 'sex '+final_image+','+final_image2+' -CATALOG_NAME '+cat2full+' -SEEING_FWHM '+seeing+' -DETECT_THRESH 1.2 -FILTER_NAME '+fname
            pt2 = ' -MAG_ZEROPOINT '+zp+' -GAIN '+gain+' -PIXEL_SCALE '+pxsc+' -PHOT_APERTURES 9 '
            sx = pt1+pt2
            print(sx)
            logging.info(sx)
            os.system(sx)

 
            
    

    # print('----------------- DONE RUNNING SOURCE EXTRACTOR ---------------')
    # print(ucen,gcen,ufull,gfull)
    logging.info('-------------------- DONE EXTRACTING SOURCES -------------')            
    # Define column names based on the letter so we can match
    # print(cat1cen,cat2cen)
#     columns = ['FLUXr', 'FLUXr_e', 'r', 'r_err', 'X', 'Y', 'RA_r', 'DEC_r', 'CLASS_r']
#     table1 = pd.read_csv(cat1cen, delim_whitespace=True, comment='#', names=columns)
    
    # ucen = 'U_cen.{}.{}.{}-{}.cat'.format(xxx,yyy,ra_sky,dec_sky)
#     ufull = 'U_full.{}.{}.{}-{}.cat'.format(xxx,yyy,ra_sky,dec_sky)
#     columns = ['FLUXu', 'FLUXu_e', 'u', 'u_err', 'Xu', 'Yu', 'RA_u', 'DEC_u', 'CLASS_u']
#     # table2 = pd.read_csv(ucen, delim_whitespace=True, comment='#', names=columns)
    
#     # gcen = 'G_cen.{}.{}.{}-{}.cat'.format(xxx,yyy,ra_sky,dec_sky)
#     gfull = 'G_full.{}.{}.{}-{}.cat'.format(xxx,yyy,ra_sky,dec_sky)    
    
#     columns = ['FLUXg', 'FLUXg_e', 'g', 'g_err', 'Xg', 'Yg', 'RA_g', 'DEC_g', 'CLASS_g']
#     # table3 = pd.read_csv(gcen, delim_whitespace=True, comment='#', names=columns)
#     table3['g'] = table3['g'] - 2.7
#     total = pd.concat((table1,table2,table3),axis=1)
#     # print(total)
    
#     # Print the stacked table
#     if total is not None:
#         print(total)
#     cat_nm = 'central-{}-{}.cat'.format(ra_sky,dec_sky)
#     total.to_csv(cat_nm)  
    
    # print(cat1full,cat2full)

    columns = ['FLUXr', 'FLUXr_e', 'r', 'r_err', 'X', 'Y', 'RA_r', 'DEC_r', 'CLASS_r']
    table1 = pd.read_csv(cat1full, delim_whitespace=True, comment='#', names=columns)
    
    
    columns = ['FLUXu', 'FLUXu_e', 'u', 'u_err', 'X', 'Y', 'RA_u', 'DEC_u', 'CLASS_u']
    table2 = pd.read_csv(ufull, delim_whitespace=True, comment='#', names=columns)
    
    columns = ['FLUXg', 'FLUXg_e', 'g', 'g_err', 'Xg', 'Yg', 'RA_g', 'DEC_g', 'CLASS_g']
    table3 = pd.read_csv(gfull, delim_whitespace=True, comment='#', names=columns)
    

    total = pd.concat((table1,table2,table3),axis=1)    
    
    g_img = '{}/G_fin_{}.{}.{}-{}'.format(folder,xxx_input,yyy_input,coords_ra,coords_dec)
    u_img = '{}/U_fin_{}.{}.{}-{}'.format(folder,xxx_input,yyy_input,coords_ra,coords_dec)
    r_img = '{}/R_fin_{}.{}.{}-{}'.format(folder,xxx_input,yyy_input,coords_ra,coords_dec)
    
#     dp = 'dolphot dphot -punions.param img0_file={} img1_file={} img2_file={} img3_file={}'.format(r_img,r_img,g_img,u_img)
#     os.system(dp)
    # print(dp)
    
    # Print the stacked table
    if total is not None:
        print(total)
    cat_nm = 'total-{}-{}.cat'.format(ra_sky,dec_sky)
    total.to_csv(cat_nm)  
    n = len(total)

    print('----------------- DONE MATCHING ---------------')
    logging.info('-------------------- DONE MATCHING -------------')            
    
    if total.empty:
        shutil.rmtree(folder)
        print('-------------- DELETING FOLDER SINCE THERE ARE NO OBJECTS --------------')
    logging.info(f'Galaxy {gal} at {ra_sky,ra_dec} has {n} objects')


    
def plotting(folder,xxx_input,yyy_input,coords_ra,coords_dec,dist):
    
    d = dist
    xxx_input = xxx_input
    yyy_input = yyy_input
        

    img_og = '{}/CFIS.{}.{}.r.fits'.format(folder,xxx_input,yyy_input)
    img_crop = '{}/R_crop_{}.{}.{}-{}.fits'.format(folder,xxx_input,yyy_input,coords_ra,coords_dec)
    img_fin = '{}/R_fin_{}.{}.{}-{}.fits'.format(folder,xxx_input,yyy_input,coords_ra,coords_dec)
    
    hdu2 = fits.open(img_fin)
    og = fits.open(img_og)
    w = WCS(og[0].header)
    og.close()
    
    (xcen), (ycen) = w.all_world2pix([coords_ra], [coords_dec], 0)

    gal_ra,gal_dec = (float(coords_ra),float(coords_dec))
    rad = np.arcsin((60*np.power(10,3))/(d*np.power(10,6)))*206265
    rad = rad/3600/2
    dec_ll, ra_ll = gal_dec-rad, gal_ra+rad
    dec_ur, ra_ur = gal_dec+rad, gal_ra-rad

    (xmin, xmax), (ymin, ymax) = w.all_world2pix([ra_ll, ra_ur], [dec_ll, dec_ur], 0)
    
    cat_nm = '{}/total-{}-{}.cat'.format(folder,coords_ra,coords_dec)
    cat = pd.read_csv(cat_nm)


    cat['rad'] = np.sqrt(np.power((cat['X']-xcen),2)+np.power((cat['Y']-ycen),2))
    cat['g'] = cat['g'] - 2.7
    cat = cat[(cat['CLASS_r']>0.8)&(cat['u']<30)&(cat['r']<30)&(cat['r_err']<0.1)&(cat['g_err']<0.5)&(cat['rad']<4000)]
    
    hdu = fits.open(img_crop)

    PGC = int(hdu[0].header['PGC'])
    file_nm = '../plots/galaxy_{}.png'.format(PGC)
    dcm_nm = '../plots/dcm_{}.png'.format(PGC)
    print(file_nm,dcm_nm)

    
#     fig = plt.figure(figsize=(15,10))
#     w = WCS(hdu[0].header)
#     ax1 = fig.add_subplot(121, projection=w)
#     tr_fk5 = ax1.get_transform("fk5")
#     plt.title(PGC)
#     plt.scatter(cat['RA_r'],cat['DEC_r'],s=10,transform=tr_fk5,c=cat['CLASS_r'],cmap='Oranges')
#     plt.xlabel('RA')
#     plt.ylabel('DEC')
#     plt.xlim(xmin,xmax)
#     plt.ylim(ymin,ymax)
#     ax2 = fig.add_subplot(122, projection=w)
#     tr_fk5 = ax2.get_transform("fk5")

#     plt.imshow(hdu[0].data,vmax=10,vmin=0)
#     plt.xlabel('RA')
#     plt.ylabel('DEC')
 
#     plt.savefig(file_nm,bbox_inches='tight')
#     plt.close

    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(8,6))
    plt.title(PGC)
    ax1.scatter(cat['u']-cat['r'],cat['g']-cat['r'],s=5,c='orange')
    ax1.set_xlabel('u-r')
    ax1.set_ylabel('g-r')
    ax1.set_ylim(-1,6)
    ax1.set_xlim(-1,6)
    ax2.scatter(cat['r'],cat['CLASS_r'],s=5,c='orange')
    ax2.set_xlabel('R')
    ax2.set_ylabel('CLASS_STAR')
    plt.savefig(dcm_nm,bbox_inches='tight')
    plt.close() 
    


    flum = '../plots/analysis_{}.png'.format(PGC)
    print(flum)

#     fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(19,6))
#     sel = cat[(cat['r_err']<0.5)&((cat['u']-cat['r'])>0)&((cat['u']-cat['r'])<3)&((cat['g']-cat['r'])<1)&((cat['g']-cat['r'])>0.5)]

#     ax1.hist(sel['r'],bins=50,ec='orange',fc='None',linewidth=3)
#     ax1.set_xlabel('R')
#     ax1.set_ylabel('N')
#     ax1.set_xlim(20,25)
#     rad_sec = sel['rad'] * 0.185
    
#     ax2.hist(rad_sec,bins=30,ec='orange',fc='None',linewidth=3)
#     ax2.set_xlabel('rad [arcsec]')
#     ax2.set_ylabel('N')
    
#     cat_nm2 = '{}/central-{}-{}.cat'.format(folder,coords_ra,coords_dec)
#     cat2 = pd.read_csv(cat_nm2)
#     ycen2 = len(hdu[0].data)/2
#     xcen2 = len(hdu[0].data[0])/2
#     cat2['rad'] = np.sqrt(np.power((cat2['X']-xcen2),2)+np.power((cat2['Y']-ycen2),2))
#     cat2 = cat2[(cat2['CLASS_r']>0.5)&(cat2['u']<30)&(cat2['r']<30)]
#     sel2 = cat2[((cat2['u']-cat2['r'])>0)&((cat2['u']-cat2['r'])<4)]
    

#     sel2 = cat2[((cat2['u']-cat2['r'])>0)&((cat2['u']-cat2['r'])<4)]
#     rad_sec2 = sel2['rad'] * 0.185
#     ax3.hist(rad_sec2,bins=15,ec='orange',fc='None',linewidth=3)
#     ax3.set_xlabel('rad [arcsec]')
#     ax3.set_ylabel('N')


    fig, ax = plt.subplots(figsize=(10,3))
    sel = cat[(cat['r_err']<0.1)&((cat['u']-cat['r'])>0)&((cat['u']-cat['r'])<3)&((cat['g']-cat['r'])<1)&((cat['g']-cat['r'])>0.5)]
    # loctx = np.max(np.log(np.histogram(sel['r'],bins=18)[0]))
    ax.hist(sel['r'],bins=18,ec='orange',fc='None',linewidth=3)
    ax.set_xlabel('r')
    ax.set_ylabel('N')
    ax.set_xlim(18,24.7)
    ax.set_yscale('log')
    ng = 'PGC {}'.format(PGC)
    # ax.text(18.5,loctx,ng,fontsize=20)
#     rad_sec = sel['rad'] * 0.185
    
#     ax2.hist(rad_sec,bins=30,ec='orange',fc='None',linewidth=3)
#     ax2.set_xlabel('rad [arcsec]')
#     ax2.set_ylabel('N')
    
#     cat_nm2 = '{}/central-{}-{}.cat'.format(folder,coords_ra,coords_dec)
#     cat2 = pd.read_csv(cat_nm2)
#     ycen2 = len(hdu[0].data)/2
#     xcen2 = len(hdu[0].data[0])/2
#     cat2['rad'] = np.sqrt(np.power((cat2['X']-xcen2),2)+np.power((cat2['Y']-ycen2),2))
#     cat2 = cat2[(cat2['CLASS_r']>0.5)&(cat2['u']<30)&(cat2['r']<30)]
#     sel2 = cat2[((cat2['u']-cat2['r'])>0)&((cat2['u']-cat2['r'])<4)]
    

#     sel2 = cat2[((cat2['u']-cat2['r'])>0)&((cat2['u']-cat2['r'])<4)]
#     rad_sec2 = sel2['rad'] * 0.185
#     ax3.hist(rad_sec2,bins=15,ec='orange',fc='None',linewidth=3)
#     ax3.set_xlabel('rad [arcsec]')
#     ax3.set_ylabel('N')

    
    plt.savefig(flum,bbox_inches='tight')
    plt.close()
        



def mag_aper(folder,file):
    os.chdir(folder)
    data = pd.read_csv(file)

    # Order the DataFrame based on the "MAG" column
    ordered_data = data.sort_values(by='MAG_r')

    # Create a new table with the first 50 lines
    processed_table = Table.from_pandas(ordered_data.head(50))

    x_values = processed_table['x'].data
    y_values = processed_table['y'].data

 

if __name__ == "__main__":
#    parser = argparse.ArgumentParser(description="Execute specified functions.")
    cat = input('Enter galaxy catalogue\n')
    test = np.loadtxt(cat)
    start_time = time.time()
    result = subprocess.check_output(['wc', '-l', cat]).decode('utf-8')
    num_lines = int(result.split()[0])
    with open('output.txt', 'w') as file:
        for i in range(int(num_lines)):

            coords_ra = test[i][0]
            coords_dec = test[i][1]
            gal = test[i][2]
            check = in_UNIONS(coords_ra,coords_dec,listing='UNIONS_CFIS_LSB_r_DR5.txt')

            if check[0]==False:
                check = in_UNIONS(coords_ra,coords_dec,listing='tiles_DR2.list.txt')

            xxx_input = int(check[1])
            yyy_input = int(check[2])
            dist = test[i][3]


            xxx_input = '{:03d}'.format(xxx_input)
            yyy_input = '{:03d}'.format(yyy_input)

            if yyy_input==0:
                yyy_input = '000'
            folder = '../CFIS.{}.{}'.format(xxx_input,yyy_input)
#             try:
#                 download(coords_ra,coords_dec)
#             except:
#                 print('contained in',folder) 
            
#             try:
#                 hdul = detect(folder,coords_ra,coords_dec,gal,dist)  
#                 print('done with',coords_ra,coords_dec,'in folder',folder)
                
#             except:
#                 print('could not detect in',coords_ra,coords_dec,'in folder',folder)

#             os.chdir('/arc/home/aennis/scripts/unions_gcs/')
#             print(folder,xxx_input,yyy_input,coords_ra,coords_dec,dist)
            try:
                plotting(folder,xxx_input,yyy_input,coords_ra,coords_dec,dist)
                file.write(f'{gal} {coords_ra} {coords_dec} {dist} {xxx_input}{yyy_input}\n')
            except:
                print('no')

    end_time = time.time()
    elapsed_time = end_time - start_time
    print('Elapsed time downloading',elapsed_time)  

    
    
 