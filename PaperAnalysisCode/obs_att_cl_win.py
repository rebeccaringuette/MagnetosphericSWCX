# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 20:33:31 2020

@author: rebec

- run gti_index_cl.py first
- get ECEF position and velocity, RA, DEC, and times from att file for 
  promising MSWCX flank observations
- output into a fits file 
"""
from ruf_2_Loopd import tai_time_conv
import sys, datetime
import numpy as np
import astropy.io.fits as pyfits
from csv import DictReader

#set directories
prefix_dir = 'C:/Users/rebec/Documents/UIowa/'
code_dir = prefix_dir+'HaloSat_AnalysisSoft/CL_Code/'
att_dir= prefix_dir+'HaloSat_AnalysisSoft/MSWCX/'
sys.path.append(code_dir)
out_dir = prefix_dir+'/MSWCX_2020Analysis/Calculations/'

#get list of target RA and DECs
t = pyfits.open(code_dir+'LHB_nH_ecl3.fits')  
hsIDs = t[1].data['Target_ID']
TargetName = t[1].data['Target_Name']
RA = t[1].data['RA']
DEC = t[1].data['Dec']
t.close()

#get list of cleaned GTIs
cl_gtis = {}
gti_file = DictReader(open(out_dir+'MSWCX2021_GTI_cl.csv','r'))
for row in gti_file:
    for column, value in row.iteritems():
        cl_gtis.setdefault(column, []).append(value)
for key in cl_gtis: 
    if ((key=='Object') or (key=='Directory')): cl_gtis[key]=np.array(cl_gtis[key],dtype=str)
    elif ((key=='Target') or (key=='Obs#')): cl_gtis[key]=np.array(cl_gtis[key],dtype=int)
    else: cl_gtis[key] = np.array(cl_gtis[key],dtype=float)
  #print key

#output fits file with desired data for each target-obs pair
target_obs = ['77-1']#['70-4']#['71-3','75-5','76-9','77-10']
#['67-1','67-2','67-7','67-8','67-9','67-11']  #HSWCX2
#['66-9','66-10','369-3'] #HSWCX1
#['74-11','76-9']#['71-3','72-5','75-5','76-7','77-1']  #promising MSWCX flank observations
for tobs in target_obs:
    #get target number, ID, and nominal pointing (RA, DEC)
    target, obs = tobs.split('-')
    hsname=str('{:04d}01'.format(int(target)))
    idx = np.where(hsIDs==int(target))[0]
    Tname, target_RA, target_DEC = TargetName[idx][0],RA[idx][0],DEC[idx][0]
    print Tname, target_RA, target_DEC
    
    #pull gtis for target, obs
    idx=np.where((cl_gtis['Target']==int(target))&(cl_gtis['Obs#']==int(obs)))[0]
    start_tai, stop_tai = cl_gtis['TAI_Start'][idx], cl_gtis['TAI_Stop'][idx]
    date1 = tai_time_conv(int(min(start_tai)),refout='')  #UTC
    date2 = tai_time_conv(int(max(stop_tai)),refout='')    
    to_idx = idx[0]  #save first index for dir retreival
    print 'Target ',target,' Observation ',obs,'\nStart (UTC):',date1,'\nStop (UTC):',date2,'\n'
    
    #read in att data
    att_file = att_dir+hsname+'/unfiltered/hs'+hsname+'.att'
    f = pyfits.open(att_file) # read in the FITS file with att data
    att_data = f[1].data # put all of the data into the FITS record array d
    f.close() # close the FITS file    
    
    #collect att data for desired times
    idx_list,date_times=[],[]
    for time in range(len(start_tai)):
        idx = np.where((att_data['TIME']>=start_tai[time])&(att_data['TIME']<=stop_tai[time]))[0]
        idx_list.extend(idx)
        for i in idx: date_times.append(tai_time_conv(int(att_data['TIME'][i]),refout=''))  #UTC
    
    #output att_out to fits file
    file_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')

    # make the primary header (standard in all files)
    prihdr = pyfits.Header()
    prihdr['TELESCOP'] = ('HALOSAT', 'Telescope (mission) name')
    prihdr['INSTRUME'] = ('SDD', 'Instrument name') 
    prihdr['OBS_ID'] = (hsname, 'Observation ID')    
    prihdr['OBJECT'] = (Tname, 'Object/Target Name')
    prihdr['OBJTYPE'] = ('SWCX', 'Object/Target type')  
    prihdr['DATE-OBS'] = (date1, 'Start date of observations in UTC')
    prihdr['DATE-END'] = (date2, 'End date of observations in UTC')
    prihdr['DATE'] = (file_time, 'File creation date')  #checksum and datasum are added with fverify etc
    prihdu = pyfits.PrimaryHDU(header=prihdr) # make the primary header    
    
    #make data table
    hk_col=[]
    hk_col.append(pyfits.Column(name='TIME', format='1D', unit='s', array=att_data['TIME'][idx_list])) 
    hk_col.append(pyfits.Column(name='DATE_TIME', format='24A', array=date_times))
    hk_col.append(pyfits.Column(name='POSITION', format='3E', unit='km', array=att_data['POSITION'][idx_list]))   
    hk_col.append(pyfits.Column(name='VELOCITY', format='3E', unit='km/s', array=att_data['VELOCITY'][idx_list]))
    hk_col.append(pyfits.Column(name='RA', format='1D', unit='deg', array=att_data['RA'][idx_list])) 
    hk_col.append(pyfits.Column(name='DEC', format='1D', unit='deg', array=att_data['DEC'][idx_list])) 
    cols = pyfits.ColDefs(hk_col) # create a ColDefs (column-definitions) object for all columns:
    
    #make data table header
    hkhdu = pyfits.BinTableHDU.from_columns(cols) #, [], 'SPECTRUM') # create a new binary table HDU object  
    hkhdu.header.comments['TTYPE1'] = 'Time (TAI)'   
    hkhdu.header.comments['TFORM1'] = 'data format of field'
    hkhdu.header.comments['TUNIT1'] = 'physical unit of field'    
    hkhdu.header.comments['TTYPE2'] = 'Date Time (UTC) (YY-MM-DD HH:MM:SS)'   
    hkhdu.header.comments['TFORM2'] = 'data format of field' 
    hkhdu.header.comments['TTYPE3'] = 'ECEF position of satellite [X,Y,Z]'   
    hkhdu.header.comments['TFORM3'] = 'data format of field'
    hkhdu.header.comments['TUNIT3'] = 'physical unit of field'  
    hkhdu.header.comments['TTYPE4'] = 'ECEF velocity of satellite [vX,vY,vZ]'   
    hkhdu.header.comments['TFORM4'] = 'data format of field'
    hkhdu.header.comments['TUNIT4'] = 'physical unit of field'      
    hkhdu.header.comments['TTYPE5'] = 'Right Ascension of the pointing' 
    hkhdu.header.comments['TFORM5'] = 'data format of field'
    hkhdu.header.comments['TUNIT5'] = 'physical unit of field'
    hkhdu.header.comments['TTYPE6'] = 'Declination of the pointing' 
    hkhdu.header.comments['TFORM6'] = 'data format of field'
    hkhdu.header.comments['TUNIT6'] = 'physical unit of field'    

    #set standard extension header keys, updated Dec12 2019 RR
    hkhdu.header['EXTNAME'] = ('HK', 'Binary table extension name')  
    hkhdu.header['HDUCLASS'] = ('OGIP', 'format conforms to OGIP/GSFC standards')
    hkhdu.header['HDUCLAS1'] = ('TEMPORALDATA', 'First class level')  
    hkhdu.header['HDUCLAS2'] = ('HK', 'Second class level')
    hkhdu.header['TELESCOP'] = ('HALOSAT', 'Telescope (mission) name')
    hkhdu.header['INSTRUME'] = ('SDD', 'Instrument name')   
    hkhdu.header['OBSERVER'] = ('PHILIP KAARET', 'Principal Investigator')
    hkhdu.header['OBS_ID'] = (hsname, 'Observation ID')
    hkhdu.header['OBJECT'] = (Tname, 'Object/Target Name')
    hkhdu.header['OBJTYPE'] = ('SWCX', 'Object/Target type')    
    hkhdu.header['EQUINOX'] = (2000, '[yr] Equinox of celestial coord system')
    hkhdu.header['RADECSYS'] = (target_RA, '[deg] R.A. of nominal aspect point [J2000]')  
    hkhdu.header['DEC_NOM'] = (target_DEC, '[deg] Dec. of nominal aspect point [J2000]')  
    hkhdu.header['RA_OBJ'] = (target_RA, '[deg] Object Right ascension [J2000]')
    hkhdu.header['DEC_OBJ'] = (target_DEC, '[deg] Object Declination [J2000]')
    hkhdu.header['TIMESYS'] = ('TT', 'Reference Time System')
    hkhdu.header['MJDREFI'] = (51544, '[d] MJD reference day (2000-01-01T00:00:00)')
    hkhdu.header['MJDREFF'] = (0.00074287037037037, '[d] MJD reference (fraction of day)')
    hkhdu.header['TIMEREF'] = ('LOCAL', 'Reference Frame')
    hkhdu.header['TASSIGN'] = ('SATELLITE', 'Time assigned by clock')
    hkhdu.header['TIMEUNIT'] = ('s', 'Time unit for timing header keyword')
    hkhdu.header['TIMEDEL'] = (8.0,'[s] Data time resolution.')   
    hkhdu.header['TIMEZERO'] = (0.0, '[s] Time Zero')
    hkhdu.header['TIMEPIXR'] = (1, 'Bin time beginning=1 middle=0.5 end=1')
    hkhdu.header['TIERRELA'] = (4.0E-06, '[s/s] relative errors expressed as rate')
    hkhdu.header['TIERABSO'] = (1.0, '[s] timing precision in seconds')
    hkhdu.header['TSTART'] = (min(start_tai), '[s] Observation start time')
    hkhdu.header['TSTOP'] = (max(stop_tai), '[s] Observation stop time')
    hkhdu.header['DATE-OBS'] = (date1, 'Start date of observations')
    hkhdu.header['DATE-END'] = (date2, 'End date of observations') 
    hkhdu.header['CLOCKAPP'] = (True, 'Clock correction applied? (F/T)')
    hkhdu.header['ORIGIN'] = ('UNIVERSITY OF IOWA', 'Origin of fits file')
    hkhdu.header['PROCVER'] = ('hsuf_20200226', 'Processing script version number')   #prev. hsuf_20200131_hscl_20200131
    hkhdu.header['SOFTVER'] = ('Hea_26jun2019_V6.26.1', 'Software version')   
    hkhdu.header['CALDBVER'] = ('hs20200129', 'CALDB index version used')   
    hkhdu.header['TLM2FITS'] = ('db_20200218', 'Telemetry converter FITS number')   #location?!  
    hkhdu.header['CREATOR'] = ('db_hsuf_obs', 'Software creator of the file')
    hkhdu.header['DATE'] = (file_time, 'File creation date') 
    
    # join primary with extension, write, and check
    hdulist = pyfits.HDUList([prihdu, hkhdu])
    att_filename='/hs'+hsname+'-'+obs+'_cl.att'
    hdulist.writeto(out_dir+att_filename, overwrite=True)    

    #gzip file - don't have imcopy on computer
    #if os.path.exists(out_dir+att_filename+'.gz'): os.remove(out_dir+att_filename+'.gz')
    #os.system('imcopy '+out_dir+att_filename+' '+out_dir+att_filename+'.gz')  
    #stop 