# -*- coding: utf-8 -*-
"""
FitHalo_v1.py: PyXspec for Halo model calculation
- Initial code taken from BKG_fit4_doubleplv5.py
- Input spectra: VLE <= 2.5 cts/s, LCHARD <= x cts/s for a single halo field
- find reduced X^2 for each fit, print on each spectral fit and output to files
- plot norm and error of norm of fitted O lines in model
- model: gauss+gauss+TBabs(vapec+bknpower+bknpower) 
    - double bknpower for CXB (for now) from Smith et al 2007
    - 2 gaussians for O lines
    - same absorption for CXB as for halo?!
    - APEC with O set to zero in model

v1p1: loop through detectors instead of count rate bands

v1p2: Using newer recommended cuts (XLowBkgv4), and new version of APEC model 
    with O VII and O VIII lines removed in O band
    
FitHalo_v2
    - change to fit only CXB and particle bkg in 3-7 keV
    - using XLowBkg9, net high halo file
    
v2p1
    - test calculated PI and norms for pegpwrlw on MidHigh gal lat halo fields
    
v2p3
    - use calcflux attribute to calculate fixed and free BKG flux in O band
    - calculate fixed flux with E range = 0.4 to 7 keV
    - fit for parameters with E range = 3 to 7 keV, then freeze
    - calculate fitted flux with E range = 0.4 to 7 keV
    
v2p4
    - use calculated BKG parameters for pegpwrlw through diag response
    - fit for CXB in 2 to 7 keV range using pow
    - use calcFlux to determine flux in 2 to 7 keV band
    - for individual fields

v2p5
    - similar to v2p4 for high halo summed field  
    - link model to all three detectors
    
FitDark_Emily.py: 
    - loops through target numbers, assuming heasarc cl structure
    - loops through set of error commands if new fit is found
    - includes all spectral data in python plot (black, red, blue) 
    - spectra somewhat hardcoded to use all three detectors
    
BKG_FullERange.py:
    - use code to determine which parameters are free instead of hardcoding it
    - fitting full halo model over 0.4 - 7 energy range, fitting BKG PI and N 
    - avoid hard-coding det#s and parameters names (mostly)
    - dynamically set plot limits in xspec plot, add reduced X^2 in xspec plot
    - commented out python version of plot since xspec version is improved
    - tested on first halo field
    
Halo_FullERange.py
    - assign BKG PI to fixed value for halo fields, diff output names
    - Checking for new fits found during error process doesn't work with error code,
      possibly because query=yes. Comparing fitted values at beginning vs end instead.
    - this method results in looping even if now new fit is found.
      
Halo_FullERange_v1.py      
    - alternate method for error loop if new fit found.
    - does not loop if no new fit is found, does if one is
    
Halo_FullERange_v2.py
    - use new BKG_PI calculation based on HARD ctrt
    - use new nH_LAB file for nH and LAB N values       

v3: for HS0373

v4: Using updated full SWCX model excluding MgXI

HSWCX2_v1.py: alter to run HSWCX2 GOOD obs. fitting BKG_PIs and BKG_Ns, one temp apec.
    Code adapted from HSWCX1_v116.py
v2: 2 apecs
v3: fixed halo, fit SWCX

v4: same as v1, using new SWAN-based helio calcs.
v5: same as v2, using new SWAN-based helio calcs.
v6: same as v3, using v4 results
v6H2: using v4 results for GOODH2 (2-temp halo)
v7: scanning through SWCX scaling factors
v8: same as v7 but using Model 1
"""

'''
############################################################################
****************************************************************************
------------------------Start PyXspec program-------------------------------
****************************************************************************
############################################################################
'''
#from __future__ import division 

import xspec as xs
import numpy as np
import glob, itertools, os
import astropy.io.fits as pyfits
from csv import DictReader
#pyxspec only works with python 2.7 commands since that is the version installed
#pyxspec also requires starting spyder from the terminal
#pyxspec won't output fitting details to screen if I print to the screen (sometimes)
xs.Xset.addModelString("APECROOT","/home/ringuette/Software/heasoft-6.25/spectral/modelData/apec_v3.0.9_201") #use new apec file

#define directories
spec_dir = '/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/Rebinned_spectra_cl1/'
code_dir='/home/ringuette/halosat_sourceanalysis/AnalysisCode/HSWCXCode/'
out_dir='/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/XspecOutput_cl1_HSWCX2/'
prefix = 'HSWCX2v8' #default is SWAN-based HSWCX calcs (model 3)
HSWCXmod = '1'   #3 for SWAN-based calcs
        
#define options
#OVII_list = np.array(['2.00','2.50','3.00','3.50','3.60','3.70','3.75','3.80','3.85','3.90','4.00','4.50','5.00','5.50'],dtype=str) #a list
#OVII_list = np.array(['1.00','1.50','2.30','2.40','2.60','2.70','2.80','2.90'],dtype=str) #b list
OVII_list = np.array(['1.60','1.70','2.55','2.56','2.57','2.58','2.59','2.61','2.62','2.63','2.64','2.65'],dtype=str) #c list
nH_type = '_nHKip' #set type of nH used options: '_nHLAB', '_nHHI4pi', '_nHSFDW', '_nHPZ'
model_note = nH_type+'OVII_c'
CXB_PI, CXB_N = 1.45, 0.382   #Capellutti et al. 2017

#initialize output format
file_header='Target,LHB_N,'+nH_type.split('_')[1]+',OVII_scale' #all other headers not hard-coded
file_format='\n{0[0]},{0[1]:.6f},{0[2]:.3f},{0[3]}'
csv_num=4  #next number in []

#pull LHB N from fits file, store in dict 
t = pyfits.open(code_dir+'LHB_nH_ecl3.fits') 
hsIDs = t[1].data['Target_ID']
RA = t[1].data['RA']
DEC = t[1].data['Dec']
glat = t[1].data['Gal_Lat']
glon = t[1].data['Gal_Lon']
LHB = t[1].data['LHB_norm']
nh_HI4Pi = t[1].data['nH_HI4pi']
nh_SFDW = t[1].data['nH_SFDW']
nh_PZ = t[1].data['nH_PlanckZhu']
nh_RT = t[1].data['nH_RT']  #Rachford_Tau
t.close()

if nH_type=='_nHRT': nh = nh_RT
elif nH_type=='_nHHI4pi': nh = nh_HI4Pi
elif nH_type=='_nHSFDW': nh = nh_SFDW
elif nH_type=='_nHPZ': nh = nh_PZ
elif nH_type=='_nHKip': nh = np.repeat(0.0661,len(hsIDs))
nh_target, LHB_target = nh[66], LHB[66]  #all are for target 67

#read helioSWCX calcs per spec
SWCX_LU = {}
SWCX_file = DictReader(open(code_dir+'helioSWCX.csv','r'))
for row in SWCX_file:
    for column, value in row.iteritems():
        SWCX_LU.setdefault(column, []).append(value)
#print SWCX_LU
for key in SWCX_LU: 
  SWCX_LU[key] = np.array(SWCX_LU[key])
  #print key, SWCX_LU[key][0], type(SWCX_LU[key][0])
#SWCX line ratios   
OVII_ratio = '0.126'
OVIII_ratio = '0.549'  #choose He ratio since through He-cone
NeIX_ratio = '0.100'
LowE_ratio = '0.925'   #also 100% He

key_list = ['GOODH2']#'Dec2018','Dec2019a','Dec2019bc','Dec2019d','Feb2020']  #HSWCX2 observations

for OVII_Oratio in OVII_list:

    #loop through target numbers
    for key in key_list:
        spec_name='HSWCX2_'+key+'_cl1_d'
        spec_list=[spec_dir+spec_name+'14_rebin40.pi',spec_dir+spec_name+'38_rebin40.pi',spec_dir+spec_name+'54_rebin40.pi']
    
        #collect parameters for observation
        idx = np.where(SWCX_LU['\xef\xbb\xbfObsName']==key)[0]
        OVII_N = float(SWCX_LU['OVII_Mod'+HSWCXmod][idx][0])*0.035*float(OVII_Oratio)
        OVIII_N = float(SWCX_LU['OVIII_Mod'+HSWCXmod][idx][0])*0.035*float(OVII_Oratio)
        NeIX_N = float(SWCX_LU['NeIX_Mod'+HSWCXmod][idx][0])*0.035*float(OVII_Oratio)
        #print key, OVII_N/0.035, OVIII_N/0.035, NeIX_N/0.035
        
        #setup xspec
        xs.AllData.clear()
        xs.AllModels.clear()
        xs.Fit.statMethod='chi' #statistic chi
        xs.Fit.statTest='chi'  #statistic test chi --> different from Phil?!
        xs.Xset.abund='wilm' #abund wilm
        xs.Xset.xsect='vern' #xsect vern
        xs.Xset.cosmo='70 0 0.73'  #cosmo 70 0 0.73
        xs.Fit.method='leven  10  0.01'  #method leven 10 0.01
        xs.Fit.query = "yes"
        
        #initialize variables
        bkg_par, xray_par, SWCX_par, TTS_par = {}, {}, {}, {}
        bmodel, xmodel, smodel, tmodel, bkg_spec = [], [], [], [], []
        halo_arr=[key,LHB_target,nh_target,OVII_Oratio]
        
        #loop through detector spectra
        print 'Initializing '+key
        for i in range(len(spec_list)): #fit model to each spectrum linked between detectors    
            bkg_spec.append(0)
            bmodel.append(0)
            xmodel.append(0)
            smodel.append(0)
            #tmodel.append(0)
            
            #set up files    
            xs.AllData(str(i+1)+':'+str(i+1)+' '+spec_list[i]) #assign each spectrum to diff data groups
            bkg_spec[i]=xs.AllData(i+1)    
            bkg_spec[i].multiresponse[0] = code_dir+'halosat_diag_20190423.rmf' #only rmf for particle bkg fit
            bkg_spec[i].multiresponse[1] = code_dir+'halosat_avgenoise_20190423.rmf'
            bkg_spec[i].multiresponse[1].arf = code_dir+'halosat_20190211.arf'
            bkg_spec[i].multiresponse[2] = code_dir+'halosat_avgenoise_20190423.rmf'
            bkg_spec[i].multiresponse[2].arf = code_dir+'halosat_20190211.arf'  
                       
            #add exposure time and bkg_pi to output
            dpu = int(spec_list[i].split('.')[0].split('_d')[1].split('_rebin40')[0]) #get dpu num
            halo_arr.append(bkg_spec[i].exposure)
            if (OVII_Oratio==OVII_list[0]) and (key==key_list[0]): 
              file_header+=',d'+str(dpu)+' Exp' #all other headers not hard-coded
              file_format+=',{0['+str(csv_num)+']:.3f}'
              csv_num+=1  #next number in [] 
                      
            #set up model
            bkg_par[i]={1:'0.75,0.2,0.4,0.6,1.0,1.3',2:'0.02,0.02,0.,0.01,0.5,1e6'}  #BKG PI fixed and N free per det
        
            #using Liu values for LHB (fixed apec), 
            #absorbed CXB (fixed to Capelluti values) and absorbed linked halo APEC
            #SWCX lines from Alpha file from Dimitra
            if (i==0):    
                xray_par[i]={1:'0.097,-1',2:'1,-1',3:'0,-1',4:str(LHB_target)+',-1',\
                              5:str(nh_target)+',-1',\
                              6:'0.22,0.05,0.1,0.15,0.3,0.35',7:'0.3,-1',8:'0,-1',9:'1.0,0.1,0,0.1,1e20,1e24',\
                              10:str(CXB_PI)+',-1',11:str(CXB_N)+',-1',\
                              12:'0.75,0.1,0.6,0.7,0.9,1.5',13:'0.3,-1',14:'0,-1',15:'1.0,0.1,0,0.1,1e20,1e24'}
                SWCX_par[i] = {1:'0.44339,-1',2:'0.001,-1',3:'=SWCX:p6*'+LowE_ratio,\
                              4:'0.56325,-1',5:'0.001,-1',6:str(OVII_N)+',-1',\
                              7:'0.67916,-1',8:'0.001,-1',9:'=SWCX:p6*'+OVII_ratio,\
                              10:'0.65310,-1',11:'0.001,-1',12:str(OVIII_N)+',-1',\
                              13:'0.80308,-1',14:'0.001,-1',15:'=SWCX:p12*'+OVIII_ratio,\
                              16:'0.90873,-1',17:'0.001,-1',18:str(NeIX_N)+',-1',\
                              19:'1.10043,-1',20:'0.001,-1',21:'=SWCX:p18*'+NeIX_ratio} 
    
            else: 
                xray_par[i]={1:'=Xray:p1',2:'=Xray:p2',3:'=Xray:p3',4:'=Xray:p4',\
                              5:'=Xray:p5',\
                              6:'=Xray:p6',7:'=Xray:p7',8:'=Xray:p8',9:'=Xray:p9',\
                              10:'=Xray:p10',11:'=Xray:p11',\
                              12:'=Xray:p12',13:'=Xray:p13',14:'=Xray:p14',15:'=Xray:p15'}
                SWCX_par[i]={1:'0.44339,-1',2:'0.001,-1',3:'=SWCX:p3',\
                              4:'0.56325,-1',5:'0.001,-1',6:'=SWCX:p6',\
                              7:'0.67916,-1',8:'0.001,-1',9:'=SWCX:p9',\
                              10:'0.65310,-1',11:'0.001,-1',12:'=SWCX:p12',\
                              13:'0.80308,-1',14:'0.001,-1',15:'=SWCX:p15',\
                              16:'0.90873,-1',17:'0.001,-1',18:'=SWCX:p18',\
                              19:'1.10043,-1',20:'0.001,-1',21:'=SWCX:p21'}
                              
        #set dictionary of component names for free parameters. change if model changes
        xray_names={6:'Halo_kT1',9:'Halo_N1',12:'Halo_kT2',15:'Halo_N2'}
        SWCX_names={3:'LowE',6:'O_VII',12:'O_VIII',18:'Ne_IX'}
        
        #set Xspec model and parameters
        bkg_model=xs.Model('pow','BKG',1)  #particle background
        xray_model=xs.Model('apec+tbabs*(apec+pow+apec)','Xray',2)   #LHB+halo+CXB
        SWCX_model=xs.Model('gauss+gauss+gauss+gauss+gauss+gauss+gauss','SWCX',3)
        #TTS_model=xs.Model('tbabs*apec','TTS',4)  #hot thermal component from TTS
        for i in range(len(spec_list)):    
            bmodel[i]=xs.AllModels(i+1,'BKG')
            xmodel[i]=xs.AllModels(i+1,'Xray')
            smodel[i]=xs.AllModels(i+1,'SWCX')
            #tmodel[i]=xs.AllModels(i+1,'TTS')
            bmodel[i].setPars(bkg_par[i])
            xmodel[i].setPars(xray_par[i])
            smodel[i].setPars(SWCX_par[i])
            #tmodel[i].setPars(TTS_par[i])
            
            xs.AllData.ignore('0.0-0.4,7.-**')
            xs.AllModels.setEnergies('reset')  #set to use energies instead of channels?
            xs.AllModels.systematic = 0  #systematic 0
        
        # display information for all loaded spectra, perform fit
        xs.AllModels.show()
        xs.AllData.show()
        xs.Fit.nIterations = 1000 # use lots of iterations
        xs.Fit.renorm()
        xs.Fit.perform() #execute fit
        
        #run steppar on all free bkg parameters
        for i, j in itertools.product(range(len(spec_list)), range(bmodel[i].nParameters)):
            if not ((bmodel[i](j+1).frozen)or(len(bmodel[i](j+1).link)>0)):
                parnum=bmodel[i].startParIndex+bmodel[i](j+1).index-1
                if (bmodel[i](j+1).index==1):  #BKG PI
                    print 'steppar: BKG', parnum
                    xs.Fit.steppar('best BKG:'+str(parnum)+' 0.4 1.3 10')                
                if (bmodel[i](j+1).index==2):  #BKG N
                    print 'steppar: BKG', parnum
                    xs.Fit.steppar('best BKG:'+str(parnum)+' 0. 0.1 10')
        #run steppar on all free xray parameters
        for i, j in itertools.product(range(len(spec_list)), range(xmodel[i].nParameters)):
            if not ((xmodel[i](j+1).frozen)or(len(xmodel[i](j+1).link)>0)):
                parnum=xmodel[i].startParIndex+xmodel[i](j+1).index-1
                if (xmodel[i](j+1).index==6):  #halo kT1
                    print 'steppar: Xray', parnum
                    xs.Fit.steppar('best Xray:'+str(parnum)+' 0.15 0.3 15')
                elif (xmodel[i](j+1).index==9):  #halo N1
                    print 'steppar: Xray', parnum
                    xs.Fit.steppar('best Xray:'+str(parnum)+' 0. 1.5 15')    
                elif (xmodel[i](j+1).index==12):  #halo kT2
                    print 'steppar: Xray', parnum
                    xs.Fit.steppar('best Xray:'+str(parnum)+' 0.7 1.3 10')
                elif (xmodel[i](j+1).index==15):  #halo N2
                    print 'steppar: Xray', parnum
                    xs.Fit.steppar('best Xray:'+str(parnum)+' 0. 1.5 15')                    
        #run steppar on all free SWCX parameters
        for i, j in itertools.product(range(len(spec_list)), range(smodel[i].nParameters)):
            if not ((smodel[i](j+1).frozen)or(len(smodel[i](j+1).link)>0)):
                parnum=smodel[i].startParIndex+smodel[i](j+1).index-1
                print 'steppar: SWCX', parnum    #line Ns
                xs.Fit.steppar('best SWCX:'+str(parnum)+' 0. 0.1 25') 

        #save current parameter values for later comparison after error runs
        test=[xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof]
                    
        #open log file, print info to file, including error runs
        logFile = xs.Xset.openLog(out_dir+prefix+model_note+str(np.where(OVII_list==OVII_Oratio)[0][0])+'_'+key+'_xspeclog.txt')
        xs.Xset.show()
        xs.AllData.show()  
        if (xs.Fit.statistic/xs.Fit.dof<2.):
            #run error on all free BKG parameters
            for i, j in itertools.product(range(len(spec_list)), range(bmodel[i].nParameters)):
                if not ((bmodel[i](j+1).frozen)or(len(bmodel[i](j+1).link)>0)):
                    parnum=bmodel[i].startParIndex+bmodel[i](j+1).index-1     
                    print 'error: BKG', parnum
                    xs.Fit.error('BKG:'+str(parnum))
            #run error on all free Xray parameters
            for i, j in itertools.product(range(len(spec_list)), range(xmodel[i].nParameters)):
                if not ((xmodel[i](j+1).frozen)or(len(xmodel[i](j+1).link)>0)):
                    parnum=xmodel[i].startParIndex+xmodel[i](j+1).index-1
                    print 'error: Xray', parnum
                    xs.Fit.error('Xray:'+str(parnum))            
            #run error on all free SWCX parameters
            for i, j in itertools.product(range(len(spec_list)), range(smodel[i].nParameters)):
                if not ((smodel[i](j+1).frozen)or(len(smodel[i](j+1).link)>0)):
                    parnum=smodel[i].startParIndex+smodel[i](j+1).index-1
                    print 'error: SWCX', parnum
                    xs.Fit.error('SWCX:'+str(parnum)) 
        xs.AllModels.show()
        xs.Fit.show()   
        xs.Xset.closeLog()  
    
        #save fit statistics for comparison with before error runs
        teste, error_check, error_loop = [xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof], 0, 0
                
        #if new parameters values differ by more than 1%, rerun errors
        if max(abs(np.array(test)-np.array(teste))/np.array(test))>0.001: error_check+=1 
        while error_check>0:
            #redo errors with log file
            print 'New fit found. Rerunning error process.'
            error_loop+=1
            logFile = xs.Xset.openLog(out_dir+prefix+model_note+str(np.where(OVII_list==OVII_Oratio)[0][0])+'_'+key+'_xspeclog'+str(error_loop)+'.txt')
            xs.Xset.show()
            xs.AllData.show()  
            if (xs.Fit.statistic/xs.Fit.dof<2.):
                #run error on all free BKG parameters
                for i, j in itertools.product(range(len(spec_list)), range(bmodel[i].nParameters)):
                    if not ((bmodel[i](j+1).frozen)or(len(bmodel[i](j+1).link)>0)):
                        parnum=bmodel[i].startParIndex+bmodel[i](j+1).index-1     
                        print 'error: BKG', parnum
                        xs.Fit.error('BKG:'+str(parnum))
                #run error on all free Xray parameters
                for i, j in itertools.product(range(len(spec_list)), range(xmodel[i].nParameters)):
                    if not ((xmodel[i](j+1).frozen)or(len(xmodel[i](j+1).link)>0)):
                        parnum=xmodel[i].startParIndex+xmodel[i](j+1).index-1
                        print 'error: Xray', parnum
                        xs.Fit.error('Xray:'+str(parnum))      
                #run error on all free SWCX parameters
                for i, j in itertools.product(range(len(spec_list)), range(smodel[i].nParameters)):
                    if not ((smodel[i](j+1).frozen)or(len(smodel[i](j+1).link)>0)):
                        parnum=smodel[i].startParIndex+smodel[i](j+1).index-1
                        print 'error: SWCX', parnum
                        xs.Fit.error('SWCX:'+str(parnum))   
            xs.AllModels.show()
            xs.Fit.show()   
            xs.Xset.closeLog()  
            
            #check new error runs for new fits, loop again until no new fits are found
            test=teste  #last values become 'original' values to compare to
            teste, error_check = [xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof], 0  #reset values     
            #if new parameters values differ by more than 1%, rerun errors
            if max(abs(np.array(test)-np.array(teste))/np.array(test))>0.001: error_check+=1      
    
        #save BKG output to array
        for i, j in itertools.product(range(len(spec_list)), range(bmodel[i].nParameters)):
            if not ((bmodel[i](j+1).frozen)or(len(bmodel[i](j+1).link)>0)):
                idxnum = bmodel[i](j+1).index   
                det = spec_list[i].split('/')[-1].split('_d')[1].split('_')[0]
                halo_arr.extend([bmodel[i](idxnum).values[0],bmodel[i](idxnum).values[0]-bmodel[i](idxnum).error[0],\
                                 bmodel[i](idxnum).error[1]-bmodel[i](idxnum).values[0],bmodel[i](idxnum).error[2]])  #save par values, LL, UL, and code
                if (OVII_Oratio==OVII_list[0]) and (key==key_list[0]):  #first key
                    if bmodel[i](idxnum).name=='PhoIndex': bname='PI'
                    elif bmodel[i](idxnum).name=='norm': bname='N'
                    file_header+=',d'+str(det)+' BKG_'+bname+',ErrL,ErrU,ErrCode'
                    file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']:.6f},{0['\
                        +str(csv_num+2)+']:.6f},{0['+str(csv_num+3)+']}'
                    csv_num+=4                  
                    
        #save Xray output to array
        for i, j in itertools.product(range(len(spec_list)), range(xmodel[i].nParameters)):
            if not ((xmodel[i](j+1).frozen)or(len(xmodel[i](j+1).link)>0)):
                idxnum = xmodel[i](j+1).index  #det# not important since all free pars are linked   
                det = spec_list[i].split('/')[-1].split('_d')[1].split('_')[0]
                halo_arr.extend([xmodel[i](idxnum).values[0],xmodel[i](idxnum).values[0]-xmodel[i](idxnum).error[0],\
                                 xmodel[i](idxnum).error[1]-xmodel[i](idxnum).values[0],xmodel[i](idxnum).error[2]])  #Line N
                if (OVII_Oratio==OVII_list[0]) and (key==key_list[0]):  #first key
                    file_header+=','+xray_names[idxnum]+',ErrL,ErrU,ErrCode'   
                    file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']:.6f},{0['\
                        +str(csv_num+2)+']:.6f},{0['+str(csv_num+3)+']}'
                    csv_num+=4   
    
        #save SWCX output to array
        for i, j in itertools.product(range(len(spec_list)), range(smodel[i].nParameters)):
            if not ((smodel[i](j+1).frozen)or(len(smodel[i](j+1).link)>0)):
                idxnum = smodel[i](j+1).index  #det# not important since all free pars are linked   
                halo_arr.extend([smodel[i](idxnum).values[0],smodel[i](idxnum).values[0]-smodel[i](idxnum).error[0],\
                                 smodel[i](idxnum).error[1]-smodel[i](idxnum).values[0],smodel[i](idxnum).error[2]])  #Line N
                if (OVII_Oratio==OVII_list[0]) and (key==key_list[0]):  #first key
                    file_header+=','+SWCX_names[idxnum]+',ErrL,ErrU,ErrCode'   
                    file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']:.6f},{0['\
                        +str(csv_num+2)+']:.6f},{0['+str(csv_num+3)+']}'
                    csv_num+=4   
                
        halo_arr.extend([xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof])  #add fit statistics to list
        if (OVII_Oratio==OVII_list[0]) and (key==key_list[0]):  #initialize csv file header
            file_header+=',X^2,dof,red_X^2'
            file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']},{0['\
                +str(csv_num+2)+']:.10f}'   
            chi_file=open(out_dir+prefix+model_note+'.csv','w') #initialize output file
            chi_file.write(file_header)
            chi_file.close()        
            
        
        #output results to csv file spec_dir+run_type
        with open(out_dir+prefix+model_note+'.csv','a') as chi_file: 
            chi_file.write(file_format.format(halo_arr))
        
        '''
        #-------------------------------------------------------------------
        #--- Plotting Section   
        #-------------------------------------------------------------------
        '''
        #create python version of plot
        xs.Plot.commands = ()    
        xs.Plot.device = '/xs' # '/cps' # 
        xs.Plot.xLog = True # True
        xs.Plot.yLog = True
        xs.Plot.splashPage = False
        xs.Plot.add = True
        xs.Plot.xAxis = 'keV'
        xs.Plot("data delchi")   #change for each error type
        redchisq_text=str('Chi^2/dof = {:.3f}'.format(xs.Fit.statistic/xs.Fit.dof))  #add reduced X^2 to plot
        
        #collect data
        psym={0:'k',1:'r',2:'b'}
        xg, rates, xErrs, yErrs, folded, x, y = {}, {}, {}, {}, {}, {}, {}
        #chisq, chisq_unc, chisq_err = {}, {}, {}
        delchi, delchi_unc, delchi_err = {}, {}, {}
        #resid, resid_unc, resid_err = {}, {}, {}
        xmax, xmin, ymax, ymin, y2min = 0., 10., 0., 10., 10.
        for i in range(len(spec_list)):
            xg[i] = np.array(xs.Plot.x(i+1,1),dtype=float)
            rates[i] = np.array(xs.Plot.y(i+1,1),dtype=float)
            xErrs[i] = np.array(xs.Plot.xErr(i+1,1),dtype=float)
            yErrs[i] = np.array(xs.Plot.yErr(i+1,1),dtype=float)
            #only plotting one spectral fit (not for each dpu)
    
            #find good maxima for spectral plots
            xmax = max([xmax,max(xg[i]+xErrs[i])])
            xmin = min([xmin,min(xg[i]-xErrs[i])]) 
            ridx = np.where(rates[i]>0.)[0]  #avoid using zero bins to determine mins
            ymax = max([ymax,max(rates[i][ridx]+yErrs[i][ridx])])
            ymin = min([ymin,min(rates[i][ridx]-yErrs[i][ridx])])  
            y2min = min([y2min,min(rates[i][ridx])])
            
        #check/adjust maxima values
        if xmin<0: xmin=0
        if ymin<=0: ymin=y2min
    
        #Xspec plot
        xs.Plot.addCommand('r y '+str(ymin/3.)+' '+str(ymax))
        xs.Plot.addCommand('LA Top '+key+' '+redchisq_text)
        xs.Plot.device = '/xs'    
        xs.Plot("data delchi")
        xs.Plot.device = out_dir+prefix+model_note+str(np.where(OVII_list==OVII_Oratio)[0][0])+'_'+key+'_delchi.ps/cps' # save xspec version of plots
        xs.Plot("data delchi") 
        xs.Plot.device = '/null'
        #stop  #testing on one observation