# -*- coding: utf-8 -*-
"""      
v101: No TTS, using new SWCX from Dimitra (fixed), only halo.   
v102: same as 101, add second high temp apec (absorbed)  
v103: low temp halo=0 (removed), high temp halo only.   
v104: same as 102, low T halo-0.224 keV from target 135
v105: same as 102, high T comp kT = 0.783 keV and N = 0.308, low kT and N free 
v106: same as v105, halo_kT1 = 0.224 keV
v107: one apec for hot emission, kT = 0.783 keV and N = 0.308, SWCX primary lines free

v109: halo=0, new SWCX data, only unabsorbed TTS (0.5 and 2.07 keV temps, equal Ns)
v110: halo kT and N free, same TTS as 109
v111: halo kT = 0.224 keV, halo_N free, same TTS as 109

v112: two apecs for halo fixed to measured values, SWCX primary lines free.
      halo_kT1 = 0.272 keV, halo_N1 = 0.074 (apec units)
      halo_kT2 = 0.776 keV, halo_N2 = 0.308 (apec units)
      
v113: loop through low_E Ns (0, 0.5*value, 1*value, 1.5*value, 2*value)

v114: halo temps and SWCX fixed, BKG_PIs and BKG_Ns free
v115: all halo and bkg pars free (trying something crazy)
v116: same as v115, but including LowE line.
v117: fitting for SWCX (including LowE) with halo pars fixed to fitted values from GOOD spectrum
      BKG_PIs and BKG_Ns all free.
v117a: using v116a halo pars (100% He ratios for all, almost identical)
v117b: 100% H ratios

v118: same as 117, but scanning fixed line ratios of LowE N relative to fitted O VII lines (0%, 50%, 100%, 150%, 200%)
v119: same as 117, but fixing line ratio of LowE N relative to calculated O VII lines (100% only)
v120: excess near 0.6 keV might be lower kT halo component. trying kT3 = 0.21 (from MSWCX average)
v121: try as a CME, unabsorbed with kT > 1.5 keV

v122: same as 116a, using new SWAN-based helio calcs.
v123: same as v117b, using halo results from v122
v124: fix kT1=0.2, fit for halo
v125: fix kT1=0.2, fit for halo and SWCX primary lines.
v125b: same, kT=0.21
v126: same as v122, but with N2=0 (to show residuals)

v127: loop through a range of OVII values for fake H1 spectrum (only SWCX and halo)
v128: like v123, but for simulation spectra, v127 halo pars
v129: same as v127, all SWCX scales with OVII scaling factor
v130: same as v128, using halo pars from v129

v131: using real ECL high stat spec, loop through scaling factors until residuals near 0.6 keV are improved (SWCX too high?)
v132: same as v131, but using Model 1
v133: same as v132, using Model 3, N2=0 for comparison
v134: same as v123, using halo results from best fit scaling factor (v131_0h)
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
import itertools
import astropy.io.fits as pyfits
from csv import DictReader
#pyxspec only works with python 2.7 commands since that is the version installed
#pyxspec also requires starting spyder from the terminal
#xspec won't output fitting details to screen if I print to the screen (sometimes)

xs.Xset.addModelString("APECROOT","/home/ringuette/Software/heasoft-6.25/spectral/modelData/apec_v3.0.9_201") #use new apec file

clean_type='cl1'

#define directories
spec_dir = '/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/Rebinned_spectra_'+clean_type+'/'
code_dir='/home/ringuette/halosat_sourceanalysis/AnalysisCode/HSWCXCode/'
out_dir='/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/XspecOutput_'+clean_type+'_v4/'
prefix = 'HSWCX1v134'

#define constants
halo_kT1, halo_N1 = 0.255, 0.79  
halo_kT2, halo_N2 = 1.01, 0.239  
CXB_PI, CXB_N = 1.45, 0.382
nH_type = '_nHKip' #set type of nH used options: '_nHRT', '_nHHI4pi', '_nHSFDW', '_nHPZ'
model_note = nH_type+'OVII_0h'  #store choices

#initialize output format
file_header='Target,LHB_N,'+nH_type.split('_')[1] #all other headers not hard-coded
file_format='\n{0[0]},{0[1]:.6f},{0[2]:.6f}'
csv_num=3  #next number in []

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
elif nH_type=='_nHKip': nh = np.repeat(0.248,len(hsIDs))

#SWCX line ratios
OVII_ratio = '0.126'
OVIII_ratio = '0.549'  #choose He ratio since through He-cone  (100% He - 0.549, 100% H - 0.320)
NeIX_ratio = '0.100'

#define dictionary of fixed LHB norm, nH from LAB, and last directory
#keep hard-coded becase Oct2018+ and Nov2018+ combine two different fields,
#and are thus weighted by the exposure times (total of dets 38 and 54)
'''             
#record exposure times to calculate LHB norms and nH for Oct2018+ and Nov2018+
Oct2018+: HSID#     det14   det38   det54     Total     w/o det14
          367       5376.0  6912.5  5312.0    17600.5   12224.5
          368       14144.0 14406.0 14528.0   43078.0   28934.0

Nov2018+: HSID#     det14   det38   det54     Total     w/o det14   
          369       9600.0  9798.5  9863.5    29262.0   19662.0
          371       11584.0 11788.0 11840.5   35212.5   23628.5
          
Since this analysis ignores det 14, the weights are the values in the last column
index numbers count from 0, target numbers count from 1
''' 
fixed_pars={'Oct2018+':[(LHB[366]*12224.5+LHB[367]*28934.0)/(12224.5+28934.0),\
                        (nh[366]*12224.5+nh[367]*28934.0)/(12224.5+28934.0)],
            'Oct2018a+':[LHB[366],nh[366]],'Oct2018b+':[LHB[367],nh[367]],\
            'Nov2018+':[(LHB[368]*19662.0+LHB[370]*23628.5)/(19662.0+23628.5),\
                        (nh[368]*19662.0+nh[370]*23628.5)/(19662.0+23628.5)],\
            'Nov2018a+':[LHB[368],nh[368]],'Nov2018b+':[LHB[370],nh[370]],\
            'Dec2018':[LHB[65],nh[65]],'Feb2019':[LHB[65],nh[65]],\
            'Sep2019':[LHB[65],nh[65]],'Oct2019':[LHB[65],nh[65]],\
             'Nov2019+':[LHB[368],nh[368]],'Nov2019':[LHB[65],nh[65]],\
             'Dec2019+':[LHB[368],nh[368]],'Dec2019':[LHB[65],nh[65]],\
             'Jan2020+':[LHB[368],nh[368]],'Feb2020+':[LHB[368],nh[368]],\
             'Mar2020+':[LHB[368],nh[368]],'GOOD':[LHB[65],nh[65]],\
             'GOOD2':[LHB[65],nh[65]],'Mar2019':[LHB[65],nh[65]]}
key_list = ['Oct2018a+','Oct2018b+','Nov2018a+','Nov2018b+','Dec2018','Feb2019','Mar2019','Sep2019','Oct2019',\
            'Nov2019+','Nov2019','Dec2019+','Dec2019','Jan2020+','Feb2020+','Mar2020+']  #total list
#key_list = ['Oct2018b+','Nov2018b+','Oct2019','Nov2019+','Nov2019','Dec2019+','Dec2019','GOOD2']  #good list only
#calculated line ratios are only done for good list

#loop through observation names
for key in key_list:
    #choose spectra that are rebinned using grppha
    if key[-1]=='+': 
        spec_name='HSWCX1+_'+key[0:-1]+'_'+clean_type+'_d'
    else: 
        spec_name='HSWCX1_'+key+'_'+clean_type+'_d'
    spec_list=[spec_dir+spec_name+'14_rebin40.pi',spec_dir+spec_name+'38_rebin40.pi',spec_dir+spec_name+'54_rebin40.pi']
    
    #collect parameters for observation
    LHB_target, nh_target = fixed_pars[key][0], fixed_pars[key][1]
    
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
    halo_arr=[key,LHB_target,nh_target]
    
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
        bkg_spec[i].multiresponse[1].arf = code_dir+'hs_sdd_all20180701v001.arf'
        bkg_spec[i].multiresponse[2] = code_dir+'halosat_avgenoise_20190423.rmf'
        bkg_spec[i].multiresponse[2].arf = code_dir+'hs_sdd_all20180701v001.arf'  
        #bkg_spec[i].multiresponse[3] = code_dir+'halosat_avgenoise_20190423.rmf'
        #bkg_spec[i].multiresponse[3].arf = code_dir+'hs_sdd_all20180701v001.arf'   
        
        #add exposure time and bkg_pi to output
        dpu = int(spec_list[i].split('.')[0].split('_d')[1].split('_rebin40')[0]) #get dpu num
        halo_arr.append(bkg_spec[i].exposure)
        if key==key_list[0]:
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
                          6:str(halo_kT1)+',-1',7:'0.3,-1',8:'0,-1',9:str(halo_N1)+',-1',\
                          10:str(CXB_PI)+',-1',11:str(CXB_N)+',-1',
                          12:str(halo_kT2)+',-1',13:'0.3,-1',14:'0,-1',15:str(halo_N2)+',-1'}
            SWCX_par[i] = {1:'0.44339,-1',2:'0.001,-1',3:'0.035,0.01,0.,0.,0.06,10.',\
                          4:'0.56325,-1',5:'0.001,-1',6:'0.035,0.01,0.,0.,0.06,10.',\
                          7:'0.67916,-1',8:'0.001,-1',9:'=SWCX:p6*'+OVII_ratio,\
                          10:'0.65310,-1',11:'0.001,-1',12:'0.004,0.01,0.,0.,0.06,10.',\
                          13:'0.80308,-1',14:'0.001,-1',15:'=SWCX:p12*'+OVIII_ratio,\
                          16:'0.90873,-1',17:'0.001,-1',18:'0.002,0.01,0.,0.,0.06,10.',\
                          19:'1.10043,-1',20:'0.001,-1',21:'=SWCX:p18*'+NeIX_ratio} 

        else: 
            xray_par[i]={1:'=Xray:p1',2:'=Xray:p2',3:'=Xray:p3',4:'=Xray:p4',\
                          5:'=Xray:p5',\
                          6:'=Xray:p6',7:'=Xray:p7',8:'=Xray:p8',9:'=Xray:p9',\
                          10:'=Xray:p10',11:'=Xray:p11',
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
    #TTS_names={1:'TTS_nH',2:'TTS_kT',34:'TTS_N'}
    
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
            if (xmodel[i](j+1).index==6):  #halo kT
                print 'steppar: Xray', parnum
                xs.Fit.steppar('best Xray:'+str(parnum)+' 0.7 0.9 20')
            elif (xmodel[i](j+1).index==9):  #halo N
                print 'steppar: Xray', parnum
                xs.Fit.steppar('best Xray:'+str(parnum)+' 0. 1.5 15')                      
    #run steppar on all free SWCX parameters
    for i, j in itertools.product(range(len(spec_list)), range(smodel[i].nParameters)):
        if not ((smodel[i](j+1).frozen)or(len(smodel[i](j+1).link)>0)):
            parnum=smodel[i].startParIndex+smodel[i](j+1).index-1
            print 'steppar: SWCX', parnum    #line Ns
            xs.Fit.steppar('best SWCX:'+str(parnum)+' 0. 0.1 10')   #should be less than 2 LU
    '''#run steppar on all free TTS parameters
    for i, j in itertools.product(range(len(spec_list)), range(tmodel[i].nParameters)):
        if not ((tmodel[i](j+1).frozen)or(len(tmodel[i](j+1).link)>0)):
            parnum=tmodel[i].startParIndex+tmodel[i](j+1).index-1
            if (tmodel[i](j+1).index==1):
                print 'steppar: TTS', parnum
                xs.Fit.steppar('best TTS:'+str(parnum)+' 0.002 0.248 25')  #abs
            elif (tmodel[i](j+1).index==2):
                print 'steppar: TTS', parnum
                xs.Fit.steppar('best TTS:'+str(parnum)+' 0.5 5. 25')  #hot temp
            elif (tmodel[i](j+1).index==34):
                print 'steppar: TTS', parnum
                xs.Fit.steppar('best TTS:'+str(parnum)+' 0.0001 0.5 25')  #hot norm 
    '''
    #save current parameter values for later comparison after error runs
    test=[xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof]
                
    #open log file, print info to file, including error runs
    logFile = xs.Xset.openLog(out_dir+prefix+model_note+'_'+key+'_xspeclog.txt')
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
        '''#run error on all free TTS parameters
        for i, j in itertools.product(range(len(spec_list)), range(tmodel[i].nParameters)):
            if not ((tmodel[i](j+1).frozen)or(len(tmodel[i](j+1).link)>0)):
                parnum=tmodel[i].startParIndex+tmodel[i](j+1).index-1
                print 'error: TTS', parnum
                xs.Fit.error('TTS:'+str(parnum))                 
        '''
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
        logFile = xs.Xset.openLog(out_dir+prefix+model_note+'_'+key+'_xspeclog'+str(error_loop)+'.txt')
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
            '''#run error on all free TTS parameters
            for i, j in itertools.product(range(len(spec_list)), range(tmodel[i].nParameters)):
                if not ((tmodel[i](j+1).frozen)or(len(tmodel[i](j+1).link)>0)):
                    parnum=tmodel[i].startParIndex+tmodel[i](j+1).index-1
                    print 'error: TTS', parnum
                    xs.Fit.error('TTS:'+str(parnum))                         
            '''
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
            if key==key_list[0]:
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
            if key==key_list[0]:  #first key
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
            if key==key_list[0]:  #first key
                file_header+=','+SWCX_names[idxnum]+',ErrL,ErrU,ErrCode'   
                file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']:.6f},{0['\
                    +str(csv_num+2)+']:.6f},{0['+str(csv_num+3)+']}'
                csv_num+=4   

    '''#save TTS output to array
    for i, j in itertools.product(range(len(spec_list)), range(tmodel[i].nParameters)):
        if not ((tmodel[i](j+1).frozen)or(len(tmodel[i](j+1).link)>0)):
            idxnum = tmodel[i](j+1).index  #det# not important since all free pars are linked   
            halo_arr.extend([tmodel[i](idxnum).values[0],tmodel[i](idxnum).values[0]-tmodel[i](idxnum).error[0],\
                             tmodel[i](idxnum).error[1]-tmodel[i](idxnum).values[0],tmodel[i](idxnum).error[2]])  
            if key==key_list[0]:  #first key
                file_header+=','+TTS_names[idxnum]+',ErrL,ErrU,ErrCode'   
                file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']:.6f},{0['\
                    +str(csv_num+2)+']:.6f},{0['+str(csv_num+3)+']}'
                csv_num+=4   
    '''
    #add fluxes to output
    xs.AllModels.calcFlux('0.4 7.0')
    for i in range(len(spec_list)):
        halo_arr.extend([bkg_spec[i].flux[6], bkg_spec[i].flux[9]])
        if key==key_list[0]:  #first key
            det = spec_list[i].split('/')[-1].split('_d')[1].split('_')[0]
            file_header+=',d'+str(det)+'_Flux(erg),d'+str(det)+'_Flux(photons)'   
            file_format+=',{0['+str(csv_num)+']:.15f},{0['+str(csv_num+1)+']:.6f}'
            csv_num+=2
            
    halo_arr.extend([xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof])  #add fit statistics to list
    if key==key_list[0]: #initialize csv file header
        file_header+=',X^2,dof,red_X^2'
        file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']},{0['\
            +str(csv_num+2)+']:.6f}'   
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
    xs.Plot.add = False
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
    xs.Plot.device = out_dir+prefix+model_note+'_'+key+'_delchi.ps/cps' # save xspec version of plots
    xs.Plot("data delchi") 
    xs.Plot.device = '/null'
    #stop  #testing on one observation