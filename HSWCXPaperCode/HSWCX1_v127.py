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
v116a: 100% He for all lines
v116aG2: new 'good' list

v117: fitting for SWCX (including LowE) with halo pars fixed to fitted values from GOOD spectrum
      BKG_PIs and BKG_Ns all free.
v117a: using v116a halo pars (100% He ratios for all, almost identical)
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

v127: loop through a range of O VII values for fake H1 spectrum (only SWCX and halo)
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
spec_dir = '/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/HSWCX1_Sim/'
code_dir='/home/ringuette/halosat_sourceanalysis/AnalysisCode/HSWCXCode/'
out_dir=spec_dir
prefix = 'HSWCX1v127' #default is SWAN-based HSWCX calcs (model 3)
HSWCXmod = '3'   #3 for SWAN-based calcs

nH_type = '_nHKip' #set type of nH used options: '_nHRT', '_nHHI4pi', '_nHSFDW', '_nHPZ'
model_note = nH_type+'_fake2'
nh = 0.248
CXB_PI, CXB_N = 1.45, 0.382
LHB_T, LHB_N = 0.097, 0.186515

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
OVII_ratio = '0.126'
OVIII_ratio = '0.549'  #choose He ratio since through He-cone
NeIX_ratio = '0.100'
LowE_ratio = '0.925'   #also 100% He
OVII_list = np.array(['1.000','0.500','1.500'],dtype=str)

key_list = ['GOODH1']#,'Oct2018b+','Nov2018b+','Oct2019','Nov2019+','Nov2019','Dec2019+','Dec2019']  #for running with SWCX fixed

for OVII_Oratio in OVII_list:
    #initialize output format
    file_header='Target,'+nH_type.split('_')[1]+',OVII_scale' #all other headers not hard-coded
    file_format='\n{0[0]},{0[1]:.6f},{0[2]:.3f}'
    csv_num=3  #next number in []
    model_note = nH_type+'OVII_'+str(np.where(OVII_list==OVII_Oratio)[0][0])
    
    #loop through observation names
    for key in key_list:
        #choose spectra that are rebinned using grppha
        spec_name='HSWCX1_'+key+'_fake2.pi'
        spec_list=[spec_dir+spec_name]
        
        #collect parameters for observation
        idx = np.where(SWCX_LU['\xef\xbb\xbfObsName']==key)[0]
        OVII_N, OVIII_N, NeIX_N = float(SWCX_LU['OVII_Mod'+HSWCXmod][idx][0])*0.035, float(SWCX_LU['OVIII_Mod'+HSWCXmod][idx][0])*0.035, \
              float(SWCX_LU['NeIX_Mod'+HSWCXmod][idx][0])*0.035
        OVII_N2 = OVII_N*float(OVII_Oratio)  #adjust predicted O VII by scaling factor
        LowE_N = OVII_N*float(LowE_ratio)
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
        xray_par, SWCX_par = {}, {}
        xmodel, smodel, bkg_spec = [], [], []
        halo_arr=[key,nh,float(OVII_Oratio)]
        
        #loop through detector spectra
        print 'Initializing '+key
        for i in range(len(spec_list)): #fit model to each spectrum linked between detectors    
            bkg_spec.append(0)
            xmodel.append(0)
            smodel.append(0)
            #tmodel.append(0)
            
            #set up files    
            xs.AllData(str(i+1)+':'+str(i+1)+' '+spec_list[i]) #assign each spectrum to diff data groups
            bkg_spec[i]=xs.AllData(i+1)    
            bkg_spec[i].multiresponse[0] = code_dir+'halosat_avgenoise_20190423.rmf'
            bkg_spec[i].multiresponse[0].arf = code_dir+'halosat_20190211.arf'
            bkg_spec[i].multiresponse[1] = code_dir+'halosat_avgenoise_20190423.rmf'
            bkg_spec[i].multiresponse[1].arf = code_dir+'halosat_20190211.arf'  
    
            #using Liu values for LHB (fixed apec), 
            #absorbed CXB (fixed to Capelluti values) and absorbed linked halo APEC
            #SWCX lines from Alpha file from Dimitra
            if (i==0):    
                xray_par[i]={1:str(nh)+',-1',\
                              2:'0.22,0.05,0.1,0.15,0.3,0.35',3:'0.3,-1',4:'0,-1',5:'1.0,0.1,0,0.1,1e20,1e24',\
                              6:'0.75,0.1,0.6,0.7,0.9,1.5',7:'0.3,-1',8:'0,-1',9:'1.0,0.1,0,0.1,1e20,1e24',\
                              10:str(CXB_PI)+',-1',11:str(CXB_N)+',-1',\
                              12:str(LHB_T)+',-1',13:'1.0,-1',14:'0.,-1',15:str(LHB_N)+',-1'}
                SWCX_par[i] = {1:'0.44339,-1',2:'0.001,-1',3:str(LowE_N)+',-1',\
                              4:'0.56325,-1',5:'0.001,-1',6:str(OVII_N2)+',-1',\
                              7:'0.67916,-1',8:'0.001,-1',9:'=SWCX:p6*'+OVII_ratio,\
                              10:'0.65310,-1',11:'0.001,-1',12:str(OVIII_N)+',-1',\
                              13:'0.80308,-1',14:'0.001,-1',15:'=SWCX:p12*'+OVIII_ratio,\
                              16:'0.90873,-1',17:'0.001,-1',18:str(NeIX_N)+',-1',\
                              19:'1.10043,-1',20:'0.001,-1',21:'=SWCX:p18*'+NeIX_ratio} 
    
            else: 
                xray_par[i]={1:'=Xray:p1',\
                              2:'=Xray:p2',3:'=Xray:p3',4:'=Xray:p4',5:'=Xray:p5',\
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
        xray_names={2:'Halo_kT1',5:'Halo_N1',6:'Halo_kT2',9:'Halo_N2'}
        SWCX_names={3:'LowE',6:'O_VII',12:'O_VIII',18:'Ne_IX'}
        
        #set Xspec model and parameters
        xray_model=xs.Model('tbabs*(apec+apec+pow)+apec','Xray',1)   #halo only
        SWCX_model=xs.Model('gauss+gauss+gauss+gauss+gauss+gauss+gauss','SWCX',2)
        for i in range(len(spec_list)):    
            xmodel[i]=xs.AllModels(i+1,'Xray')
            smodel[i]=xs.AllModels(i+1,'SWCX')
            xmodel[i].setPars(xray_par[i])
            smodel[i].setPars(SWCX_par[i])
            
            xs.AllData.ignore('0.0-0.4,7.-**')
            xs.AllModels.setEnergies('reset')  #set to use energies instead of channels?
            xs.AllModels.systematic = 0  #systematic 0
        
        # display information for all loaded spectra, perform fit
        xs.AllModels.show()
        xs.AllData.show()
        xs.Fit.nIterations = 1000 # use lots of iterations
        xs.Fit.renorm()
        xs.Fit.perform() #execute fit
        
        #run steppar on all free xray parameters
        for i, j in itertools.product(range(len(spec_list)), range(xmodel[i].nParameters)):
            if not ((xmodel[i](j+1).frozen)or(len(xmodel[i](j+1).link)>0)):
                parnum=xmodel[i].startParIndex+xmodel[i](j+1).index-1
                if (xmodel[i](j+1).index==2):  #halo kT1
                    print 'steppar: Xray', parnum
                    xs.Fit.steppar('best Xray:'+str(parnum)+' 0.15 0.3 15')
                elif (xmodel[i](j+1).index==5):  #halo N1
                    print 'steppar: Xray', parnum
                    xs.Fit.steppar('best Xray:'+str(parnum)+' 0. 1.5 15')    
                elif (xmodel[i](j+1).index==6):  #halo kT2
                    print 'steppar: Xray', parnum
                    xs.Fit.steppar('best Xray:'+str(parnum)+' 0.7 0.9 20')
                elif (xmodel[i](j+1).index==9):  #halo N2
                    print 'steppar: Xray', parnum
                    xs.Fit.steppar('best Xray:'+str(parnum)+' 0. 1.5 15')                    
    
        #save current parameter values for later comparison after error runs
        test=[xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof]
                    
        #open log file, print info to file, including error runs
        logFile = xs.Xset.openLog(out_dir+prefix+model_note+'_'+key+'_xspeclog.txt')
        xs.Xset.show()
        xs.AllData.show()  
        if (xs.Fit.statistic/xs.Fit.dof<2.):
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
            logFile = xs.Xset.openLog(out_dir+prefix+model_note+'_'+key+'_xspeclog'+str(error_loop)+'.txt')
            xs.Xset.show()
            xs.AllData.show()  
            if (xs.Fit.statistic/xs.Fit.dof<2.):
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
                
        #save Xray output to array
        for i, j in itertools.product(range(len(spec_list)), range(xmodel[i].nParameters)):
            if not ((xmodel[i](j+1).frozen)or(len(xmodel[i](j+1).link)>0)):
                idxnum = xmodel[i](j+1).index  #det# not important since all free pars are linked   
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
        xs.Plot.device = out_dir+prefix+model_note+'_'+key+'_delchi.ps/cps' # save xspec version of plots
        xs.Plot("data delchi") 
        xs.Plot.device = '/null'
        #stop  #testing on one observation