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

TailFit_v1.py
    - fit all MSWCX flank spectra
    - fixed CXB, LHB, nH
    - SWCX O7, halo kt, halo N, bkg N all free
    - bkg pi fixed to calculated values from SH 

v2: updated MgXIb/MgXIa ratio, using file now
v3: using only highest ratio lines (NVI, CVI, OVII, OVIII, NeIX); not MgXI b/c ratio looks way too high
v4: use only detected lines (OVII, OVIII, NeIX)
v5: use only OVII and OVIII as done for halo
v6: all except MgXI
v8: change halo metallicity to 0.3, remove Mg XI, use one line for lower SWCX (free)
v8b: LowE_N=0
v9: all SWCX=0, halo metallicity=0.3
v10: fix SWCX LU to values from Dimitra, adjusted for exposure time. secondary lines included, 100% H line ratios
v10b: 100% He line ratios
v10a: updated SWCX calculations, 100% H ratios
v10c: updated SWCX calculations, 100% H ratios, restricted halo temp to 0.15-0.35
v11: no secondary lines, same as v10 otherwise.
v12: same as v10a, BKG_PIs free, 50% H, 50% He ratios
v13: only clean tail observations combined, one each available in spec_dir

v14: working with new HSWCX calculations from Dimitra, tail files filtered by net hard ctrt
     - run through all combined MSWCX tail spectra with no HSWCX scaling factor
v15: scan through SWCX scaling factors, diff targets in diff files
v16: same scan as v15, with LowE line free

v17: same procedure as HSWCX paper, mod2 scaling=1, single apec
Hot: add 2nd apec for all
v18: new SWCX file
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
#xspec won't output fitting details to screen if I print to the screen (sometimes)
#xs.Xset.addModelString("APECROOT","/home/ringuette/Software/heasoft-6.25/spectral/modelData/apec_v3.0.9_201") #use new apec file


#for folder, choose tail spectra from folder 
source_list=np.array(['HS0068','HS0069','HS0070','HS0071','HS0072','HS0073','HS0074','HS0075','HS0076','HS0077','HS0078','HS0079'],dtype=str) 
# need HSWCX calcs for: 'MSWCX1','MSWCX2','MSWCX3','MSWCX11','MSWCX12'    
# mod2 done: 'MSWCX4','MSWCX5','MSWCX7','MSWCX8','MSWCX9','MSWCX10'
spec_dir='/home/ringuette/MSWCX/RebinnedFiles/'
code_dir='/home/ringuette/MSWCX/MSWCXCode/'
out_dir ='/home/ringuette/MSWCX/Results/'
if not os.path.exists(out_dir): os.makedirs(out_dir)
nH_type = '_nHRT' #set type of nH used options: '_nHHI4pi', '_nHSFDW', '_nHPZ', '_nHRT'
model_note = 'MSWCXTailsHot_v18'+nH_type+'_mod2'  #store choices, fit type
 
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

#read helioSWCX calcs per spec
SWCX_LU = {}
SWCX_file = DictReader(open(code_dir+'NewHelioSWCX.csv','r'))
for row in SWCX_file:
    for column, value in row.iteritems():
        SWCX_LU.setdefault(column, []).append(value)
for key in SWCX_LU: 
  SWCX_LU[key] = np.array(SWCX_LU[key])
  #print key, SWCX_LU[key][3], type(SWCX_LU[key][3]) 
#print(SWCX_LU.keys())
    
#choose SWCX line ratios
OVII_ratio = '0.126'  #same for 100% H and 100% He
OVIII_ratio = '0.4345'  #100% H -> 0.320, 100% He -> 0.549, 50/50 of each -> 0.4345
NeIX_ratio = '0.100'  #same for 100% H and 100% He
LowE_ratio = '0.946'   #100% H = 0.966, 100% He = 0.925, 50% each = 0.946, not sensitive

#initialize output format
file_header='Target,RA,DEC,GLon,Glat,LHB_N,'+nH_type.split('_')[1] #all other headers not hard-coded
file_format='\n{0[0]},{0[1]:.6f},{0[2]:.6f},{0[3]:.6f},{0[4]:.6f},{0[5]:.6f},{0[6]:.6f}'
csv_num=7  #next number in []
#print source_name, SWCX_scale[source_name]

#loop through MSWCX#s (one tail spec per target)
for source_name in source_list:
    #choose nH and LHB EM for target
	source_idx = np.where(hsIDs==int(source_name[2:]))[0]  #'HS00XX -> Target_ID position
	nh_target, LHB_target = str(nh[source_idx][0]), str(LHB[source_idx][0])  #nh and LHB_EM
	spec_list=np.sort(glob.glob(spec_dir+source_name+'_tail_cl1_d*_rebin40.pi'))      #get file names
	obsname = source_name+'_tail'  #set obsname

	#retrieve SWCX_LU values for current combined tail obs   
	idx = np.where(SWCX_LU['\xef\xbb\xbfObsName']==source_name+'_tail')[0]
	if SWCX_LU['OVII_2'][idx][0]=='': continue   #skip if no SWCX calculations available
	OVII_N = float(SWCX_LU['OVII_2'][idx][0])*0.035
	OVIII_N = float(SWCX_LU['OVIII_2'][idx][0])*0.035
	NeIX_N = float(SWCX_LU['NeIX_2'][idx][0])*0.035  
	#print spec_list[0].split('/')[-1], OVII_N/0.035, OVIII_N/0.035, NeIX_N/0.035

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
	bkg_par, xray_par, SWCX_par = {}, {}, {}
	bmodel, xmodel, smodel, bkg_spec = [], [], [], []
	
	#store values
	halo_arr=[source_name,RA[source_idx][0],DEC[source_idx][0],glon[source_idx][0],glat[source_idx][0],LHB[source_idx][0],nh[source_idx][0]]
	
	#loop through detector spectra
	print 'Initializing '+source_name
	for i in range(len(spec_list)): #fit model to each spectrum linked between detectors        
		#add placeholders in lists
		bkg_spec.append(0)  #spectrum objects
		bmodel.append(0)    #background models
		xmodel.append(0)    #xray models
		smodel.append(0)    #swcx models
		
		#set up files    
		xs.AllData(str(i+1)+':'+str(i+1)+' '+spec_list[i]) #assign each spectrum to diff data groups
		bkg_spec[i]=xs.AllData(i+1)    
		bkg_spec[i].multiresponse[0] = code_dir+'halosat_diag_20190423.rmf'   #only rmf for particle bkg fit
		bkg_spec[i].multiresponse[1] = code_dir+'halosat_avgenoise_20190423.rmf'
		bkg_spec[i].multiresponse[1].arf = code_dir+'hs_sdd_all20180701v001.arf'
		bkg_spec[i].multiresponse[2] = code_dir+'halosat_avgenoise_20190423.rmf'
		bkg_spec[i].multiresponse[2].arf = code_dir+'hs_sdd_all20180701v001.arf' 
		
		#save exposure time to array
		halo_arr.append(bkg_spec[i].exposure/1000.)  #assume same order for each halo target
		name_arr = spec_list[i].split('/')[-1].split('_')
		dpu = int(name_arr[3].split('d')[1])
		if source_name==source_list[0]:
			file_header+=',d'+str(dpu)+'_Exp (ks)' #all exposure time to header
			file_format+=',{0['+str(csv_num)+']:.3f}'
			csv_num+=1        
			
		#set up model
		bkg_par[i]={1:'0.75,0.2,0.4,0.6,1.0,1.3',2:'0.02,0.01,0.,0.,1e6,1e6'}  #BKG PI free and N free per det
		if (i==0):    
			xray_par[i]={1:'0.097,-1',2:'1,-1',3:'0,-1',4:LHB_target+',-1',\
						  5:nh_target+',-1',\
						  6:'0.22,0.1,0.1,0.15,0.3,0.35',7:'0.3,-1',8:'0,-1',9:'1.0,0.1,0,0.1,1e20,1e24',\
						  10:'1.45,-1',11:'0.382,-1',
						  12:'0.75,0.1,0.6,0.7,0.9,1.5',13:'0.3,-1',14:'0,-1',15:'1.0,0.1,0,0.1,1e20,1e24'}
			SWCX_par[i] = {1:'0.44339,-1',2:'0.001,-1',3:'=SWCX:p6*'+LowE_ratio,\
						  4:'0.56325,-1',5:'0.001,-1',6:str(OVII_N)+',-1',\
						  7:'0.67916,-1',8:'0.001,-1',9:'=SWCX:p6*'+OVII_ratio,\
						  10:'0.65310,-1',11:'0.001,-1',12:str(OVIII_N)+',-1',\
						  13:'0.80308,-1',14:'0.001,-1',15:'=SWCX:p12*'+OVIII_ratio,\
						  16:'0.90873,-1',17:'0.001,-1',18:str(NeIX_N)+',-1',\
						  19:'1.10043,-1',20:'0.001,-1',21:'=SWCX:p18*'+NeIX_ratio} 
		else: 
			xray_par[i]={1:'0.097,-1',2:'1,-1',3:'0,-1',4:LHB_target+',-1',\
						  5:nh_target+',-1',\
						  6:'=Xray:p6',7:'=Xray:p7',8:'0,-1',9:'=Xray:p9',\
						  10:'1.45,-1',11:'0.382,-1',\
						  12:'=Xray:p12',13:'=Xray:p13',14:'=Xray:p14',15:'=Xray:p15'}      
			SWCX_par[i] = {1:'0.44339,-1',2:'0.001,-1',3:'=SWCX:p3',\
						  4:'0.56325,-1',5:'0.001,-1',6:'=SWCX:p6',\
						  7:'0.67916,-1',8:'0.001,-1',9:'=SWCX:p9',\
						  10:'0.65310,-1',11:'0.001,-1',12:'=SWCX:p12',\
						  13:'0.80308,-1',14:'0.001,-1',15:'=SWCX:p15',\
						  16:'0.90873,-1',17:'0.001,-1',18:'=SWCX:p18',\
						  19:'1.10043,-1',20:'0.001,-1',21:'=SWCX:p21'} 
	#set names for parameters
	xray_names={6:'Halo_kT1',9:'Halo_N1',12:'Halo_kT2',15:'Halo_N2'}
	SWCX_names={3:'LowE',6:'O_VII',12:'O_VIII',18:'Ne_IX'}    
	
	#set models and parameters
	bkg_model=xs.Model('pow','BKG',1)  #particle background
	xray_model=xs.Model('apec+tbabs*(apec+pow+apec)','Xray',2)   #LHB+absorbed halo+absorbed CXB 
	SWCX_model=xs.Model('gauss+gauss+gauss+gauss+gauss+gauss+gauss','SWCX',3)
	for i in range(len(spec_list)):    
		bmodel[i]=xs.AllModels(i+1,'BKG')
		xmodel[i]=xs.AllModels(i+1,'Xray')
		smodel[i]=xs.AllModels(i+1,'SWCX')
		bmodel[i].setPars(bkg_par[i])
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
	
	#run steppar on all free bkg parameters
	for i, j in itertools.product(range(len(spec_list)), range(bmodel[i].nParameters)):
		if not ((bmodel[i](j+1).frozen)or(len(bmodel[i](j+1).link)>0)):
			parnum=bmodel[i].startParIndex+bmodel[i](j+1).index-1
			if (bmodel[i](j+1).index==1):  #BKG PI
					print 'steppar: BKG', parnum
					xs.Fit.steppar('best BKG:'+str(parnum)+' 0.4 1.3 10')
			if (bmodel[i](j+1).index==2):  #BKG N
				print 'steppar: BKG', parnum
				xs.Fit.steppar('best BKG:'+str(parnum)+' 0. 0.05 15')                
	#run steppar on all free xray parameters
	for i, j in itertools.product(range(len(spec_list)), range(xmodel[i].nParameters)):
		if not ((xmodel[i](j+1).frozen)or(len(xmodel[i](j+1).link)>0)):
			parnum=xmodel[i].startParIndex+xmodel[i](j+1).index-1
			if (xmodel[i](j+1).index==6):  #halo kT1
				print 'steppar: Xray', parnum
				xs.Fit.steppar('best Xray:'+str(parnum)+' 0.15 .3 15')
			if (xmodel[i](j+1).index==9) or (xmodel[i](j+1).index==15):  #halo N
				print 'steppar: Xray', parnum
				xs.Fit.steppar('best Xray:'+str(parnum)+' 0. 2. 20')  
			if (xmodel[i](j+1).index==12):  #halo kT2
				print 'steppar: Xray', parnum
				xs.Fit.steppar('best Xray:'+str(parnum)+' 0.7 0.9 20')				
	#run steppar on all free SWCX parameters
	for i, j in itertools.product(range(len(spec_list)), range(smodel[i].nParameters)):
		if not ((smodel[i](j+1).frozen)or(len(smodel[i](j+1).link)>0)):
			parnum=smodel[i].startParIndex+smodel[i](j+1).index-1
			print 'steppar: SWCX', parnum
			xs.Fit.steppar('best SWCX:'+str(parnum)+' 0.0 0.2 20')                          

	#save current parameter values for later comparison after error runs
	test=[xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof]
				
	#open log file, print info to file, including error runs
	logFile = xs.Xset.openLog(out_dir+model_note+'_'+source_name+'_xspeclog.txt')
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
	
	#save current parameter values for comparison with before error runs
	teste, error_check, error_loop = [xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof], 0, 0
			
	#if new parameters values differ by more than 1%, rerun errors
	if max(abs(np.array(test)-np.array(teste))/np.array(test))>0.001: error_check+=1 
	while error_check>0:
		#redo errors with log file
		print 'New fit found. Rerunning error process.'
		error_loop+=1
		logFile = xs.Xset.openLog(out_dir+model_note+'_'+source_name+'_xspeclog'+str(error_loop)+'.txt')
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
			det = spec_list[i].split('_d')[-1].split('_')[0]
			halo_arr.extend([bmodel[i](idxnum).values[0],bmodel[i](idxnum).values[0]-bmodel[i](idxnum).error[0],\
							 bmodel[i](idxnum).error[1]-bmodel[i](idxnum).values[0],bmodel[i](idxnum).error[2]])  #save par values, LL, UL, and code
			if source_name==source_list[0]:
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
			halo_arr.extend([xmodel[i](idxnum).values[0],xmodel[i](idxnum).values[0]-xmodel[i](idxnum).error[0],\
							 xmodel[i](idxnum).error[1]-xmodel[i](idxnum).values[0],xmodel[i](idxnum).error[2]])  #Line N
			if source_name==source_list[0]:
				if xmodel[i](idxnum).name=='kT': xname='kT'
				elif xmodel[i](idxnum).name=='norm': xname='N'
				file_header+=',halo_'+xname+',ErrL,ErrU,ErrCode'   
				file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']:.6f},{0['\
					+str(csv_num+2)+']:.6f},{0['+str(csv_num+3)+']}'
				csv_num+=4                 

	#save SWCX output to array
	for i, j in itertools.product(range(len(spec_list)), range(smodel[i].nParameters)):
		if not ((smodel[i](j+1).frozen)or(len(smodel[i](j+1).link)>0)):
			idxnum = smodel[i](j+1).index  #det# not important since all free pars are linked   
			halo_arr.extend([smodel[i](idxnum).values[0],smodel[i](idxnum).values[0]-smodel[i](idxnum).error[0],\
							 smodel[i](idxnum).error[1]-smodel[i](idxnum).values[0],smodel[i](idxnum).error[2]])  
			if source_name==source_list[0]:  #first key
				file_header+=','+SWCX_names[idxnum]+',ErrL,ErrU,ErrCode'   
				file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']:.6f},{0['\
					+str(csv_num+2)+']:.6f},{0['+str(csv_num+3)+']}'
				csv_num+=4 
	
	halo_arr.extend([xs.Fit.statistic,xs.Fit.dof,xs.Fit.statistic/xs.Fit.dof])  #add fit statistics to list
	if source_name==source_list[0]: #initialize csv file header
		file_header+=',X^2,dof,red_X^2'
		file_format+=',{0['+str(csv_num)+']:.6f},{0['+str(csv_num+1)+']},{0['\
			+str(csv_num+2)+']:.6f}'        
		chi_file=open(out_dir+model_note+'.csv','w') #initialize output file
		chi_file.write(file_header)
		chi_file.close()        
	
	#output results to csv file spec_dir
	with open(out_dir+model_note+'.csv','a') as chi_file: 
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
	xs.Plot.add = False  #only want net models plotting, not individual pieces
	xs.Plot.xAxis = 'keV'
	xs.Plot("data delchi")   #change for each error type
	redchisq_text=str('Chi^2/dof = {:.3f}'.format(xs.Fit.statistic/xs.Fit.dof))  #add reduced X^2 to plot
	
	#collect data
	xg, rates, xErrs, yErrs, x, y = {}, {}, {}, {}, {}, {}
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
	xs.Plot.addCommand('r y '+str(ymin/5.)+' '+str(ymax))
	xs.Plot.addCommand('LA Top '+source_name+'_Tails '+redchisq_text)
	xs.Plot.device = '/xs'    
	xs.Plot("data delchi")
	xs.Plot.device = out_dir+model_note+'_'+source_name+'_delchi.ps/cps' # save xspec version of plots
	xs.Plot("data delchi") 
	xs.Plot.device = '/null'
	#stop
	