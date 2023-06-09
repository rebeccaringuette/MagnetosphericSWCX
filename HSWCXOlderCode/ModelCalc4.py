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
    - Checking for new fits found during error process doesn't work with error code,
      possibly because query=yes. Comparing fitted values at beginning vs end instead. 
    - Also reports fit statistics for 0.4 to 0.55 and 2 - 7 keV to judge bkg performance
    
Halo_FullERange.py
    - assign BKG PI to fixed value for halo fields, diff output names
    - Checking for new fits found during error process doesn't work with error code,
      possibly because query=yes. Comparing fitted values at beginning vs end instead.
    - this method results in looping even if now new fit is found.
      
Halo_FullERange_v1.py      
    - alternate method for error loop if new fit found.
    - does not loop if no new fit is found, does if one is

BKG_FullERange_v1.py:
    - applying new method for error loop to this program
    
BKG_FullERange_v2.py:
    - selecting by galactic latitude and exposure time    
    - skipping targets with less then 15ks total exposure time
    - file format now hardcoded

ModelCalc.py
    - calculate the flux in photons/cm^2/s and ergs/cm^s/s of an apec component with fixed parameters
    - output results in a csv file
    
ModelCalc3.py: double absorbed apecs (based on Gudel 2007 and Masui 2009)

ModelCalc4.py: calculate flux in 0.4-7 and 0.3-10 keV ranges for fitted TTS components
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
xs.Xset.addModelString('APECROOT', '/home/ringuette/Software/heasoft-6.25/spectral/modelData/apec_v3.0.9_201')
#pyxspec only works with python 2.7 commands since that is the version installed
#pyxspec also requires starting spyder from the terminal
#xspec won't output fitting details to screen if I print to the screen (sometimes)

#calc using longest exposure tile, only matters if using calculated LHB numbers
spec_dir = '/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/Rebinned_spectra_cl1/'
code_dir='/home/ringuette/halosat_sourceanalysis/AnalysisCode/'
out_dir ='/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/'
spec_file = spec_dir+'HSWCX1+_Dec2019_cl1_d14_rebin40.pi'  #Dec2019+ HSWCX1 observation, target 369
if not os.path.exists(out_dir): os.makedirs(out_dir)
file_root = 'ModelCalc4'

#energies for flux calculation
LowE1 = 0.4
HighE1 = 7.0
LowE2 = 0.3
HighE2 = 10.0

#initialize output format
file_header='ObsName,kt1 (keV),kt2 (keV),N (apec units),ErrL,ErrU' #0-5
file_header+=',Low E1 (keV),High E1 (keV),Flux1 (ergs/cm^2/s),ErrL,ErrU'  #6-10
file_header+=',Low E2 (keV),High E2 (keV),Flux2 (ergs/cm^2/s),ErrL,ErrU'  #11-15
file_format='\n{0[0]},{0[1]:.3f},{0[2]:.3f},{0[3]:.6f},{0[4]:.6f},{0[5]:.6f}'
file_format+=',{0[6]:.1f},{0[7]:.1f},{0[8]:.16f},{0[9]:.16f},{0[10]:.16f}'
file_format+=',{0[11]:.1f},{0[12]:.1f},{0[13]:.16f},{0[14]:.16f},{0[15]:.16f}'
chi_file=open(out_dir+file_root+'.csv','w') #initialize output file
chi_file.write(file_header)
chi_file.close()   

#dictionary of TTS_N values
TTS_Nd={}  #TTS_N,ErrL,ErrU 
TTS_Nd['Oct2018b+']=np.array([0.043967,0.016364,0.016364],dtype=float)
TTS_Nd['Nov2018+']=np.array([0.053897,0.014108,0.014108],dtype=float)
TTS_Nd['Oct2019']=np.array([0.04443,0.010617,0.010617],dtype=float)
TTS_Nd['Nov2019+']=np.array([0.068387,0.018573,0.018573],dtype=float)
TTS_Nd['Nov2019']=np.array([0.044666,0.009737,0.009737],dtype=float)
TTS_Nd['Dec2019+']=np.array([0.049949,0.006981,0.006981],dtype=float)
TTS_Nd['Dec2019']=np.array([0.054988,0.01527,0.01527],dtype=float)
TTS_T1, TTS_T2 = 0.808, 2.07   #keV

#setup xspec
xs.AllData.clear()
xs.AllModels.clear()
xs.Xset.abund='wilm' #abund wilm
xs.Xset.xsect='vern' #xsect vern
xs.Xset.cosmo='70 0 0.73'  #cosmo 70 0 0.73

#set up files    
xs.AllData('1:1 '+spec_file) #assign each spectrum to diff data groups
bkg_spec=xs.AllData(1)    
    
#loop through parameter sets
for key in TTS_Nd: 
    n1, ErrL, ErrU = TTS_Nd[key]
    out_arr = [key,TTS_T1,TTS_T2,n1,ErrL,ErrU,LowE1,HighE1,0.,0.,0.,LowE2,HighE2,0.,0.,0.]
    xs.AllModels.clear()  #clear previous model 
    
    #set up model for value
    xray_par={1:str('{:.3f}'.format(TTS_T1))+',-1',2:'1,-1',3:'0,-1',4:str('{:.3f}'.format(n1))+',-1',\
            5:str('{:.3f}'.format(TTS_T2))+',-1',6:'1,-1',7:'0,-1',8:str('{:.3f}'.format(n1))+',-1'}
    xray_model=xs.Model('apec+apec','Xray',1)   #TTS component
    xmodel=xs.AllModels(1,'Xray')
    xmodel.setPars(xray_par)
    xs.AllModels.setEnergies('reset')  #set to use energies instead of channels?
    xs.AllModels.systematic = 0  #systematic 0
    #store fluxes              
    xs.AllModels.calcFlux(str(LowE1)+' '+str(HighE1))
    out_arr[8] = xmodel.flux[0]
    xs.AllModels.calcFlux(str(LowE2)+' '+str(HighE2))
    out_arr[13] = xmodel.flux[0]      
    
    #set up model for ErrL
    xray_par={1:str('{:.3f}'.format(TTS_T1))+',-1',2:'1,-1',3:'0,-1',4:str('{:.3f}'.format(n1-ErrL))+',-1',\
            5:str('{:.3f}'.format(TTS_T2))+',-1',6:'1,-1',7:'0,-1',8:str('{:.3f}'.format(n1-ErrL))+',-1'}
    xray_model=xs.Model('apec+apec','Xray',1)   #TTS component
    xmodel=xs.AllModels(1,'Xray')
    xmodel.setPars(xray_par)
    xs.AllModels.setEnergies('reset')  #set to use energies instead of channels?
    xs.AllModels.systematic = 0  #systematic 0
    #store fluxes              
    xs.AllModels.calcFlux(str(LowE1)+' '+str(HighE1))
    out_arr[9] = out_arr[8]-xmodel.flux[0]
    xs.AllModels.calcFlux(str(LowE2)+' '+str(HighE2))
    out_arr[14] = out_arr[13]-xmodel.flux[0]       
    
    #set up model for ErrU
    xray_par={1:str('{:.3f}'.format(TTS_T1))+',-1',2:'1,-1',3:'0,-1',4:str('{:.3f}'.format(n1+ErrU))+',-1',\
            5:str('{:.3f}'.format(TTS_T2))+',-1',6:'1,-1',7:'0,-1',8:str('{:.3f}'.format(n1+ErrU))+',-1'}
    xray_model=xs.Model('apec+apec','Xray',1)   #TTS component
    xmodel=xs.AllModels(1,'Xray')
    xmodel.setPars(xray_par)
    xs.AllModels.setEnergies('reset')  #set to use energies instead of channels?
    xs.AllModels.systematic = 0  #systematic 0
    #store fluxes              
    xs.AllModels.calcFlux(str(LowE1)+' '+str(HighE1))
    out_arr[10] = xmodel.flux[0]-out_arr[8]
    xs.AllModels.calcFlux(str(LowE2)+' '+str(HighE2))
    out_arr[15] = xmodel.flux[0]-out_arr[13]      
    #print out_arr

    #output results to csv file spec_dir+run_type
    with open(out_dir+file_root+'.csv','a') as chi_file: 
        chi_file.write(file_format.format(out_arr))
    #stop
 