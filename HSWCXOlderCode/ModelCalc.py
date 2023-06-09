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
nH_type='_nHSFDW' #set type of nH used options: '_nHLAB', '_nHHI4pi', '_nHSFDW', '_nHPZ'
file_root = 'ModelCalc'+nH_type

#energies for flux calculation
LowE1 = 0.56341 - 0.04
HighE1 = 0.56341 + 0.04
LowE2 = 0.65310 - 0.04
HighE2 = 0.65310 + 0.04
LowE3 = 0.90881 - 0.04
HighE3 = 0.90881 + 0.04

#initialize output format
file_header='nH (10^22),kt (keV),N (apec units),EM (pc/cm^6),Total Flux (ph/cm2/s/sr),Total Flux (ergs/cm2/s)'  #0-5
file_header+=',Low E1 (keV),High E1 (keV),Flux1 (ph/cm2/s/sr)'  #6-8
file_header+=',Low E2 (keV),High E2 (keV),Flux2 (ph/cm2/s/sr)'  #9-11
file_header+=',Low E3 (keV),High E3 (keV),Flux3 (ph/cm2/s/sr)'  #12-14
file_format='\n{0[0]:.5f},{0[1]:.5f},{0[2]:.5f},{0[3]:.10f},{0[4]:.6f},{0[5]:.16f}'
file_format+=',{0[6]:.6f},{0[7]:.6f},{0[8]:.16f}'
file_format+=',{0[9]:.6f},{0[10]:.6f},{0[11]:.16f}'
file_format+=',{0[12]:.6f},{0[13]:.6f},{0[14]:.16f}'
chi_file=open(out_dir+file_root+'.csv','w') #initialize output file
chi_file.write(file_header)
chi_file.close()   

#pull LHB N from fits file, store in dict 
t = pyfits.open(code_dir+'LHB_nH_ecl2.fits')  #use new file
hsIDs = t[1].data['Target_ID']
RA = t[1].data['RA']
DEC = t[1].data['Dec']
Target_Name = t[1].data['Target_Name']
glat = t[1].data['Gal_Lat']
glon = t[1].data['Gal_Lon']
LHB_kt = t[1].data['LHB_kT'] 
LHB_N = t[1].data['LHB_norm']
nh_LAB = t[1].data['nH_LAB']  #nH from LAB. Try using Daniel's total nH
nh_HI4Pi = t[1].data['nH_HI4pi']
nh_SFDW = t[1].data['nH_SFDW']
nh_PZ = t[1].data['nH_PlanckZhu']
t.close()

if nH_type=='_nHLAB': nh = nh_LAB
elif nH_type=='_nHHI4pi': nh = nh_HI4Pi
elif nH_type=='_nHSFDW': nh = nh_SFDW
elif nH_type=='_nHPZ': nh = nh_PZ

#list of nH, kt, N for apec in each run
par_list = [[0.,0.097,LHB_N[368]],\
            [nh[368],0.222,0.72],[nh[368],0.238,0.72],[nh[368],0.202,0.72],[nh[368],0.222,0.88],[nh[368],0.222,0.61],\
            [nh[368],0.222,0.48],[nh[368],0.238,0.48],[nh[368],0.202,0.48],[nh[368],0.222,0.587],[nh[368],0.222,0.407],\
            [nh[368],0.222,0.36],[nh[368],0.238,0.36],[nh[368],0.202,0.36],[nh[368],0.222,0.44],[nh[368],0.222,0.305],\
            [nh[368],0.222,0.24],[nh[368],0.238,0.24],[nh[368],0.202,0.24],[nh[368],0.222,0.293],[nh[368],0.222,0.203],\
            [nh[368],0.222,0.18],[nh[368],0.238,0.18],[nh[368],0.202,0.18],[nh[368],0.222,0.22],[nh[368],0.222,0.1525],\
            [nh[368],0.222,0.144],[nh[368],0.238,0.144],[nh[368],0.202,0.144],[nh[368],0.222,0.176],[nh[368],0.222,0.122],\
            [nh[368],0.222,0.090],[nh[368],0.238,0.090],[nh[368],0.202,0.090],[nh[368],0.222,0.110],[nh[368],0.222,0.07625],\
            [nh[368],0.222,0.072],[nh[368],0.238,0.072],[nh[368],0.202,0.072],[nh[368],0.222,0.088],[nh[368],0.222,0.061]]
#par_list = [[0.,0.1852737775,0.266424141285]]  #testing against Phil's
#setup xspec
xs.AllData.clear()
xs.AllModels.clear()
xs.Xset.abund='wilm' #abund wilm
xs.Xset.xsect='vern' #xsect vern
xs.Xset.cosmo='70 0 0.73'  #cosmo 70 0 0.73

#set up files    
xs.AllData('1:1 '+spec_file) #assign each spectrum to diff data groups
bkg_spec=xs.AllData(1)    
bkg_spec.multiresponse[0] = code_dir+'halosat_avgenoise_20190423.rmf'
bkg_spec.multiresponse[0].arf = code_dir+'halosat_20190211.arf'
bkg_spec.multiresponse[1] = code_dir+'halosat_diag_20190423.rmf'   #only diag rmf for particle bkg fit 
    
#loop through parameter sets
for trio in par_list: 
    nh, kt, N = trio
    EM = (1E14/3.0857E+18)*4*np.pi*N/0.035
    out_arr = [nh,kt,N,EM,0.,0.,LowE1,HighE1,0.,LowE2,HighE2,0.,LowE3,HighE3,0.]
    xs.AllModels.clear()  #clear previous model 
    
    #set up model
    xray_par={1:str(nh)+',-1',2:str(kt)+',-1',3:'1,-1',4:'0,-1',5:str(N)+',-1'}
    xray_model=xs.Model('phabs*apec','Xray',1)   #absorbed apec
    xmodel=xs.AllModels(1,'Xray')
    xmodel.setPars(xray_par)
    xs.AllData.ignore('0.0-0.4,7.-**')
    xs.AllModels.setEnergies('reset')  #set to use energies instead of channels?
    xs.AllModels.systematic = 0  #systematic 0
                  
    #add fluxes to output
    xs.AllModels.calcFlux('0.4 7.0')
    out_arr[4], out_arr[5] = bkg_spec.flux[3]/0.035, bkg_spec.flux[0]    
    xs.AllModels.calcFlux(str(LowE1)+' '+str(HighE1))
    out_arr[8] = bkg_spec.flux[3]/0.035  
    xs.AllModels.calcFlux(str(LowE2)+' '+str(HighE2))
    out_arr[11] = bkg_spec.flux[3]/0.035    
    xs.AllModels.calcFlux(str(LowE3)+' '+str(HighE3))
    out_arr[14] = bkg_spec.flux[3]/0.035       
    #print out_arr

    #output results to csv file spec_dir+run_type
    with open(out_dir+file_root+'.csv','a') as chi_file: 
        chi_file.write(file_format.format(out_arr))
    #stop
 