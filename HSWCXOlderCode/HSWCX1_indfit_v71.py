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
    
HSWCX1_indfit_v1.py
    - adjusted for HSWCX1 fits
    
HSWCX1_indfit_v2.py
    - oxygen norms set to zero. meant for testing oct2019 obs specifically
    - output names have 'HSWCX1v2' prefix
    
HSWCX1_indfit_v3.py
    - halo parameters fixed to kT =  and N = 
    - output names have 'HSWCX1v3' prefix    
    
HSWCX1_indfit_v3b.py
    - halo parameters fixed to kT =  and N = 
    - output names have 'HSWCX1v3b' prefix    
    - testing alternate method for error loop if new fit found
    
HSWCX1_indfit_v4.py
    - halo parameters allowed to vary with 1 and 2 sigma (soft and hard limits)
    - starting values are set by user
    - output names have 'HSWCX1v4' prefix    
    
HSWCX1_indfit_v5.py
    - halo parameters limited to region loosely defined by nearby halo field values
    - starting values are set by user
    - output names have 'HSWCX1v5' prefix

HSWCX1_indfit_v6.py
    - halo N limited to region loosely defined by nearby halo field values, kT fixed
    - starting values are set by user
    - output names have 'HSWCX1v6' prefix
    
HSWCX1_indfit_v7.py
    - halo kT and N fixed, using vapec for halo instead,O fixed at 1, others free
    - starting values are set by user
    - output names have 'HSWCX1v7' prefix    
    
HSWCX1_indfit_v8.py
    - halo kT and N fixed, using vapec for halo instead, all fixed at 1, Ne free
    - starting values are set by user
    - output names have 'HSWCX1v8' prefix     
    - performed well!
    
HSWCX1_indfit_v9.py
    - halo kT and N fixed, using vapec for halo instead, all fixed at 1, Ne free
    - starting values are set by user
    - background fit with absorbed CXB, LHB, and fixed absorbed halo
    - above 2 keV with below 0.5 keV included, then BKG Ns fixed.
    - output names have 'HSWCX1v9' prefix      
    - not performing any better below 0.5 keV in residuals than fitting simultaneously
    
HSWCX1_indfit_v10.py
    - halo kT and N fixed, using vapec for halo instead, all fixed at 1, Ne free
    - starting values are set by user
    - background fit with absorbed CXB, LHB, and fixed absorbed halo
    - fitting PI and N for BKG
    - above 2 keV with below 0.5 keV included, then BKG PIs and Ns fixed.
    - output names have 'HSWCX1v10' prefix      
    - not performing any better below 0.5 keV in residuals than fitting simultaneously    

HSWCX1_indfit_v11.py
    - halo kT and N fixed, using vapec for halo instead, all fixed at 1, Ne and Mg free
    - starting values are set by user
    - output names have 'HSWCX1v11' prefix     
    - performed well!    

HSWCX1_indfit_v12.py
    - halo parameters fixed to kT =  and N = 
    - nH allowed to be variable, starting at value type chosen
    - output names have 'HSWCX1v12' prefix      

HSWCX1_indfit_v13.py
    - halo parameters fixed to kT =  and N = 
    - nH allowed to be variable, starting at value type chosen
    - Ne allowed to be variable
    - output names have 'HSWCX1v13' prefix 

HSWCX1_indfit_v14.py
    - halo parameters fixed to kT =  and N = 
    - BKG PI and N fit simultaneously
    - output names have 'HSWCX1v14' prefix    

HSWCX1_indfit_v15.py
    - halo parameters fixed to kT =  and N = 
    - BKG PI and N fit simultaneously
    - no O lines or vapec. Only BKG pars free.
    - output names have 'HSWCX1v15' prefix    
    
HSWCX1_indfit_v16.py
    - halo parameters fixed to kT =  and N = 
    - BKG N fit over 3-7 keV, then fixed
    - no O lines or vapec. Only BKG pars free for 3-7 keV, then fixed.
    - output names have 'HSWCX1v16' prefix   

HSWCX1_indfit_v17.py
    - halo parameters fixed to kT =  and N = 
    - BKG PI and N fit over 3-7 keV, then fixed
    - no O lines or vapec. Only BKG pars free for 3-7 keV, then fixed.
    - output names have 'HSWCX1v17' prefix   

HSWCX1_indfit_v18.py
    - halo parameters fixed to kT =  and N = 
    - BKG PI and N fit over 2-7 keV, then fixed
    - no O lines or vapec. Only BKG pars free for 3-7 keV, then fixed.
    - output names have 'HSWCX1v18' prefix 

HSWCX1_indfit_v19.py
    - halo parameters fixed to kT =  and N = 
    - second apec added with higher temp, free N
    - output names have 'HSWCX1v19' prefix     
    
HSWCX1_indfit_v20.py
    - halo parameters fixed to kT =  and N = 
    - second apec added with higher temp, free N
    - added gaussian for Ne IX line (905 eV from Fujimoto 2007)
    - added gaussian for C IV line (459 eV from Fujimoto 2007)
    - output names have 'HSWCX1v20' prefix       
    
HSWCX1_indfit_v21.py
    - halo parameters fixed to kT =  and N = 
    - powerlaw instead of second apec, free PI and N
    - added gaussian for Ne IX line (905 eV from Fujimoto 2007)
    - added gaussian for C IV line (459 eV from Fujimoto 2007)
    - output names have 'HSWCX1v21' prefix      
    
HSWCX1_indfit_v22.py
    - halo parameters limited to region loosely defined by nearby halo field values
    - powerlaw instead of second apec, free PI and N
    - added gaussian for Ne IX line (905 eV from Fujimoto 2007)
    - added gaussian for C IV line (459 eV from Fujimoto 2007)
    - output names have 'HSWCX1v22' prefix   

HSWCX1_indfit_v23.py
    - halo parameters fixed to kT =  and N = 
    - powerlaw instead of second apec, free PI and N
    - removed gaussian for Ne IX line (905 eV from Fujimoto 2007)
    - removed gaussian for C IV line (459 eV from Fujimoto 2007)
    - output names have 'HSWCX1v23' prefix       
    
HSWCX1_indfit_v24.py
    - halo parameters fixed to kT =  and N = 
    - BKG PI and N fit over entire energy range
    - output names have 'HSWCX1v24' prefix  
    
HSWCX1_indfit_v25.py
    - halo parameters fixed to kT =  and N = 
    - BKG PI and N fit over entire energy range
    - attempting bknpower instad of pow for bkg
    - output names have 'HSWCX1v25' prefix 
    - froze
    
HSWCX1_indfit_v26.py
    - halo parameters fixed to kT =  and N = 
    - BKG PI and N fit over entire energy range, upper fixed to HH values
    - attempting bknpower instead of pow for bkg
    - output names have 'HSWCX1v26' prefix
    
HSWCX1_indfit_v27.py
    - halo N limited to region loosely defined by nearby halo field values, kT fixed
    - BKG PI and N fit over entire energy range
    - output names have 'HSWCX1v27' prefix    
    
HSWCX1_indfit_v28.py
    - halo N and kT fixed
    - BKG PI and N fit over entire energy range
    - CXB N slightly higher and fixed.
    - output names have 'HSWCX1v28' prefix    
    
HSWCX1_indfit_v29.py
    - halo kT and N fixed, using vapec for halo instead, all fixed at 1, Ne free
    - starting values are set by user
    - BKG PI and N fit simultaneously
    - output names have 'HSWCX1v29' prefix     

HSWCX1_indfit_v30.py
    - halo kT and N fixed
    - starting values are set by user
    - BKG N fit simultaneously
    - output names have 'HSWCX1v30' prefix  
    - trying nei for contribution from plasma sheet?

HSWCX1_indfit_v30.py
    - halo kT and N fixed
    - starting values are set by user
    - BKG N fit simultaneously
    - output names have 'HSWCX1v31' prefix  
    - trying npshock for contribution from plasma sheet?

HSWCX1_indfit_v31.py
    - halo kT and N fixed
    - starting values are set by user
    - BKG N fit simultaneously
    - output names have 'HSWCX1v31' prefix 
    
HSWCX1_indfit_v32.py
    - halo kT and N fixed
    - BKG N fit simultaneously
    - Adding gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - output names have HSWCX1v32 prefix    

HSWCX1_indfit_v33.py
    - halo kT and N fixed
    - BKG PI and N fit simultaneously
    - Adding gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - output names have HSWCX1v33 prefix       
    
HSWCX1_indfit_v34.py
    - halo kT fixed, N free, loosely restricted
    - BKG PI and N fit simultaneously
    - Adding gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - output names have HSWCX1v34 prefix  
    - froze  

HSWCX1_indfit_v35.py
    - halo kT and N fixed
    - BKG PI and N fit over 3-7 keV, then fixed
    - Same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - output names have HSWCX1v35 prefix  

HSWCX1_indfit_v36.py
    - halo kT and N fixed
    - BKG PI and N fit simultaneously
    - Same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - Added another linked pow with harder PI and free N, fit simultaneously
    - output names have HSWCX1v36 prefix   
    - froze
    
HSWCX1_indfit_v37.py
    - halo kT and N fixed
    - BKG PI and N fit simultaneously
    - Same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - Added another non-linked pow with harder PI and free N, fit simultaneously
    - output names have HSWCX1v36 prefix  

HSWCX1_indfit_v38.py
    - halo kT and N fixed
    - BKG PI and N fit simultaneously
    - Same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - Added another pow with harder PI and free N, fit simultaneously
    - output names have HSWCX1v38 prefix   
    - additional pow through diagonal response, not linked (particles?)    
    
HSWCX1_indfit_v39.py
    - halo kT and N fixed
    - BKG PI and N fit simultaneously
    - same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - O7 line seems slightly off. Freeing but restricted.
    - output names have HSWCX1v39 prefix 

HSWCX1_indfit_v40.py
    - halo kT and N fixed
    - BKG N fit simultaneously, PI fixed
    - Same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - Added another pow with harder PI and free N, fit simultaneously
    - output names have HSWCX1v40 prefix   
    - additional pow through diagonal response, not linked (particles?)    
    - froze    
    
B: negative power law for second one, froze
C: fixed PI to zero

HSWCX1_indfit_v41.py
    - halo kT and N fixed
    - BKG PI and N fit simultaneously
    - same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - O7 line seems slightly off. Fixing to slightly higher value
    - output names have HSWCX1v41 prefix 
    
HSWCX1_indfit_v42.py
    - halo kT and N fixed
    - BKG N fit simultaneously, PI fixed
    - Same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - Added another pow with harder PI and free N, fit simultaneously
    - output names have HSWCX1v42 prefix   
    - additional pow through full response, (electron xray spectra?)
    - froze, fixing PI of additional pow to zero (B)      
    
HSWCX1_indfit_v43.py
    - halo kT and N fixed
    - BKG PI and N fit simultaneously
    - Same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - Added another pow with PI=0 and free N, fit simultaneously (Xray)
    - output names have HSWCX1v43 prefix   
    
HSWCX1_indfit_v44.py
    - halo kT and N fixed
    - BKG PI and N fit simultaneously
    - Same gaussians for CVI, NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - output names have HSWCX1v44 prefix   
    - additional pow through diagonal response with PI=0, not linked (BKG)   
    
HSWCX1_indfit_v45.py
    - halo kT and N fixed
    - BKG PI and N fit simultaneously
    - Adding gaussians for NeIX, and Mg XI
    - same error loop method as HSWCX1v3b
    - output names have HSWCX1v45 prefix      
    - removed gaussians for O VIII and C VI    
    
HSWCX1_indfit_v46.py
    - adapted from v33
    - relative line intensities fixed (gaussians) relative to O VII
    - using AlphaHe data only, ignoring lines with less than 10% rel int compared to 506.9 eV line
    - attempt to fit halo parameters simultaneously with background
    - (later: try again with bkg PI fixed?)    

HSWCX1_indfit_v47.py
    - adapted from v33
    - relative line intensities fixed (gaussians) relative to O VII
    - using AlphaHe data only, ignoring lines with less than 10% rel int compared to 506.9 eV line
    - attempt to fit halo parameters simultaneously with background with bkg PI fixed  
    
HSWCX1_indfit_v48.py
    - adapted from v33
    - relative line intensities fixed (gaussians) relative to O VII
    - using AlphaHe data only, ignoring lines with less than 10% rel int compared to 506.9 eV line
    - attempt to fit halo parameters simultaneously with background with bkg PI fixed
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt  
    
HSWCX1_indfit_v49.py
    - using acx instead of gaussians, model 7
    - halo temp and N fixed
    - now using calculated BKG_PI from Phil's southern halo analysis, based on hard ctrt
    
HSWCX1_indfit_v50.py
    - adapted from v33
    - relative line intensities fixed (gaussians) relative to O VII
    - using AlphaHe data only, ignoring lines with less than 10% rel int compared to 506.9 eV line
    - halo fixed, to look at residuals
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt        
    
HSWCX1_indfit_v51.py
    - adapted from v33
    - relative line intensities fixed (gaussians) relative to O VII
    - using AlphaH data only, ignoring lines with less than 10% rel int compared to 506.9 eV line
    - halo fixed, to look at residuals
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 

HSWCX1_indfit_v52.py
    - adapted from v33
    - relative line intensities free, linked between detectors
    - ignoring lines with less than 10% rel int compared to 506.9 eV line
    - halo fixed, to look at residuals
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt               

HSWCX1_indfit_v53.py
    - adapted from v33
    - relative line intensities fixed (gaussians) relative to O VII
    - using AlphaHe data only, ignoring lines with less than 10% rel int compared to 506.9 eV line
    - halo temp fixed, halo N free
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 

HSWCX1_indfit_v54.py
    - adapted from v33
    - relative line intensities free (gaussians) relative to O VII
    - halo temp fixed, halo N set to zero
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    
HSWCX1_indfit_v55.py
    - adapted from v33
    - relative line intensities fixed (gaussians) relative to O VII
    - using AlphaHe data only, ignoring lines with less than 10% rel int compared to 506.9 eV line
    - halo temp fixed, halo N = 0
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt   
    
HSWCX1_indfit_v56.py
    - adapted from v33
    - relative line intensities fixed (gaussians) relative to O VII
    - using AlphaHe data only, ignoring lines with less than 10% rel int compared to 506.9 eV line
    - halo temp fixed, halo N free (larger range than v53)
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt     

HSWCX1_indfit_v57.py
    - adapted from v33
    - relative line intensities fixed (gaussians) relative to O VII
    - using AlphaH data only, ignoring lines with less than 10% rel int compared to 506.9 eV line
    - halo temp fixed, halo N free (larger range than v53)
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt  
    
HSWCX1_indfit_v58.py
    - adapted from v33
    - relative line intensities free (gaussians) relative to O VII
    - halo temp fixed, halo N set to zero
    - bkg PI free (compare X^2 to v54), started at calc1 values     
    - froze for most observations  
    
HSWCX1_indfit_v59.py
    - halo kT and N fixed
    - BKG N fit simultaneously
    - Adding gaussians for CVI, NeIX, and Mg XI (same lines as v33)
    - same error loop method as HSWCX1v3b
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt    
    
HSWCX1_indfit_v60.py
    - halo kT and N fixed
    - BKG N fit simultaneously
    - Gaussians for O VII, NeIX, and Mg XI (same lines as v45)
    - same error loop method as HSWCX1v3b
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt   
    
HSWCX1_indfit_v61.py
    - adapted from v54
    - relative line intensities free (gaussians) relative to O VII
    - halo fixed (not zero)
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt          

HSWCX1_indfit_v62.py
    - adapted from v54
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero or fixed
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 

HSWCX1_indfit_v63.py
    - adapted from v54
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero or fixed
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - halo abundance set to 0.3 as in Phil's halo analysis

HSWCX1_indfit_v64.py
    - adapted from v54
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero or fixed
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - halo abundance set to 0.3 as in Phil's halo analysis
    - using more thorough calculations for all SWCX lines, ions linked if same
    
HSWCX1_indfit_v65.py
    - adapted from v54
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg over higher energy range (2-7 keV)

HSWCX1_indfit_v66.py
    - adapted from v54
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg over higher energy range (3-7 keV)
    
HSWCX1_indfit_v67.py
    - adapted from v54
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg over entire energy range simultaneously (choose best method from 65, 66, 67)       

HSWCX1_indfit_v68.py
    - adapted from v54
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg N over entire energy range simultaneously 
    - free nH

HSWCX1_indfit_v69.py
    - adapted from v67
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg over entire energy range simultaneously
    - adding second Ne IX line at 1.08 keV

HSWCX1_indfit_v70.py
    - adapted from v67
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg over entire energy range simultaneously
    - adding second Ne IX line at 1.08 keV
    - adding low T Fe XVII line at 0.73164 keV

HSWCX1_indfit_v71.py
    - adapted from v67
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg over entire energy range simultaneously
    - adding second Ne IX line at 1.08 keV
    - adding low T Fe XVII line at 0.73164 keV
    - adding Fe XVI line at 0.9945 keV 
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

#xs.Xset.addModelString("APECROOT","/home/ringuette/Software/heasoft-6.25/spectral/modelData/apec_v3.0.9_halo") #use custom apec file

clean_type='cl1'

#define directories
spec_dir = '/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/Rebinned_spectra_'+clean_type+'/'
code_dir='/home/ringuette/halosat_sourceanalysis/AnalysisCode/HSWCXCode/'
out_dir='/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/XspecOutput_'+clean_type+'/'
prefix = 'HSWCX1v71'

#define bkg type and values
bkg_pi, bkg_type = {}, '_CALC1'
bkgpi_file = DictReader(open(code_dir+'HSWCX1_ObsData_cl1.csv','r'))
for row in bkgpi_file:
    #print row.keys(), row['ObsName']
    bkg_pi[row['ObsName']]={14:float(row['Det14_PI']),\
        38:float(row['Det38_PI']),54:float(row['Det54_PI'])}
#stop
#bkg_pi, bkg_type = {38:0.82032, 54:0.84247}, '_SH'  #Phil's numbers using Southern Halo, 0.4-7 keV
#bkg_pi, bkg_type = {38:0.689, 54:0.635}, '_HH'  #RR numbers using High Halo, 3-7 keV           
#for key in bkg_pi.keys():  #use a fixed PI instead
    #bkg_pi[key] = {38:0.82032, 54:0.84247}
    #bkg_pi[key] = {38:0.689, 54:0.635} #likely only for high hard count rates?
#bkg_type = '_SH'
#bkg_type = '_HH'


#halo_kT, halo_N =  0.216, 0.62  #set halo kT and N, halo1
#halo_kT, halo_N =  0.216, 0.46  #set halo kT and N, halo2
#halo_kT, halo_N =  0.242, 0.442  #set halo kT and N, halo3 same as halo0128 field
#halo_kT, halo_N = 0.216, 0.2   #set halo kT, start halo N low but free, halo4
#halo_kT, halo_N =  0.216, 0.06  #halo5
halo_kT, halo_N = 0.216, 0.
#halo_kT, halo_N =  0.216, 2.07  #set halo kT and N, halo6 (N/0.3 b/c halo abund=0.3)
#halo_kT, halo_N =  0.216, 1.12  #set halo kT and N, halo7 (N higher b/c halo abund=0.3)
#halo_kT, halo_N =  0.216, 1.62  #set halo kT and N, halo8 (N higher b/c halo abund=0.3)
#prefix+='_halo8'   #for kT = 0.216 and N = 0.62
CXB_PI, CXB_N = 1.45, 0.38
#CXB_PI, CXB_N = 1.45, 0.382  #CXB0
#CXB_PI, CXB_N = 1.45, 0.376  #CXB1
#CXB_PI, CXB_N = 1.45, 0.387  #CXB2
#CXB_PI, CXB_N = 1.47, 0.382  #CXB3
#CXB_PI, CXB_N = 1.43, 0.382  #CXB4
#CXB_PI, CXB_N = 1.45, 0.393  #CXB5
#CXB_PI, CXB_N = 1.45, 0.399  #CXB6
#CXB_PI, CXB_N = 1.45, 0.475  #CXB0p475
#prefix+='_CXB0p475'

nH_type = '_nHSFDW' #set type of nH used options: '_nHLAB', '_nHHI4pi', '_nHSFDW', '_nHPZ'
model_note = nH_type+bkg_type  #store choices

#define ion ratios per obs (only calculated for good observations)
line_ints = {}
line_file = DictReader(open(code_dir+'HSWCX1_EO7R.csv','r'))
for row in line_file:
    line_ints[row['\xef\xbb\xbfObsName']]=row

#initialize output format
file_header='Target,LHB_N,'+nH_type.split('_')[1] #all other headers not hard-coded
file_format='\n{0[0]},{0[1]:.6f},{0[2]:.6f}'
csv_num=3  #next number in []

#pull LHB N from fits file, store in dict 
t = pyfits.open(code_dir+'LHB_nH_ecl2.fits')  #use new file
hsIDs = t[1].data['Target_ID']
RA = t[1].data['RA']
DEC = t[1].data['Dec']
LHB = t[1].data['LHB_norm']
nh_LAB = t[1].data['nH_LAB']  #nH from LAB. Try using Daniel's total nH
nh_HI4Pi = t[1].data['nH_HI4pi']
nh_SFDW = t[1].data['nH_SFDW']
nh_PZ = t[1].data['nH_PlanckZhu']
t.close()

if nH_type=='_nHLAB': nh = nh_LAB
elif nH_type=='_nHHI4pi': nh = nh_HI4Pi
elif nH_type=='_nHSFDW': nh = nh_SFDW
elif nH_type=='_nHPZ': nh = nh_PZ
elif nH_type=='_nHtest4': nh = np.repeat(0.4,len(nh_LAB))

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
            'Oct2018b+':[LHB[367],nh[367]],\
            'Nov2018+':[(LHB[368]*19662.0+LHB[370]*23628.5)/(19662.0+23628.5),\
                        (nh[368]*19662.0+nh[370]*23628.5)/(19662.0+23628.5)],\
            'Dec2018':[LHB[65],nh[65]],'Feb2019':[LHB[65],nh[65]],\
            'Sep2019':[LHB[65],nh[65]],'Oct2019':[LHB[65],nh[65]],\
             'Nov2019+':[LHB[368],nh[368]],'Nov2019':[LHB[65],nh[65]],\
             'Dec2019+':[LHB[368],nh[368]],'Dec2019':[LHB[65],nh[65]],\
             'Jan2020+':[LHB[368],nh[368]],'Feb2020+':[LHB[368],nh[368]],\
             'Mar2020+':[LHB[368],nh[368]]}
#key_list = ['Oct2018+','Oct2018b+','Nov2018+','Dec2018','Feb2019','Sep2019','Oct2019',\
#            'Nov2019+','Nov2019','Dec2019+','Dec2019','Jan2020+','Feb2020+','Mar2020+']
key_list = ['Oct2018b+','Nov2018+','Oct2019','Nov2019+','Nov2019','Dec2019+','Dec2019']  #good list only
#key_list = ['Oct2019','Dec2019+']
#calculated line ratios are only done for good list, average line ratios used otherwise

#loop through observation names
for key in key_list:
    #choose spectra that are rebinned using grppha
    if key[-1]=='+': 
        spec_name='HSWCX1+_'+key[0:-1]+'_'+clean_type+'_d'
    else: 
        spec_name='HSWCX1_'+key+'_'+clean_type+'_d'
    spec_list=[spec_dir+spec_name+'38_rebin40.pi',spec_dir+spec_name+'54_rebin40.pi']
    
    #collect parameters for observation
    #LHB_target, nh_target = fixed_pars[key][0], fixed_pars[key][1]
    LHB_target, nh_target = 0., fixed_pars[key][1]   #set norm of LHB to zero
    NVI_ratio = str('{:.5f}'.format(float(line_ints[key]['N VIb'])/float(line_ints[key]['N VIa'])))
    NVII_ratio = str('{:.5f}'.format(float(line_ints[key]['N VIIb'])/float(line_ints[key]['N VIIa'])))
    OVIII_ratio = str('{:.5f}'.format(float(line_ints[key]['O VIIIb'])/float(line_ints[key]['O VIIIa'])))
    MgXI_ratio = str('{:.5f}'.format(float(line_ints[key]['Mg XIb'])/float(line_ints[key]['Mg XIa'])))
    print NVI_ratio, NVII_ratio, OVIII_ratio, MgXI_ratio

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
    bkg_par, xray_par = {}, {}
    bmodel, xmodel, bkg_spec = [], [], []
    halo_arr=[key,LHB_target,nh_target]
    
    #loop through detector spectra
    print 'Initializing '+key
    for i in range(len(spec_list)): #fit model to each spectrum linked between detectors    
        bkg_spec.append(0)
        bmodel.append(0)
        xmodel.append(0)
        
        #set up files    
        xs.AllData(str(i+1)+':'+str(i+1)+' '+spec_list[i]) #assign each spectrum to diff data groups
        bkg_spec[i]=xs.AllData(i+1)    
        bkg_spec[i].multiresponse[0] = code_dir+'halosat_avgenoise_20190423.rmf'
        bkg_spec[i].multiresponse[0].arf = code_dir+'halosat_20190211.arf'
        bkg_spec[i].multiresponse[1] = code_dir+'halosat_diag_20190423.rmf'   #only rmf for particle bkg fit
        
        #add exposure time and bkg_pi to output
        dpu = int(spec_list[i].split('.')[0].split('_d')[1].split('_rebin40')[0]) #get dpu num
        halo_arr.append(bkg_spec[i].exposure)
        if key==key_list[0]:
          file_header+=',d'+str(dpu)+' Exp' #all other headers not hard-coded
          file_format+=',{0['+str(csv_num)+']:.3f}'
          csv_num+=1  #next number in [] 
        if len(bkg_pi)>3:  #if using calculated BKG_PIs
          halo_arr.append(bkg_pi[key][dpu])       
          if key==key_list[0]:
            file_header+=',d'+str(dpu)+' BKG_PI' #all other headers not hard-coded
            file_format+=',{0['+str(csv_num)+']:.3f}'
            csv_num+=1  #next number in [] 
                  
        #set up model
        bkg_par[i]={1:str(bkg_pi[key][dpu])+',-1',2:'0.5,0.1,0.,0.,1e6,1e6'}  #BKG PI fixed and N free per det
    
        #using O VII, O VIII, Liu values for LHB (fixed apec), 
        #absorbed CXB (fixed to Capelluti values) and absorbed linked halo APEC
        #SWCX lines from Alpha file from Dimitra, AlphaHe section
        if (i==0):    
            xray_par[i]={1:str(nh_target)+',-1',\
                          2:str(halo_kT)+',-1',3:'1.,-1',4:'0,-1',5:str(halo_N)+',-1',\
                          6:str(CXB_PI)+',-1',7:str(CXB_N)+',-1',\
                          8:'0.42171,-1',9:'0.001,-1',10:'0.1,0.1,0,0,1e+20,1e+20',\
                          11:'0.42162,-1',12:'0.001,-1',13:'0.1,0.1,0,0,1e+20,1e+20',\
                          14:'0.44696,-1',15:'0.001,-1',16:'0.1,0.1,0,0,1e+20,1e+20',\
                          17:'0.49790,-1',18:'0.001,-1',19:'=Xray:p13*'+NVI_ratio,\
                          20:'0.50000,-1',21:'0.001,-1',22:'0.1,0.1,0,0,1e+20,1e+20',\
                          23:'0.56341,-1',24:'0.001,-1',25:'0.1,0.1,0,0,1e+20,1e+20',\
                          26:'0.61835,-1',27:'0.001,-1',28:'=Xray:p22*'+NVII_ratio,\
                          29:'0.65310,-1',30:'0.001,-1',31:'0.1,0.1,0,0,1e+20,1e+20',\
                          32:'0.80329,-1',33:'0.001,-1',34:'=Xray:p31*'+OVIII_ratio,\
                          35:'0.90881,-1',36:'0.001,-1',37:'0.1,0.1,0,0,1e+20,1e+20',\
                          38:'1.33498,-1',39:'0.001,-1',40:'0.1,0.1,0,0,1e+20,1e+20',\
                          41:'1.57000,-1',42:'0.001,-1',43:'=Xray:p40*'+MgXI_ratio,\
                          44:'1.08623,-1',45:'0.001,-1',46:'0.1,0.1,0,0,1e+20,1e+20',\
                          47:'0.73164,-1',48:'0.001,-1',49:'0.1,0.1,0,0,1e+20,1e+20',\
                          50:'0.9945,-1',51:'0.001,-1',52:'0.1,0.1,0,0,1e+20,1e+20'}
        else: 
            xray_par[i]={1:'=Xray:p1',2:'=Xray:p2',3:'=Xray:p3',4:'0,-1',5:'=Xray:p5',\
                          6:'=Xray:p6',7:'=Xray:p7',\
                          8:'=Xray:p8',9:'0.001,-1',10:'=Xray:p10',\
                          11:'=Xray:p11',12:'0.001,-1',13:'=Xray:p13',\
                          14:'=Xray:p14',15:'0.001,-1',16:'=Xray:p16',\
                          17:'=Xray:p17',18:'0.001,-1',19:'=Xray:p19',\
                          20:'=Xray:p20',21:'0.001,-1',22:'=Xray:p22',\
                          23:'=Xray:p23',24:'0.001,-1',25:'=Xray:p25',\
                          26:'=Xray:p26',27:'0.001,-1',28:'=Xray:p28',\
                          29:'=Xray:p29',30:'0.001,-1',31:'=Xray:p31',\
                          32:'=Xray:p32',33:'0.001,-1',34:'=Xray:p34',\
                          35:'=Xray:p35',36:'0.001,-1',37:'=Xray:p37',\
                          38:'=Xray:p38',39:'0.001,-1',40:'=Xray:p40',\
                          41:'=Xray:p41',42:'0.001,-1',43:'=Xray:p43',\
                          44:'=Xray:p44',45:'0.001,-1',46:'=Xray:p46',\
                          47:'=Xray:p47',48:'0.001,-1',49:'=Xray:p49',\
                          50:'=Xray:p50',51:'0.001,-1',52:'=Xray:p52'}
                          
    #set dictionary of component names for free parameters. change if model changes
    xray_names={2:'Halo_kT',5:'Halo_N',10:'SIX_N',13:'NVI_N',16:'CVI_N',22:'NVII_N',\
                  25:'OVII_N',31:'OVIII_N',37:'NeIX_N',40:'MgXI_N',46:'NeIXb_N',\
                  49:'FeXVIIlowT_N',52:'FeXVI_N'}  #for O lines, det# is added later
    
    #set Xspec model and parameters
    bkg_model=xs.Model('pow','BKG',2)  #particle background
    xray_model=xs.Model('phabs*(apec+pow)+gauss+gauss+gauss+gauss+gauss+gauss+gauss+gauss+gauss+gauss+gauss+gauss+gauss+gauss+gauss','Xray',1)   #OVII+OVIII+LHB+absorbed halo+absorbed CXB+lines
    for i in range(len(spec_list)):    
        bmodel[i]=xs.AllModels(i+1,'BKG')
        xmodel[i]=xs.AllModels(i+1,'Xray')
        bmodel[i].setPars(bkg_par[i])
        xmodel[i].setPars(xray_par[i])
        
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
                xs.Fit.steppar('best BKG:'+str(parnum)+' 0.4 2. 25')                
            if (bmodel[i](j+1).index==2):  #BKG N
                print 'steppar: BKG', parnum
                xs.Fit.steppar('best BKG:'+str(parnum)+' 0. 0.2 25')
    #run steppar on all free xray parameters
    for i, j in itertools.product(range(len(spec_list)), range(xmodel[i].nParameters)):
        if not ((xmodel[i](j+1).frozen)or(len(xmodel[i](j+1).link)>0)):
            parnum=xmodel[i].startParIndex+xmodel[i](j+1).index-1
            if (xmodel[i](j+1).index==2):  #halo kT
                print 'steppar: Xray', parnum
                xs.Fit.steppar('best Xray:'+str(parnum)+' 0.01 2.0 25')
            elif (xmodel[i](j+1).index==5):  #halo N
                print 'steppar: Xray', parnum
                xs.Fit.steppar('best Xray:'+str(parnum)+' 0. 0.5 25')    
            else:  #line Ns
                print 'steppar: Xray', parnum
                xs.Fit.steppar('best Xray:'+str(parnum)+' 0. 0.5 25')                   

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
        folded[i] = np.array(xs.Plot.model(i+1,1),dtype=float)
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