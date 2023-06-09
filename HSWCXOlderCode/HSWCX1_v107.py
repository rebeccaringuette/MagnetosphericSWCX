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
    
HSWCX1_indfit_v71b.py
    - LHB and halo = 0
    - use new SWCX ratios with upper lines (csv file v2)
    - separate SWCX model from Xray model
    - put results in a new folder (new SWCX model)
    - also, new nH from Kip

HSWCX1_indfit_v71c.py
    - LHB fixed (not zero)    
    
HSWCX1_indfit_v72.py
    - adapted from v71
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg over entire energy range simultaneously
    - adding second Ne IX line at 1.08 keV
    - adding low T Fe XVII line at 0.73164 keV
    - adding Fe XVI line at 0.9945 keV
    - fitting CXB_N    (froze)

HSWCX1_indfit_v73.py
    - adapted from v67
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg over entire energy range simultaneously
    - adding second Ne IX line at 1.08 keV
    - adding low T Fe XVII line at 0.73164 keV
    - adding Fe XVI line at 0.9945 keV 
    - looping through CXB Ns

HSWCX1_indfit_v74.py
    - adapted from v71
    - relative line intensities free (gaussians) relative to O VII
    - LHB zero, halo zero 
    - now using calculated PI from Phil's southern halo analysis, based on hard ctrt 
    - using more thorough calculations for all SWCX lines, ions linked if same
    - fitting bkg over entire energy range simultaneously
    - adding second Ne IX line at 1.08 keV
    - adding low T Fe XVII line at 0.73164 keV
    - adding Fe XVI line at 0.9945 keV 
    - T Tauri model: absorbed 2-temp component with special abundances
    - switching to anders and grevesse abundances for ease of TTS model
    
HSWCX1_indfit_v75.py
    - v74 + fixed LHB    

HSWCX1_indfit_v76.py
    - v75 + fitted halo  (froze)
    
HSWCX1_indfit_v77.py
    - v73, looping through nH instead of CXB_N
    
HSWCX1_indfit_v78.py
    - v74 + Ne abundance free in TTS model (linked)       (froze)
    
HSWCX1_indfit_v79.py
    - v74 + scan of Ne abundances (linked)
    
HSWCX1_indfit_v80.py
    - v67, LHB and Halo fixed
    - bkg Ns fixed to v71HH values, using HH 
    
HSWCX1_indfit_v81.py
    - v76, fixed CXB, nH (SFDW)
    - fixed LHB, halo temp fixed to 0.222 from target 135
    - calc bkg PI, free bkg N
    - removed thermal lines
    - SWCX ratios fixed
    - approx TTS as hot vvapec, special abund based on wilms    
    - try free halo and TTS Ns, free TTS kt, small nH on TTS?
    
HSWCX1_indfit_v82.py
    - v81 and nH, halo N both free
    
HSWCX1_indfit_v83.py
    - v82 with halo kT free    
    
HSWCX1_indfit_v84.py
    - v83, also free gal nh to be larger than SFDW (more absorption due to other molecules in cloud?)
    - fix TTS nh at SFDW for now
    - allow high halo_Ns
    
HSWCX1_indfit_v85.py
    - v84, add low temp vvapec to TTS model (more physical)
    
HSWCX1_indfit_v86.py
    - v84, all nH more free
    
HSWCX1_indfit_v87.py
    - v85, all nH more free, halo=0 and fixed  
    
HSWCX1_indfit_v88.py
    - v85, second TTS comp with higher temp as in Gudel et al 2007    
    
HSWCX1_indfit_v89.py
    - v88 copied, trying lower temp for LHB
    
HSWCX1_indfit_v90.py
    - v86 with lower temp for LHB   
    
HSWCX1_indfit_v91.py
    - v86 with TTS nH = 0.    
    
HSWCX1_indfit_v92.py
    - v71c new SWCX model, nH = 0.248 from Kip
    - halo kT and N free, absorbed
    - TTS nH = 0, kT and N free

v92a - halo kT fixed to 0.224 (okay, TTS kT wants to be higher than 10 keV)
v92b - TTS N = 0, halo kT still fixed to 0.224 keV  (not good near O VII and Ne IX
v92c - TTS nH free, limited by halo nH, halo kT fixed
v92d - TTS nH fixed to ROSAT value, halo kT fixed

v93: read OVII LU from file and fix.

v94: test TTS model, solar abundances, fixed temps, and linked Ns;
      all lines below O VII combined, N VIIb ignored for simplicity (LU<0.01 always),
      MgXI also removed. Changed to tbabs and 0.3 metallicity for halo.
v94nH: TTS completely absorbed.
v94a: TTS nH free (froze)

v95: TTS (absorbed AND unabsorbed), from v94 (froze)

v96: scan through nH for TTS component. removed lowE since always zero. Halo N free
v96a: Halo_N = 0

v97: SWCX OVII free, TTS_N free, halo_N=0, LowE_N =0.
v98: all primary SWCX free, TTS_N free, halo_N=0, lowE_N=0.

v99: LowE=0, SWCX fixed, halo_N, TTS_N, and lowT on TTS free
v99a: halo_N=0
v99b: halo_N=0, lowT fixed to 0.831 keV, TTS_N free
v99c: halo_N=0, lowT fixed to 0.808 keV, TTS_N free
v99dalla: halo_N=0, all detectors included
v99dallc: halo_N=0, all detectors included, TTS_T1 = 0.805

v100: pull TTS_N from v99c output, lowT fixed to 0.808 keV, halo_N=0, 
      SWCX primary lines free, LowE_N=0 
v100dall: pull TTS_N from v99dallc output, lowT fixed to 0.805 keV, halo_N=0, 
      SWCX primary lines free, LowE_N=0   
v100dalla: TTS_nH=0, halo_kT=0.224 keV, halo_N=0.239, TTS_T1=0.805 keV, pull
      TTS_N from v99dalle, LowE_N=0    
v100dallb:  TTS_nH=0, halo_kT=0.224 keV, halo_N=0.238, TTS_T1=0.807 keV, pull
      TTS_N from v99dalle3, LowE_N free with other SWCX lines.    
      
v101: No TTS, using new SWCX from Dimitra (fixed), only halo.   
v102: same as 101, add second high temp apec (absorbed)  
v103: low temp halo=0 (removed), high temp halo only.   
v104: same as 102, low T halo-0.224 keV from target 135
v105: same as 102, high T comp kT = 0.783 keV and N = 0.308, low kT and N free 
v106: same as v105, halo_kT1 = 0.224 keV
v107: one apec for hot emission, kT = 0.783 keV and N = 0.308, SWCX primary lines free
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
out_dir='/home/ringuette/halosat_sourceanalysis/HSWCX_2020Analysis/XspecOutput_'+clean_type+'_v3/'
prefix = 'HSWCX1v107_GOOD'

#define bkg type and values
bkg_pi, bkg_type = {}, '_CALC1'
bkgpi_file = DictReader(open(code_dir+'HSWCX1_ObsData_cl1.csv','r'))
for row in bkgpi_file:
    #print row.keys(), row['ObsName']
    bkg_pi[row['\xef\xbb\xbfObsName']]={14:float(row['Det14_PI']),\
        38:float(row['Det38_PI']),54:float(row['Det54_PI'])}

halo_kT, halo_N = 0.776, 0.308
CXB_PI, CXB_N = 1.45, 0.382
nH_type = '_nHKip' #set type of nH used options: '_nHRT', '_nHHI4pi', '_nHSFDW', '_nHPZ'
model_note = nH_type+bkg_type  #store choices

#define ion ratios per obs (only calculated for good observations)
line_ints = {}
line_file = DictReader(open(code_dir+'HSWCX1_EO7R_v2.csv','r'))
for row in line_file:
    line_ints[row['\xef\xbb\xbfObsName']]=row

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

#read helioSWCX calcs per spec
SWCX_LU = {}
SWCX_file = DictReader(open(code_dir+'HSWCX1_helioSWCXavg.csv','r'))
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
             'Mar2020+':[LHB[368],nh[368]],'GOOD':[LHB[65],nh[65]]}
#key_list = ['Oct2018+','Oct2018b+','Nov2018+','Dec2018','Feb2019','Sep2019','Oct2019',\
#            'Nov2019+','Nov2019','Dec2019+','Dec2019','Jan2020+','Feb2020+','Mar2020+']
#key_list = ['Oct2018b+','Nov2018+','Oct2019','Nov2019+','Nov2019','Dec2019+','Dec2019']  #good list only
#key_list = ['Dec2019+']
key_list=['GOOD']
#calculated line ratios are only done for good list, average line ratios used otherwise

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
    idx = np.where(SWCX_LU['\xef\xbb\xbfMonth']==key)[0]
    OVII_N, OVIII_N, NeIX_N = float(SWCX_LU['OVII'][idx][0])*0.035, float(SWCX_LU['OVIII'][idx][0])*0.035, \
          float(SWCX_LU['NeIX'][idx][0])*0.035
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
        bkg_spec[i].multiresponse[1].arf = code_dir+'halosat_20190211.arf'
        bkg_spec[i].multiresponse[2] = code_dir+'halosat_avgenoise_20190423.rmf'
        bkg_spec[i].multiresponse[2].arf = code_dir+'halosat_20190211.arf'  
        #bkg_spec[i].multiresponse[3] = code_dir+'halosat_avgenoise_20190423.rmf'
        #bkg_spec[i].multiresponse[3].arf = code_dir+'halosat_20190211.arf'  
        
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
        bkg_par[i]={1:str(bkg_pi[key][dpu])+',-1',2:'0.02,0.02,0.,0.01,0.5,1e6'}  #BKG PI fixed and N free per det
    
        #using Liu values for LHB (fixed apec), 
        #absorbed CXB (fixed to Capelluti values) and absorbed linked halo APEC
        #SWCX lines from Alpha file from Dimitra
        if (i==0):    
            xray_par[i]={1:'0.097,-1',2:'1,-1',3:'0,-1',4:str(LHB_target)+',-1',\
                          5:str(nh_target)+',-1',\
                          6:str(halo_kT)+',-1',7:'0.3,-1',8:'0,-1',9:str(halo_N)+',-1',\
                          10:str(CXB_PI)+',-1',11:str(CXB_N)+',-1'}
            SWCX_par[i] = {1:line_ints['Energy']['LowE']+',-1',2:'0.001,-1',3:'0.,-1',\
                          4:line_ints['Energy']['O VIIa']+',-1',5:'0.001,-1',6:'0.035,0.01,0.,0.,0.06,10.',\
                          7:line_ints['Energy']['O VIIb']+',-1',8:'0.001,-1',9:'=SWCX:p6*'+OVII_ratio,\
                          10:line_ints['Energy']['O VIIIa']+',-1',11:'0.001,-1',12:'0.004,0.01,0.,0.,0.06,10.',\
                          13:line_ints['Energy']['O VIIIb']+',-1',14:'0.001,-1',15:'=SWCX:p12*'+OVIII_ratio,\
                          16:line_ints['Energy']['Ne IXa']+',-1',17:'0.001,-1',18:'0.002,0.01,0.,0.,0.06,10.',\
                          19:line_ints['Energy']['Ne IXb']+',-1',20:'0.001,-1',21:'=SWCX:p18*'+NeIX_ratio} 

        else: 
            xray_par[i]={1:'=Xray:p1',2:'=Xray:p2',3:'=Xray:p3',4:'=Xray:p4',\
                          5:'=Xray:p5',\
                          6:'=Xray:p6',7:'=Xray:p7',8:'=Xray:p8',9:'=Xray:p9',\
                          10:'=Xray:p10',11:'=Xray:p11'}
            SWCX_par[i]={1:line_ints['Energy']['LowE']+',-1',2:'0.001,-1',3:'=SWCX:p3',\
                          4:line_ints['Energy']['O VIIa']+',-1',5:'0.001,-1',6:'=SWCX:p6',\
                          7:line_ints['Energy']['O VIIb']+',-1',8:'0.001,-1',9:'=SWCX:p9',\
                          10:line_ints['Energy']['O VIIIa']+',-1',11:'0.001,-1',12:'=SWCX:p12',\
                          13:line_ints['Energy']['O VIIIb']+',-1',14:'0.001,-1',15:'=SWCX:p15',\
                          16:line_ints['Energy']['Ne IXa']+',-1',17:'0.001,-1',18:'=SWCX:p18',\
                          19:line_ints['Energy']['Ne IXb']+',-1',20:'0.001,-1',21:'=SWCX:p21'}
                          
    #set dictionary of component names for free parameters. change if model changes
    xray_names={6:'Halo_kT',9:'Halo_N'}
    SWCX_names={3:'LowE',6:'O_VII',12:'O_VIII',18:'Ne_IX'}
    #TTS_names={1:'TTS_nH',2:'TTS_kT',34:'TTS_N'}
    
    #set Xspec model and parameters
    bkg_model=xs.Model('pow','BKG',1)  #particle background
    xray_model=xs.Model('apec+tbabs*(apec+pow)','Xray',2)   #LHB+halo+CXB
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
                xs.Fit.steppar('best BKG:'+str(parnum)+' 0.4 2. 25')                
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
            xs.Fit.steppar('best SWCX:'+str(parnum)+' 0. 0.07 14')   #should be less than 2 LU
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