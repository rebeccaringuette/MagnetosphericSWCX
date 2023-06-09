This folder contains the python scripts (pyXSPEC in Python 2.6 or 2.7), the arf, rmfs and other input files used in the analysis. Updates to the code will be needed to accommodate more recent versions of XSPEC and Python before adapting to an updated analysis. The scripts listed here belong to four categories:

- 'Tail' scripts analyze the cleaned spectra obtained from lines of sight down the tail of the magnetosphere,
- 'Flank' scripts analyze the cleaned spectra obtained from lines of sight through the flank of the magnetosphere, and
- 'Hot' scripts analyze the cleaned spectra from the indicated type of line of sight, but with a second absorbed thermal component for the Galactic Halo.
- If 'Hot' is not in the script name, then the script analyzes the cleaned spectra from the indicated type of line of sight with the default spectral analysis, namely excluding the second absorbed thermal component for the Galactic Halo.  

Other files include the rmfs used, the XSPEC scripts to sum the composite spectra, the input files for the heliospheric SWCX calculations per observation and the Hydrogen column absorptions. The 'obs_att_cl_win.py' script produced the fits file containing spacecraft position and velocity information used to determine the O VII line intensity predicted by the MHD model for the magnetospheric SWCX contribution. 
