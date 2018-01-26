pyplanet
========

planetary atmosphere code

README file for pyPlanet

Before you start:
    make sure the Jupiter/config.par file is accurate 
    make sure the gas/cloud files exist
        make sure that Jupiter/JupiterTweak.py is what you think it is
    use.py to confirm which constituents are being used 
        use.py --nh3 to toggle e.g. nh3

start ipython --pylab (you set the alias ipy for that)

import planet
j = planet.planet('jupiter')
j.run(freqs=..., b=...)

time-stamped data file is written to Output and log file to Logs

options for freqs:
    freqs = 1.42    ==> single frequency at 1.42 GHz
    freqs = [1,4,1]  ==> (i.e. three values within []) interprets as start,stop,step
    freqs = [1.4,2.5,5.5,8.4,12.1]  ==> (i.e. a list within [] of length other than 3) runs at those frequencies
    freqs = 'freq.dat'   ==> (string within '') reads in those frequencies, one per line


options for b:
        Output is "Image"
        b = 0.1  ==> generates a full image at that resolution (see blocks)
        b = 'stamp' ==> generates a small image (queries for extents)
        Output is "Spectrum"
        b = [[0.0,0.0],[0.1,0.0],...]  ==> generates at listed points
        b = [0.0,0.0] ==> same as above at that one point
        b = 'disc' ==> disc-averaged brightness temperature
    Output is "Profile"
        b = [45.0,0.0,0.1,0.2,0.3,0.4,0.5,0.9,0.95] ==> generates a line at angle of first term (45deg) at magnitude of rest
        b = [45.0,0.0,1.0,0.02]  ==> if length is 4, it assumes start, stop, step for last three
==================================================================================================



*****12/7/17 David DeBoer - created README file
This README file explains the structure and execution of the planetary atmosphere code.  It comprises the following directories:

.(pyModel) : contains the high-level python modules
	planet.py - executive module for automatically running the pipeline
	atmosphere.py - module to compute the atmospheric structure
	microwave.py - module to compute the opacity at each layer
	brightness.py - module to compute the overall brightness temperature
[planetName]: contains planet-specific data input files [Jupiter, Saturn, Uranus, Neptune]
constituents:  sub-directories contain the modules and files to compute the microwave opacity of the constituents
modelSupportInfo:  back-up information/data etc
=================================================================================================================

*****12/12/7 David DeBoer 
Initial posting of pyPlanet (renamed from pyModel) to googleCode
need to set:
	git remote https://code.google.com/p/planetary-atmospheres/
need to add to .netrc
	machine code.google.com/p/planetary-atmospheres/ login username@gmail.com password [currently Vc9aw3ZE5qR6]
then you can:
	git push https://code.google.com/p/planetary-atmospheres/ 
or
	git clone "
Additional files/directories:
	raypath.py - used to calculate raypaths in atmosphere
	regrid.py - regrids data to common array
[planetName]:  added tweakFile.py to contain planet/run dependent atmospheric tweaks
Logs/:  directory where log files are written
Output/: directory where output files are written

*****13/3/20 David DeBoer
git commit -a -m "This skips the adding/staging part..."

