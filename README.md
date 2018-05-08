pyplanet
========

planetary atmosphere code

README file for pyPlanet

Before you start:

1. Set up your **config file**
    1. Nearly all of the parameters are set within this configuration file
    2. The defaults are set in config.json, which also sets up which parameters are contained within the config
    3. Default name is `config.par` in each of the planet sub-directories
    4. You may use different files, just call planet with config='filename', which must reside within <planet_name>
    5. The tweakmodule filename gets set within the config, which adjusts the read in atmosphere - make sure it is what you want or make it so.
    6. For mcmc, the scalemodule also gets set here (and the scalefilename)
2. When you make a planet instance you may also set state_variables, they are initialized in state_variables.py
3. When you make a planet instance you may specify a mode, which sets the state_variables to various configurations (see state_variables.py)
4. The planet call variables and defaults are shown below.  kwargs may be one of the state_variables
    planet:  name, mode='normal', config='config.par', \*\*kwargs
5. By convention, the instance is set to the lower-case first letter of the planet (not required)


Within a python environment here is an example:

\> import planet
\> j = planet.Planet('jupiter',config='config_testing.par')
\> catch_data = j.run(freqs='1:100:5', b='disc')

time-stamped data file is written to Output and log file to Logs

options for freqs:
    freqs = 1.42    ==> single frequency at 1.42 GHz
    freqs = [1.4,2.5,5.5,8.4,12.1]  ==> uses these frequencies (can be a list, csv string or numpy array)
    freqs = '1:10:1' ==> this will generate a range as start:stop:step (stop is always included)
    freqs = '1:100:20' ==> this will generate a log range as start;stop;nvalues
    freqs = 'freq.dat'   ==> (string within '') reads in those frequencies, one per line


options for b:
        b = 0.1  ==> generates a full image at that resolution (see blocks)
        b = 'stamp' ==> generates a small image (queries for extents)
        b = [[0.0,0.0],[0.1,0.0],...]  ==> generates at listed points
        b = [0.0,0.0] ==> same as above at that one point
        b = 'disc' ==> disc-averaged brightness temperature
        b = '0.0,0.2,0.4,0.6,0.8,1.0<45' ==> will use those values as an angle of 45deg, (default angle is 0.0)
        b = '0.0:1.0:0.1<45' ==> range start:stop:step<angle                                                  >
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

