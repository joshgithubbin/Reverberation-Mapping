###################################

REVERBERATION MAPPING PROGRAMMING

###################################

INTRO	.	.	.	RMP1

FILES	.	.	.	RMP2

CONTACT	.	.	.	RMP3


###

INTRO - RPM1

These programs were written in Python 3.8 between October 2019 and March 2020
for my MPhys with Astronomy project. The project involved cross-correlation analysis
between K and V band emission data for a number of low-redshift Active Galactic Nuclei
to determine a relationship between AGN size and luminosity.

The data involved for each AGN included a set of flux recordings and ten simulated flux recordings
in an alternate band. While the programming is geared towards analysing this data
minor alterations could fit it to another purpose.

###

FILES - RMP2

Report.pdf : Contains information on the purpose of the code, its use in context and an explanation of
the mathematical functions involved. See here if you're unsure about anything.

correlationfunctions.py : Contains functions relevant to cross-correlation analysis and obtaining the lag in data.
Requires scipy, numpy.

ensemblesamplerfunctions.py : Contains functions relevant to ensemble sampling and examining the lag-luminosity
relation. For further information on emcee see https://emcee.readthedocs.io/en/stable/user/sampler/ , 
https://emcee.readthedocs.io/en/stable/tutorials/line/ .
Requires numpy, pylab, emcee, corner.

project-correlation.py : Correlation functions as used in the project. Left as example.

project-ensemblesampler.py : Ensemble sampler functions as used in the project. Left as example.

###

CONTACT - RMP3

I hope this repository provides useful or interesting reading to someone. In the event that anything is confusing,
I'm happy to be reached via jwestonstem@gmail.com . 





"On this cable system we receive over one million channels from the furthest reaches of the galaxy."

"You get H.B.O.?"
