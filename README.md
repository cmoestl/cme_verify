# cme_verify

This code is concerned with forecast verification of solar storms (coronal mass ejections, CMEs) and derivation of statistics of CME parameters from in situ magnetic field and plasma observations by spacecraft such as Wind, STEREO, Venus Express, MESSENGER, and MAVEN.

It is distributed with an MIT license.

If this code is used for peer-reviewed scientific publications, the paper below needs to cited and in any case, please also contact me on twitter https://twitter.com/chrisoutofspace or at christian.moestl@oeaw.ac.at.

The code cme_verify.py was used for producing the figures and results in the publication:
Möstl et al. (2017) https://arxiv.org/abs/1703.00705 
http://cdsads.u-strasbg.fr/abs/2017SpWea..15..955M

On CME statistics, the code cme_stats.py is used, and a paper is in preparation (early 2019).

## Dependencies
* The code is written in python, and I run it over the command line in MAC OS X (High Sierra) in ipython.

* Install python anaconda (v3 or up, I use 3.5.5).

https://www.anaconda.com/download/#macos

* Add the packages sunpy and seaborn. 

http://docs.sunpy.org/en/stable/guide/installation/index.html

    $ conda config --add channels conda-forge
     
    $ conda install sunpy

    
https://seaborn.pydata.org/installing.html

    $ conda install seaborn    
    

## Running the forecast verification code
* Download the repository, start an OS X command line and go to the "cme_verify" directory
* This code has been made for CME forecast verification with an predicted arrival catalog for the years 2007-2014 in the EU funded HELCATS project. Here this verification is done for CME modeling with the SSEF30 method for heliospheric imagers.

Start ipython and run the code file:

      $ ipython
      
      $ run cme_verify

or in the terminal:

      $ python cme_verify.py     
  
## Running the CME statistics code 

* needs sunpy, astropy, seaborn, pickle, urllib, json
* for a paper in preparation:

ipython:
      $ run cme_stats

terminal:      
      $ python cme_stats.py

* for the poster on CME observations with the Parker Solar Probe, shown at EGU 2018
* https://doi.org/10.6084/m9.figshare.6110636.v1

      $ run cme_stats_parker.py
