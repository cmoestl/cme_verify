# 3DCORE

This code is concerned with forecast verification of solar storms. It is currently (February 2018) work in progress.

It is distribued with an MIT license.

If this code is used for peer-reviewed scientific publications, the paper below needs to cited and in any case, please also contact me on twitter https://twitter.com/chrisoutofspace or at christian.moestl@oeaw.ac.at.

This code was used for producing the figures and results in the publication:
Möstl et al. (2017) https://arxiv.org/abs/1703.00705
http://cdsads.u-strasbg.fr/abs/2017SpWea..15..955M

## Dependencies
* The code is written in python, and I run it over the command line in MAC OS X (High Sierra) in ipython.

* Install python anaconda (v3 or up, I use 3.5.2).

https://www.anaconda.com/download/#macos

* Add the packages sunpy and seaborn. 

http://docs.sunpy.org/en/stable/guide/installation/index.html

    $ conda config --add channels conda-forge
     
    $ conda install sunpy

    
https://seaborn.pydata.org/installing.html

    $ conda install seaborn    
    

## Running the code
* Download the repository, start an OS X command line and go to the "cme_verify" directory
* Start ipython and run the only code file:

      $ ipython
      
      $ run cme_verify_v1.py
  
* ... 
* ...



