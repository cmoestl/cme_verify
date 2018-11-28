# cme_verify_v1.py
#
# analyses HELCATS ARRCAT and ICMECAT data for deli 4.2 and subsequent publication:
# "Prediction of solar coronal mass ejections verified with the Heliophysics System Observatory"
# Author: C. Moestl, SRI and University of Graz, Austria
# started 11 May 2015
# last update: May 2017
#

# MIT License
# 
# Copyright (c) 2018 Christian Möstl
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


from scipy import stats
import scipy.io
from matplotlib import cm
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import sunpy.time
import time
import pickle
from functools import reduce
import seaborn as sns


def getcat(filename):
  print('reading CAT')
  cat=scipy.io.readsav(filename, verbose='true')  
  print('done CAT')
  return cat  
    
def decode_array(bytearrin):
 #for decoding the strings from the IDL .sav file to a list of python strings, not bytes 
 #make list of python lists with arbitrary length
 bytearrout= ['' for x in range(len(bytearrin))]
 for i in range(0,len(bytearrin)-1):
  bytearrout[i]=bytearrin[i].decode()
 #has to be np array so to be used with numpy "where"
 bytearrout=np.array(bytearrout)
 return bytearrout  
  
def time_to_num_cat(time_in):  

  #for time conversion from catalogue .sav to numerical time
  #this for 1-minute data or lower time resolution

  #for all catalogues
  #time_in is the time in format: 2007-11-17T07:20:00 or 2007-11-17T07:20Z
  #for times help see: 
  #http://docs.sunpy.org/en/latest/guide/time.html
  #http://matplotlib.org/examples/pylab_examples/date_demo2.html
  
  j=0
  #time_str=np.empty(np.size(time_in),dtype='S19')
  time_str= ['' for x in range(len(time_in))]
  #=np.chararray(np.size(time_in),itemsize=19)
  time_num=np.zeros(np.size(time_in))
  
  for i in time_in:

   #convert from bytes (output of scipy.readsav) to string
   time_str[j]=time_in[j][0:16].decode()+':00'
   year=int(time_str[j][0:4])
   time_str[j]
   #convert time to sunpy friendly time and to matplotlibdatetime
   #only for valid times so 9999 in year is not converted
   #pdb.set_trace()
   if year < 2100:
    	  time_num[j]=mdates.date2num(sunpy.time.parse_time(time_str[j]))
   j=j+1  
   #the date format in matplotlib is e.g. 735202.67569444
   #this is time in days since 0001-01-01 UTC, plus 1.
   
   #return time_num which is already an array and convert the list of strings to an array
  return time_num, np.array(time_str)








######################################################
#main program

plt.close('all')
print('Start catpy main program. Analyses and plots for HELCATS deli 4.2.')




latitude_exclusion=0




#check where target_arrival_num +/- time_window has a hit in ICMECAT

timewindow=1.5 # in days



sns.set_context("talk")     
#sns.set_style("darkgrid")  
sns.set_style("ticks")


#-------------------------------------------------------- get cats
filename_arrcat='ALLCATS/HELCATS_ARRCAT_v6.sav'
a=getcat(filename_arrcat)

filename_icmecat='ALLCATS/HELCATS_ICMECAT_v10_SCEQ.sav'
i=getcat(filename_icmecat)

#filename_linkcat='ALLCATS/HELCATS_LINKCAT_v10.sav'
#l=getcat(filename_linkcat)

filename_higeocat='ALLCATS/HELCATS_HICAT_v03.sav'
h=getcat(filename_higeocat)


#now these are structured arrays  
#access each element of the array see http://docs.scipy.org/doc/numpy/user/basics.rec.html
#access variables
#i.icmecat['id']
#look at contained variables
#print(a.arrcat.dtype)
#print(i.icmecat.dtype)

#get spacecraft and planet positions
pos=getcat('DATACAT/positions_2007_2018_HEEQ_6hours.sav')
pos_time_num=time_to_num_cat(pos.time)[0]


#---------------------------- get all parameters from ICMECAT
iid=i.icmecat['id']
#need to decode all strings
iid=decode_array(iid)

isc=i.icmecat['sc_insitu'] #string
isc=decode_array(isc)

icme_start_time=i.icmecat['ICME_START_TIME']
[icme_start_time_num,icme_start_time_str]=time_to_num_cat(icme_start_time)

mo_start_time=i.icmecat['MO_START_TIME']
[mo_start_time_num,mo_start_time_str]=time_to_num_cat(mo_start_time)

mo_end_time=i.icmecat['MO_END_TIME']
[mo_end_time_num,mo_end_time_str]=time_to_num_cat(mo_end_time)

icme_end_time=i.icmecat['ICME_END_TIME']
[icme_end_time_num,icme_end_time_str]=time_to_num_cat(icme_end_time)

sc_heliodistance=i.icmecat['SC_HELIODISTANCE']
sc_long_heeq=i.icmecat['SC_LONG_HEEQ']
sc_lat_heeq=i.icmecat['SC_LAT_HEEQ']
mo_bmax=i.icmecat['MO_BMAX']
mo_bmean=i.icmecat['MO_BMEAN']
mo_bstd=i.icmecat['MO_BSTD']
mo_bzmean=i.icmecat['MO_BZMEAN']
mo_bzmin=i.icmecat['MO_BZMIN']
mo_duration=i.icmecat['MO_DURATION']
mo_mva_axis_long=i.icmecat['MO_MVA_AXIS_LONG']
mo_mva_axis_lat=i.icmecat['MO_MVA_AXIS_LAT']
mo_mva_ratio=i.icmecat['MO_MVA_RATIO']
sheath_speed=i.icmecat['SHEATH_SPEED']
sheath_speed_std=i.icmecat['SHEATH_SPEED_STD']
mo_speed=i.icmecat['MO_SPEED']
mo_speed_std=i.icmecat['MO_SPEED_STD']
sheath_density=i.icmecat['SHEATH_DENSITY']
sheath_density_std=i.icmecat['SHEATH_DENSITY_STD']
mo_density=i.icmecat['MO_DENSITY']
mo_density_std=i.icmecat['MO_DENSITY_STD']
sheath_temperature=i.icmecat['SHEATH_TEMPERATURE']
sheath_temperature_std=i.icmecat['SHEATH_TEMPERATURE_STD']
mo_temperature=i.icmecat['MO_TEMPERATURE']
mo_temperature_std=i.icmecat['MO_TEMPERATURE_STD']

#get indices of events in different spacecraft
ivexind=np.where(isc == 'VEX')[0]
istaind=np.where(isc == 'STEREO-A')[0]
istbind=np.where(isc == 'STEREO-B')[0]
iwinind=np.where(isc == 'Wind')[0]
imesind=np.where(isc == 'MESSENGER')[0]
iulyind=np.where(isc == 'ULYSSES')[0]



#------------------- get all parameters from ARRCAT
aid=a.arrcat['id']
aid=decode_array(aid)
asc=a.arrcat['sc']
asc=decode_array(asc)

sse_launch=a.arrcat['sse_launch']
#times need to be converted to matplotlib date and string
[sse_launch_num,sse_launch_str]=time_to_num_cat(sse_launch)

target_name=a.arrcat['TARGET_NAME']
#need to decode all strings
target_name=decode_array(target_name)
sse_heeq_long=a.arrcat['sse_heeq_long']

#not included latitude of sse fit
#sse_heeq_lat=a.arrcat['sse_heeq_lat']

target_delta=a.arrcat['TARGET_DELTA']
sse_speed=a.arrcat['SSE_SPEED']
target_speed=a.arrcat['target_speed']
target_arrival=a.arrcat['target_arrival']
#times need to be converted to matplotlib date and string
[target_arrival_num,target_arrival_str]=time_to_num_cat(target_arrival)

target_distance=a.arrcat['target_distance']
target_heeq_lat=a.arrcat['target_heeq_lat']
target_heeq_long=a.arrcat['target_heeq_long']
target_pa=a.arrcat['target_pa']
pa_fit=a.arrcat['pa_fit']
pa_n=a.arrcat['pa_n']
pa_s=a.arrcat['pa_s']
pa_center=a.arrcat['pa_center']


######################

#get hicat parameters

hicat_sc=decode_array(h.hicat['sc'])
hicat_sc_a=np.where(hicat_sc =='A')
hicat_sc_b=np.where(hicat_sc =='B')

######## now get all the indices that match targets and HI spacecraft needed for the calculations

# so lists of indices in ARRCAT of HIA and Venus, HIB and Wind etc. are created

#first get a list of all separated arrivals
arr_hia_obs_ind=np.where(asc == 'A')[0]
arr_hib_obs_ind=np.where(asc == 'B')[0]






############## ** do all this next steps with intersect1d

#get indices of the arrivals of for HIA and HIB each in situ spacecraft separately

#select those events which have Wind as target
arrwinind_all=np.where(target_name == 'EARTH_L1')[0]
#merge both arrays for wind und hia into single array
merge=np.concatenate((arr_hia_obs_ind,arrwinind_all),axis=0)
#bincount checks how many times each value is represented in this array
counts=np.bincount(merge)
#if a value comes up twice the event is both arriving at Wind and observed by HIA
arrwinind_a=np.where(counts > 1)[0]

#one can check this with
#asc[arrwinind_a]
#target_name[arrwinind_a]

#same for HIB
merge=np.concatenate((arr_hib_obs_ind,arrwinind_all),axis=0)
counts=np.bincount(merge)
arrwinind_b=np.where(counts > 1)[0]

#Venus
arrvexind_all=np.where(target_name  == 'VENUS')[0]
merge=np.concatenate((arr_hia_obs_ind,arrvexind_all),axis=0)
counts=np.bincount(merge)
arrvexind_a=np.where(counts > 1)[0]
merge=np.concatenate((arr_hib_obs_ind,arrvexind_all),axis=0)
counts=np.bincount(merge)
arrvexind_b=np.where(counts > 1)[0]

#MESSENGER      
arrmesind_all=np.where(target_name  == 'MESSENGER')[0]
merge=np.concatenate((arr_hia_obs_ind,arrmesind_all),axis=0)
counts=np.bincount(merge)
arrmesind_a=np.where(counts > 1)[0]
merge=np.concatenate((arr_hib_obs_ind,arrmesind_all),axis=0)
counts=np.bincount(merge)
arrmesind_b=np.where(counts > 1)[0]

#STEREO-B
arrstbind_all=np.where(target_name  == 'STEREO-B')[0]
merge=np.concatenate((arr_hia_obs_ind,arrstbind_all),axis=0)
counts=np.bincount(merge)
arrstbind_a=np.where(counts > 1)[0]
merge=np.concatenate((arr_hib_obs_ind,arrstbind_all),axis=0)
counts=np.bincount(merge)
arrstbind_b=np.where(counts > 1)[0]

#STEREO-A
arrstaind_all=np.where(target_name  == 'STEREO-A')[0]
merge=np.concatenate((arr_hia_obs_ind,arrstaind_all),axis=0)
counts=np.bincount(merge)
arrstaind_a=np.where(counts > 1)[0]
merge=np.concatenate((arr_hib_obs_ind,arrstaind_all),axis=0)
counts=np.bincount(merge)
arrstaind_b=np.where(counts > 1)[0]


print()
print()
print()
print()
print()
print()



############ latitude exclusion on/off


if latitude_exclusion == 1: 

  print('latitude exclusion on:')
  
  ####################### latitude
  
  #degree +/-270 and 90
  pa_window=20
  
  print('PA window of +/- ', pa_window, ' around 90 deg (HIA) and 270 deg (HIB)' )


  #get all CMEs with a central PA around 90 +/- 10 for Ahead
  #***************** is this correct? one should take the pa_center of all hia_obs_ind

  lat_a=np.where(np.logical_and(pa_center < 90+pa_window,pa_center > 90-pa_window))[0]
  
  #lat_fit_a=np.where(np.logical_and(pa_fit[arr_hia_obs_ind] < 90+pa_window,pa_fit[arr_hia_obs_ind] > 90-pa_window))
  #lat_fit_b=np.where(np.logical_and(pa_fit[arr_hib_obs_ind] < 270+pa_window,pa_fit[arr_hib_obs_ind] > 270-pa_window))


   #same for B with PA of 270
  lat_b=np.where(np.logical_and(pa_center < 270+pa_window,pa_center > 270-pa_window))[0]

  #HIA
  arrwinind_a=np.intersect1d(arrwinind_a,lat_a)
  arrvexind_a=np.intersect1d(arrvexind_a,lat_a)
  arrmesind_a=np.intersect1d(arrmesind_a,lat_a)
  arrstaind_a=np.intersect1d(arrstaind_a,lat_a)
  arrstbind_a=np.intersect1d(arrstbind_a,lat_a)

  #HIB
  arrwinind_b=np.intersect1d(arrwinind_b,lat_b)
  arrvexind_b=np.intersect1d(arrvexind_b,lat_b)
  arrmesind_b=np.intersect1d(arrmesind_b,lat_b)
  arrstaind_b=np.intersect1d(arrstaind_b,lat_b)
  arrstbind_b=np.intersect1d(arrstbind_b,lat_b)




print()



######################## OVERALL STATISTICS ###############################################

print()
print('------------------------------- STATS ----------------------')
print()

print('predicted arrivals HIA vs. in situ ICME detections' )
print('Wind ',np.size(arrwinind_a),'   ',np.size(iwinind))
print('VEX  ',np.size(arrvexind_a),'   ',np.size(ivexind))
print('MES  ',np.size(arrmesind_a),'   ',np.size(imesind))
print('STB  ',np.size(arrstbind_a),'   ',np.size(istbind))
print('STA  ',np.size(arrstaind_a),'   ',np.size(istaind))

print('predicted arrivals HIB vs. in situ ICME detections' )
print('Wind ',np.size(arrwinind_b),'   ',np.size(iwinind))
print('VEX  ',np.size(arrvexind_b),'   ',np.size(ivexind))
print('MES  ',np.size(arrmesind_b),'   ',np.size(imesind))
print('STB  ',np.size(arrstbind_b),'   ',np.size(istbind))
print('STA  ',np.size(arrstaind_b),'   ',np.size(istaind))

print('However, note that VEX/MES go in and out of the HI FoVs.')
print('Also there are less predictions than arrivals at both STEREO - check confusion matrix')

print('ARRCAT total number of arrivals at all spacecraft from HIA',np.size(arr_hia_obs_ind))
print('ARRCAT total number of arrivals at all spacecraft from HIB',np.size(arr_hib_obs_ind))

numberat5_a=np.size(arrwinind_a)+np.size(arrvexind_a)+np.size(arrmesind_a)+np.size(arrstbind_a)+np.size(arrstaind_a)
numberat5_b=np.size(arrwinind_b)+np.size(arrvexind_b)+np.size(arrmesind_b)+np.size(arrstbind_b)+np.size(arrstaind_b)

print()

print('ARRCAT total number of arrivals at 5 in situ spacecraft from HIA',numberat5_a)
print('ARRCAT total number of arrivals at 5 in situ spacecraft from HIB',numberat5_b)










#################### (1) HIT AND MISS PREDICTIONS ######################################## 


print()
print()
print()
print('------------------------------------ Hit and miss calculation -------------------')
print()



print('Search for correct hits predicted by STEREO-A HI for all 5 in situ spacecraft:')

#this contains the number of correctly predicted impacts; 0 for false alarm
winhits_a=np.zeros(np.size(arrwinind_a))
vexhits_a=np.zeros(np.size(arrvexind_a))
stahits_a=np.zeros(np.size(arrstaind_a))
stbhits_a=np.zeros(np.size(arrstbind_a))
meshits_a=np.zeros(np.size(arrmesind_a))

#these are arrays for the arrival time difference for each (single) hit,0 if no hit
winhits_time_diff_a=np.zeros(np.size(winhits_a))
vexhits_time_diff_a=np.zeros(np.size(vexhits_a))
stahits_time_diff_a=np.zeros(np.size(stahits_a))
stbhits_time_diff_a=np.zeros(np.size(stbhits_a))
meshits_time_diff_a=np.zeros(np.size(meshits_a))

#these are arrays for the arrival speed difference for each (single) hit,0 if no hit
winhits_speed_diff_a=np.zeros(np.size(winhits_a))
stahits_speed_diff_a=np.zeros(np.size(stahits_a))
stbhits_speed_diff_a=np.zeros(np.size(stbhits_a))



#make an array with pairs of ARRCAT and ICMECAT for which the hit is correct
overlap_a_hi=np.zeros(np.size(iid))
overlap_a_icme=np.zeros(np.size(iid))
overcount=0

kwin=0
kvex=0
ksta=0
kstb=0
kmes=0


#if latitude exclusion is on, need to select only those arrcat events with correct PA
#instead of all HIA entries in ARRCAT
if latitude_exclusion == 1: 
  arr_hia_obs_ind=lat_a

#same for B for the loop after the next
if latitude_exclusion == 1: 
  arr_hib_obs_ind=lat_b


#go through all arrivals predicted by STEREO-A

for k in np.arange(0,np.size(arr_hia_obs_ind)-1):

 #check for each STA predicted arrival ...
 this_target_arrival_num=target_arrival_num[arr_hia_obs_ind][k]
 this_target_speed=target_speed[arr_hia_obs_ind][k]
 this_target_name=target_name[arr_hia_obs_ind][k]
 
 
 #search all in situ data for the corresponding arrival time for the correct spacecraft:

 if this_target_name == 'EARTH_L1':
   #calculate the time difference to each arrival at Wind from the in situ list as calculated-observed
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[iwinind]
   #same for the speed C-O
   speed_diff_to_all_insitu=this_target_speed-sheath_speed[iwinind]
   
   #search for all correctly predicted hits, defined by the time difference being less than the time window, 
   #in absolute numbers
   #kwin is the index of all Wind correctly predicted hits
   winhits_a[kwin]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))

   
   #save the time and speed difference for plot histogram 
   if winhits_a[kwin] > 0:
   
      insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
      
      #print()
      #print(time_diff_to_all_insitu[insitu_index])
      #print(time_diff_to_all_insitu[insitu_index[0][0]])
      #print(icme_start_time_str[iwinind][insitu_index])
      
      #if there is more than 1 ICME within the time window only take the closer one for arrival time / speed comparison
      index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))      
      
      winhits_time_diff_a[kwin]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
      #print(winhits_time_diff_a[kwin])

      #this is the ICME speed of the corresponding in situ arrival
      winhits_speed_diff_a[kwin]=speed_diff_to_all_insitu[insitu_index[0][index_closer]] 
      
      #make an entry for k which is index of arr_hia_obs_ind to get to ARRCAT event
      overlap_a_hi[overcount]=arr_hia_obs_ind[k]
      #and for the in situ index of the event in ICMECAT - take the index from iwinind
      overlap_a_icme[overcount]=iwinind[insitu_index[0][index_closer]]
      overcount=overcount+1                  
       
      #for checking if the above lines are correct
      #print(icme_start_time_str[iwinind][insitu_index])
      #print(target_arrival_str[arr_hia_obs_ind][k])
      #print(this_target_speed)
      #print(sheath_speed[iwinind][insitu_index])
      #print(winhits_speed_diff_a[kwin])
      #print()

   kwin=kwin+1



 if this_target_name == 'VENUS':
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[ivexind]
   vexhits_a[kvex]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))
   if vexhits_a[kvex] > 0:
      insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
      index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))      
      vexhits_time_diff_a[kvex]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
      overlap_a_hi[overcount]=arr_hia_obs_ind[k]
      overlap_a_icme[overcount]=ivexind[insitu_index[0][index_closer]]
      overcount=overcount+1             
   kvex=kvex+1

 if this_target_name == 'MESSENGER':
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[imesind]
   meshits_a[kmes]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))
   if meshits_a[kmes] > 0: 
        insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
        index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))      
        meshits_time_diff_a[kmes]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
        overlap_a_hi[overcount]=arr_hia_obs_ind[k]
        overlap_a_icme[overcount]=insitu_index[0][index_closer]
        overlap_a_icme[overcount]=imesind[insitu_index[0][index_closer]]
        overcount=overcount+1             
   kmes=kmes+1

 if this_target_name == 'STEREO-B':
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[istbind]
   speed_diff_to_all_insitu=this_target_speed-sheath_speed[istbind]   
   stbhits_a[kstb]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))
   if stbhits_a[kstb] > 0:
        insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
        index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))      
        stbhits_time_diff_a[kstb]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
        stbhits_speed_diff_a[kstb]=speed_diff_to_all_insitu[insitu_index[0][index_closer]]      
        overlap_a_hi[overcount]=arr_hia_obs_ind[k]
        overlap_a_icme[overcount]=istbind[insitu_index[0][index_closer]]
        overcount=overcount+1             
   kstb=kstb+1

 if this_target_name == 'STEREO-A':
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[istaind]
   speed_diff_to_all_insitu=this_target_speed-sheath_speed[istaind]
   stahits_a[ksta]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))
   if stahits_a[ksta] > 0:
        insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
        index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))      
        stahits_time_diff_a[ksta]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
        stahits_speed_diff_a[ksta]=speed_diff_to_all_insitu[insitu_index[0][index_closer]]      
        overlap_a_hi[overcount]=arr_hia_obs_ind[k]
        overlap_a_icme[overcount]=istaind[insitu_index[0][index_closer]]        
        overcount=overcount+1             
   ksta=ksta+1



percentage_double_hits_a=np.size(np.where(winhits_a>1))/np.size(np.where(winhits_a>0))*100
print('Percentage of double hits at Wind',percentage_double_hits_a)
   
print('time window in days',timewindow)   
print('Wind: HI predicted impacts with clear ICME percent:', round(np.size(np.where(winhits_a > 0))/np.size(winhits_a)*100))
print('VEX: HI predicted impacts percent:', round(np.size(np.where(vexhits_a > 0))/np.size(vexhits_a)*100))
print('MESSENGER: HI predicted impacts percent:', round(np.size(np.where(meshits_a > 0))/np.size(meshits_a)*100))
print('STA: HI predicted impacts percent:', round(np.size(np.where(stahits_a > 0))/np.size(stahits_a)*100))
print('STB: HI predicted impacts percent:', round(np.size(np.where(stbhits_a > 0))/np.size(stbhits_a)*100))



#calculate hit percentage overall for HIA

#how many hits in all arrays
icme_hits_a=np.size(np.where(np.concatenate((winhits_a,vexhits_a,meshits_a,stahits_a,stbhits_a), axis=0) > 0))

#total number of predicted impacts
pred_hits_a=np.size(np.concatenate((winhits_a,vexhits_a,meshits_a,stahits_a,stbhits_a), axis=0))

print('HIA: all spacecraft total number of predicted impacts: ',pred_hits_a)
print('HIA: all spacecraft HI predicted impacts with clear ICME in percent: ', round(icme_hits_a/pred_hits_a*100) )


print()
print()
print()



################################
print('Search for correct hits for STEREO-B HI for all 5 in situ spacecraft:')

#this contains the number of correctly predicted impacts; 0 for false alarm
winhits_b=np.zeros(np.size(arrwinind_b))
vexhits_b=np.zeros(np.size(arrvexind_b))
stahits_b=np.zeros(np.size(arrstaind_b))
stbhits_b=np.zeros(np.size(arrstbind_b))
meshits_b=np.zeros(np.size(arrmesind_b))


#these are arrays for the arrival time difference for each (single) hit,0 if no hit
winhits_time_diff_b=np.zeros(np.size(winhits_b))
vexhits_time_diff_b=np.zeros(np.size(vexhits_b))
stahits_time_diff_b=np.zeros(np.size(stahits_b))
stbhits_time_diff_b=np.zeros(np.size(stbhits_b))
meshits_time_diff_b=np.zeros(np.size(meshits_b))

#these are arrays for the arrival speed difference for each (single) hit,0 if no hit
winhits_speed_diff_b=np.zeros(np.size(winhits_b))
stahits_speed_diff_b=np.zeros(np.size(stahits_b))
stbhits_speed_diff_b=np.zeros(np.size(stbhits_b))

#make an array with pairs of ARRCAT and ICMECAT for which the hit is correct
overlap_b_hi=np.zeros(np.size(iid))
overlap_b_icme=np.zeros(np.size(iid))
overcount=0



kwin=0
kvex=0
ksta=0
kstb=0
kmes=0


#go through all arrivals predicted by STB
for k in np.arange(0,np.size(arr_hib_obs_ind)-1):

 this_target_arrival_num=target_arrival_num[arr_hib_obs_ind][k]
 this_target_speed=target_speed[arr_hib_obs_ind][k]
 this_target_name=target_name[arr_hib_obs_ind][k]
 
 if this_target_name == 'EARTH_L1':
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[iwinind]
   speed_diff_to_all_insitu=this_target_speed-sheath_speed[iwinind]
   winhits_b[kwin]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))
   if winhits_b[kwin] > 0:   
      insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
      index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))       
      winhits_time_diff_b[kwin]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
      winhits_speed_diff_b[kwin]=speed_diff_to_all_insitu[insitu_index[0][index_closer]]      
      overlap_b_hi[overcount]=arr_hib_obs_ind[k]
      overlap_b_icme[overcount]=iwinind[insitu_index[0][index_closer]]
      overcount=overcount+1             
   kwin=kwin+1

 if this_target_name == 'VENUS':
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[ivexind]
   vexhits_b[kvex]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))
   if vexhits_b[kvex] > 0:
      insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
      index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))       
      vexhits_time_diff_b[kvex]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
      overlap_b_hi[overcount]=arr_hib_obs_ind[k]
      overlap_b_icme[overcount]=ivexind[insitu_index[0][index_closer]]
      overcount=overcount+1             
   kvex=kvex+1

 if this_target_name == 'MESSENGER':
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[imesind]
   meshits_b[kmes]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))
   if meshits_b[kmes] > 0: 
        insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
        index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))       
        meshits_time_diff_b[kmes]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
        overlap_b_hi[overcount]=arr_hib_obs_ind[k]
        overlap_b_icme[overcount]=imesind[insitu_index[0][index_closer]]
        overcount=overcount+1             
   kmes=kmes+1

 if this_target_name == 'STEREO-B':
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[istbind]
   speed_diff_to_all_insitu=this_target_speed-sheath_speed[istbind]   
   stbhits_b[kstb]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))
   if stbhits_b[kstb] > 0:
        insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
        index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))       
        stbhits_time_diff_b[kstb]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
        stbhits_speed_diff_b[kstb]=speed_diff_to_all_insitu[insitu_index[0][index_closer]]      
        overlap_b_hi[overcount]=arr_hib_obs_ind[k]
        overlap_b_icme[overcount]=istbind[insitu_index[0][index_closer]]
        overcount=overcount+1             
   kstb=kstb+1

 if this_target_name == 'STEREO-A':
   time_diff_to_all_insitu=this_target_arrival_num-icme_start_time_num[istaind]
   speed_diff_to_all_insitu=this_target_speed-sheath_speed[istaind]
   stahits_b[ksta]=np.size(np.where(abs(time_diff_to_all_insitu) < timewindow))
   if stahits_b[ksta] > 0: 
        insitu_index=np.where(abs(time_diff_to_all_insitu) < timewindow)
        index_closer=np.argmin(abs(time_diff_to_all_insitu[insitu_index]))       
        stahits_time_diff_b[ksta]=time_diff_to_all_insitu[insitu_index[0][index_closer]]
        stahits_speed_diff_b[ksta]=speed_diff_to_all_insitu[insitu_index[0][index_closer]]      
        overlap_b_hi[overcount]=arr_hib_obs_ind[k]
        overlap_b_icme[overcount]=istaind[insitu_index[0][index_closer]]
        overcount=overcount+1             
   ksta=ksta+1


percentage_double_hits_b=np.size(np.where(winhits_b>1))/np.size(np.where(winhits_b>0))*100
print('Percentage of double hits HIB at Wind',percentage_double_hits_b)

   
print('time window in days',timewindow)   
print('Wind: HI predicted impacts with clear ICME percent:', round(np.size(np.where(winhits_b > 0))/np.size(winhits_b)*100))
print('VEX: HI predicted impacts percent:', round(np.size(np.where(vexhits_b > 0))/np.size(vexhits_b)*100))
print('MESSENGER: HI predicted impacts percent:', round(np.size(np.where(meshits_b > 0))/np.size(meshits_b)*100))
print('STA: HI predicted impacts percent:', round(np.size(np.where(stahits_b > 0))/np.size(stahits_b)*100))
print('STB: HI predicted impacts percent:', round(np.size(np.where(stbhits_b > 0))/np.size(stbhits_b)*100))


#calculate hit percentage overall for HIB

#how many hits in all arrays
icme_hits_b=np.size(np.where(np.concatenate((winhits_b,vexhits_b,meshits_b,stahits_b,stbhits_b), axis=0) > 0))

#total number of predicted impacts
pred_hits_b=np.size(np.concatenate((winhits_b,vexhits_b,meshits_b,stahits_b,stbhits_b), axis=0))

print('HIB: all spacecraft total number of predicted impacts: ',pred_hits_b)
print('HIB: all spacecraft HI predicted impacts with clear ICME in percent: ', round(icme_hits_b/pred_hits_b*100) )

print()
print()
print()



############################# PLOTS for hit and miss ############################################


##################### (1) arrival frequencies in ICMECAT and ARRCAT
fsize=10
matplotlib.rcParams.update({'font.size': fsize})

yearly_bin_edges=[mdates.date2num(sunpy.time.parse_time('2007-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2008-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2009-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2010-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2011-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2012-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2013-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2014-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2015-01-01'))]

#bin width in days         
binweite=360/6

plt.figure(1)
ax1 = plt.subplot(211) 
plt.ylabel('Number of ARRCAT impacts HIA/HIB', fontsize=8)

(histwin, bin_edgeswin) = np.histogram(target_arrival_num[arrwinind_all], yearly_bin_edges)
(histvex, bin_edgesvex) = np.histogram(target_arrival_num[arrvexind_all], yearly_bin_edges)
(histmes, bin_edgesmes) = np.histogram(target_arrival_num[arrmesind_all], yearly_bin_edges)
(histstb, bin_edgesstb) = np.histogram(target_arrival_num[arrstbind_all], yearly_bin_edges)
(histsta, bin_edgessta) = np.histogram(target_arrival_num[arrstaind_all], yearly_bin_edges)

#bin_edgeswin is cut off because there is no data in 2015 (ARRCAT stops before STEREO conjunction)
ax1.bar(bin_edgeswin[:-1]+30,histwin, width=binweite,color='mediumseagreen', alpha=0.5)
ax1.bar(bin_edgesvex[:-1]+30+binweite,histvex, width=binweite,color='orange', alpha=0.5)
ax1.bar(bin_edgesmes[:-1]+30+binweite*2,histmes, width=binweite,color='darkgrey', alpha=0.5)
ax1.bar(bin_edgesstb[:-1]+30+binweite*3,histstb, width=binweite,color='royalblue', alpha=0.5)
ax1.bar(bin_edgessta[:-1]+30+binweite*4,histsta, width=binweite,color='red', alpha=0.5)
plt.xlim(yearly_bin_edges[0],yearly_bin_edges[8])
#plt.
ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)

fsize=10 #same for all elements

#panel labels
plt.figtext(0.03,0.95,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.03,0.45,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')

plt.xticks(fontsize=fsize)
plt.yticks(fontsize=fsize)

#sets planet / spacecraft labels
xoff=0.15
yoff=0.87

plt.figtext(xoff,yoff,'Earth L1',color='mediumseagreen', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.05*1,'VEX',color='orange', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.05*2,'MESSENGER',color='dimgrey', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.05*3,'STEREO-A',color='red', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.05*4,'STEREO-B',color='royalblue', fontsize=fsize, ha='left')

ax2 = plt.subplot(212) 
plt.ylabel('Number of ICMECAT events', fontsize=8)
plt.xlabel('Year', fontsize=fsize)
(histwin, bin_edgeswin) = np.histogram(icme_start_time_num[iwinind], yearly_bin_edges)
(histvex, bin_edgesvex) = np.histogram(icme_start_time_num[ivexind], yearly_bin_edges)
(histmes, bin_edgesmes) = np.histogram(icme_start_time_num[imesind], yearly_bin_edges)
(histstb, bin_edgesstb) = np.histogram(icme_start_time_num[istbind], yearly_bin_edges)
(histsta, bin_edgessta) = np.histogram(icme_start_time_num[istaind], yearly_bin_edges)

ax2.bar(bin_edgeswin[:-1]+30,histwin, width=binweite,color='mediumseagreen', alpha=0.5)
ax2.bar(bin_edgesvex[:-1]+30+binweite,histvex, width=binweite,color='orange', alpha=0.5)
ax2.bar(bin_edgesmes[:-1]+30+ binweite*2,histmes, width=binweite,color='darkgrey', alpha=0.5)
ax2.bar(bin_edgesstb[:-1]+30+binweite*3,histstb, width=binweite,color='royalblue', alpha=0.5)
ax2.bar(bin_edgessta[:-1]+30+binweite*4,histsta, width=binweite,color='red', alpha=0.5)
plt.xlim(yearly_bin_edges[0],yearly_bin_edges[8])
ax2.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax2.xaxis.set_major_formatter(myformat)

plt.xticks(fontsize=fsize)
plt.yticks(fontsize=fsize)


sns.despine()

plt.savefig('plots/arrcat_icmecat_frequency_per_year.pdf', dpi=300)
plt.savefig('plots/arrcat_icmecat_frequency_per_year.png', dpi=300)








################### (2) correct hits and false alarms ratio per planet per year


#get the times of all hits for the "ghost bars"


#these are all predicted times by HIA
win_all_pred_times_a=target_arrival_num[arrwinind_a]
#get all the times of the correct hits
win_hit_times_a=(winhits_a >0)*target_arrival_num[arrwinind_a]
#clean out all the zeros in this array
win_hit_times_clean_a=win_hit_times_a[np.where(win_hit_times_a != 0)]
#percentage of correct predicted hits, rest is false alarms
win_perc_hit_a=np.size(win_hit_times_clean_a)/np.size(win_hit_times_a)*100
#this is done with times so we can plot them in the yearly histogram below

vex_all_pred_times_a=target_arrival_num[arrvexind_a]
vex_hit_times_a=(vexhits_a >0)*target_arrival_num[arrvexind_a]
vex_hit_times_clean_a=vex_hit_times_a[np.where(vex_hit_times_a != 0)]
vex_perc_hit_a=np.size(vex_hit_times_clean_a)/np.size(vex_hit_times_a)*100

mes_all_pred_times_a=target_arrival_num[arrmesind_a]
mes_hit_times_a=(meshits_a >0)*target_arrival_num[arrmesind_a]
mes_hit_times_clean_a=mes_hit_times_a[np.where(mes_hit_times_a != 0)]
mes_perc_hit_a=np.size(mes_hit_times_clean_a)/np.size(mes_hit_times_a)*100

stb_all_pred_times_a=target_arrival_num[arrstbind_a]
stb_hit_times_a=(stbhits_a >0)*target_arrival_num[arrstbind_a]
stb_hit_times_clean_a=stb_hit_times_a[np.where(stb_hit_times_a != 0)]
stb_perc_hit_a=np.size(stb_hit_times_clean_a)/np.size(stb_hit_times_a)*100

sta_all_pred_times_a=target_arrival_num[arrstaind_a]
sta_hit_times_a=(stahits_a >0)*target_arrival_num[arrstaind_a]
sta_hit_times_clean_a=sta_hit_times_a[np.where(sta_hit_times_a != 0)]
sta_perc_hit_a=np.size(sta_hit_times_clean_a)/np.size(sta_hit_times_a)*100




binweite=360/6

#make histogram for each year

plt.figure(2)
ax1 = plt.subplot(211) 
plt.ylabel('HIA number of hits',fontsize=fsize)
#these are the correct hits
(histwin, bin_edgeswin) = np.histogram(win_hit_times_clean_a, yearly_bin_edges)
#these are all predicted arrivals
(histwin_all, bin_edgeswin) = np.histogram(win_all_pred_times_a, yearly_bin_edges)
(histvex, bin_edgesvex) = np.histogram(vex_hit_times_clean_a, yearly_bin_edges)
(histvex_all, bin_edgeswin) = np.histogram(vex_all_pred_times_a, yearly_bin_edges)
(histmes, bin_edgesmes) = np.histogram(mes_hit_times_clean_a, yearly_bin_edges)
(histmes_all, bin_edgeswin) = np.histogram(mes_all_pred_times_a, yearly_bin_edges)
(histstb, bin_edgesstb) = np.histogram(stb_hit_times_clean_a, yearly_bin_edges)
(histstb_all, bin_edgeswin) = np.histogram(stb_all_pred_times_a, yearly_bin_edges)
(histsta, bin_edgessta) = np.histogram(sta_hit_times_clean_a,yearly_bin_edges)
(histsta_all, bin_edgeswin) = np.histogram(sta_all_pred_times_a, yearly_bin_edges)


ax1.bar(bin_edgeswin[:-1]+30,histwin_all, width=binweite,color='mediumseagreen', alpha=0.3)
ax1.bar(bin_edgeswin[:-1]+30,histwin, width=binweite,color='mediumseagreen', alpha=0.7)

ax1.bar(bin_edgesvex[:-1]+30+binweite,histvex_all, width=binweite,color='orange', alpha=0.3)
ax1.bar(bin_edgesvex[:-1]+30+binweite,histvex, width=binweite,color='orange', alpha=0.7)

ax1.bar(bin_edgesmes[:-1]+30+ binweite*2,histmes_all, width=binweite,color='darkgrey', alpha=0.3)
ax1.bar(bin_edgesmes[:-1]+30+ binweite*2,histmes, width=binweite,color='darkgrey', alpha=0.7)

ax1.bar(bin_edgesstb[:-1]+30+binweite*3,histstb_all, width=binweite,color='royalblue', alpha=0.3)
ax1.bar(bin_edgesstb[:-1]+30+binweite*3,histstb, width=binweite,color='royalblue', alpha=0.7)

ax1.bar(bin_edgessta[:-1]+30+binweite*4,histsta_all, width=binweite,color='red', alpha=0.3)
ax1.bar(bin_edgessta[:-1]+30+binweite*4,histsta, width=binweite,color='red', alpha=0.7)

plt.xlim(yearly_bin_edges[0],yearly_bin_edges[8])

ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)


plt.xticks(fontsize=fsize)
plt.yticks(fontsize=fsize)

fsize=10
xoff=0.15
yoff=0.87

plt.figtext(xoff,yoff,'Earth L1',color='mediumseagreen', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.05*1,'VEX',color='orange', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.05*2,'MESSENGER',color='dimgrey', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.05*3,'STEREO-A',color='red', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.05*4,'STEREO-B',color='royalblue', fontsize=fsize, ha='left')

plt.figtext(0.03,0.90,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.03,0.45,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')


##########same for HIB

win_all_pred_times_b=target_arrival_num[arrwinind_b]
win_hit_times_b=(winhits_b >0)*target_arrival_num[arrwinind_b]
win_hit_times_clean_b=win_hit_times_b[np.where(win_hit_times_b != 0)]
win_perc_hit_b=np.size(win_hit_times_clean_b)/np.size(win_hit_times_b)*100

vex_all_pred_times_b=target_arrival_num[arrvexind_b]
vex_hit_times_b=(vexhits_b >0)*target_arrival_num[arrvexind_b]
vex_hit_times_clean_b=vex_hit_times_b[np.where(vex_hit_times_b != 0)]
vex_perc_hit_b=np.size(vex_hit_times_clean_b)/np.size(vex_hit_times_b)*100

mes_all_pred_times_b=target_arrival_num[arrmesind_b]
mes_hit_times_b=(meshits_b >0)*target_arrival_num[arrmesind_b]
mes_hit_times_clean_b=mes_hit_times_b[np.where(mes_hit_times_b != 0)]
mes_perc_hit_b=np.size(mes_hit_times_clean_b)/np.size(mes_hit_times_b)*100

stb_all_pred_times_b=target_arrival_num[arrstbind_b]
stb_hit_times_b=(stbhits_b >0)*target_arrival_num[arrstbind_b]
stb_hit_times_clean_b=stb_hit_times_b[np.where(stb_hit_times_b != 0)]
stb_perc_hit_b=np.size(stb_hit_times_clean_b)/np.size(stb_hit_times_b)*100

sta_all_pred_times_b=target_arrival_num[arrstaind_b]
sta_hit_times_b=(stahits_b >0)*target_arrival_num[arrstaind_b]
sta_hit_times_clean_b=sta_hit_times_b[np.where(sta_hit_times_b != 0)]
sta_perc_hit_b=np.size(sta_hit_times_clean_b)/np.size(sta_hit_times_b)*100




ax1 = plt.subplot(212) 
plt.ylabel('HIB number of hits', fontsize=fsize)
plt.xlabel('Year', fontsize=fsize)
(histwin, bin_edgeswin) = np.histogram(win_hit_times_clean_b, yearly_bin_edges)
(histwin_all, bin_edgeswin) = np.histogram(win_all_pred_times_b, yearly_bin_edges)
(histvex, bin_edgesvex) = np.histogram(vex_hit_times_clean_b, yearly_bin_edges)
(histvex_all, bin_edgeswin) = np.histogram(vex_all_pred_times_b, yearly_bin_edges)
(histmes, bin_edgesmes) = np.histogram(mes_hit_times_clean_b, yearly_bin_edges)
(histmes_all, bin_edgeswin) = np.histogram(mes_all_pred_times_b, yearly_bin_edges)
(histstb, bin_edgesstb) = np.histogram(stb_hit_times_clean_b, yearly_bin_edges)
(histstb_all, bin_edgeswin) = np.histogram(stb_all_pred_times_b, yearly_bin_edges)
(histsta, bin_edgessta) = np.histogram(sta_hit_times_clean_b,yearly_bin_edges)
(histsta_all, bin_edgeswin) = np.histogram(sta_all_pred_times_b, yearly_bin_edges)


ax1.bar(bin_edgeswin[:-1]+30,histwin_all, width=binweite,color='mediumseagreen', alpha=0.3)
ax1.bar(bin_edgeswin[:-1]+30,histwin, width=binweite,color='mediumseagreen', alpha=0.7)

ax1.bar(bin_edgesvex[:-1]+30+binweite,histvex_all, width=binweite,color='orange', alpha=0.3)
ax1.bar(bin_edgesvex[:-1]+30+binweite,histvex, width=binweite,color='orange', alpha=0.7)

ax1.bar(bin_edgesmes[:-1]+30+ binweite*2,histmes_all, width=binweite,color='darkgrey', alpha=0.3)
ax1.bar(bin_edgesmes[:-1]+30+ binweite*2,histmes, width=binweite,color='darkgrey', alpha=0.7)

ax1.bar(bin_edgesstb[:-1]+30+binweite*3,histstb_all, width=binweite,color='royalblue', alpha=0.3)
ax1.bar(bin_edgesstb[:-1]+30+binweite*3,histstb, width=binweite,color='royalblue', alpha=0.7)

ax1.bar(bin_edgessta[:-1]+30+binweite*4,histsta_all, width=binweite,color='red', alpha=0.3)
ax1.bar(bin_edgessta[:-1]+30+binweite*4,histsta, width=binweite,color='red', alpha=0.7)

plt.xlim(yearly_bin_edges[0],yearly_bin_edges[8])

ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)

fsize=10
xoff=0.15
yoff=0.87

sns.despine()

plt.savefig('plots/arrcat_icmecat_correct_hits.pdf',  dpi=300)
plt.savefig('plots/arrcat_icmecat_correct_hits.png', dpi=300)









################### (3) ratio of correct hits at Earth as function of STEREO separation from Earth
# STEREO separation to Earth is not in ARRCAT, 
# so plot against time and positions from position file

fsize=8
#histogram for each year of the correct hits
(histwin_hits, bin_edgeswin) = np.histogram(win_hit_times_clean_a, yearly_bin_edges)
(histwin_preds, bin_edgeswin) = np.histogram(win_all_pred_times_a, yearly_bin_edges)
yearly_ratio_hia=histwin_hits/histwin_preds*100

(histwin_hits, bin_edgeswin) = np.histogram(win_hit_times_clean_b, yearly_bin_edges)
(histwin_preds, bin_edgeswin) = np.histogram(win_all_pred_times_b, yearly_bin_edges)
yearly_ratio_hib=histwin_hits/histwin_preds*100

plt.figure(3)
ax1 = plt.subplot(211) 
plt.ylabel('Correct hits at Earth L1 [%]',fontsize=fsize)
plt.bar(bin_edgeswin[:-1]+30,yearly_ratio_hib, width=100,color='royalblue', alpha=0.7,label='HIB')
plt.bar(bin_edgeswin[:-1]+160,yearly_ratio_hia, width=100,color='red', alpha=0.7,label='HIA')



ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)
ax1.legend(loc=2, fontsize=fsize)
plt.xlim(yearly_bin_edges[0],yearly_bin_edges[8])
plt.xticks(fontsize=fsize)
plt.yticks(fontsize=fsize)

ax2 = plt.subplot(212) 
plt.ylabel('HI separation from Earth [deg]',fontsize=fsize)
plt.xlabel('Year',fontsize=fsize)
plt.plot_date(pos_time_num,abs(pos['stb'][1]*180/np.pi),'-',color='royalblue',label='HIB',alpha=0.7)
plt.plot_date(pos_time_num,abs(pos['sta'][1]*180/np.pi),'-',color='red',label='HIA',lw=2,alpha=0.7)
ax2.legend(loc=2,fontsize=fsize)
plt.xlim(yearly_bin_edges[0],yearly_bin_edges[8])

plt.figtext(0.01,0.90,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.01,0.45,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')


#check at which index STEREO-B was at closest to L5 at 60° HEEQ longitude from Earth
l5_index=np.argmin(abs(pos['stb'][1]*180/np.pi+60))
plt.axvline(x=pos_time_num[l5_index], ymin=0, ymax = 180, linewidth=2, color='royalblue')
plt.figtext(0.4,0.4,'L5', color='royalblue', fontsize=fsize, ha='left')

plt.xticks(fontsize=fsize)
plt.yticks(fontsize=fsize)

plt.savefig('plots/correct_hits_time_separation.pdf', dpi=300)
plt.savefig('plots/correct_hits_time_separation.png', dpi=300)



print()





################################ full contingency table and skill scores






#############################(1) Earth L1

#all predicted arrivals
np.size(arrwinind_a)

#1 part of contingency table:
#0 for false alarm, 1 or 2 for hits
np.size(winhits_a)

#HIA
false_alarms_a=np.size(np.where(winhits_a ==0))
hits_a=np.size(np.where(winhits_a > 0))

#HIB
false_alarms_b=np.size(np.where(winhits_b ==0))
hits_b=np.size(np.where(winhits_b > 0))

#all Wind in situ ICMEs
np.size(iwinind)
winstart=icme_start_time_num[iwinind]

############ define false rejection: no prediction, but in situ detection
#this means looking for all icme start times for iwinind that do not have corresponding arrival
#each arrival has a 1.5 days window is used
#these are all predicted arrival times at L1: win_all_pred_times_a
#these are all icme start times at L1: winstart

#will contain 1 for false rejection, 0 if not
false_rejection_array_a=np.zeros(np.size(iwinind))
false_rejection_array_b=np.zeros(np.size(iwinind))

#go through all icmes
for k in np.arange(0,np.size(winstart)):
   
   #for every icme, calculate array with arrival differences to all predictions
   #take absolute minimum
   min_diff_a=np.min(abs(winstart[k]-win_all_pred_times_a))
   min_diff_b=np.min(abs(winstart[k]-win_all_pred_times_b))
   
   #if this difference is > timewindow days **crucial for numbers!
   if min_diff_a > timewindow: false_rejection_array_a[k]=1
   if min_diff_b > timewindow: false_rejection_array_b[k]=1
  
false_rejections_a=np.size(np.where(false_rejection_array_a > 0))
false_rejections_b=np.size(np.where(false_rejection_array_b > 0))


##########################################################
#define correct rejections: no prediction, no in situ detection
#take number of all observed CMEs and subtract other rates 


#n_days=April 1 2007 to 19 August 2014 for HIA
n_days_a=2697
#n_days=April 1 2007 to 27 September 2014 for HIB
n_days_b=2736


correct_rejections_a=n_days_a/(timewindow*2)-false_rejections_a-hits_a-false_alarms_a
correct_rejections_b=n_days_b/(timewindow*2)-false_rejections_b-hits_b-false_alarms_b





print('------------------------- Earth L1 ---------------------------')

print('contingency table HIA:')

print('number of all predicted arrivals', np.size(arrwinind_a))
print('number of in situ ICMEs', np.size(iwinind))

print('')

print('Hits (TP)', hits_a)
print('False alarms (FP)', false_alarms_a)
print('False rejections (FN)', false_rejections_a)
print('Correct rejections (TN)', correct_rejections_a)
print('')
print('')




############################  all scores for HIA

TP=hits_a
FP=false_alarms_a
FN=false_rejections_a
TN=correct_rejections_a
n=TP+FP+FN+TN

############### simple scores
TPR=TP/(FN+TP)
FNR=FN/(TP+FN)
PPV=TP/(TP+FP)
FAR=FP/(TP+FP)

print('TPR ',round(TPR,2))
print('FNR ',round(FNR,2))
print('PPV ',round(PPV,2))
print('FAR ',round(FAR,2))

############### advanced scores

#Threat score
TS=TP/(TP+FP+FN)

#Bias
BS=(TP+FP)/(TP+FN)

#Heidke skill score HSS
#Joliffe 2006 page 48
E=((TP+FN)/n)*((TP+FP)/n)+((FP+TN)/n)*((FN+TN)/n)
PC=(TP+TN)/n
HSS=(PC-E)/(1-E)
#Martins formel:
#HSS2=2*(TP*TN-FP*FN)/((TP+FN)*(FN+TN)+(TP+FP)*(FP+TN))

#TSS, Peirce's 
TSS=(TP*TN-FP*FN)/((FP+TN)*(TP+FN))

print('TS ',round(TS,2))
print('BS ',round(BS,2))
print('HSS ',round(HSS,2))
#print('HSS2 ',round(HSS2,2))
print('TSS ',round(TSS,2))



#######################

print('')
print('')
print('')

print('-------------------------')

print('contingency table HIB:')
print('number of all predicted arrivals', np.size(arrwinind_b))
print('number of in situ ICMEs', np.size(iwinind))
print('')
print('Hits (TP)', hits_b)
print('False alarms (FP)', false_alarms_b)
print('False rejections (FN)', false_rejections_b)
print('Correct rejections (TN)', correct_rejections_b)
print('')
print('')

TP=hits_b
FP=false_alarms_b
FN=false_rejections_b
TN=correct_rejections_b
n=TP+FP+FN+TN

############### simple scores
TPR=TP/(FN+TP)
FNR=FN/(TP+FN)
PPV=TP/(TP+FP)
FAR=FP/(TP+FP)

print('TPR ',round(TPR,2))
print('FNR ',round(FNR,2))
print('PPV ',round(PPV,2))
print('FAR ',round(FAR,2))

############### advanced scores

#Threat score
TS=TP/(TP+FP+FN)

#Bias
BS=(TP+FP)/(TP+FN)

#Heidke skill score HSS
#Joliffe 2006 page 48
E=((TP+FN)/n)*((TP+FP)/n)+((FP+TN)/n)*((FN+TN)/n)
PC=(TP+TN)/n
HSS=(PC-E)/(1-E)
#Martins formel:
#HSS=2*(TP*TN-FP*FN)/[(TP+FN)(FN+TN)+(TP+FP)(FP+TN)]

#TSS, Peirce's 
TSS=(TP*TN-FP*FN)/((FP+TN)*(TP+FN))

print('TS ',round(TS,2))
print('BS ',round(BS,2))
print('HSS ',round(HSS,2))
print('TSS ',round(TSS,2))

print('')
print('')
print('')




##################################################################





















#############################(2) STEREO-A self prediction

#all predicted arrivals
np.size(arrstaind_a)

#1 part of contingency table:
#0 for false alarm, 1 or 2 for hits
np.size(stahits_a)

#HIA
false_alarms_a=np.size(np.where(stahits_a ==0))
hits_a=np.size(np.where(stahits_a > 0))

#HIB
false_alarms_b=np.size(np.where(stahits_b ==0))
hits_b=np.size(np.where(stahits_b > 0))

#all Wind in situ ICMEs
np.size(istaind)
stastart=icme_start_time_num[istaind]

############ define false rejection: no prediction, but in situ detection
#this means looking for all icme start times for iwinind that do not have corresponding arrival
#each arrival has a 1.5 days window is used
#these are all predicted arrival times at L1: win_all_pred_times_a
#these are all icme start times at L1: winstart
 
#will contain 1 for false rejection, 0 if not
false_rejection_array_a=np.zeros(np.size(istaind))
false_rejection_array_b=np.zeros(np.size(istaind))

#go through all icmes
for k in np.arange(0,np.size(stastart)):
   
   #for every icme, calculate array with arrival differences to all predictions
   #take absolute minimum
   min_diff_a=np.min(abs(stastart[k]-sta_all_pred_times_a))
   
   #if this difference is > timewindow days **crucial for numbers!
   if min_diff_a > timewindow: false_rejection_array_a[k]=1
  
false_rejections_a=np.size(np.where(false_rejection_array_a > 0))


##########################################################
#define correct rejections: no prediction, no in situ detection
#take number of all observed CMEs and subtract other rates 


#n_days=April 1 2007 to 19 August 2014 for HIA
n_days_a=2697
#n_days=April 1 2007 to 27 September 2014 for HIB
n_days_b=2736

correct_rejections_a=n_days_a/(timewindow*2)-false_rejections_a-hits_a-false_alarms_a

print('------------------------- STEREO-A self prediction ---------------------------')

print('contingency table HIA self prediction:')

print('number of all predicted arrivals', np.size(arrstaind_a))
print('number of in situ ICMEs', np.size(istaind))

print('')

print('Hits (TP)', hits_a)
print('False alarms (FP)', false_alarms_a)
print('False rejections (FN)', false_rejections_a)
print('Correct rejections (TN)', correct_rejections_a)
print('')
print('')




############################  all scores for HIA

TP=hits_a
FP=false_alarms_a
FN=false_rejections_a
TN=correct_rejections_a
n=TP+FP+FN+TN

############### simple scores
TPR=TP/(FN+TP)
FNR=FN/(TP+FN)
PPV=TP/(TP+FP)
FAR=FP/(TP+FP)

print('TPR ',round(TPR,2))
print('FNR ',round(FNR,2))
print('PPV ',round(PPV,2))
print('FAR ',round(FAR,2))

############### advanced scores

#Threat score
TS=TP/(TP+FP+FN)

#Bias
BS=(TP+FP)/(TP+FN)

#Heidke skill score HSS
#Joliffe 2006 page 48
E=((TP+FN)/n)*((TP+FP)/n)+((FP+TN)/n)*((FN+TN)/n)
PC=(TP+TN)/n
HSS=(PC-E)/(1-E)
#Martins formel:
#HSS2=2*(TP*TN-FP*FN)/((TP+FN)*(FN+TN)+(TP+FP)*(FP+TN))

#TSS, Peirce's 
TSS=(TP*TN-FP*FN)/((FP+TN)*(TP+FN))

print('TS ',round(TS,2))
print('BS ',round(BS,2))
print('HSS ',round(HSS,2))
#print('HSS2 ',round(HSS2,2))
print('TSS ',round(TSS,2))





























#############################(3) STEREO-B self prediction

#all predicted arrivals
np.size(arrstbind_b)

#1 part of contingency table:
#0 for false alarm, 1 or 2 for hits
np.size(stbhits_b)

#HIB
false_alarms_b=np.size(np.where(stbhits_b ==0))
hits_b=np.size(np.where(stbhits_b > 0))

#all Wind in situ ICMEs
np.size(istbind)
stbstart=icme_start_time_num[istbind]

false_rejection_array_b=np.zeros(np.size(istbind))

#go through all icmes
for k in np.arange(0,np.size(stbstart)):
   
   #for every icme, calculate array with arrival differences to all predictions
   #take absolute minimum
   min_diff_b=np.min(abs(stbstart[k]-stb_all_pred_times_b))
   
   #if this difference is > timewindow days **crucial for numbers!
   if min_diff_b > timewindow: false_rejection_array_b[k]=1
  
false_rejections_b=np.size(np.where(false_rejection_array_b > 0))


##########################################################
#define correct rejections: no prediction, no in situ detection
#take number of all observed CMEs and subtract other rates 

#n_days=April 1 2007 to 27 September 2014 for HIB
n_days_b=2736


correct_rejections_b=n_days_b/(timewindow*2)-false_rejections_b-hits_b-false_alarms_b

print('------------------------- STEREO-B self prediction ---------------------------')

print('contingency table HIB self prediction:')

print('number of all predicted arrivals', np.size(arrstbind_b))
print('number of in situ ICMEs', np.size(istbind))

print('')

print('Hits (TP)', hits_b)
print('False alarms (FP)', false_alarms_b)
print('False rejections (FN)', false_rejections_b)
print('Correct rejections (TN)', correct_rejections_b)
print('')
print('')




############################  all scores for HIB

TP=hits_b
FP=false_alarms_b
FN=false_rejections_b
TN=correct_rejections_b
n=TP+FP+FN+TN

############### simple scores
TPR=TP/(FN+TP)
FNR=FN/(TP+FN)
PPV=TP/(TP+FP)
FAR=FP/(TP+FP)

print('TPR ',round(TPR,2))
print('FNR ',round(FNR,2))
print('PPV ',round(PPV,2))
print('FAR ',round(FAR,2))

############### advanced scores

#Threat score
TS=TP/(TP+FP+FN)

#Bias
BS=(TP+FP)/(TP+FN)

#Heidke skill score HSS
#Joliffe 2006 page 48
E=((TP+FN)/n)*((TP+FP)/n)+((FP+TN)/n)*((FN+TN)/n)
PC=(TP+TN)/n
HSS=(PC-E)/(1-E)
#TSS, Peirce's 
TSS=(TP*TN-FP*FN)/((FP+TN)*(TP+FN))

print('TS ',round(TS,2))
print('BS ',round(BS,2))
print('HSS ',round(HSS,2))
print('TSS ',round(TSS,2))










#############################(4) Venus

#all predicted arrivals
np.size(arrvexind_a)

#1 part of contingency table:
#0 for false alarm, 1 or 2 for hits
np.size(vexhits_a)

#HIA
false_alarms_a=np.size(np.where(vexhits_a ==0))
hits_a=np.size(np.where(vexhits_a > 0))

#HIB
false_alarms_b=np.size(np.where(vexhits_b ==0))
hits_b=np.size(np.where(vexhits_b > 0))

#all Wind in situ ICMEs
np.size(ivexind)
vexstart=icme_start_time_num[ivexind]

############ define false rejection: no prediction, but in situ detection
#this means looking for all icme start times for iwinind that do not have corresponding arrival
#each arrival has a 1.5 days window is used
#these are all predicted arrival times at L1: win_all_pred_times_a
#these are all icme start times at L1: winstart

#will contain 1 for false rejection, 0 if not
false_rejection_array_a=np.zeros(np.size(ivexind))
false_rejection_array_b=np.zeros(np.size(ivexind))

#go through all icmes
for k in np.arange(0,np.size(vexstart)):
   
   #for every icme, calculate array with arrival differences to all predictions
   #take absolute minimum
   min_diff_a=np.min(abs(vexstart[k]-vex_all_pred_times_a))
   min_diff_b=np.min(abs(vexstart[k]-vex_all_pred_times_b))
   
   #if this difference is > timewindow days **crucial for numbers!
   if min_diff_a > timewindow: false_rejection_array_a[k]=1
   if min_diff_b > timewindow: false_rejection_array_b[k]=1
  
false_rejections_a=np.size(np.where(false_rejection_array_a > 0))
false_rejections_b=np.size(np.where(false_rejection_array_b > 0))


##########################################################
#define correct rejections: no prediction, no in situ detection
#take number of all observed CMEs and subtract other rates 

#n_days=April 1 2007 to 19 August 2014 for HIA
#n_days_a=2697
#n_days=April 1 2007 to 27 September 2014 for HIB
#n_days_b=2736

#correct number of all observations days for Venus orbit - Venus must be inside HIA and HIB FoVs

#time resolution is 6 hours of positions, so we divide by 4 when counting days

vexlong=pos.venus[1]*180/np.pi
stalong=pos.sta[1]*180/np.pi
vexok_a_alltime=np.where(np.logical_and(vexlong < stalong+30,vexlong-stalong > -180))
n_days_a=np.size(np.where(np.logical_and(pos_time_num[vexok_a_alltime] > mdates.date2num(sunpy.time.parse_time('2007-Apr-1')),pos_time_num[vexok_a_alltime] < mdates.date2num(sunpy.time.parse_time('2014-Aug-19')) )))/4

stblong=pos.stb[1]*180/np.pi
vexok_b_alltime=np.where(np.logical_and(vexlong > stblong-30,vexlong-stblong < 180))
n_days_b=np.size(np.where(np.logical_and(pos_time_num[vexok_b_alltime] > mdates.date2num(sunpy.time.parse_time('2007-Apr-1')),pos_time_num[vexok_b_alltime] < mdates.date2num(sunpy.time.parse_time('2014-Sep-27')) )))/4


correct_rejections_a=n_days_a/(timewindow*2)-false_rejections_a-hits_a-false_alarms_a
correct_rejections_b=n_days_b/(timewindow*2)-false_rejections_b-hits_b-false_alarms_b

print('------------------------- Venus --------------------------')
print('contingency table HIA:')
print('number of all predicted arrivals', np.size(arrvexind_a))
print('number of in situ ICMEs', np.size(ivexind))
print('')
print('Hits (TP)', hits_a)
print('False alarms (FP)', false_alarms_a)
print('False rejections (FN)', false_rejections_a)
print('Correct rejections (TN)', correct_rejections_a)
print('')
print('')
############################  all scores for HIA
TP=hits_a
FP=false_alarms_a
FN=false_rejections_a
TN=correct_rejections_a
n=TP+FP+FN+TN
############### simple scores
TPR=TP/(FN+TP)
FNR=FN/(TP+FN)
PPV=TP/(TP+FP)
FAR=FP/(TP+FP)
print('TPR ',round(TPR,2))
print('FNR ',round(FNR,2))
print('PPV ',round(PPV,2))
print('FAR ',round(FAR,2))
############### advanced scores
#Threat score
TS=TP/(TP+FP+FN)
#Bias
BS=(TP+FP)/(TP+FN)
#Heidke skill score HSS
#Joliffe 2006 page 48
E=((TP+FN)/n)*((TP+FP)/n)+((FP+TN)/n)*((FN+TN)/n)
PC=(TP+TN)/n
HSS=(PC-E)/(1-E)
#TSS, Peirce's 
TSS=(TP*TN-FP*FN)/((FP+TN)*(TP+FN))
print('TS ',round(TS,2))
print('BS ',round(BS,2))
print('HSS ',round(HSS,2))
print('TSS ',round(TSS,2))
#######################
print('')
print('')
print('')
print('-------------------------')
print('contingency table HIB:')
print('number of all predicted arrivals', np.size(arrvexind_b))
print('number of in situ ICMEs', np.size(ivexind))
print('')
print('Hits (TP)', hits_b)
print('False alarms (FP)', false_alarms_b)
print('False rejections (FN)', false_rejections_b)
print('Correct rejections (TN)', correct_rejections_b)
print('')
print('')
TP=hits_b
FP=false_alarms_b
FN=false_rejections_b
TN=correct_rejections_b
n=TP+FP+FN+TN
############### simple scores
TPR=TP/(FN+TP)
FNR=FN/(TP+FN)
PPV=TP/(TP+FP)
FAR=FP/(TP+FP)
print('TPR ',round(TPR,2))
print('FNR ',round(FNR,2))
print('PPV ',round(PPV,2))
print('FAR ',round(FAR,2))
############### advanced scores
#Threat score
TS=TP/(TP+FP+FN)
#Bias
BS=(TP+FP)/(TP+FN)
#Heidke skill score HSS
#Joliffe 2006 page 48
E=((TP+FN)/n)*((TP+FP)/n)+((FP+TN)/n)*((FN+TN)/n)
PC=(TP+TN)/n
HSS=(PC-E)/(1-E)
#Martins formel:
#HSS=2*(TP*TN-FP*FN)/[(TP+FN)(FN+TN)+(TP+FP)(FP+TN)]
#TSS, Peirce's 
TSS=(TP*TN-FP*FN)/((FP+TN)*(TP+FN))
print('TS ',round(TS,2))
print('BS ',round(BS,2))
print('HSS ',round(HSS,2))
print('TSS ',round(TSS,2))
print('')
print('')
print('')














####################################### (2) ARRIVAL TIMES and SPEEDS ###########################


print('--------------------------- Arrival time and speed comparison -------------------')


############# 1 ARRIVAL TIME 


#all STEREO-A predictions compared to observed in situ time
winhits_time_diff_clean_a=winhits_time_diff_a[np.where(winhits_time_diff_a != 0)]*24 # in hours
vexhits_time_diff_clean_a=vexhits_time_diff_a[np.where(vexhits_time_diff_a != 0)]*24 # in hours
meshits_time_diff_clean_a=meshits_time_diff_a[np.where(meshits_time_diff_a != 0)]*24 # in hours
stbhits_time_diff_clean_a=stbhits_time_diff_a[np.where(stbhits_time_diff_a != 0)]*24 # in hours
stahits_time_diff_clean_a=stahits_time_diff_a[np.where(stahits_time_diff_a != 0)]*24 # in hours


#all STEREO-B predictions compared to observed in situ time
winhits_time_diff_clean_b=winhits_time_diff_b[np.where(winhits_time_diff_b != 0)]*24 # in hours
vexhits_time_diff_clean_b=vexhits_time_diff_b[np.where(vexhits_time_diff_b != 0)]*24 # in hours
meshits_time_diff_clean_b=meshits_time_diff_b[np.where(meshits_time_diff_b != 0)]*24 # in hours
stbhits_time_diff_clean_b=stbhits_time_diff_b[np.where(stbhits_time_diff_b != 0)]*24 # in hours
stahits_time_diff_clean_b=stahits_time_diff_b[np.where(stahits_time_diff_b != 0)]*24 # in hours


#stats:

print()

print('time window used in days  ', timewindow)

print('mean/std of arrival time C-O:')
print()
print('HIA: Earth L1 {:.1f} +/- {:.1f} hours'.format(np.mean(winhits_time_diff_clean_a),np.std(winhits_time_diff_clean_a))) 
print('HIB: Earth L1 {:.1f} +/- {:.1f} hours'.format(np.mean(winhits_time_diff_clean_b),np.std(winhits_time_diff_clean_b))) 
print()
print('HIA: VEX {:.1f} +/- {:.1f} hours'.format(np.mean(vexhits_time_diff_clean_a),np.std(vexhits_time_diff_clean_a))) 
print('HIB: VEX {:.1f} +/- {:.1f} hours'.format(np.mean(vexhits_time_diff_clean_b),np.std(vexhits_time_diff_clean_b))) 
print()
print('HIA: MES {:.1f} +/- {:.1f} hours'.format(np.mean(meshits_time_diff_clean_a),np.std(meshits_time_diff_clean_a))) 
print('HIB: MES {:.1f} +/- {:.1f} hours'.format(np.mean(meshits_time_diff_clean_b),np.std(meshits_time_diff_clean_b))) 
print()
print('HIA: STB {:.1f} +/- {:.1f} hours'.format(np.mean(stbhits_time_diff_clean_a),np.std(stbhits_time_diff_clean_a))) 
print('HIB: STB {:.1f} +/- {:.1f} hours'.format(np.mean(stbhits_time_diff_clean_b),np.std(stbhits_time_diff_clean_b))) 
print()
print('HIA: STA {:.1f} +/- {:.1f} hours'.format(np.mean(stahits_time_diff_clean_a),np.std(stahits_time_diff_clean_a))) 
print('HIB: STA {:.1f} +/- {:.1f} hours'.format(np.mean(stahits_time_diff_clean_b),np.std(stahits_time_diff_clean_b))) 
print()

#make an array containing all arrival time differences o-c

all_arrival_diff_a=np.concatenate((winhits_time_diff_clean_a,vexhits_time_diff_clean_a,meshits_time_diff_clean_a,stbhits_time_diff_clean_a,stahits_time_diff_clean_a),axis=0)
all_arrival_diff_b=np.concatenate((winhits_time_diff_clean_b,vexhits_time_diff_clean_b,meshits_time_diff_clean_b,stbhits_time_diff_clean_b,stahits_time_diff_clean_b),axis=0)
print()
print('np.mean absolute errors ALL: ')
print('HIA:', np.mean(abs(all_arrival_diff_a)))
print('HIB:', np.mean(abs(all_arrival_diff_b)))

#just for 1 AU
AU1_arrival_diff_a=np.concatenate((winhits_time_diff_clean_a,stbhits_time_diff_clean_a,stahits_time_diff_clean_a),axis=0)
AU1_arrival_diff_b=np.concatenate((winhits_time_diff_clean_b,stbhits_time_diff_clean_b,stahits_time_diff_clean_b),axis=0)

print('np.mean absolute errors 1 AU: ')
print('HIA:', np.mean(abs(AU1_arrival_diff_a)))
print('HIB:', np.mean(abs(AU1_arrival_diff_b)))

print('np.mean absolute errors Earth: ')
print('HIA:', np.mean(abs(winhits_time_diff_clean_a)))
print('HIB:', np.mean(abs(winhits_time_diff_clean_b)))



print('total number of comparisons for HIA at 1 AU:',len(AU1_arrival_diff_a))
print('HIA: {:.1f} +/- {:.1f} hours'.format(np.mean(AU1_arrival_diff_a),np.std(AU1_arrival_diff_a))) 

print('total number of comparisons for HIB at 1 AU:',len(AU1_arrival_diff_b))
print('HIB: {:.1f} +/- {:.1f} hours'.format(np.mean(AU1_arrival_diff_b),np.std(AU1_arrival_diff_b))) 


print()
print()

print('total number of comparisons for HIA',len(all_arrival_diff_a))
print('HIA: all spacecraft {:.1f} +/- {:.1f} hours'.format(np.mean(all_arrival_diff_a),np.std(all_arrival_diff_a))) 

print('total number of comparisons for HIB',len(all_arrival_diff_b))
print('HIB: all spacecraft {:.1f} +/- {:.1f} hours'.format(np.mean(all_arrival_diff_b),np.std(all_arrival_diff_b))) 


print()
print()

###################### plot histogram for arrival time differences:

fsize=10

#for labels
xoff=0.4
yoff=0.98

fig=plt.figure(4,figsize=(5,9))

print('in the arrival time histogram the left plot is for STEREO-A HI and the right for STEREO-B HI.')

plt.figtext(0.05,0.98,'HIA',color='red', fontsize=fsize+3, ha='center')
plt.figtext(0.95,0.98,'HIB',color='royalblue', fontsize=fsize+3, ha='center')

plt.figtext(0.03,0.96,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.03,0.96-0.185*1,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.03,0.96-0.185*2,'c',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.03,0.96-0.185*3,'d',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.03,0.96-0.185*4,'e',color='black', fontsize=fsize, ha='left',fontweight='bold')


plt.figtext(0.97,0.96,'f',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.97,0.96-0.185*1,'g',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.97,0.96-0.185*2,'h',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.97,0.96-0.185*3,'i',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.97,0.96-0.185*4,'j',color='black', fontsize=fsize, ha='left',fontweight='bold')


#Wind - STEREO-A HI predictions
ax1 = plt.subplot2grid((5,2), (0, 0))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(winhits_time_diff_clean_a, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='mediumseagreen', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff,'Earth L1',color='mediumseagreen', fontsize=fsize, ha='left')
plt.ylabel('Number of CMEs',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

#Wind - STEREO-B HI predictions
ax1 = plt.subplot2grid((5,2), (0, 1))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(winhits_time_diff_clean_b, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='mediumseagreen', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff,'Earth L1',color='mediumseagreen', fontsize=fsize, ha='left')
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

#VEX - STEREO-A HI predictions
ax1 = plt.subplot2grid((5,2), (1, 0))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(vexhits_time_diff_clean_a, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='orange', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff-0.19*1,'VEX',color='orange', fontsize=fsize, ha='left')
plt.ylabel('Number of CMEs',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

#VEX - STEREO-B HI predictions
ax1 = plt.subplot2grid((5,2), (1, 1))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(vexhits_time_diff_clean_b, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='orange', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff-0.19*1,'VEX',color='orange', fontsize=fsize, ha='left')
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


#MES - STEREO-A HI predictions
ax1 = plt.subplot2grid((5,2), (2, 0))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(meshits_time_diff_clean_a, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='dimgrey', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff-0.19*2,'MESSENGER',color='dimgrey', fontsize=fsize, ha='left')
plt.ylabel('Number of CMEs',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


#MES - STEREO-B HI predictions
ax1 = plt.subplot2grid((5,2), (2, 1))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(meshits_time_diff_clean_b, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='dimgrey', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff-0.19*2,'MESSENGER',color='dimgrey', fontsize=fsize, ha='left')
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

#STA - STEREO-A HI predictions
ax1 = plt.subplot2grid((5,2), (3, 0))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(stahits_time_diff_clean_a, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='red', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff-0.19*3,'STEREO-A',color='red', fontsize=fsize, ha='left')
plt.ylabel('Number of CMEs',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

#STA - STEREO-B HI predictions
ax1 = plt.subplot2grid((5,2), (3, 1))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(stahits_time_diff_clean_b, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='red', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff-0.19*3,'STEREO-A',color='red', fontsize=fsize, ha='left')
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

#STB - STEREO-A HI predictions
ax1 = plt.subplot2grid((5,2), (4, 0))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(stbhits_time_diff_clean_a, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='royalblue', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff-0.19*4,'STEREO-B',color='royalblue', fontsize=fsize, ha='left')
plt.ylabel('Number of CMEs',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 
plt.xlabel('$\mathregular{\Delta}$  C-O [hours]', fontsize=fsize)

#STB - STEREO-B HI predictions
ax1 = plt.subplot2grid((5,2), (4, 1))
hours_bin_edges=np.linspace(-48,48,num=13)
(hist, bin_edges) = np.histogram(stbhits_time_diff_clean_b, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='royalblue', alpha=0.5)
plt.xlim(-48,48)
plt.figtext(xoff,yoff-0.19*4,'STEREO-B',color='royalblue', fontsize=fsize, ha='left')
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 
plt.xlabel('$\mathregular{\Delta}$  C-O [hours]', fontsize=fsize)

sns.despine()
plt.tight_layout()

plt.savefig('plots/histogram_arrival_times.pdf', dpi=300)
plt.savefig('plots/histogram_arrival_times.png', dpi=300)







##################### 2 ARRIVAL SPEED


#all STEREO-A predictions compared to observed in situ speed
winhits_speed_diff_clean_a=winhits_speed_diff_a[np.where(winhits_speed_diff_a != 0)]
stbhits_speed_diff_clean_a=stbhits_speed_diff_a[np.where(stbhits_speed_diff_a != 0)]
stahits_speed_diff_clean_a=stahits_speed_diff_a[np.where(stahits_speed_diff_a != 0)]

#get rid of nan's
winhits_speed_diff_clean_a=winhits_speed_diff_clean_a[np.where(np.isfinite(winhits_speed_diff_clean_a))]
stbhits_speed_diff_clean_a=stbhits_speed_diff_clean_a[np.where(np.isfinite(stbhits_speed_diff_clean_a))]
stahits_speed_diff_clean_a=stahits_speed_diff_clean_a[np.where(np.isfinite(stahits_speed_diff_clean_a))]


#all STEREO-B predictions compared to observed in situ speed
winhits_speed_diff_clean_b=winhits_speed_diff_b[np.where(winhits_speed_diff_b != 0)]
stbhits_speed_diff_clean_b=stbhits_speed_diff_b[np.where(stbhits_speed_diff_b != 0)]
stahits_speed_diff_clean_b=stahits_speed_diff_b[np.where(stahits_speed_diff_b != 0)]

#get rid of nan's
winhits_speed_diff_clean_b=winhits_speed_diff_clean_b[np.where(np.isfinite(winhits_speed_diff_clean_b))]
stbhits_speed_diff_clean_b=stbhits_speed_diff_clean_b[np.where(np.isfinite(stbhits_speed_diff_clean_b))]
stahits_speed_diff_clean_b=stahits_speed_diff_clean_b[np.where(np.isfinite(stahits_speed_diff_clean_b))]


#stats:

print()

print('mean/std of arrival speed to icme sheath speed C-O:')
print()
print('HIA: Earth L1 {:.1f} +/- {:.1f} km/s'.format(np.nanmean(winhits_speed_diff_clean_a),np.nanstd(winhits_speed_diff_clean_a))) 
print('HIB: Earth L1 {:.1f} +/- {:.1f} km/s'.format(np.nanmean(winhits_speed_diff_clean_b),np.nanstd(winhits_speed_diff_clean_b))) 
print()
#print('HIA: STB {:.1f} +/- {:.1f} km/s'.format(np.nanmean(stbhits_speed_diff_clean_a),np.nanstd(stbhits_speed_diff_clean_a))) 
#print('HIB: STB {:.1f} +/- {:.1f} km/s'.format(np.nanmean(stbhits_speed_diff_clean_b),np.nanstd(stbhits_speed_diff_clean_b))) 
#print()
#print('HIA: STA L1 {:.1f} +/- {:.1f} km/s'.format(np.nanmean(stahits_speed_diff_clean_a),np.nanstd(stahits_speed_diff_clean_a))) 
#print('HIB: STA L1 {:.1f} +/- {:.1f} km/s'.format(np.nanmean(stahits_speed_diff_clean_b),np.nanstd(stahits_speed_diff_clean_b))) 
#print()

print('note that the statistics are too low for STEREO-A and STEREO-B in situ, so they are omitted in the plot')
print('The reason is that the in situ sheath speeds are compared, but not many CMEs had sheahts during minimum when STEREO was still close together')

print()
print()

#all_speed_diff_a=np.concatenate((winhits_speed_diff_clean_a,stbhits_speed_diff_clean_a,stahits_speed_diff_clean_a),axis=0)
#all_speed_diff_b=np.concatenate((winhits_speed_diff_clean_b,stbhits_speed_diff_clean_b,stahits_speed_diff_clean_b),axis=0)


print('total number of comparisons for HIA',len(winhits_speed_diff_clean_a))
print('total number of comparisons for HIB',len(winhits_speed_diff_clean_b))

print()
print()


###################### plot histogram for arrival speed differences at Earth ###############


sns.set_context("talk")     
#sns.set_style("darkgrid")  
sns.set_style("ticks")

fsize=10

#for labels
xoff=0.3

yoff=0.6

fig=plt.figure(5,figsize=(10,5))

print('in the arrival speed histogram the left plot is for STEREO-A HI and the right for STEREO-B HI.')


plt.figtext(0.33,0.8,'HIA',color='red', fontsize=fsize+3, ha='center')
plt.figtext(0.8,0.8,'HIB',color='royalblue', fontsize=fsize+3, ha='center')


#Wind - STEREO-A HI predictions
ax1 = plt.subplot2grid((1,2), (0, 0))
hours_bin_edges=np.linspace(-300,1000,num=13)
(hist, bin_edges) = np.histogram(winhits_speed_diff_clean_a, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='mediumseagreen', alpha=0.5)
plt.xlim(-300,1000)
plt.figtext(xoff,yoff,'Earth L1',color='mediumseagreen', fontsize=fsize, ha='left')
plt.ylabel('Number of CMEs',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 
plt.xlabel('$\mathregular{\Delta V_{IPo} - V_{sheath}}$ $\mathregular{ [km \\ s^{-1}]}}$', fontsize=fsize)


#Wind - STEREO-B HI predictions
ax1 = plt.subplot2grid((1,2), (0, 1))
hours_bin_edges=np.linspace(-300,1000,num=13)
(hist, bin_edges) = np.histogram(winhits_speed_diff_clean_b, hours_bin_edges)
width = bin_edges[1] - bin_edges[0]
ax1.bar(bin_edges[:-1],hist, width=width, color='mediumseagreen', alpha=0.5)
plt.xlim(-300,1000)
plt.figtext(xoff,yoff,'Earth L1',color='mediumseagreen', fontsize=fsize, ha='left')
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 
plt.xlabel('$\mathregular{\Delta V_{IPo} - V_{sheath}}$ $\mathregular{ [km \\ s^{-1}]}}$', fontsize=fsize)

#panel labels
plt.figtext(0.02,0.93,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.50,0.93,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')

plt.tight_layout()


 
plt.savefig('plots/histogram_speeds.pdf', dpi=300)
plt.savefig('plots/histogram_speeds.png', dpi=300)














################################## plots relating HI to B




#stack overlap arrays

overlap_a=np.asarray(np.column_stack((np.trim_zeros(overlap_a_hi),np.trim_zeros(overlap_a_icme))),int)
overlap_b=np.asarray(np.column_stack((np.trim_zeros(overlap_b_hi),np.trim_zeros(overlap_b_icme))),int)

#events that come together can be found by, e.g. for event 10
#aid[int(overlap_a[10,0])]
#Out[225]: 'HCME_A__20100203_01'
#iid[int(overlap_a[10,1])]
#Out[224]: 'ICME_Wind_NASA_20100207_01'

#####################HIA

ai_a=aid[overlap_a[:,0]]
ii_a=iid[overlap_a[:,1]]

#define parameters

#asarray converts the index in overlap_a to integer
sl_a=sse_launch_num[overlap_a[:,0]]
ta_a=target_arrival_num[overlap_a[:,0]]
#TARGET_SPEED - corrected for shape
tspeed_a=target_speed[overlap_a[:,0]]
#Target_delta
tdelta=target_delta[overlap_a[:,0]]
#target_distance
tdistance_a=target_distance[overlap_a[:,0]]
tname_a=target_name[overlap_a[:,0]]


pan_a=pa_n[overlap_a[:,0]]
pas_a=pa_s[overlap_a[:,0]]
paext_a=abs(pan_a-pas_a)
pacenter_a=pa_center[overlap_a[:,0]]



#get travel time from HI target_arrival-sse_launch, convert from days to hours
TT_a=(ta_a-sl_a)*24

#alternatively get travel time  from icme_start_time-sse_launch, any difference in plot?
icme_start_a=icme_start_time_num[overlap_a[:,1]] 
TTalt_a=(icme_start_a-sl_a)*24



#TTalt contains a few < values 0! for time window 1.5 days


#insitu b in magnetic obstacle
bmean_a=mo_bmean[overlap_a[:,1]] 
#error bar:
bstd_a=mo_bstd[overlap_a[:,1]] 
bmax_a=mo_bmax[overlap_a[:,1]] 

bzmin_a=mo_bzmean[overlap_a[:,1]] 
bzmean_a=mo_bzmin[overlap_a[:,1]] 


##################### HIB

ai_b=aid[overlap_b[:,0]]
ii_b=iid[overlap_b[:,1]]

sl_b=sse_launch_num[overlap_b[:,0]]
ta_b=target_arrival_num[overlap_b[:,0]]
tspeed_b=target_speed[overlap_b[:,0]]
tdelta=target_delta[overlap_b[:,0]]
tdistance_b=target_distance[overlap_b[:,0]]
tname_b=target_name[overlap_b[:,0]]

pan_b=pa_n[overlap_b[:,0]]
pas_b=pa_s[overlap_b[:,0]]
paext_b=abs(pan_b-pas_b)
pacenter_b=pa_center[overlap_b[:,0]]



TT_b=(ta_b-sl_b)*24

icme_start_b=icme_start_time_num[overlap_b[:,1]] 

TTalt_b=(icme_start_b-sl_b)*24

bmean_b=mo_bmean[overlap_b[:,1]] 
bstd_b=mo_bstd[overlap_b[:,1]] 
bmax_b=mo_bmax[overlap_b[:,1]] 

bzmin_b=mo_bzmean[overlap_b[:,1]] 
bzmean_b=mo_bzmin[overlap_b[:,1]] 








win_ind_over_a=np.where(tname_a=='EARTH_L1')
vex_ind_over_a=np.where(tname_a=='VENUS')
mes_ind_over_a=np.where(tname_a=='MESSENGER')
sta_ind_over_a=np.where(tname_a=='STEREO-A')
stb_ind_over_a=np.where(tname_a=='STEREO-B')

win_ind_over_b=np.where(tname_b=='EARTH_L1')
vex_ind_over_b=np.where(tname_b=='VENUS')
mes_ind_over_b=np.where(tname_b=='MESSENGER')
sta_ind_over_b=np.where(tname_b=='STEREO-A')
stb_ind_over_b=np.where(tname_b=='STEREO-B')










################################################ Figure with travel time vs bmean, bmax


#s bmean vs target_speed


#errorbars bstd for bmean, 10% on TT
fig=plt.figure(6,figsize=(15,8))
ax1 = plt.subplot2grid((1,2), (0, 0))

#errorbars 10% on speed taken from Moestl et al. 2014, and bstd from ICMECAT

plt.errorbar(tspeed_a[win_ind_over_a],bmean_a[win_ind_over_a],xerr=tspeed_a[win_ind_over_a]*0.1,yerr=bstd_a[win_ind_over_a],color='mediumseagreen', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(tspeed_a[vex_ind_over_a],bmean_a[vex_ind_over_a],xerr=tspeed_a[vex_ind_over_a]*0.1,yerr=bstd_a[vex_ind_over_a],color='orange', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(tspeed_a[mes_ind_over_a],bmean_a[mes_ind_over_a],xerr=tspeed_a[mes_ind_over_a]*0.1,yerr=bstd_a[mes_ind_over_a],color='dimgrey', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(tspeed_a[sta_ind_over_a],bmean_a[sta_ind_over_a],xerr=tspeed_a[sta_ind_over_a]*0.1,yerr=bstd_a[sta_ind_over_a],color='red', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(tspeed_a[stb_ind_over_a],bmean_a[stb_ind_over_a],xerr=tspeed_a[stb_ind_over_a]*0.1,yerr=bstd_a[stb_ind_over_a],color='royalblue', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)

plt.errorbar(tspeed_b[win_ind_over_b],bmean_b[win_ind_over_b],xerr=tspeed_b[win_ind_over_b]*0.1,yerr=bstd_b[win_ind_over_b],color='mediumseagreen', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9,label='Earth L1')
plt.errorbar(tspeed_b[sta_ind_over_b],bmean_b[sta_ind_over_b],xerr=tspeed_b[sta_ind_over_b]*0.1,yerr=bstd_b[sta_ind_over_b],color='red', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9, label='STEREO-A')
plt.errorbar(tspeed_b[stb_ind_over_b],bmean_b[stb_ind_over_b],xerr=tspeed_b[stb_ind_over_b]*0.1,yerr=bstd_b[stb_ind_over_b],color='royalblue', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9, label='STEREO-B')
plt.errorbar(tspeed_b[vex_ind_over_b],bmean_b[vex_ind_over_b],xerr=tspeed_b[vex_ind_over_b]*0.1,yerr=bstd_b[vex_ind_over_b],color='orange', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9,label='VEX')
plt.errorbar(tspeed_b[mes_ind_over_b],bmean_b[mes_ind_over_b],xerr=tspeed_b[mes_ind_over_b]*0.1,yerr=bstd_b[mes_ind_over_b],color='dimgrey', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9, label='MESSENGER')



fsize=10
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

plt.legend(loc=1,fontsize=fsize)

#plt.xlabel('Predicted speed at target $\mathregular{ V_{IPo} [km \\ s^{-1}]}}$',fontsize=fsize)
plt.xlabel('Predicted speed at target $\mathregular{[km \\ s^{-1}]}$',fontsize=fsize)

plt.ylabel('mean B in magnetic obstacle [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.tight_layout()







#same for bmean


#errorbars bstd for bmean, 10% on TT
ax2 = plt.subplot2grid((1,2), (0, 1))

#errorbars 10% on TT taken from Moestl et al. 2014, and bstd from ICMECAT

plt.errorbar(TT_a[win_ind_over_a],bmean_a[win_ind_over_a],xerr=TT_a[win_ind_over_a]*0.1,yerr=bstd_a[win_ind_over_a],color='mediumseagreen', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(TT_a[vex_ind_over_a],bmean_a[vex_ind_over_a],xerr=TT_a[vex_ind_over_a]*0.1,yerr=bstd_a[vex_ind_over_a],color='orange', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(TT_a[mes_ind_over_a],bmean_a[mes_ind_over_a],xerr=TT_a[mes_ind_over_a]*0.1,yerr=bstd_a[mes_ind_over_a],color='dimgrey', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(TT_a[sta_ind_over_a],bmean_a[sta_ind_over_a],xerr=TT_a[sta_ind_over_a]*0.1,yerr=bstd_a[sta_ind_over_a],color='red', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(TT_a[stb_ind_over_a],bmean_a[stb_ind_over_a],xerr=TT_a[stb_ind_over_a]*0.1,yerr=bstd_a[stb_ind_over_a],color='royalblue', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)

plt.errorbar(TT_b[win_ind_over_b],bmean_b[win_ind_over_b],xerr=TT_b[win_ind_over_b]*0.1,yerr=bstd_b[win_ind_over_b],color='mediumseagreen', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9,label='Earth L1')
plt.errorbar(TT_b[sta_ind_over_b],bmean_b[sta_ind_over_b],xerr=TT_b[sta_ind_over_b]*0.1,yerr=bstd_b[sta_ind_over_b],color='red', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9, label='STEREO-A')
plt.errorbar(TT_b[stb_ind_over_b],bmean_b[stb_ind_over_b],xerr=TT_b[stb_ind_over_b]*0.1,yerr=bstd_b[stb_ind_over_b],color='royalblue', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9, label='STEREO-B')
plt.errorbar(TT_b[vex_ind_over_b],bmean_b[vex_ind_over_b],xerr=TT_b[vex_ind_over_b]*0.1,yerr=bstd_b[vex_ind_over_b],color='orange', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9,label='VEX')
plt.errorbar(TT_b[mes_ind_over_b],bmean_b[mes_ind_over_b],xerr=TT_b[mes_ind_over_b]*0.1,yerr=bstd_b[mes_ind_over_b],color='dimgrey', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9, label='MESSENGER')


fsize=10
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

plt.legend(loc=1,fontsize=fsize)

plt.xlabel('Predicted travel time [hours]',fontsize=fsize)
#plt.ylabel('mean B in magnetic obstacle [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


#panel labels
plt.figtext(0.02,0.93,'a',color='black', fontsize=fsize+3, ha='left',fontweight='bold')
plt.figtext(0.51,0.93,'b',color='black', fontsize=fsize+3, ha='left',fontweight='bold')

plt.tight_layout()


plt.savefig('plots/bmean_traveltime_bmean_speed.png', dpi=300)
plt.savefig('plots/bmean_traveltime_bmean_speed.pdf', dpi=300)


################################# results for Wind, VEX, MESSENGER

bmean_wind_both=np.mean(np.append(bmean_a[win_ind_over_a], bmean_b[win_ind_over_b])) 
bstd_wind_both=np.std(np.append(bmean_a[win_ind_over_a], bmean_b[win_ind_over_b])) 

print('Wind Bmean ', bmean_wind_both, ' +/- ',bstd_wind_both)

bmean_mes_both=np.mean(np.append(bmean_a[mes_ind_over_a], bmean_b[mes_ind_over_b])) 
bstd_mes_both=np.std(np.append(bmean_a[mes_ind_over_a], bmean_b[mes_ind_over_b])) 

print('MESSENGER Bmean ', bmean_mes_both, ' +/- ',bstd_mes_both)


bmean_vex_both=np.mean(np.append(bmean_a[vex_ind_over_a], bmean_b[vex_ind_over_b])) 
bstd_vex_both=np.std(np.append(bmean_a[vex_ind_over_a], bmean_b[vex_ind_over_b])) 

print('VEX Bmean ', bmean_vex_both, ' +/- ',bstd_vex_both)














sys.exit()





corr_win=np.corrcoef(tspeed_a[win_ind_over_a],bmean_a[win_ind_over_a])[0][1]
corr_mes=np.corrcoef(tspeed_a[mes_ind_over_a],bmean_a[mes_ind_over_a])[0][1]
corr_vex=np.corrcoef(tspeed_a[vex_ind_over_a],bmean_a[vex_ind_over_a])[0][1]
corr_sta=np.corrcoef(tspeed_a[sta_ind_over_a],bmean_a[sta_ind_over_a])[0][1]
corr_stb=np.corrcoef(tspeed_a[stb_ind_over_a],bmean_a[stb_ind_over_a])[0][1]


print(corr_win, corr_mes,corr_vex,corr_sta, corr_stb)








################################### Experiment for ARRCAT systematic errors in longitude

print()
print('Experiment for ARRCAT systematic errors in longitude')
print()
##############VENUS
diff_a=np.zeros(np.size(arrvexind_a))

#HIA
for k in np.arange(0,np.size(arrvexind_a)):
  predlong=target_heeq_long[arrvexind_a[k]]
  time_ind=np.max(np.where(pos_time_num < target_arrival_num[arrvexind_a[k]]))
  reallong=pos.venus[1][time_ind]*180/np.pi  
  #print(target_arrival_num[arrvexind_a[k]])
  #print(pos_time_num[time_ind])
  #print(predlong)
  #print(reallong)
  #print('diff in deg', predlong-reallong)
  #print()
  diff_a[k]=abs(predlong-reallong)

diff_b=np.zeros(np.size(arrvexind_a))


#HIB
for k in np.arange(0,np.size(arrvexind_b)):
  predlong=target_heeq_long[arrvexind_b[k]]
  time_ind=np.max(np.where(pos_time_num < target_arrival_num[arrvexind_b[k]]))
  reallong=pos.venus[1][time_ind]*180/np.pi  
  #print(target_arrival_num[arrvexind_a[k]])
  #print(pos_time_num[time_ind])
  #print(predlong)
  #print(reallong)
  #print('diff in deg', predlong-reallong)
  #print()
  diff_b[k]=abs(predlong-reallong)



print('Venus predicted arrival minus real longitude at arrival time')

print('HIA: {:.1f} +/- {:.1f}, max {:.1f} degree'.format(np.mean(diff_a),np.std(diff_a),np.max(diff_a))) 
print('HIB: {:.1f} +/- {:.1f}, max {:.1f} degree'.format(np.mean(diff_b),np.std(diff_b),np.max(diff_b))) 
both=np.concatenate((diff_a,diff_b), axis=0)
print('both: {:.1f} +/- {:.1f}, max {:.1f} degree'.format(np.mean(both),np.std(both),np.max(both))) 


#calculate difference in arrival time

di=0.72 #AU
vsse=513 #km/s
AU=149500000 #km
lamda=30 #deg


delta=15 #deg
risse=di*(1+np.sin(np.deg2rad(lamda)))/(np.cos(np.deg2rad(delta))+ np.sqrt((np.sin(np.deg2rad(lamda)))**2-(np.sin(np.deg2rad(delta)))**2)  ) #AU
tcsse=(risse-di)*AU/vsse/3600 #seconds to hours

print('tcsse 1, hours')
print(delta)
print(tcsse)

delta=15-1.8 #deg
risse=di*(1+np.sin(np.deg2rad(lamda)))/(np.cos(np.deg2rad(delta))+ np.sqrt((np.sin(np.deg2rad(lamda)))**2-(np.sin(np.deg2rad(delta)))**2)  ) #AU
tcsse=(risse-di)*AU/vsse/3600 #seconds to hours

print('tcsse 2, hours')
print(delta)
print(tcsse)


########## MESSENGER
diff_a=np.zeros(np.size(arrmesind_a))
for k in np.arange(0,np.size(arrmesind_a)):
  predlong=target_heeq_long[arrmesind_a[k]]
  time_ind=np.max(np.where(pos_time_num < target_arrival_num[arrmesind_a[k]]))
  reallong=pos.messenger[1][time_ind]*180/np.pi  
  diff_a[k]=abs(predlong-reallong)

diff_b=np.zeros(np.size(arrmesind_b))
for k in np.arange(0,np.size(arrmesind_b)):
  predlong=target_heeq_long[arrmesind_b[k]]
  time_ind=np.max(np.where(pos_time_num < target_arrival_num[arrmesind_b[k]]))
  reallong=pos.messenger[1][time_ind]*180/np.pi  
  diff_b[k]=abs(predlong-reallong)

print()
print('MESSENGER predicted arrival minus real longitude at arrival time')
print('HIA: {:.1f} +/- {:.1f}, max {:.1f} degree'.format(np.mean(diff_a),np.std(diff_a),np.max(diff_a))) 
print('HIB: {:.1f} +/- {:.1f}, max {:.1f} degree'.format(np.mean(diff_b),np.std(diff_b),np.max(diff_b))) 
both=np.concatenate((diff_a,diff_b), axis=0)
print('both: {:.1f} +/- {:.1f}, max {:.1f} degree'.format(np.mean(both),np.std(both),np.max(both))) 



#calculate difference in arrival time

di=0.35 #AU
vsse=527.5 #km/s
AU=149500000 #km
lamda=30 #deg


delta=15 #deg
risse=di*(1+np.sin(np.deg2rad(lamda)))/(np.cos(np.deg2rad(delta))+ np.sqrt((np.sin(np.deg2rad(lamda)))**2-(np.sin(np.deg2rad(delta)))**2)  ) #AU
tcsse=(risse-di)*AU/vsse/3600 #seconds to hours

print('tcsse 1, hours')
print(delta)
print(tcsse)

delta=15+4.2 #deg
risse=di*(1+np.sin(np.deg2rad(lamda)))/(np.cos(np.deg2rad(delta))+ np.sqrt((np.sin(np.deg2rad(lamda)))**2-(np.sin(np.deg2rad(delta)))**2)  ) #AU
tcsse=(risse-di)*AU/vsse/3600 #seconds to hours

print('tcsse 2, hours')
print(delta)
print(tcsse)


#################### Mars


#arrivals at Mars HIA
arrmarsind_a=np.intersect1d(np.where(target_name=='MARS'),arr_hia_obs_ind)


#arrivals at Mars HIB
arrmarsind_b=np.intersect1d(np.where(target_name=='MARS'),arr_hib_obs_ind)


########## MESSENGER
diff_a=np.zeros(np.size(arrmarsind_a))
for k in np.arange(0,np.size(arrmarsind_a)):
  predlong=target_heeq_long[arrmarsind_a[k]]
  time_ind=np.max(np.where(pos_time_num < target_arrival_num[arrmarsind_a[k]]))
  reallong=pos.mars[1][time_ind]*180/np.pi  
  diff_a[k]=abs(predlong-reallong)

diff_b=np.zeros(np.size(arrmarsind_b))
for k in np.arange(0,np.size(arrmarsind_b)):
  predlong=target_heeq_long[arrmarsind_b[k]]
  time_ind=np.max(np.where(pos_time_num < target_arrival_num[arrmarsind_b[k]]))
  reallong=pos.mars[1][time_ind]*180/np.pi  
  diff_b[k]=abs(predlong-reallong)


print()
print('MARS predicted arrival minus real longitude at arrival time')
print('HIA: {:.1f} +/- {:.1f}, max {:.1f} degree'.format(np.mean(diff_a),np.std(diff_a),np.max(diff_a))) 
print('HIB: {:.1f} +/- {:.1f}, max {:.1f} degree'.format(np.mean(diff_b),np.std(diff_b),np.max(diff_b))) 
both=np.concatenate((diff_a,diff_b), axis=0)
print('both: {:.1f} +/- {:.1f}, max {:.1f} degree'.format(np.mean(both),np.std(both),np.max(both))) 




#calculate difference in arrival time

di=1.45 #AU
vsse=527.5 #km/s
AU=149500000 #km
lamda=30 #deg


delta=15 #deg
risse=di*(1+np.sin(np.deg2rad(lamda)))/(np.cos(np.deg2rad(delta))+ np.sqrt((np.sin(np.deg2rad(lamda)))**2-(np.sin(np.deg2rad(delta)))**2)  ) #AU
tcsse=(risse-di)*AU/vsse/3600 #seconds to hours

print('tcsse 1, hours')
print(delta)
print(tcsse)

delta=15+4.2 #deg
risse=di*(1+np.sin(np.deg2rad(lamda)))/(np.cos(np.deg2rad(delta))+ np.sqrt((np.sin(np.deg2rad(lamda)))**2-(np.sin(np.deg2rad(delta)))**2)  ) #AU
tcsse=(risse-di)*AU/vsse/3600 #seconds to hours

print('tcsse 2, hours')
print(delta)
print(tcsse)



















sys.exit()
















############## NOT USED FOR PAPER FROM HERE, future study HI parameters -> Dst influence



################################################ Figure with travel time vs bmax


#Travel time here launch to in situ




#errorbars 10% on TT taken from Moestl et al. 2014
plt.figure(10)
plt.errorbar(TT_a[win_ind_over_a],bmax_a[win_ind_over_a],xerr=TT_a[win_ind_over_a]*0.1,yerr=0,color='mediumseagreen', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(TT_a[vex_ind_over_a],bmax_a[vex_ind_over_a],xerr=TT_a[vex_ind_over_a]*0.1,yerr=0,color='orange', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(TT_a[mes_ind_over_a],bmax_a[mes_ind_over_a],xerr=TT_a[mes_ind_over_a]*0.1,yerr=0,color='dimgrey', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(TT_a[sta_ind_over_a],bmax_a[sta_ind_over_a],xerr=TT_a[sta_ind_over_a]*0.1,yerr=0,color='red', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(TT_a[stb_ind_over_a],bmax_a[stb_ind_over_a],xerr=TT_a[stb_ind_over_a]*0.1,yerr=0,color='royalblue', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)

plt.errorbar(TT_b[win_ind_over_b],bmax_b[win_ind_over_b],xerr=TT_b[win_ind_over_b]*0.1,yerr=0,color='mediumseagreen', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9,label='Earth L1')
plt.errorbar(TT_b[sta_ind_over_b],bmax_b[sta_ind_over_b],xerr=TT_b[sta_ind_over_b]*0.1,yerr=0,color='red', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9, label='STEREO-A')
plt.errorbar(TT_b[stb_ind_over_b],bmax_b[stb_ind_over_b],xerr=TT_b[stb_ind_over_b]*0.1,yerr=0,color='royalblue', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9, label='STEREO-B')
plt.errorbar(TT_b[vex_ind_over_b],bmax_b[vex_ind_over_b],xerr=TT_b[vex_ind_over_b]*0.1,yerr=0,color='orange', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9,label='VEX')
plt.errorbar(TT_b[mes_ind_over_b],bmax_b[mes_ind_over_b],xerr=TT_b[mes_ind_over_b]*0.1,yerr=0,color='dimgrey', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9, label='MESSENGER')


fsize=10
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

plt.legend(loc=1,fontsize=fsize)

plt.xlabel('Travel time [hours]',fontsize=fsize)
plt.ylabel('maximum B in magnetic obstacle [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.tight_layout()


plt.savefig('plots/bmax_traveltime.png', dpi=300)
plt.savefig('plots/bmax_traveltime.pdf', dpi=300)


sys.exit()



################################### PA EXTENT and BZ


####################### for Bzmin

plt.figure(13)
plt.errorbar(bzmin_a,paext_a,xerr=0,yerr=0,color='r', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(bzmin_b,paext_b,xerr=0,yerr=0,color='b', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)

fsize=10
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

plt.legend(loc=1,fontsize=fsize)

plt.ylabel('Absolute CME extent in latitude PA [deg]',fontsize=fsize)
plt.xlabel('minimum Bz in magnetic obstacle [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.tight_layout()


plt.savefig('plots/bzmin_paextent.png', dpi=300)
plt.savefig('plots/bzmin_paextent.pdf', dpi=300)


################################### same for Bzmean


plt.figure(14)

plt.errorbar(bzmean_a/bmean_a,paext_a,xerr=0,yerr=0,color='r', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(bzmean_b/bmean_b,paext_b,xerr=0,yerr=0,color='b', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)


fsize=10
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

plt.legend(loc=1,fontsize=fsize)

plt.ylabel('Absolute CME extent in latitude PA [deg]',fontsize=fsize)
plt.xlabel('<Bz> / <B> in magnetic obstacle [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 



plt.tight_layout()


plt.savefig('plots/bzmean_paextent.png', dpi=300)
plt.savefig('plots/bzmean_paextent.pdf', dpi=300)







################################### PA center vs. bmean (latitude- middle B higher?


plt.figure(15)

plt.errorbar(bmean_a,pacenter_a,xerr=0,yerr=0,color='r', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)
plt.errorbar(bmean_b,pacenter_b-180,xerr=0,yerr=0,color='b', lw=1,capthick=1,markersize=8,fmt='o',alpha=0.9)


fsize=10
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

plt.legend(loc=1,fontsize=fsize)

plt.ylabel('PA center of CME in latitude PA [deg]',fontsize=fsize)
plt.xlabel('<B> in magnetic obstacle [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 



plt.tight_layout()


plt.savefig('plots/bmean_pacenter.png', dpi=300)
plt.savefig('plots/bmean_pacenter.pdf', dpi=300)






















































#################################### redo Moestl et al. 2014 Figures 


#Figure 12, only with Wind

wind_ind_over_a=np.where(tname_a=='EARTH_L1')
wind_ind_over_b=np.where(tname_b=='EARTH_L1')

#Travel time here launch to in situ
plt.figure(10)
plt.plot(tspeed_a[wind_ind_over_a],TTalt_a[wind_ind_over_a],'go')  
plt.plot(tspeed_b[wind_ind_over_b],TTalt_b[wind_ind_over_b],'go')  

#relationship from 2014 paper
ttdash=9537*(np.arange(200,3000))**-0.76

plt.plot(np.arange(200,3000),ttdash,'k--')  


corr_=np.corrcoef(tspeed_a[wind_ind_over_a],TTalt_a[wind_ind_over_a])



#Figure 13 only with Wind


#Travel time here launch to in situ
plt.figure(11)
plt.plot(tspeed_a[wind_ind_over_a],bmax_a[wind_ind_over_a],'go')  
plt.plot(tspeed_b[wind_ind_over_b],bmax_b[wind_ind_over_b],'go')  

#relationship from 2014 paper
bdash=0.0189*(np.arange(200,3000))+6.73

plt.plot(np.arange(200,3000),bdash,'k--')  









############################## check on new relationships for next paper

#delta vs. bmax shows barely a relationship - direction error is too large


#look at Bz too! #check extent in PA vs. Bz -> does high extent correlate with high Bz?
#mean Bz
#maybe this works instead of MVA results
#check mo_duration

# -> eigenes paper: 
#How well are factors influencing Dst predictable from HI? Speed, density, Bz
#



#latitude and bmean, bmax

plt.figure(10)
plt.plot(tspeed_a,bmax_a,'ro')  
plt.plot(tspeed_b,bmax_b,'bo')  

plt.title('Tspeed Bmax')

plt.figure(11)
plt.plot(tspeed_a,bmean_a,'ro')  
plt.plot(tspeed_b,bmean_b,'bo')  

plt.title('Tspeed Bmean')

plt.figure(12)
plt.plot(bmean_a,TT_a,'ro')  
plt.plot(bmean_b,TT_b,'bo')  

plt.title('Bmean TTalt')

plt.figure(13)
plt.plot(bmax_a,TTalt_a,'ro')  
plt.plot(bmax_b,TTalt_b,'ro')  
plt.title('Bmax TT')

plt.figure(14)
plt.plot(bmean_a,TTalt_a,'ro')  
plt.plot(bmean_b,TTalt_b,'ro')  
plt.title('Bmean TTa')

plt.figure(15)
plt.plot(bmax_a,TTalt_a,'ro')  
plt.plot(bmax_b,TTalt_b,'bo')  
plt.title('Bmax TTa')







# MIT License
# 
# Copyright (c) 2018 Christian Möstl
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.













