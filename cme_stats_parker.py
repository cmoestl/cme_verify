# cme_stats_parker.py
#
# analyses HELCATS ICMECAT data for EGU 2018 poster on CME statistics
# Author: C. Moestl, Space Research Institute IWF Graz, Austria
# started May 2015
# last update: February 2018


from scipy import stats
import scipy.io
from matplotlib import cm
import sys
import matplotlib
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import sunpy.time
import time
import pickle
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


def IDL_time_to_num(time_in):  
 #convert IDL time to matplotlib datetime
 time_num=np.zeros(np.size(time_in))
 for ii in np.arange(0,np.size(time_in)):
   time_num[ii]=mdates.date2num(sunpy.time.parse_time(time_in[ii]))   
 return time_num 
  


def gaussian(x, amp, mu, sig):
     return amp * exp(-(x-cen)**2 /wid)



#define global variables from OMNI2 dataset
#see http://omniweb.gsfc.nasa.gov/html/ow_data.html

dataset=473376;
#global Variables
spot=np.zeros(dataset) 
btot=np.zeros(dataset) #floating points
bx=np.zeros(dataset) #floating points
by=np.zeros(dataset) #floating points
bz=np.zeros(dataset) #floating points
bzgsm=np.zeros(dataset) #floating points
bygsm=np.zeros(dataset) #floating points

speed=np.zeros(dataset) #floating points
speedx=np.zeros(dataset) #floating points
speed_phi=np.zeros(dataset) #floating points
speed_theta=np.zeros(dataset) #floating points

dst=np.zeros(dataset) #float
kp=np.zeros(dataset) #float

den=np.zeros(dataset) #float
pdyn=np.zeros(dataset) #float
year=np.zeros(dataset)
day=np.zeros(dataset)
hour=np.zeros(dataset)
t=np.zeros(dataset) #index time
times1=np.zeros(dataset) #datetime time


def convertomnitime():
 #http://docs.sunpy.org/en/latest/guide/time.html
 #http://matplotlib.org/examples/pylab_examples/date_demo2.html

 print('convert time start')
 for index in range(0,dataset):
      #first to datetimeobject 
      timedum=datetime.datetime(int(year[index]), 1, 1) + datetime.timedelta(day[index] - 1) +datetime.timedelta(hours=hour[index])
      #then to matlibplot dateformat:
      times1[index] = matplotlib.dates.date2num(timedum)
      #print time
      #print year[index], day[index], hour[index]
 print('convert time done')   #for time conversion


def getomnidata():

 #statt NaN waere besser linear interpolieren
 #lese file ein:
 
 #FORMAT(2I4,I3,I5,2I3,2I4,14F6.1,F9.0,F6.1,F6.0,2F6.1,F6.3,F6.2, F9.0,F6.1,F6.0,2F6.1,F6.3,2F7.2,F6.1,I3,I4,I6,I5,F10.2,5F9.2,I3,I4,2F6.1,2I6,F5.1)
 #1963   1  0 1771 99 99 999 999 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 9999999. 999.9 9999. 999.9 999.9 9.999 99.99 9999999. 999.9 9999. 999.9 999.9 9.999 999.99 999.99 999.9  7  23    -6  119 999999.99 99999.99 99999.99 99999.99 99999.99 99999.99  0   3 999.9 999.9 99999 99999 99.9

 
 j=0
 print('start reading variables from file')
 with open('/Users/chris/python/data/omni_data/omni2_all_years.dat') as f:
  for line in f:
   line = line.split() # to deal with blank 
   #print line #41 is Dst index, in nT
   dst[j]=line[40]
   kp[j]=line[38]
   
   if dst[j] == 99999: dst[j]=np.NaN
   #40 is sunspot number
   spot[j]=line[39]
   #if spot[j] == 999: spot[j]=NaN

   #25 is bulkspeed F6.0, in km/s
   speed[j]=line[24]
   if speed[j] == 9999: speed[j]=np.NaN
 
   #get speed angles F6.1
   speed_phi[j]=line[25]
   if speed_phi[j] == 999.9: speed_phi[j]=np.NaN

   speed_theta[j]=line[26]
   if speed_theta[j] == 999.9: speed_theta[j]=np.NaN
   #convert speed to GSE x see OMNI website footnote
   speedx[j] = - speed[j] * np.cos(np.radians(speed_theta[j])) * np.cos(np.radians(speed_phi[j]))



   #9 is total B  F6.1 also fill ist 999.9, in nT
   btot[j]=line[9]
   if btot[j] == 999.9: btot[j]=np.NaN

   #GSE components from 13 to 15, so 12 to 14 index, in nT
   bx[j]=line[12]
   if bx[j] == 999.9: bx[j]=np.NaN
   by[j]=line[13]
   if by[j] == 999.9: by[j]=np.NaN
   bz[j]=line[14]
   if bz[j] == 999.9: bz[j]=np.NaN
 
   #GSM
   bygsm[j]=line[15]
   if bygsm[j] == 999.9: bygsm[j]=np.NaN
 
   bzgsm[j]=line[16]
   if bzgsm[j] == 999.9: bzgsm[j]=np.NaN 	
 
 
   #24 in file, index 23 proton density /ccm
   den[j]=line[23]
   if den[j] == 999.9: den[j]=np.NaN
 
   #29 in file, index 28 Pdyn, F6.2, fill values sind 99.99, in nPa
   pdyn[j]=line[28]
   if pdyn[j] == 99.99: pdyn[j]=np.NaN 		
 
   year[j]=line[0]
   day[j]=line[1]
   hour[j]=line[2]
   j=j+1     
 
 print('done reading variables from file')
 print(j, ' datapoints')   #for reading data from OMNI file



######################################################
#main program

plt.close('all')
print('Start catpy main program. Analyses and plots for ICME duration and planetary (Mars!) impacts')

#getomnidata()
#convertomnitime()

#-------------------------------------------------------- get cats
#filename_arrcat='ALLCATS/HELCATS_ARRCAT_v6.sav'
#a=getcat(filename_arrcat)

filename_icmecat='ALLCATS/HELCATS_ICMECAT_v11_SCEQ.sav'
i=getcat(filename_icmecat)

#filename_linkcat='ALLCATS/HELCATS_LINKCAT_v10.sav'
#l=getcat(filename_linkcat)


#now this is a structured array  
#access each element of the array see http://docs.scipy.org/doc/numpy/user/basics.rec.html
#access variables
#i.icmecat['id']
#look at contained variables
#print(a.arrcat.dtype)
#print(i.icmecat.dtype)


#get spacecraft and planet positions
pos=getcat('../catpy/DATACAT/positions_2007_2018_HEEQ_6hours.sav')
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
mo_speed_st=i.icmecat['MO_SPEED_STD']
sheath_density=i.icmecat['SHEATH_DENSITY']
sheath_density_std=i.icmecat['SHEATH_DENSITY_STD']
mo_density=i.icmecat['MO_DENSITY']
mo_density_std=i.icmecat['MO_DENSITY_STD']
sheath_temperature=i.icmecat['SHEATH_TEMPERATURE']
sheath_temperature_std=i.icmecat['SHEATH_TEMPERATURE_STD']
mo_temperature=i.icmecat['MO_TEMPERATURE']
mo_temperature_std=i.icmecat['MO_TEMPERATURE_STD']


#get indices of events in different spacecraft
ivexind=np.where(isc == 'VEX')
istaind=np.where(isc == 'STEREO-A')
istbind=np.where(isc == 'STEREO-B')
iwinind=np.where(isc == 'Wind')
imesind=np.where(isc == 'MESSENGER')
iulyind=np.where(isc == 'ULYSSES')
imavind=np.where(isc == 'MAVEN')


#take MESSENGER only at Mercury, only events after orbit insertion
imercind=np.where(np.logical_and(isc =='MESSENGER',icme_start_time_num > mdates.date2num(sunpy.time.parse_time('2011-03-18'))))

#limits of solar minimum, rising phase and solar maximum

minstart=mdates.date2num(sunpy.time.parse_time('2007-01-01'))
minend=mdates.date2num(sunpy.time.parse_time('2009-12-31'))

risestart=mdates.date2num(sunpy.time.parse_time('2010-01-01'))
riseend=mdates.date2num(sunpy.time.parse_time('2011-06-30'))

maxstart=mdates.date2num(sunpy.time.parse_time('2011-07-01'))
maxend=mdates.date2num(sunpy.time.parse_time('2014-12-31'))

#extract events by limits of solar min, rising, max, too few events for MAVEN and Ulysses

iallind_min=np.where(np.logical_and(icme_start_time_num > minstart,icme_start_time_num < minend))[0]
iallind_rise=np.where(np.logical_and(icme_start_time_num > risestart,icme_start_time_num < riseend))[0]
iallind_max=np.where(np.logical_and(icme_start_time_num > maxstart,icme_start_time_num < maxend))[0]

iwinind_min=iallind_min[np.where(isc[iallind_min]=='Wind')]
iwinind_rise=iallind_rise[np.where(isc[iallind_rise]=='Wind')]
iwinind_max=iallind_max[np.where(isc[iallind_max]=='Wind')]

ivexind_min=iallind_min[np.where(isc[iallind_min]=='VEX')]
ivexind_rise=iallind_rise[np.where(isc[iallind_rise]=='VEX')]
ivexind_max=iallind_max[np.where(isc[iallind_max]=='VEX')]

imesind_min=iallind_min[np.where(isc[iallind_min]=='MESSENGER')]
imesind_rise=iallind_rise[np.where(isc[iallind_rise]=='MESSENGER')]
imesind_max=iallind_max[np.where(isc[iallind_max]=='MESSENGER')]

istaind_min=iallind_min[np.where(isc[iallind_min]=='STEREO-A')]
istaind_rise=iallind_rise[np.where(isc[iallind_rise]=='STEREO-A')]
istaind_max=iallind_max[np.where(isc[iallind_max]=='STEREO-A')]

istbind_min=iallind_min[np.where(isc[iallind_min]=='STEREO-B')]
istbind_rise=iallind_rise[np.where(isc[iallind_rise]=='STEREO-B')]
istbind_max=iallind_max[np.where(isc[iallind_max]=='STEREO-B')]


#imercind_min=iallind_min[np.where(isc[iallind_min]=='VEX')]
#imercind_rise=iallind_min[np.where(isc[iallind_rise]=='VEX')]
#imercind_max=iallind_min[np.where(isc[iallind_max]=='VEX')]
#old: use these indices like  mo_bmean[imercind][imercind_rise] to get all MO_BMEAN at Mercury in the rising phase



Rs_in_AU=7e5/149.5e6























###################################################################################

##################### (1) DURATION PLOT and linear fit  ############################




sns.set_context("talk")     
#sns.set_style("darkgrid")  
sns.set_style("ticks",{'grid.linestyle': '--'})

fig=plt.figure(1,figsize=(12,11	))
fsize=15
ax1 = plt.subplot2grid((2,1), (0, 0))


xfit=np.linspace(0,2,1000)
icme_durations=(mo_end_time_num-icme_start_time_num)*24 #hours

#make linear fits - no forcing through origin
durfit=np.polyfit(sc_heliodistance,icme_durations,1)
durfitmin=np.polyfit(sc_heliodistance[iallind_min],icme_durations[iallind_min],1)
durfitrise=np.polyfit(sc_heliodistance[iallind_rise],icme_durations[iallind_rise],1)
durfitmax=np.polyfit(sc_heliodistance[iallind_max],icme_durations[iallind_max],1)


#force through origin, fit with y=kx
scx=sc_heliodistance[:,np.newaxis]
durfit_f, _, _, _ =np.linalg.lstsq(scx,icme_durations)
scx=sc_heliodistance[iallind_min][:,np.newaxis]
durfitmin_f, _, _, _ =np.linalg.lstsq(scx,icme_durations[iallind_min])
scx=sc_heliodistance[iallind_rise][:,np.newaxis]
durfitrise_f, _, _, _ =np.linalg.lstsq(scx,icme_durations[iallind_rise])
scx=sc_heliodistance[iallind_max][:,np.newaxis]
durfitmax_f, _, _, _ =np.linalg.lstsq(scx,icme_durations[iallind_max])


#this is similar to D=durfit[0]*xfit+durfit[1]
durfitall=np.poly1d(durfit)
durfitmin=np.poly1d(durfitmin)
durfitrise=np.poly1d(durfitrise)
durfitmax=np.poly1d(durfitmax)


#make the y axis for the fits forced through the origin
ydurfitall_f=durfit_f*xfit
ydurfitmin_f=durfitmin_f*xfit
ydurfitrise_f=durfitrise_f*xfit
ydurfitmax_f=durfitmax_f*xfit

#for fit plotting
print('ICME duration linear function: D[hours]={:.2f}r[AU]+{:.2f}'.format(durfit[0],durfit[1]))



plt.plot(sc_heliodistance,icme_durations,'o',color='blue',markersize=5, alpha=0.3,label='D')
#plt.plot(sc_heliodistance[iallind_min],icme_durations[iallind_min],'o',color='dimgrey',markersize=3, alpha=0.4,label='D min')
#plt.plot(sc_heliodistance[iallind_rise],icme_durations[iallind_rise],'o',color='grey',markersize=3, alpha=0.7,label='D rise')
#plt.plot(sc_heliodistance[iallind_max],icme_durations[iallind_max],'o',color='black',markersize=3, alpha=0.8,label='D max')



#plot fits

plt.plot(xfit,ydurfitall_f,'-',color='blue', lw=2.5, alpha=0.9,label='fit')
plt.plot(xfit,ydurfitmin_f,'--',color='black', lw=2, alpha=0.9,label='min fit')
plt.plot(xfit,ydurfitrise_f,'-.',color='black', lw=2, alpha=0.9,label='rise fit')
plt.plot(xfit,ydurfitmax_f,'-',color='black', lw=2, alpha=0.9,label='max fit')


#these don't go through the origin
#plt.plot(xfit,durfitall(xfit),'-',color='blue', lw=2.5, alpha=0.9,label='fit')
#plt.plot(xfit,durfitmin(xfit),'--',color='black', lw=2, alpha=0.9,label='min fit')
#plt.plot(xfit,durfitrise(xfit),'-.',color='black', lw=2, alpha=0.9,label='rise fit')
#plt.plot(xfit,durfitmax(xfit),'-',color='black', lw=2, alpha=0.9,label='max fit')


plt.annotate('overall: D[h]={:.2f} R[AU] '.format(durfit_f[0]),xy=(0.1,60),fontsize=11)
plt.annotate('minimum: D[h]={:.2f} R[AU] '.format(durfitmin_f[0]),xy=(0.1,55),fontsize=11)
plt.annotate('rising phase: D[h]={:.2f} R[AU]'.format(durfitrise_f[0]),xy=(0.1,50),fontsize=11)
plt.annotate('maximum: D[h]={:.2f} R[AU]'.format(durfitmax_f[0]),xy=(0.1,45),fontsize=11)



#plt.annotate('overall: D[h]={:.2f} R[AU] + {:.2f}'.format(durfitall[0],durfitall[1]),xy=(0.1,120),fontsize=12)
#plt.annotate('minimum: D[h]={:.2f} R[AU] + {:.2f}'.format(durfitmin[0],durfitmin[1]),xy=(0.1,100),fontsize=12)
#plt.annotate('rising phase: D[h]={:.2f} R[AU] + {:.2f}'.format(durfitrise[0],durfitrise[1]),xy=(0.1,80),fontsize=12)
#plt.annotate('maximum: D[h]={:.2f} R[AU] + {:.2f}'.format(durfitmax[0],durfitmax[1]),xy=(0.1,60),fontsize=12)


#planet limits
plt.axvspan(np.min(pos.mars[0]),np.max(pos.mars[0]), color='orangered', alpha=0.2)
plt.axvspan(np.min(pos.mercury[0]),np.max(pos.mercury[0]), color='darkgrey', alpha=0.2)
plt.axvspan(np.min(pos.venus[0]),np.max(pos.venus[0]), color='orange', alpha=0.2)
plt.axvspan(np.min(pos.earth[0]),np.max(pos.earth[0]), color='mediumseagreen', alpha=0.2)
plt.axvspan(0.044,0.3,color='magenta', alpha=0.2)

plt.annotate('Mars', xy=(1.5,65), ha='center',fontsize=fsize)
plt.annotate('PSP', xy=(0.11,65), ha='center',fontsize=fsize)
plt.annotate('Mercury', xy=(0.38,65), ha='center',fontsize=fsize)
plt.annotate('Venus', xy=(0.72,65), ha='center',fontsize=fsize)
plt.annotate('Earth', xy=(1,65), ha='center',fontsize=fsize)

ax1.set_xticks(np.arange(0,2,0.2))


plt.xlim(0,max(sc_heliodistance)+0.3)
#plt.ylim(0,max(icme_durations)+30)
plt.ylim(0,70)


plt.legend(loc=4,fontsize=fsize-1)

plt.xlabel('Heliocentric distance R [AU]',fontsize=fsize)
plt.ylabel('ICME duration D [hours]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


#plt.grid()


print('results for durations and distance for PSP from 0.04 to 0.3 AU')

print('all fits')
psp_dist=np.where(np.logical_and(xfit < 0.3, xfit > 0.04))
#psp_durs=durfitall(xfit[psp_dist])

psp_durs=durfit_f*xfit[psp_dist]
pspdur_mean=np.mean(psp_durs)
pspdur_std=np.std(psp_durs)
pspdur_min=np.min(psp_durs)
pspdur_max=np.max(psp_durs)

print('Parker predicted durations from 0.044 to 0.3 AU: mean +/ std, min, max')
print(pspdur_mean, ' +/- ',pspdur_std)
print(pspdur_min, ' to ',pspdur_max)

print('min fits')

#psp_dursmin=durfitmin(xfit[psp_dist])

psp_dursmin=durfitmin_f*xfit[psp_dist]
np.mean(psp_dursmin)
np.min(psp_dursmin)
np.max(psp_dursmin)

print('rise fits')
#psp_dursrise=durfitrise(xfit[psp_dist])

psp_dursrise=durfitrise_f*xfit[psp_dist]

np.mean(psp_dursrise)
np.min(psp_dursrise)
np.max(psp_dursrise)

print('max fits')
psp_dursmax=durfitmax_f*xfit[psp_dist]
np.mean(psp_dursmax)
np.min(psp_dursmax)
np.max(psp_dursmax)




############# plot 2 vs time



#exclude STEREO for better visibility
#plt.plot_date(icme_start_time_num[istbind],np.log10(mo_bmean[istbind]),'o',color='royalblue',markersize=markers,linestyle='-',linewidth=linew)
#plt.plot_date(icme_start_time_num[istaind],np.log10(mo_bmean[istaind]),'o',color='red',markersize=markers,linestyle='-',linewidth=linew)

#Wind
tfit=mdates.date2num(sunpy.time.parse_time('2009-04-01'))+np.arange(0,365*10)
t0=mdates.date2num(sunpy.time.parse_time('2009-01-01'))

#is a gaussian better?
#sigma=1000
#bfitmax=30
#mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
#ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')

#or is this better, like sunspot cycle?
#Hathaway 2015 equation 6 page 40
#average cycle sunspot number 
A=100 #amplitude ##195 for sunspot
b=100*12 #56*12 for months to days
c=0.8

#4 free parameters A, b, c, t0

Fwind=A*(((tfit-t0)/b)**3) * 1/(np.exp((((tfit-t0)/b)**2))-c)
#plt.plot_date(tfit, Fwind,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')

#xaxis: 10 years, daily data point
xfit2=mdates.date2num(sunpy.time.parse_time('2007-01-01'))+np.arange(0,365*10)
#MESSENGER
sigma=1000
bfitmax=10
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))

ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit2-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='darkgrey',linestyle='-',markersize=0, label='Mercury fit')

#VEX
#inital guess
sigma=1000
bfitmax=20
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit2-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit2, ygauss*1/max(ygauss)*bfitmax,'o',color='orange',linestyle='-',markersize=0, label='Venus fit')

#for Mars: reconstruct likely parameters if sigma is quite similar for all fits, take mean of those sigmas and adjust bfitmax as function of distance with power law)
#plot reconstructed function for Mars
bfitmax=40
#plt.plot_date(xfit2, Fwind,'o',color='steelblue',linestyle='--',markersize=0, label='Mars reconstr.')




#####plot 

ax2 = plt.subplot2grid((2,1), (1, 0))


markers=6
linew=0




ax2.plot_date(icme_start_time_num[imesind],icme_durations[imesind],'o',color='darkgrey',markersize=markers,linestyle='-',linewidth=linew,label='MESSENGER')
ax2.plot_date(icme_start_time_num[ivexind],icme_durations[ivexind],'o',color='orange',markersize=markers,linestyle='-',linewidth=linew, label='Venus')
ax2.plot_date(icme_start_time_num[iwinind],icme_durations[iwinind],'o',color='mediumseagreen',markersize=markers, linestyle='-', linewidth=linew, label='Earth')
ax2.plot_date(icme_start_time_num[imavind],icme_durations[imavind],'o',color='steelblue',markersize=markers,linestyle='-',linewidth=linew, label='Mars')




#ax3 = ax2.twinx()
#ax3.plot_date(times1, spot, '-', color='black', alpha=0.5)
#ax3.set_ylabel('Sunspot number')



#limits solar min/rise/max, and means as horizontal lines for each sub interval

vlevel=60

plt.axvspan(minstart,minend, color='green', alpha=0.1)
plt.annotate('solar minimum',xy=(minstart+(minend-minstart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(minstart+10,vlevel),ha='left')
plt.annotate('>',xy=(minend-10,vlevel),ha='right')

plt.plot_date( [minstart,minend], [np.mean(icme_durations[iwinind_min]),np.mean(icme_durations[iwinind_min])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [minstart,minend], [np.mean(icme_durations[ivexind_min]),np.mean(icme_durations[ivexind_min])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [minstart,minend], [np.mean(icme_durations[imesind_min]),np.mean(icme_durations[imesind_min])], color='darkgrey', linestyle='-', markersize=0) 



plt.axvspan(risestart,riseend, color='yellow', alpha=0.1)
plt.annotate('rising phase',xy=(risestart+(riseend-risestart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(risestart+10,vlevel),ha='left')
plt.annotate('>',xy=(riseend-10,vlevel),ha='right')


plt.plot_date( [risestart,riseend], [np.mean(icme_durations[iwinind_rise]),np.mean(icme_durations[iwinind_rise])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [risestart,riseend], [np.mean(icme_durations[ivexind_rise]),np.mean(icme_durations[ivexind_rise])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [risestart,riseend], [np.mean(icme_durations[imesind_rise]),np.mean(icme_durations[imesind_rise])], color='darkgrey', linestyle='-', markersize=0) 


plt.axvspan(maxstart,maxend, color='red', alpha=0.1)
plt.annotate('solar maximum',xy=(maxstart+(maxend-maxstart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(maxstart+10,vlevel),ha='left')
plt.annotate('>',xy=(maxend,vlevel),ha='right')


plt.plot_date( [maxstart,maxend], [np.mean(icme_durations[iwinind_max]),np.mean(icme_durations[iwinind_max])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [maxstart,maxend], [np.mean(icme_durations[ivexind_max]),np.mean(icme_durations[ivexind_max])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [maxstart,maxend], [np.mean(icme_durations[imesind_max]),np.mean(icme_durations[imesind_max])], color='darkgrey', linestyle='-', markersize=0) 



plt.ylim(0,70)
plt.xlim(mdates.date2num(sunpy.time.parse_time('2007-01-01')), mdates.date2num(sunpy.time.parse_time('2016-12-31')))

plt.ylabel('ICME duration D [hours]',fontsize=fsize)
plt.xlabel('year',fontsize=fsize)

plt.tight_layout()

plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.legend(loc=4,fontsize=fsize-1)


#panel labels
plt.figtext(0.01,0.98,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.01,0.485,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')

plt.show()
plt.savefig('plots_psp/icme_durations_distance_time_parker.pdf', dpi=300)
plt.savefig('plots_psp/icme_durations_distance_time_parker.png', dpi=300)




#results on durations




#################################################

#D

# 
# print()
# print()
# print('--------------------------------------------------')
# print()
# print()
# 
# print('DURATION ')
# 
# print()
# print('MESSENGER +/-')
# print(round(np.mean(icme_durations[imercind]),1))
# print(round(np.std(icme_durations[imercind]),1))
# 
# #print('min')
# #np.mean(icme_durations[imercind][imercind_min])
# #np.std(icme_durations[imercind][imercind_min])
# print('rise')
# print(round(np.mean(icme_durations[imesind][imercind_rise]),1))
# print(round(np.std(icme_durations[imercind][imercind_rise]),1))
# print('max')
# print(round(np.mean(icme_durations[imercind][imercind_max]),1))
# print(round(np.std(icme_durations[imercind][imercind_max]),1))
# 
# 
# print()
# print('Venus')
# print(round(np.mean(icme_durations[ivexind]),1))
# print(round(np.std(icme_durations[ivexind]),1))
# print('min')
# print(round(np.mean(icme_durations[ivexind][ivexind_min]),1))
# print(round(np.std(icme_durations[ivexind][ivexind_min]),1))
# print('rise')
# print(round(np.mean(icme_durations[ivexind][ivexind_rise]),1))
# print(round(np.std(icme_durations[ivexind][ivexind_rise]),1))
# print('max')
# print(round(np.mean(icme_durations[ivexind][ivexind_max]),1))
# print(round(np.std(icme_durations[ivexind][ivexind_max]),1))
# 
# print()
# print('Earth')
# print(round(np.mean(icme_durations[iwinind]),1))
# print(round(np.std(icme_durations[iwinind]),1))
# print('min')
# print(round(np.mean(icme_durations[iwinind][iwinind_min]),1))
# print(round(np.std(icme_durations[iwinind][iwinind_min]),1))
# print('rise')
# print(round(np.mean(icme_durations[iwinind][iwinind_rise]),1))
# print(round(np.std(icme_durations[iwinind][iwinind_rise]),1))
# print('max')
# print(round(np.mean(icme_durations[iwinind][iwinind_max]),1))
# print(round(np.std(icme_durations[iwinind][iwinind_max]),1))
# 
# print()
# 
# 
# #only declining phase
# print('MAVEN')
# print(round(np.mean(icme_durations[imavind]),1))
# print(round(np.std(icme_durations[imavind]),1))
# 
# 
# 
# 
# 
# 
# 
# 

















###################################################################################

##################### (2) Bfield plot ICMECAT  ############################



#-------------------------------------------------------------- Bfield plot


sns.set_context("talk")     
#sns.set_style("darkgrid")  
sns.set_style("ticks",{'grid.linestyle': '--'})

fig=plt.figure(2,figsize=(12,12	))
#fig=plt.figure(2,figsize=(12,6	))
fsize=15

ax1 = plt.subplot2grid((2,2), (0, 0))
#ax1 = plt.subplot2grid((1,2), (0, 0))


xfit=np.linspace(0,2,1000)

#power law fits for all events
bmaxfit=np.polyfit(np.log10(sc_heliodistance),np.log10(mo_bmax),1)
b=10**bmaxfit[1]
bmaxfitfun=b*(xfit**bmaxfit[0])
print('exponent for bmax fit:', bmaxfit[0])

bmeanfit=np.polyfit(np.log10(sc_heliodistance),np.log10(mo_bmean),1)
b=10**bmeanfit[1]
bmeanfitfun=b*(xfit**bmeanfit[0])
print('exponent for bmean fit:', bmeanfit[0])



################ dont take next 3 fits too seriously -> too few events for other distances during min for example (only VEX/Wind)
##fit with only minimum events
bmeanfit_min=np.polyfit(np.log10(sc_heliodistance[iallind_min]),np.log10(mo_bmean[iallind_min]),1)
bmeanfitfun_min=(10**bmeanfit_min[1])*(xfit**bmeanfit_min[0])
print('exponent for bmean_min fit:', bmeanfit_min[0])

##fit with only rising events
bmeanfit_rise=np.polyfit(np.log10(sc_heliodistance[iallind_rise]),np.log10(mo_bmean[iallind_rise]),1)
bmeanfitfun_rise=(10**bmeanfit_rise[1])*(xfit**bmeanfit_rise[0])
print('exponent for bmean_rise fit:', bmeanfit_rise[0])

##fit with only maximum events
bmeanfit_max=np.polyfit(np.log10(sc_heliodistance[iallind_max]),np.log10(mo_bmean[iallind_max]),1)
bmeanfitfun_max=(10**bmeanfit_max[1])*(xfit**bmeanfit_max[0])
print('exponent for bmean_max fit:', bmeanfit_max[0])






plt.plot(sc_heliodistance,mo_bmean,'o',color='black',markersize=5, alpha=0.7,label='$\mathregular{<B>}$')
plt.plot(xfit,bmeanfitfun,'-',color='black', lw=2, alpha=0.7,label='$\mathregular{<B> \\ fit}$')

plt.plot(sc_heliodistance,mo_bmax,'o',color='dodgerblue',markersize=5, alpha=0.7,label='$\mathregular{B_{max}}$')
plt.plot(xfit,bmaxfitfun,'-',color='dodgerblue', lw=2, alpha=0.7,label='$\mathregular{B_{max} \\ fit}$')


plt.text(1.1,120,'$\mathregular{<B> [nT]= 8.9 R[AU]^{-1.68}}$', fontsize=10)
plt.text(1.1,100,'$\mathregular{B_{max} [nT]= 12.3 R[AU]^{-1.73}}$', fontsize=10)


#plt.annotate('$Bmax[nT]={:.1f} R[AU]^{:.2f} $'.format(10**bmeanfit[1],bmeanfit[0]),xy=(0.9,80),fontsize=9)
#plt.annotate('$\mathregular{<B>}$: B[nT]={:.1f} R[AU]^{:.2f} '.format(bmeanfit[0],10**bmeanfit[1]),xy=(0.1,60),fontsize=11)
#plt.annotate('$\mathregular{B_{max}}$: B[nT]={:.1f} R[AU]^{:.2f} '.format(bmaxfit[0],10**bmaxfit[1]),xy=(0.1,60),fontsize=11)




#mars limits
plt.axvspan(np.min(pos.mars[0]),np.max(pos.mars[0]), color='orangered', alpha=0.2)
#plt.figtext(0.8,0.8,'Mars',color='orangered')

plt.axvspan(np.min(pos.mercury[0]),np.max(pos.mercury[0]), color='darkgrey', alpha=0.2)
#plt.figtext(0.25,0.8,'Mercury',color='darkgrey')

plt.axvspan(np.min(pos.venus[0]),np.max(pos.venus[0]), color='orange', alpha=0.2)
#plt.figtext(0.42,0.8,'Venus',color='orange')

plt.axvspan(np.min(pos.earth[0]),np.max(pos.earth[0]), color='mediumseagreen', alpha=0.2)
#plt.figtext(0.6,0.8,'Earth',color='mediumseagreen')

#solar probe plus 10 to 36 Rs close approaches
Rs_in_AU=7e5/149.5e6
plt.axvspan(Rs_in_AU*10,Rs_in_AU*36,color='magenta', alpha=0.2)


#plt.figtext(0.65,0.2,' D[h]={:.2f} R[AU] + {:.2f}'.format(durfit[0],durfit[1]))
plt.xlim(0,1.8)
plt.ylim(0,max(mo_bmax)+20)



plt.legend(loc=1,fontsize=fsize)

plt.xlabel('Heliocentric distance R [AU]',fontsize=fsize)
plt.ylabel('Magnetic field in MO B [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 
#plt.grid()



######################## logarithmic plot with Sun


#for the bmean fit, append one value for the coronal field at 0.007 AU for 1.5 Rs with 1 Gauss or 10^5 nT
#or better
#patsourakos georgoulis 2016: 0.03 G for 10 Rs #10^5 nT is 1 Gauss
mo_bmean_sun=np.append(mo_bmean,10**5*0.03) 
mo_bmax_sun=np.append(mo_bmax,10**5*0.03) 


sc_heliodistance_sun=np.append(sc_heliodistance,10*Rs_in_AU)

ax3 = plt.subplot2grid((2,2), (0, 1))

bmeanfit_sun=np.polyfit(np.log10(sc_heliodistance_sun),np.log10(mo_bmean_sun),1)
b=10**bmeanfit_sun[1]
bmeanfitfun_sun=b*(xfit**bmeanfit_sun[0])
print('exponent for bmean fit sun:', bmeanfit_sun[0])


bmaxfit_sun=np.polyfit(np.log10(sc_heliodistance_sun),np.log10(mo_bmax_sun),1)
b=10**bmaxfit_sun[1]
bmaxfitfun_sun=b*(xfit**bmaxfit_sun[0])
print('exponent for bmean fit sun:', bmaxfit_sun[0])


plt.plot(sc_heliodistance_sun,np.log10(mo_bmean_sun),'o',color='black',markersize=5, alpha=0.7,label='$\mathregular{<B>}$')
plt.plot(xfit,np.log10(bmeanfitfun_sun),'-',color='black', lw=2, alpha=0.7,label='$\mathregular{<B> fit}$')
plt.plot(xfit,np.log10(bmaxfitfun_sun),'-',color='dodgerblue', lw=2, alpha=0.7,label='$\mathregular{B_{max} fit}$')

plt.ylim(0,6)



plt.text(1.1,3,'$\mathregular{<B> [nT]= 8.9 R[AU]^{-1.70}}$', fontsize=10)
plt.text(1.1,2.5,'$\mathregular{B_{max} [nT]= 12.3 R[AU]^{-1.79}}$', fontsize=10)


ax3.annotate('Mars', xy=(1.5,4), ha='center',fontsize=fsize-2)
ax3.annotate('PSP', xy=(0.12,4), ha='center',fontsize=fsize-2)
ax3.annotate('Mercury', xy=(0.38,4), ha='center',fontsize=fsize-2)
ax3.annotate('Venus', xy=(0.72,4), ha='center',fontsize=fsize-2)
ax3.annotate('Earth', xy=(1,4), ha='center',fontsize=fsize-2)


plt.legend(loc=1,fontsize=fsize)

plt.xlabel('Heliocentric distance R [AU]',fontsize=fsize)
plt.ylabel('Magnetic field in MO log(B) [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


#mars limits
plt.axvspan(np.min(pos.mars[0]),np.max(pos.mars[0]), color='orangered', alpha=0.2)
#plt.figtext(0.8,0.8,'Mars',color='orangered')
plt.axvspan(np.min(pos.mercury[0]),np.max(pos.mercury[0]), color='darkgrey', alpha=0.2)
#plt.figtext(0.25,0.8,'Mercury',color='darkgrey')
plt.axvspan(np.min(pos.venus[0]),np.max(pos.venus[0]), color='orange', alpha=0.2)
#plt.figtext(0.42,0.8,'Venus',color='orange')
plt.axvspan(np.min(pos.earth[0]),np.max(pos.earth[0]), color='mediumseagreen', alpha=0.2)
#plt.figtext(0.6,0.8,'Earth',color='mediumseagreen')
plt.xlim(0,1.8)

#PSP from 0.044 AU to 0.3 AU unknown territory
plt.axvspan(0.044,0.3, color='magenta', alpha=0.2)



#panel labels
plt.figtext(0.03,0.96,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.515,0.96,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.03,0.49,'c',color='black', fontsize=fsize, ha='left',fontweight='bold')



#sns.despine()
#plt.tight_layout()



print('results for bmean in MO and distance for PSP from 0.04 to 0.3 AU')
print('all fits')
psp_dist=np.where(np.logical_and(xfit < 0.3, xfit > 0.04))
#calculate function only with psp_dist values
psp_b=10**bmeanfit_sun[1]*(xfit[psp_dist]**bmeanfit_sun[0])
pspbmean_mean=np.mean(psp_b)
pspbmean_std=np.std(psp_b)
pspbmean_min=np.min(psp_b)
pspbmean_max=np.max(psp_b)
print('Parker predicted mean field in MO from 0.044 to 0.3 AU: mean +/ std, min, max')
print(pspbmean_mean, ' +/- ',pspbmean_std)
print(pspbmean_min, ' to ',pspbmean_max)



print('results for bmean in MO and distance for PSP from 0.04 to 0.3 AU')
print('all fits')
psp_dist=np.where(np.logical_and(xfit < 0.3, xfit > 0.04))
#calculate function only with psp_dist values
psp_b=10**bmaxfit_sun[1]*(xfit[psp_dist]**bmaxfit_sun[0])
pspbmax_mean=np.mean(psp_b)
pspbmax_std=np.std(psp_b)
pspbmax_min=np.min(psp_b)
pspbmax_max=np.max(psp_b)
print('Parker predicted max field in MO from 0.044 to 0.3 AU: mean +/ std, min, max')
print(pspbmax_mean, ' +/- ',pspbmax_std)
print(pspbmax_min, ' to ',pspbmax_max)









################################# B vs. time


ax2 = plt.subplot2grid((2,2), (1, 0), colspan=2)

markers=6
linew=0


plt.plot_date(icme_start_time_num[imesind],mo_bmean[imesind],'o',color='darkgrey',markersize=markers,linestyle='-',linewidth=linew,label='MESSENGER')
plt.plot_date(icme_start_time_num[ivexind],mo_bmean[ivexind],'o',color='orange',markersize=markers,linestyle='-',linewidth=linew, label='Venus')
plt.plot_date(icme_start_time_num[iwinind],mo_bmean[iwinind],'o',color='mediumseagreen',markersize=markers, linestyle='-', linewidth=linew, label='Earth')
plt.plot_date(icme_start_time_num[imavind],mo_bmean[imavind],'o',color='steelblue',markersize=markers,linestyle='-',linewidth=linew, label='Mars')


#exclude STEREO for better visibility
#plt.plot_date(icme_start_time_num[istbind],np.log10(mo_bmean[istbind]),'o',color='royalblue',markersize=markers,linestyle='-',linewidth=linew)
#plt.plot_date(icme_start_time_num[istaind],np.log10(mo_bmean[istaind]),'o',color='red',markersize=markers,linestyle='-',linewidth=linew)





#add gaussian fits for MESSENGER, VEX, Wind (MAVEN too few data points)


##########just gaussian, no fit yet


#instead of gaussian, fit solar cycle functions in Hathaway 2015 solar cycle living reviews equation 6

#xaxis: 10 years, daily data point
xfit=mdates.date2num(sunpy.time.parse_time('2007-01-01'))+np.arange(0,365*10)

#MESSENGER
sigma=1000
bfitmax=80
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))

ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='darkgrey',linestyle='-',markersize=0, label='Mercury fit')


#VEX
#inital guess
sigma=1000
bfitmax=40
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='orange',linestyle='-',markersize=0, label='Venus fit')

#Wind
sigma=1000
bfitmax=10
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')


#for Mars: reconstruct likely parameters if sigma is quite similar for all fits, take mean of those sigmas and adjust bfitmax as function of distance with power law)
#plot reconstructed function for Mars
bfitmax=6
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='steelblue',linestyle='--',markersize=0, label='Mars reconstr.')


plt.legend(loc=1,fontsize=fsize-2)

#limits solar min/rise/max

vlevel=150

plt.axvspan(minstart,minend, color='green', alpha=0.1)
plt.annotate('solar minimum',xy=(minstart+(minend-minstart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(minstart+10,vlevel),ha='left')
plt.annotate('>',xy=(minend-10,vlevel),ha='right')


plt.axvspan(risestart,riseend, color='yellow', alpha=0.1)
plt.annotate('rising phase',xy=(risestart+(riseend-risestart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(risestart+10,vlevel),ha='left')
plt.annotate('>',xy=(riseend-10,vlevel),ha='right')

plt.axvspan(maxstart,maxend, color='red', alpha=0.1)
plt.annotate('solar maximum',xy=(maxstart+(maxend-maxstart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(maxstart+10,vlevel),ha='left')
plt.annotate('>',xy=(maxend,vlevel),ha='right')



plt.plot_date( [minstart,minend], [np.mean(mo_bmean[iwinind_min]),np.mean(mo_bmean[iwinind_min])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [minstart,minend], [np.mean(mo_bmean[ivexind_min]),np.mean(mo_bmean[ivexind_min])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [minstart,minend], [np.mean(mo_bmean[imesind_min]),np.mean(mo_bmean[imesind_min])], color='darkgrey', linestyle='-', markersize=0) 


plt.plot_date( [risestart,riseend], [np.mean(mo_bmean[iwinind_rise]),np.mean(mo_bmean[iwinind_rise])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [risestart,riseend], [np.mean(mo_bmean[ivexind_rise]),np.mean(mo_bmean[ivexind_rise])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [risestart,riseend], [np.mean(mo_bmean[imesind_rise]),np.mean(mo_bmean[imesind_rise])], color='darkgrey', linestyle='-', markersize=0) 


plt.plot_date( [maxstart,maxend], [np.mean(mo_bmean[iwinind_max]),np.mean(mo_bmean[iwinind_max])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [maxstart,maxend], [np.mean(mo_bmean[ivexind_max]),np.mean(mo_bmean[ivexind_max])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [maxstart,maxend], [np.mean(mo_bmean[imesind_max]),np.mean(mo_bmean[imesind_max])], color='darkgrey', linestyle='-', markersize=0) 



plt.ylabel('Magnetic field in MO [nT]', fontsize=fsize)

plt.xlabel('Year', fontsize=fsize)

#sets planet / spacecraft labels
xoff=0.15
yoff=0.8
fsize=14

#plt.figtext(xoff,yoff,'Earth L1',color='mediumseagreen', fontsize=fsize, ha='left')
#plt.figtext(xoff,yoff-0.03*1,'VEX',color='orange', fontsize=fsize, ha='left')
#plt.figtext(xoff,yoff-0.03*2,'MESSENGER',color='dimgrey', fontsize=fsize, ha='left')
#plt.figtext(xoff,yoff-0.03*3,'STEREO-A',color='red', fontsize=fsize, ha='left')
#plt.figtext(xoff,yoff-0.03*4,'STEREO-B',color='royalblue', fontsize=fsize, ha='left')
#plt.figtext(xoff,yoff-0.03*5,'MAVEN',color='steelblue', fontsize=fsize, ha='left')


#plt.ylim(0,45)
#plt.xlim(yearly_start_times[0],yearly_end_times[9])

#plt.grid()
plt.tight_layout()

#plt.show()
plt.savefig('plots_psp/icme_total_field_distance_time.pdf', dpi=300)
plt.savefig('plots_psp/icme_total_field_distance_time.png', dpi=300)











































###############################################################
##########time spent inside ICMEs, in %
##########################################################


icme_durations=(mo_end_time_num-icme_start_time_num)*24 #hours





sns.set_context("talk")     
#sns.set_style("darkgrid")  
sns.set_style("ticks",{'grid.linestyle': '--'})


yearly_start_times=[mdates.date2num(sunpy.time.parse_time('2007-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2008-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2009-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2010-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2011-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2012-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2013-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2014-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2015-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2016-01-01'))]
                  
yearly_end_times=[mdates.date2num(sunpy.time.parse_time('2007-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2008-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2009-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2010-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2011-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2012-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2013-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2014-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2015-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2016-12-31'))]

yearly_mid_times=[mdates.date2num(sunpy.time.parse_time('2007-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2008-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2009-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2010-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2011-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2012-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2013-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2014-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2015-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2016-07-01'))]



#***check further : for each year, calculate how much data is NaN for each spacecraft
#and use in calculation

#converted times are here:
[vex_time,wind_time,sta_time,stb_time,mav_time,mes_time]=pickle.load( open( "../catpy/DATACAT/insitu_times_mdates_maven_interp.p", "rb" ) )
sta= pickle.load( open( "../catpy/DATACAT/STA_2007to2015_SCEQ.p", "rb" ) )
stb= pickle.load( open( "../catpy/DATACAT/STB_2007to2014_SCEQ.p", "rb" ) )
wind=pickle.load( open( "../catpy/DATACAT/WIND_2007to2016_HEEQ.p", "rb" ) )

total_data_days_sta=np.zeros(np.size(yearly_mid_times))
total_data_days_sta.fill(np.nan)


total_data_days_stb=np.zeros(np.size(yearly_mid_times))
total_data_days_stb.fill(np.nan)


total_data_days_wind=np.zeros(np.size(yearly_mid_times))
total_data_days_wind.fill(np.nan)



#go through each year and search for data gaps, ok for solar wind missions 

for i in range(np.size(yearly_mid_times)):
 
  #Wind
  thisyear=np.where(np.logical_and((wind_time > yearly_start_times[i]),(wind_time < yearly_end_times[i])))
  nan=np.isnan(wind.btot[thisyear]) 
  notnan=np.where(nan == False)
  if np.size(notnan) >0: total_data_days_wind[i]=365
  if np.size(nan) > 0: total_data_days_wind[i]=np.size(notnan)/np.size(nan)*365
  
  #STA
  thisyear=np.where(np.logical_and((sta_time > yearly_start_times[i]),(sta_time < yearly_end_times[i])))
  nan=np.isnan(sta.btot[thisyear]) 
  notnan=np.where(nan == False)
  if np.size(notnan) >0: total_data_days_sta[i]=365
  if np.size(nan) > 0: total_data_days_sta[i]=np.size(notnan)/np.size(nan)*365

  #STB
  thisyear=np.where(np.logical_and((stb_time > yearly_start_times[i]),(stb_time < yearly_end_times[i])))
  nan=np.isnan(stb.btot[thisyear]) 
  notnan=np.where(nan == False)
  if np.size(notnan) >0: total_data_days_stb[i]=365
  if np.size(nan) > 0: total_data_days_stb[i]=np.size(notnan)/np.size(nan)*365

#for MESSENGER; VEX, MAVEN this is not correct because the nans during orbits are counted; 
#thats why we search manually for longer data gaps, and manually set the total_data_days_vex for each year




##################longer data gaps for MESSENGER and Mercury
total_data_days_mes=np.zeros(np.size(yearly_mid_times))
total_data_days_mes.fill(np.nan)

total_data_days_merc=np.zeros(np.size(yearly_mid_times))
total_data_days_merc.fill(np.nan)
#total_data_days_mes.fill(365)

#2007

jump1beg=mdates.date2num(sunpy.time.parse_time('2007-Jan-1'))
jump1end=mdates.date2num(sunpy.time.parse_time('2007-Mar-31'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2007-Jul-1'))
jump2end=mdates.date2num(sunpy.time.parse_time('2007-Jul-21'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2007-Aug-25'))
jump3end=mdates.date2num(sunpy.time.parse_time('2007-Dec-21'))

total_data_days_mes[0]=yearly_end_times[0]-yearly_start_times[0]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)

#2008
jump1beg=mdates.date2num(sunpy.time.parse_time('2008-Feb-24'))
jump1end=mdates.date2num(sunpy.time.parse_time('2008-Dec-31'))
#data gap too long - set nan
total_data_days_mes[1]=np.nan


#2009
jump1beg=mdates.date2num(sunpy.time.parse_time('2009-Jan-01'))
jump1end=mdates.date2num(sunpy.time.parse_time('2009-Jan-12'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2009-Jul-08'))
jump2end=mdates.date2num(sunpy.time.parse_time('2009-Jul-22'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2009-Aug-18'))
jump3end=mdates.date2num(sunpy.time.parse_time('2009-Aug-24'))

jump4beg=mdates.date2num(sunpy.time.parse_time('2009-Sep-01'))
jump4end=mdates.date2num(sunpy.time.parse_time('2009-Sep-04'))

jump5beg=mdates.date2num(sunpy.time.parse_time('2009-Sep-30'))
jump5end=mdates.date2num(sunpy.time.parse_time('2009-Oct-02'))

jump6beg=mdates.date2num(sunpy.time.parse_time('2009-Oct-13'))
jump6end=mdates.date2num(sunpy.time.parse_time('2009-Dec-10'))

total_data_days_mes[2]=yearly_end_times[2]-yearly_start_times[2]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)-(jump4end-jump4beg)-(jump5end-jump5beg)-(jump6end-jump6beg)

#2010
total_data_days_mes[3]=yearly_end_times[3]-yearly_start_times[3]

#2011
jump1beg=mdates.date2num(sunpy.time.parse_time('2011-Mar-16'))
jump1end=mdates.date2num(sunpy.time.parse_time('2011-Mar-24'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2011-Jan-01'))
jump2end=mdates.date2num(sunpy.time.parse_time('2011-Mar-24'))

total_data_days_merc[4]=yearly_end_times[4]-yearly_start_times[4]-(jump2end-jump2beg)
total_data_days_mes[4]=yearly_end_times[4]-yearly_start_times[4]-(jump1end-jump1beg)


#2012
jump1beg=mdates.date2num(sunpy.time.parse_time('2012-Jun-06'))
jump1end=mdates.date2num(sunpy.time.parse_time('2012-Jun-13'))

total_data_days_merc[5]=yearly_end_times[5]-yearly_start_times[5]-(jump1end-jump1beg)
total_data_days_mes[5]=yearly_end_times[5]-yearly_start_times[5]-(jump1end-jump1beg)

#2013
jump1beg=mdates.date2num(sunpy.time.parse_time('2013-Apr-01'))
jump1end=mdates.date2num(sunpy.time.parse_time('2013-May-01'))


total_data_days_merc[6]=yearly_end_times[6]-yearly_start_times[6]-(jump1end-jump1beg)
total_data_days_mes[6]=yearly_end_times[6]-yearly_start_times[6]-(jump1end-jump1beg)


#2014
jump1beg=mdates.date2num(sunpy.time.parse_time('2014-Feb-25'))
jump1end=mdates.date2num(sunpy.time.parse_time('2014-Mar-26'))

total_data_days_merc[7]=yearly_end_times[7]-yearly_start_times[7]-(jump1end-jump1beg)
total_data_days_mes[7]=yearly_end_times[7]-yearly_start_times[7]-(jump1end-jump1beg)

#2015
jump1beg=mdates.date2num(sunpy.time.parse_time('2015-Apr-30'))
jump1end=mdates.date2num(sunpy.time.parse_time('2015-Dec-31'))

total_data_days_merc[8]=yearly_end_times[8]-yearly_start_times[8]-(jump1end-jump1beg)

total_data_days_mes[8]=yearly_end_times[8]-yearly_start_times[8]-(jump1end-jump1beg)





############################ MAVEN

startdata=mdates.date2num(sunpy.time.parse_time('2014-Nov-27'))

#2015
jump1beg=mdates.date2num(sunpy.time.parse_time('2015-Mar-8'))
jump1end=mdates.date2num(sunpy.time.parse_time('2015-Jun-17'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2015-Oct-3'))
jump2end=mdates.date2num(sunpy.time.parse_time('2015-Dec-18'))

#2016
jump3beg=mdates.date2num(sunpy.time.parse_time('2016-Mar-27'))
jump3end=mdates.date2num(sunpy.time.parse_time('2016-Jun-3'))

jump4beg=mdates.date2num(sunpy.time.parse_time('2016-Sep-30'))
jump4end=mdates.date2num(sunpy.time.parse_time('2016-Dec-6'))

enddata=mdates.date2num(sunpy.time.parse_time('2016-Dec-31'))

total_data_days_mav=np.zeros(np.size(yearly_mid_times))
total_data_days_mav.fill(np.nan)

#this is 2014 - too few data points MAVEN (only 1 month)
#total_data_days_mav[7]=yearly_end_times[7]-startdata
total_data_days_mav[7]=np.nan


#this is 2015
total_data_days_mav[8]=yearly_end_times[8]-yearly_start_times[8]-(jump1end-jump1beg)-(jump2end-jump2beg)


#this is 2016
total_data_days_mav[9]=yearly_end_times[9]-yearly_start_times[9]-(jump3end-jump3beg)-(jump4end-jump4beg)






#################VEX 
total_data_days_vex=np.zeros(np.size(yearly_mid_times))
total_data_days_vex.fill(np.nan)


#times of longer data gaps

#2007
#from Jan 1 - Apr 1 there is not data gap

jump1beg=mdates.date2num(sunpy.time.parse_time('2007-Jul-5'))
jump1end=mdates.date2num(sunpy.time.parse_time('2007-Jul-12'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2007-Aug-23'))
jump2end=mdates.date2num(sunpy.time.parse_time('2007-Aug-28'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2007-Sep-18'))
jump3end=mdates.date2num(sunpy.time.parse_time('2007-Sep-20'))

total_data_days_vex[0]=yearly_end_times[0]-yearly_start_times[0]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)


#2008
jump1beg=mdates.date2num(sunpy.time.parse_time('2008-May-28'))
jump1end=mdates.date2num(sunpy.time.parse_time('2008-Jun-21'))

total_data_days_vex[1]=yearly_end_times[1]-yearly_start_times[1]-(jump1end-jump1beg)

#2009
jump1beg=mdates.date2num(sunpy.time.parse_time('2009-Apr-17'))
jump1end=mdates.date2num(sunpy.time.parse_time('2009-Apr-28'))
jump2beg=mdates.date2num(sunpy.time.parse_time('2009-Dec-28'))
jump2end=mdates.date2num(sunpy.time.parse_time('2009-Dec-31'))

total_data_days_vex[2]=yearly_end_times[2]-yearly_start_times[2]-(jump1end-jump1beg)-(jump2end-jump2beg)



#2010
jump1beg=mdates.date2num(sunpy.time.parse_time('2010-Jan-01'))
jump1end=mdates.date2num(sunpy.time.parse_time('2010-Jan-23'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2010-Apr-12'))
jump2end=mdates.date2num(sunpy.time.parse_time('2010-Apr-17'))

total_data_days_vex[3]=yearly_end_times[3]-yearly_start_times[3]-(jump1end-jump1beg)-(jump2end-jump2beg)


#2011
jump1beg=mdates.date2num(sunpy.time.parse_time('2011-Jan-24'))
jump1end=mdates.date2num(sunpy.time.parse_time('2011-Jan-27'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2011-Aug-05'))
jump2end=mdates.date2num(sunpy.time.parse_time('2011-Sep-01'))
total_data_days_vex[4]=yearly_end_times[4]-yearly_start_times[4]-(jump1end-jump1beg)-(jump2end-jump2beg)


#2012
jump1beg=mdates.date2num(sunpy.time.parse_time('2012-Mar-10'))
jump1end=mdates.date2num(sunpy.time.parse_time('2012-Mar-12'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2012-Jun-04'))
jump2end=mdates.date2num(sunpy.time.parse_time('2012-Jun-07'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2012-Jul-13'))
jump3end=mdates.date2num(sunpy.time.parse_time('2012-Jul-15'))

jump4beg=mdates.date2num(sunpy.time.parse_time('2012-Dec-29'))
jump4end=mdates.date2num(sunpy.time.parse_time('2012-Dec-31'))

total_data_days_vex[5]=yearly_end_times[5]-yearly_start_times[5]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg) -(jump4end-jump4beg)


#2013
jump1beg=mdates.date2num(sunpy.time.parse_time('2013-Mar-17'))
jump1end=mdates.date2num(sunpy.time.parse_time('2013-Apr-14'))
total_data_days_vex[6]=yearly_end_times[6]-yearly_start_times[6]-(jump1end-jump1beg)


#2014
jump1beg=mdates.date2num(sunpy.time.parse_time('2014-Feb-25'))
jump1end=mdates.date2num(sunpy.time.parse_time('2014-Mar-26'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2014-May-16'))
jump2end=mdates.date2num(sunpy.time.parse_time('2014-May-21'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2014-Jul-12'))
jump3end=mdates.date2num(sunpy.time.parse_time('2014-Jul-21'))

jump4beg=mdates.date2num(sunpy.time.parse_time('2014-Oct-13'))
jump4end=mdates.date2num(sunpy.time.parse_time('2014-Nov-11'))

jump5beg=mdates.date2num(sunpy.time.parse_time('2014-Nov-26'))
jump5end=mdates.date2num(sunpy.time.parse_time('2014-Dec-31'))

total_data_days_vex[7]=yearly_end_times[7]-yearly_start_times[7]-(jump1end-jump1beg)-(jump2end-jump2beg) -(jump3end-jump3beg) -(jump4end-jump4beg) -(jump5end-jump5beg)




#drop ulysses because too short, so 6 spacecraft

#make array with inside percentage

inside_wind_perc=np.zeros(np.size(yearly_mid_times))
inside_wind_perc.fill(np.nan)

inside_sta_perc=np.zeros(np.size(yearly_mid_times))
inside_sta_perc.fill(np.nan)

inside_stb_perc=np.zeros(np.size(yearly_mid_times))
inside_stb_perc.fill(np.nan)

inside_mes_perc=np.zeros(np.size(yearly_mid_times))
inside_mes_perc.fill(np.nan)

inside_merc_perc=np.zeros(np.size(yearly_mid_times))
inside_merc_perc.fill(np.nan)

inside_vex_perc=np.zeros(np.size(yearly_mid_times))
inside_vex_perc.fill(np.nan)

inside_mav_perc=np.zeros(np.size(yearly_mid_times))
inside_mav_perc.fill(np.nan)


#go through each year 
for i in range(np.size(yearly_mid_times)):
  
  #Wind:
  
  #select those icmes that are inside the current year
  thisyear=np.where(np.logical_and((icme_start_time_num[iwinind] > yearly_start_times[i]),(icme_start_time_num[iwinind] < yearly_end_times[i])))
  #summarize durations per year and convert to days
  total_icme_days=np.sum(icme_durations[thisyear])/24
  #get percentage
  if total_icme_days > 0:   inside_wind_perc[i]=total_icme_days/total_data_days_wind[i]*100

  thisyear=np.where(np.logical_and((icme_start_time_num[istaind] > yearly_start_times[i]),(icme_start_time_num[istaind] < yearly_end_times[i])))
  total_icme_days=np.sum(icme_durations[thisyear])/24
  if total_icme_days > 0:   inside_sta_perc[i]=total_icme_days/total_data_days_sta[i]*100
  
  thisyear=np.where(np.logical_and((icme_start_time_num[istbind] > yearly_start_times[i]),(icme_start_time_num[istbind] < yearly_end_times[i])))
  total_icme_days=np.sum(icme_durations[thisyear])/24
  if total_icme_days > 0:   inside_stb_perc[i]=total_icme_days/total_data_days_stb[i]*100

  thisyear=np.where(np.logical_and((icme_start_time_num[imesind] > yearly_start_times[i]),(icme_start_time_num[imesind] < yearly_end_times[i])))
  total_icme_days=np.sum(icme_durations[thisyear])/24
  if total_icme_days > 0:   inside_mes_perc[i]=total_icme_days/total_data_days_mes[i]*100
  
  thisyear=np.where(np.logical_and((icme_start_time_num[imercind] > yearly_start_times[i]),(icme_start_time_num[imercind] < yearly_end_times[i])))
  total_icme_days=np.sum(icme_durations[thisyear])/24
  if total_icme_days > 0:   inside_merc_perc[i]=total_icme_days/total_data_days_merc[i]*100

  
  thisyear=np.where(np.logical_and((icme_start_time_num[ivexind] > yearly_start_times[i]),(icme_start_time_num[ivexind] < yearly_end_times[i])))
  total_icme_days=np.sum(icme_durations[thisyear])/24
  if total_icme_days > 0:   inside_vex_perc[i]=total_icme_days/total_data_days_vex[i]*100

  thisyear=np.where(np.logical_and((icme_start_time_num[imavind] > yearly_start_times[i]),(icme_start_time_num[imavind] < yearly_end_times[i])))
  total_icme_days=np.sum(icme_durations[thisyear])/24
  if total_icme_days > 0:   inside_mav_perc[i]=total_icme_days/total_data_days_mav[i]*100













#######*********************** check next section again; if min/rise/max shifts, manually adjust data gaps

#####################################################################################
############ make the same thing not yearly, but for the 3 solar cycle phases


cycle_start_times=[minstart, risestart, maxstart]
cycle_end_times=[minend, riseend, maxend]


total_data_days_sta_cycle=np.zeros(np.size(cycle_start_times))
total_data_days_sta_cycle.fill(np.nan)

total_data_days_stb_cycle=np.zeros(np.size(cycle_start_times))
total_data_days_stb_cycle.fill(np.nan)

total_data_days_wind_cycle=np.zeros(np.size(cycle_start_times))
total_data_days_wind_cycle.fill(np.nan)


#define manually the data time ranges inside min, rise, max

#go through each year and search for data gaps, ok for solar wind missions 

for i in range(np.size(cycle_start_times)):
 
  #Wind
  phase=np.where(np.logical_and((wind_time > cycle_start_times[i]),(wind_time < cycle_end_times[i])))
  nan=np.isnan(wind.btot[phase]) 
  notnan=np.where(nan == False)
  if np.size(notnan) >0: total_data_days_wind_cycle[i]=cycle_end_times[i]-cycle_start_times[i]
  if np.size(nan) > 0: total_data_days_wind_cycle[i]=np.size(notnan)/np.size(nan)*(cycle_end_times[i]-cycle_start_times[i])
  
  #STA
  phase=np.where(np.logical_and((sta_time > cycle_start_times[i]),(sta_time < cycle_end_times[i])))
  nan=np.isnan(sta.btot[phase]) 
  notnan=np.where(nan == False)
  if np.size(notnan) >0: total_data_days_sta_cycle[i]=cycle_end_times[i]-cycle_start_times[i]
  if np.size(nan) > 0: total_data_days_sta_cycle[i]=np.size(notnan)/np.size(nan)*(cycle_end_times[i]-cycle_start_times[i])

  #STB
  phase=np.where(np.logical_and((stb_time > cycle_start_times[i]),(stb_time < cycle_end_times[i])))
  nan=np.isnan(stb.btot[phase]) 
  notnan=np.where(nan == False)
  if np.size(notnan) >0: total_data_days_stb_cycle[i]=cycle_end_times[i]-cycle_start_times[i]
  if np.size(nan) > 0: total_data_days_stb_cycle[i]=np.size(notnan)/np.size(nan)*(cycle_end_times[i]-cycle_start_times[i])

#for MESSENGER; VEX, MAVEN this is not correct because the nans during orbits are counted; 
#thats why we search manually for longer data gaps, and manually set the total_data_days_vex_cycle for each phase


##################longer data gaps for MESSENGER and Mercury
total_data_days_mes_cycle=np.zeros(np.size(cycle_start_times))
total_data_days_mes_cycle.fill(np.nan)

total_data_days_merc_cycle=np.zeros(np.size(cycle_start_times))
total_data_days_merc_cycle.fill(np.nan)
#total_data_days_mes.fill(365)

#min

jump1beg=mdates.date2num(sunpy.time.parse_time('2007-Jan-1'))
jump1end=mdates.date2num(sunpy.time.parse_time('2007-Mar-31'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2007-Jul-1'))
jump2end=mdates.date2num(sunpy.time.parse_time('2007-Jul-21'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2007-Aug-25'))
jump3end=mdates.date2num(sunpy.time.parse_time('2007-Dec-21'))

jump4beg=mdates.date2num(sunpy.time.parse_time('2008-Feb-24'))
jump4end=mdates.date2num(sunpy.time.parse_time('2008-Dec-31'))

jump5beg=mdates.date2num(sunpy.time.parse_time('2009-Jan-01'))
jump5end=mdates.date2num(sunpy.time.parse_time('2009-Jan-12'))

jump6beg=mdates.date2num(sunpy.time.parse_time('2009-Jul-08'))
jump6end=mdates.date2num(sunpy.time.parse_time('2009-Jul-22'))

jump7beg=mdates.date2num(sunpy.time.parse_time('2009-Aug-18'))
jump7end=mdates.date2num(sunpy.time.parse_time('2009-Aug-24'))

jump8beg=mdates.date2num(sunpy.time.parse_time('2009-Sep-01'))
jump8end=mdates.date2num(sunpy.time.parse_time('2009-Sep-04'))

jump9beg=mdates.date2num(sunpy.time.parse_time('2009-Sep-30'))
jump9end=mdates.date2num(sunpy.time.parse_time('2009-Oct-02'))

jump10beg=mdates.date2num(sunpy.time.parse_time('2009-Oct-13'))
jump10end=mdates.date2num(sunpy.time.parse_time('2009-Dec-10'))

total_data_days_mes_cycle[0]=cycle_end_times[0]-cycle_start_times[0]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)-(jump4end-jump4beg)-(jump5end-jump5beg)-(jump6end-jump6beg)-(jump7end-jump7beg)-(jump8end-jump8beg)-(jump9end-jump9beg)-(jump10end-jump10beg)



#rise 2010 anfang bis mitte 2011

jump1beg=mdates.date2num(sunpy.time.parse_time('2011-Mar-16'))
jump1end=mdates.date2num(sunpy.time.parse_time('2011-Mar-24'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2011-Jan-01'))
jump2end=mdates.date2num(sunpy.time.parse_time('2011-Mar-24'))

total_data_days_mes_cycle[1]=cycle_end_times[1]-cycle_start_times[1]-(jump1end-jump1beg)-(jump2end-jump2beg)
##*******correct?
total_data_days_merc_cycle[1]=cycle_end_times[1]-mdates.date2num(sunpy.time.parse_time('2011-Mar-24'))

#max 
#2011 mitte bis ende 2015
jump1beg=mdates.date2num(sunpy.time.parse_time('2012-Jun-06'))
jump1end=mdates.date2num(sunpy.time.parse_time('2012-Jun-13'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2013-Apr-01'))
jump2end=mdates.date2num(sunpy.time.parse_time('2013-May-01'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2014-Feb-25'))
jump3end=mdates.date2num(sunpy.time.parse_time('2014-Mar-26'))

total_data_days_mes_cycle[2]=cycle_end_times[2]-cycle_start_times[2]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)
total_data_days_merc_cycle[2]=cycle_end_times[2]-cycle_start_times[2]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)



############################ MAVEN

total_data_days_mav_cycle=np.zeros(np.size(cycle_start_times))
total_data_days_mav_cycle.fill(np.nan)

#only declining phase

jump1beg=mdates.date2num(sunpy.time.parse_time('2015-Mar-8'))
jump1end=mdates.date2num(sunpy.time.parse_time('2015-Jun-17'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2015-Oct-3'))
jump2end=mdates.date2num(sunpy.time.parse_time('2015-Dec-18'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2016-Mar-27'))
jump3end=mdates.date2num(sunpy.time.parse_time('2016-Jun-3'))

jump4beg=mdates.date2num(sunpy.time.parse_time('2016-Sep-30'))
jump4end=mdates.date2num(sunpy.time.parse_time('2016-Dec-6'))

startdata=mdates.date2num(sunpy.time.parse_time('2014-Nov-27'))
enddata=mdates.date2num(sunpy.time.parse_time('2016-Dec-31'))

total_data_days_mav_cycle[1]=enddata-startdata-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)-(jump4end-jump4beg)


#################VEX 
total_data_days_vex_cycle=np.zeros(np.size(cycle_start_times))
total_data_days_vex_cycle.fill(np.nan)

#times of longer data gaps

#min
#2007 #from Jan 1 - Apr 1 there is not data gap

jump1beg=mdates.date2num(sunpy.time.parse_time('2007-Jul-5'))
jump1end=mdates.date2num(sunpy.time.parse_time('2007-Jul-12'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2007-Aug-23'))
jump2end=mdates.date2num(sunpy.time.parse_time('2007-Aug-28'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2007-Sep-18'))
jump3end=mdates.date2num(sunpy.time.parse_time('2007-Sep-20'))

jump4beg=mdates.date2num(sunpy.time.parse_time('2008-May-28'))
jump4end=mdates.date2num(sunpy.time.parse_time('2008-Jun-21'))

jump5beg=mdates.date2num(sunpy.time.parse_time('2009-Apr-17'))
jump5end=mdates.date2num(sunpy.time.parse_time('2009-Apr-28'))

jump6beg=mdates.date2num(sunpy.time.parse_time('2009-Dec-28'))
jump6end=mdates.date2num(sunpy.time.parse_time('2009-Dec-31'))

total_data_days_vex_cycle[0]=cycle_end_times[0]-cycle_start_times[0]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)-(jump4end-jump4beg)-(jump5end-jump5beg)-(jump6end-jump6beg)

#rise
#2010 bis mitte 2011
jump1beg=mdates.date2num(sunpy.time.parse_time('2010-Jan-01'))
jump1end=mdates.date2num(sunpy.time.parse_time('2010-Jan-23'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2010-Apr-12'))
jump2end=mdates.date2num(sunpy.time.parse_time('2010-Apr-17'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2011-Jan-24'))
jump3end=mdates.date2num(sunpy.time.parse_time('2011-Jan-27'))

total_data_days_vex_cycle[1]=cycle_end_times[1]-cycle_start_times[1]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)


#max 2011 mitte bis ende 2015 (end of mission VEX)
jump1beg=mdates.date2num(sunpy.time.parse_time('2011-Aug-05'))
jump1end=mdates.date2num(sunpy.time.parse_time('2011-Sep-01'))

jump2beg=mdates.date2num(sunpy.time.parse_time('2012-Mar-10'))
jump2end=mdates.date2num(sunpy.time.parse_time('2012-Mar-12'))

jump3beg=mdates.date2num(sunpy.time.parse_time('2012-Jun-04'))
jump3end=mdates.date2num(sunpy.time.parse_time('2012-Jun-07'))

jump4beg=mdates.date2num(sunpy.time.parse_time('2012-Jul-13'))
jump4end=mdates.date2num(sunpy.time.parse_time('2012-Jul-15'))

jump5beg=mdates.date2num(sunpy.time.parse_time('2012-Dec-29'))
jump5end=mdates.date2num(sunpy.time.parse_time('2012-Dec-31'))

jump6beg=mdates.date2num(sunpy.time.parse_time('2013-Mar-17'))
jump6end=mdates.date2num(sunpy.time.parse_time('2013-Apr-14'))

jump7beg=mdates.date2num(sunpy.time.parse_time('2014-Feb-25'))
jump7end=mdates.date2num(sunpy.time.parse_time('2014-Mar-26'))

jump8beg=mdates.date2num(sunpy.time.parse_time('2014-May-16'))
jump8end=mdates.date2num(sunpy.time.parse_time('2014-May-21'))

jump9beg=mdates.date2num(sunpy.time.parse_time('2014-Jul-12'))
jump9end=mdates.date2num(sunpy.time.parse_time('2014-Jul-21'))

jump10beg=mdates.date2num(sunpy.time.parse_time('2014-Oct-13'))
jump10end=mdates.date2num(sunpy.time.parse_time('2014-Nov-11'))

jump11beg=mdates.date2num(sunpy.time.parse_time('2014-Nov-26'))
jump11end=mdates.date2num(sunpy.time.parse_time('2014-Dec-31'))

total_data_days_vex_cycle[2]=cycle_end_times[2]-cycle_start_times[2]-(jump1end-jump1beg)-(jump2end-jump2beg)-(jump3end-jump3beg)-(jump4end-jump4beg)-(jump5end-jump5beg)-(jump6end-jump6beg)-(jump7end-jump7beg)-(jump8end-jump8beg)-(jump9end-jump9beg)-(jump10end-jump10beg)-(jump11end-jump11beg)









############################################


inside_wind_perc_cycle=np.zeros(np.size(cycle_start_times))
inside_wind_perc_cycle.fill(np.nan)

inside_sta_perc_cycle=np.zeros(np.size(cycle_start_times))
inside_sta_perc_cycle.fill(np.nan)

inside_stb_perc_cycle=np.zeros(np.size(cycle_start_times))
inside_stb_perc_cycle.fill(np.nan)

inside_mes_perc_cycle=np.zeros(np.size(cycle_start_times))
inside_mes_perc_cycle.fill(np.nan)

inside_merc_perc_cycle=np.zeros(np.size(cycle_start_times))
inside_merc_perc_cycle.fill(np.nan)

inside_vex_perc_cycle=np.zeros(np.size(cycle_start_times))
inside_vex_perc_cycle.fill(np.nan)

inside_mav_perc_cycle=np.zeros(np.size(cycle_start_times))
inside_mav_perc_cycle.fill(np.nan)

#maven only declining phase
inside_mav_perc_cycle[1]=(np.sum(icme_durations[imavind]))/24/total_data_days_mav_cycle[1]*100

#go through solar cycle phases min, rise, max for Wind, VEX, MES, STA, STB
for i in range(np.size(cycle_start_times)):
  
  #Wind:
  #select those icmes that are inside min, rise, max
  phase=np.where(np.logical_and((icme_start_time_num[iwinind] > cycle_start_times[i]),(icme_start_time_num[iwinind] < cycle_end_times[i])))
  #summarize durations per phase and convert to days
  total_icme_days=np.sum(icme_durations[phase])/24
  #get percentage
  if total_icme_days > 0:   inside_wind_perc_cycle[i]=total_icme_days/total_data_days_wind_cycle[i]*100

  phase=np.where(np.logical_and((icme_start_time_num[istaind] > cycle_start_times[i]),(icme_start_time_num[istaind] < cycle_end_times[i])))
  total_icme_days=np.sum(icme_durations[phase])/24
  if total_icme_days > 0:   inside_sta_perc_cycle[i]=total_icme_days/total_data_days_sta_cycle[i]*100
  
  phase=np.where(np.logical_and((icme_start_time_num[istbind] > cycle_start_times[i]),(icme_start_time_num[istbind] < cycle_end_times[i])))
  total_icme_days=np.sum(icme_durations[phase])/24
  if total_icme_days > 0:   inside_stb_perc_cycle[i]=total_icme_days/total_data_days_stb_cycle[i]*100

  phase=np.where(np.logical_and((icme_start_time_num[imesind] > cycle_start_times[i]),(icme_start_time_num[imesind] < cycle_end_times[i])))
  total_icme_days=np.sum(icme_durations[phase])/24
  if total_icme_days > 0:   inside_mes_perc_cycle[i]=total_icme_days/total_data_days_mes_cycle[i]*100
  
  phase=np.where(np.logical_and((icme_start_time_num[imercind] > cycle_start_times[i]),(icme_start_time_num[imercind] < cycle_end_times[i])))
  total_icme_days=np.sum(icme_durations[phase])/24
  if total_icme_days > 0:   inside_merc_perc_cycle[i]=total_icme_days/total_data_days_merc_cycle[i]*100
  
  phase=np.where(np.logical_and((icme_start_time_num[ivexind] > cycle_start_times[i]),(icme_start_time_num[ivexind] < cycle_end_times[i])))
  total_icme_days=np.sum(icme_durations[phase])/24
  if total_icme_days > 0:   inside_vex_perc_cycle[i]=total_icme_days/total_data_days_vex_cycle[i]*100

  #phase=np.where(np.logical_and((icme_start_time_num[imavind] > cycle_start_times[i]),(icme_start_time_num[imavind] < cycle_end_times[i])))
  #total_icme_days=np.sum(icme_durations[phase])/24
  #if total_icme_days > 0:   inside_mav_perc_cycle[i]=total_icme_days/total_data_days_mav_cycle[i]*100





















  
  

### fix that VEX MESSENGER impact frequency is less than 1 AU by multiplying with a factor of 1.5
#check exact values with frequency plot

#inside_vex_perc=inside_vex_perc*1.5
#inside_mes_perc=inside_mes_perc*1.5


fig=plt.figure(5,figsize=(10,10	))

ax1 = plt.subplot(211) 

plt.plot_date(yearly_mid_times,inside_wind_perc,'o',color='mediumseagreen',markersize=8, linestyle='-')
plt.plot_date(yearly_mid_times,inside_merc_perc,'o',color='darkgrey',markersize=8,linestyle='-')
plt.plot_date(yearly_mid_times,inside_vex_perc,'o',color='orange',markersize=8,linestyle='-')
plt.plot_date(yearly_mid_times,inside_stb_perc,'o',color='royalblue',markersize=8,linestyle='-')
plt.plot_date(yearly_mid_times,inside_sta_perc,'o',color='red',markersize=8,linestyle='-')
plt.plot_date(yearly_mid_times,inside_mav_perc,'o',color='steelblue',markersize=8,linestyle='-')

plt.ylabel('Time inside ICME [%]')

#plt.xlim(yearly_bin_edges[0],yearly_bin_edges[10])
ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)

#sets planet / spacecraft labels
xoff=0.15
yoff=0.85
fsize=14

plt.figtext(xoff,yoff-0.03*0,'Mercury',color='darkgrey', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*1,'Venus',color='orange', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*2,'Earth',color='mediumseagreen', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*3,'Mars',color='steelblue', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*4,'STEREO-A',color='red', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*5,'STEREO-B',color='royalblue', fontsize=fsize, ha='left')
#panel labels
plt.figtext(0.02,0.98,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.02,0.48,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')



#limits solar min/rise/max

vlevel=22
fsize=11

plt.axvspan(minstart,minend, color='green', alpha=0.1)
plt.annotate('solar minimum',xy=(minstart+(minend-minstart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(minstart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(minend-10,vlevel),ha='right', fontsize=fsize)


plt.axvspan(risestart,riseend, color='yellow', alpha=0.1)
plt.annotate('rising phase',xy=(risestart+(riseend-risestart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(risestart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(riseend-10,vlevel),ha='right', fontsize=fsize)

plt.axvspan(maxstart,maxend, color='red', alpha=0.1)
plt.annotate('solar maximum',xy=(maxstart+(maxend-maxstart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(maxstart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(maxend,vlevel),ha='right', fontsize=fsize)


plt.ylim((0,25))
fsize=15
plt.ylabel('Time inside ICME [%]')
plt.xlabel('Year',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.tight_layout()










#plt.ylim(0,45)
plt.xlim(yearly_start_times[0],yearly_end_times[9])

#sns.despine()




#### plot time inside vs. heliocentric distance

pos_wind_perc=np.zeros(np.size(yearly_mid_times))
pos_wind_perc.fill(np.nan)
pos_wind_perc_std=np.zeros(np.size(yearly_mid_times))
pos_wind_perc_std.fill(np.nan)

pos_sta_perc=np.zeros(np.size(yearly_mid_times))
pos_sta_perc.fill(np.nan)
pos_sta_perc_std=np.zeros(np.size(yearly_mid_times))
pos_sta_perc_std.fill(np.nan)

pos_stb_perc=np.zeros(np.size(yearly_mid_times))
pos_stb_perc.fill(np.nan)
pos_stb_perc_std=np.zeros(np.size(yearly_mid_times))
pos_stb_perc_std.fill(np.nan)

#pos_mes_perc=np.zeros(np.size(yearly_mid_times))
#pos_mes_perc.fill(np.nan)
#pos_mes_perc_std=np.zeros(np.size(yearly_mid_times))
#pos_mes_perc_std.fill(np.nan)


pos_merc_perc=np.zeros(np.size(yearly_mid_times))
pos_merc_perc.fill(np.nan)
pos_merc_perc_std=np.zeros(np.size(yearly_mid_times))
pos_merc_perc_std.fill(np.nan)

pos_vex_perc=np.zeros(np.size(yearly_mid_times))
pos_vex_perc.fill(np.nan)
pos_vex_perc_std=np.zeros(np.size(yearly_mid_times))
pos_vex_perc_std.fill(np.nan)

pos_mav_perc=np.zeros(np.size(yearly_mid_times))
pos_mav_perc.fill(np.nan)
pos_mav_perc_std=np.zeros(np.size(yearly_mid_times))
pos_mav_perc_std.fill(np.nan)


allpositions=np.zeros([np.size(yearly_mid_times), 6])
allinside=np.zeros([np.size(yearly_mid_times), 6])

#calculate average distance +/- std for each year
#go through each year 
for i in range(np.size(yearly_mid_times)):
  
  #select those positions that are inside the current year
  thisyear=np.where(np.logical_and((pos_time_num > yearly_start_times[i]),(pos_time_num < yearly_end_times[i])))
  
  #pos_mes_perc[i]=np.mean(pos.messenger[0][thisyear])
  #pos_mes_perc_std[i]=np.std(pos.messenger[0][thisyear])
  pos_merc_perc[i]=np.mean(pos.mercury[0][thisyear])
  pos_merc_perc_std[i]=np.std(pos.mercury[0][thisyear])
  

  pos_mav_perc[i]=np.mean(pos.mars[0][thisyear])
  pos_mav_perc_std[i]=np.std(pos.mars[0][thisyear])

  pos_vex_perc[i]=np.mean(pos.venus[0][thisyear])
  pos_vex_perc_std[i]=np.std(pos.venus[0][thisyear])

  pos_wind_perc[i]=np.mean(pos.earth_l1[0][thisyear])
  pos_wind_perc_std[i]=np.std(pos.earth_l1[0][thisyear])

  pos_sta_perc[i]=np.mean(pos.sta[0][thisyear])
  pos_sta_perc_std[i]=np.std(pos.sta[0][thisyear])

  pos_stb_perc[i]=np.mean(pos.stb[0][thisyear])
  pos_stb_perc_std[i]=np.std(pos.stb[0][thisyear])
  
  allpositions[i][:]=(pos_merc_perc[i], pos_mav_perc[i], pos_vex_perc[i],pos_wind_perc[i],pos_sta_perc[i],pos_stb_perc[i])
  allinside[i][:]=(inside_merc_perc[i], inside_mav_perc[i], inside_vex_perc[i],inside_wind_perc[i],inside_sta_perc[i],inside_stb_perc[i])
  
 
  



#***make alpha variable for each year?

ax3 = plt.subplot(212) 


#for every year linear fit **check if power law works better

#for fit plotting
xfit=np.linspace(0,2,1000)

#allpositions[i] and allinside[i] are the data for each year
#no fit for 2016 as only MAVEN data is available


for i in range(np.size(yearly_mid_times)-2):
 #make linear fits ignoring NaN
 notnan=np.where(np.isfinite(allinside[i]) > 0)
 durfit=np.polyfit(allpositions[i][notnan],allinside[i][notnan],1)
 #this is similar to D=durfit[0]*xfit+durfit[1]
 durfitfun=np.poly1d(durfit)
 print('year',i+2007)
 print('time inside linear fit: D[hours]={:.2f}r[AU]+{:.2f}'.format(durfit[0],durfit[1]))
 plt.plot(xfit,durfitfun(xfit),'-',color='black', lw=2, alpha=i/10+0.2)#,label='fit')
 
 plt.errorbar(pos_merc_perc[i], inside_merc_perc[i], xerr=pos_merc_perc_std[i],yerr=0,fmt='o',color='darkgrey',markersize=8,linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_mav_perc[i], inside_mav_perc[i],xerr=pos_mav_perc_std[i],fmt='o',color='steelblue',markersize=8,linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_sta_perc[i], inside_sta_perc[i],xerr=pos_sta_perc_std[i],fmt='o',color='red',markersize=8,linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_stb_perc[i], inside_stb_perc[i],xerr=pos_stb_perc_std[i],fmt='o',color='royalblue',markersize=8,linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_wind_perc[i], inside_wind_perc[i],xerr=pos_wind_perc_std[i],fmt='o',color='mediumseagreen',markersize=8, linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_vex_perc[i], inside_vex_perc[i],xerr=pos_vex_perc_std[i],fmt='o',color='orange',markersize=8,linestyle=' ',alpha=i/10+0.2)
 
 plt.annotate(str(i+2007), xy=(0.1,5+2.5*i), alpha=i/10+0.2)
 
 
 #reconstruct Mars time inside from linear fits but not for years 2015 /2016
 if i < 8: inside_mav_perc[i]=durfitfun(pos_mav_perc[i])


#mars limits
plt.axvspan(np.min(pos.mars[0]),np.max(pos.mars[0]), color='orangered', alpha=0.2)
#plt.figtext(0.8,0.8,'Mars',color='orangered')
plt.axvspan(np.min(pos.mercury[0]),np.max(pos.mercury[0]), color='darkgrey', alpha=0.2)
#plt.figtext(0.25,0.8,'Mercury',color='darkgrey')
plt.axvspan(np.min(pos.venus[0]),np.max(pos.venus[0]), color='orange', alpha=0.2)
#plt.figtext(0.42,0.8,'Venus',color='orange')
plt.axvspan(np.min(pos.earth[0]),np.max(pos.earth[0]), color='mediumseagreen', alpha=0.2)
#plt.figtext(0.6,0.8,'Earth',color='mediumseagreen')
plt.xlim(0,1.8)

#solar probe plus 10 to 36 Rs close approaches

#plt.axvspan(Rs_in_AU*10,Rs_in_AU*36,color='magenta', alpha=0.2)

plt.ylabel('Time inside ICME [%]')
plt.xlabel('Heliocentric distance [AU]')
ax3.set_xticks(np.arange(0,2,0.2))

#add reconstructed Mars time inside on previous plot
ax1.plot_date(yearly_mid_times,inside_mav_perc,'o',color='steelblue',markersize=8,linestyle='--')


plt.ylim((0,25))

plt.tight_layout()

plt.show()
plt.savefig('plots_psp/inside.pdf', dpi=300)
plt.savefig('plots_psp/inside.png', dpi=300)
































###################################################################################

##################### (4) arrival frequencies in ICMECAT  ##############


yearly_bin_edges=[mdates.date2num(sunpy.time.parse_time('2007-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2008-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2009-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2010-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2011-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2012-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2013-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2014-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2015-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2016-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2017-01-01'))]

#bin width in days         
binweite=360/8


sns.set_context("talk")     
#sns.set_style("darkgrid")  
sns.set_style("ticks",{'grid.linestyle': '--'})

fig=plt.figure(4,figsize=(12,10	))


fsize=15

ax1 = plt.subplot(211) 

plt.plot_date(icme_start_time_num[iwinind],sc_heliodistance[iwinind],fmt='o',color='mediumseagreen',markersize=5)
plt.plot_date(icme_start_time_num[imesind],sc_heliodistance[imesind],fmt='o',color='darkgrey',markersize=5)
plt.plot_date(icme_start_time_num[ivexind],sc_heliodistance[ivexind],fmt='o',color='orange',markersize=5)
plt.plot_date(icme_start_time_num[istbind],sc_heliodistance[istbind],fmt='o',color='royalblue',markersize=5)
plt.plot_date(icme_start_time_num[istaind],sc_heliodistance[istaind],fmt='o',color='red',markersize=5)
plt.plot_date(icme_start_time_num[iulyind],sc_heliodistance[iulyind],fmt='o',color='brown',markersize=5)
plt.plot_date(icme_start_time_num[imavind],sc_heliodistance[imavind],fmt='o',color='steelblue',markersize=5)




fsize=15
plt.ylabel('Heliocentric distance R [AU]',fontsize=fsize)
plt.xlabel('Year',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.xlim(yearly_bin_edges[0],yearly_bin_edges[10])
ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)



##############

ax2 = plt.subplot(212) 

(histwin, bin_edgeswin) = np.histogram(icme_start_time_num[iwinind], yearly_bin_edges)
(histvex, bin_edgesvex) = np.histogram(icme_start_time_num[ivexind], yearly_bin_edges)
(histmes, bin_edgesmes) = np.histogram(icme_start_time_num[imesind], yearly_bin_edges)
(histstb, bin_edgesstb) = np.histogram(icme_start_time_num[istbind], yearly_bin_edges)
(histsta, bin_edgessta) = np.histogram(icme_start_time_num[istaind], yearly_bin_edges)
(histmav, bin_edgesmav) = np.histogram(icme_start_time_num[imavind], yearly_bin_edges)

#********
#recalculate number of ICMEs as events per month or day, including data gaps


cycle_bin_edges=[minstart, minend, riseend, maxend]

(histwincyc, bin_edgescyc) = np.histogram(icme_start_time_num[iwinind], cycle_bin_edges)
(histvexcyc, bin_edgescyc) = np.histogram(icme_start_time_num[ivexind], cycle_bin_edges)
(histmescyc, bin_edgescyc) = np.histogram(icme_start_time_num[imesind], cycle_bin_edges)
(histstbcyc, bin_edgescyc) = np.histogram(icme_start_time_num[istbind], cycle_bin_edges)
(histstacyc, bin_edgescyc) = np.histogram(icme_start_time_num[istaind], cycle_bin_edges)
(histmavcyc, bin_edgescyc) = np.histogram(icme_start_time_num[imavind], cycle_bin_edges)

#use total_data_days_vex etc. from previous plot 
histwincyc=histwincyc/total_data_days_wind_cycle*365
histvexcyc=histvexcyc/total_data_days_vex_cycle*365
histmescyc=histmescyc/total_data_days_mes_cycle*365
histstbcyc=histstbcyc/total_data_days_stb_cycle*365
histstacyc=histstacyc/total_data_days_sta_cycle*365
histmavcyc=histmavcyc/total_data_days_mav_cycle*365


#normalize each dataset for data gaps

histwin=histwin/total_data_days_wind*365
histvex=histvex/total_data_days_vex*365
histmes=histmes/total_data_days_mes*365
histsta=histsta/total_data_days_sta*365
histstb=histstb/total_data_days_stb*365
histmav=histmav/total_data_days_mav*365

binedges=bin_edgeswin
pickle.dump([binedges,histwin,histvex,histmes,histsta,histstb,histmav], open( "plots_psp/icme_frequency.p", "wb" ), protocol=2 )
#[binedges,histwin,histvex,histmes,histsta,histstb,histmav]=pickle.load( open( "plots/stats_psp/icme_frequency.p", "rb" ) )

#binweite=45
ax2.bar(bin_edgeswin[:-1]+30,histwin, width=binweite,color='mediumseagreen', alpha=0.5,edgecolor='black')
ax2.bar(bin_edgesvex[:-1]+30+binweite,histvex, width=binweite,color='orange', alpha=0.5,edgecolor='black')
ax2.bar(bin_edgesmes[:-1]+30+ binweite*2,histmes, width=binweite,color='darkgrey', alpha=0.5,edgecolor='black')
ax2.bar(bin_edgesstb[:-1]+30+binweite*3,histstb, width=binweite,color='royalblue', alpha=0.5,edgecolor='black')
ax2.bar(bin_edgessta[:-1]+30+binweite*4,histsta, width=binweite,color='red', alpha=0.5,edgecolor='black')
#ax2.bar(bin_edgessta[:-1]+30+binweite*5,histuly, width=binweite,color='brown', alpha=0.5)
ax2.bar(bin_edgesmav[:-1]+30+binweite*6,histmav, width=binweite,color='steelblue', alpha=0.5,edgecolor='black')

plt.xlim(yearly_bin_edges[0],yearly_bin_edges[10])
ax2.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax2.xaxis.set_major_formatter(myformat)
#sets planet / spacecraft labels
xoff=0.85
yoff=0.45
fsize=12
plt.figtext(xoff,yoff,'Earth L1',color='mediumseagreen', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*1,'VEX',color='orange', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*2,'MESSENGER',color='dimgrey', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*3,'STEREO-A',color='red', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*4,'STEREO-B',color='royalblue', fontsize=fsize, ha='left')
#plt.figtext(xoff,yoff-0.03*5,'Ulysses',color='brown', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*5,'MAVEN',color='steelblue', fontsize=fsize, ha='left')
#panel labels
plt.figtext(0.02,0.98,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.02,0.48,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')

plt.ylim(0,48)


#limits solar min/rise/max

vlevel=44
fsize=13

plt.axvspan(minstart,minend, color='green', alpha=0.1)
plt.annotate('solar minimum',xy=(minstart+(minend-minstart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(minstart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(minend-10,vlevel),ha='right', fontsize=fsize)


plt.axvspan(risestart,riseend, color='yellow', alpha=0.1)
plt.annotate('rising phase',xy=(risestart+(riseend-risestart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(risestart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(riseend-10,vlevel),ha='right', fontsize=fsize)

plt.axvspan(maxstart,maxend, color='red', alpha=0.1)
plt.annotate('solar maximum',xy=(maxstart+(maxend-maxstart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(maxstart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(maxend,vlevel),ha='right', fontsize=fsize)


fsize=15
plt.ylabel('ICMEs per year',fontsize=fsize)
plt.xlabel('Year',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.tight_layout()

#sns.despine()
plt.show()
plt.savefig('plots_psp/frequency.pdf', dpi=300)
plt.savefig('plots_psp/frequency.png', dpi=300)

print('for solar min 2007-2009 average ICME per year rate:')
mean07=np.mean([histwin[0],histvex[0],histsta[0],histstb[0],histmes[0]])
mean08=np.mean([histwin[1],histvex[1],histsta[1],histstb[1],histmes[1]])
mean09=np.mean([histwin[2],histvex[2],histsta[2],histstb[2],histmes[2]])
print(np.nanmean([mean07,mean08,mean09]))

print('for 2010 2011')
mean10=np.mean([histwin[3],histvex[3],histsta[3],histstb[3],histmes[3]])
mean11=np.mean([histwin[4],histvex[4],histsta[4],histstb[4],histmes[4]])
print(np.mean([mean10,mean11]))


print('for 2012 2013 2014')
mean12=np.mean([histwin[5],histvex[5],histsta[5],histstb[5],histmes[5]])
mean13=np.mean([histwin[6],histvex[6],histsta[6],histstb[6],histmes[6]])
mean14=np.mean([histwin[7],histvex[7],histsta[7],histstb[7],histmes[7]])

print(np.mean([mean12,mean13,mean14]))







sys.exit()













































###################################### ALL RESULTS


##########

#B 

print('--------------------------------------------------')
print('Magnetic field B MO_BMEAN')

print()
print('Mercury +/-')
print(round(np.mean(mo_bmean[imercind]),1))
print(round(np.std(mo_bmean[imercind]),1))
#print('min')
#np.mean(mo_bmean[imercind][imercind_min])
#np.std(mo_bmean[imercind][imercind_min])
print('rise')
print(round(np.mean(mo_bmean[imercind][imercind_rise]),1))
print(round(np.std(mo_bmean[imercind][imercind_rise]),1))
print('max')
print(round(np.mean(mo_bmean[imercind][imercind_max]),1))
print(round(np.std(mo_bmean[imercind][imercind_max]),1))


print()
print('Venus')
print(round(np.mean(mo_bmean[ivexind]),1))
print(round(np.std(mo_bmean[ivexind]),1))
print('min')
print(round(np.mean(mo_bmean[ivexind][ivexind_min]),1))
print(round(np.std(mo_bmean[ivexind][ivexind_min]),1))
print('rise')
print(round(np.mean(mo_bmean[ivexind][ivexind_rise]),1))
print(round(np.std(mo_bmean[ivexind][ivexind_rise]),1))
print('max')
print(round(np.mean(mo_bmean[ivexind][ivexind_max]),1))
print(round(np.std(mo_bmean[ivexind][ivexind_max]),1))

print()
print('Earth')
print(round(np.mean(mo_bmean[iwinind]),1))
print(round(np.std(mo_bmean[iwinind]),1))
print('min')
print(round(np.mean(mo_bmean[iwinind][iwinind_min]),1))
print(round(np.std(mo_bmean[iwinind][iwinind_min]),1))
print('rise')
print(round(np.mean(mo_bmean[iwinind][iwinind_rise]),1))
print(round(np.std(mo_bmean[iwinind][iwinind_rise]),1))
print('max')
print(round(np.mean(mo_bmean[iwinind][iwinind_max]),1))
print(round(np.std(mo_bmean[iwinind][iwinind_max]),1))

print()


#only declining phase
print('MAVEN')
print(round(np.mean(mo_bmean[imavind]),1))
print(round(np.std(mo_bmean[imavind]),1))


###################################################
#F


#TI

print()
print()
print('--------------------------------------------------')
print()
print()

print('Time Inside')

print()
print('Mercury +/-')
print(round(np.nanmean(inside_merc_perc_cycle),1))
print(round(np.nanstd(inside_merc_perc_cycle),1))
#print('min')
#print(round(inside_merc_perc_cycle[0],1))
print('rise')
print(round(inside_merc_perc_cycle[1],1))
print('max')
print(round(inside_merc_perc_cycle[2],1))


print()
print('Venus +/-')
print(round(np.nanmean(inside_vex_perc_cycle),1))
print(round(np.nanstd(inside_vex_perc_cycle),1))
print('min')
print(round(inside_vex_perc_cycle[0],1))
print('rise')
print(round(inside_vex_perc_cycle[1],1))
print('max')
print(round(inside_vex_perc_cycle[2],1))


print()
print('Earth +/-')
print(round(np.nanmean(inside_wind_perc_cycle),1))
print(round(np.nanstd(inside_wind_perc_cycle),1))
print('min')
print(round(inside_wind_perc_cycle[0],1))
print('rise')
print(round(inside_wind_perc_cycle[1],1))
print('max')
print(round(inside_wind_perc_cycle[2],1))


#only declining phase
print('MAVEN')
print(round(inside_mav_perc_cycle[1],1))

histwincyc






###################### MAVEN

#from processing program
#all in days
totaldays=385 
total_icme_duration_maven=np.sum(icme_durations[imavind])/24

print()
print()
print()


print('MAVEN results from 385 days of data, Dec 2014-Feb 2016, with gaps where no solar wind is available')
print('MAVEN total days of observations with solar wind data:')
print(totaldays)
print('MAVEN total ICME durations:')
print(total_icme_duration_maven)
print('Mars is in percent of time inside ICMEs, for intervals in 2014-2016 (declining phase):')
print(total_icme_duration_maven/totaldays*100)
print('on average, Mars is hit by a CME every ... days')
print(totaldays/np.size(imavind))
print('The ICME average duration is, in hours')
print(np.mean(icme_durations[imavind]))












































