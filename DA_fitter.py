import numpy as np
import fitting_scripts
import sys
import matplotlib.pyplot as plt
from scipy import optimize
import scipy
import new_err
from multiprocessing import Process, Queue
import multiprocessing

model_c='da2014'#'pier' #da2014
basedir='/Users/christophermanser/Storage/PhD_files/DESI'

infile=sys.argv[1]
spec1,spec2,spec3=np.loadtxt(infile,usecols=(0,1,2),unpack=True)
spec1=spec1[np.isnan(spec2)==False]
spec3=spec3[np.isnan(spec2)==False]
spec2=spec2[np.isnan(spec2)==False]


#load lines to fit
line_crop = np.loadtxt(basedir+'/WDFitting/line_crop.dat')


spec2=spec2[(spec1>3500)] #& (spec1<6200)]
spec3=spec3[(spec1>3500)]
spec1=spec1[(spec1>3500)] #& (spec1<6200)]

#-----------------------conversion from air to vacuum----------------------------------------------
#spec1= spec_v / (1.0 + 2.735182e-4 + 131.4182/np.power(spec_v,2.) + 2.76249e8/np.power(spec_v,4.))
#--------------------------------------------------------------------------------------------------


spectra=np.dstack((spec1,spec2,spec3))[0]
#fit entire grid to find good starting point
best=fitting_scripts.fit_line(spectra, model_in=None,  quick=True,  model=model_c)

first_T= best[2][0]
first_g= best[2][1]


all_chi=best[4]
all_TL=best[3]

#------find starting point for secondary solution
if first_T < 13000.:
    other_TL=all_TL[all_TL[:,0]>13000.]
    other_chi=all_chi[all_TL[:,0]>13000.]
    other_sol= other_TL[other_chi==np.min(other_chi)]


if first_T > 13000.:
    other_TL=all_TL[all_TL[:,0]<13000.]
    other_chi=all_chi[all_TL[:,0]<13000.]
    other_sol= other_TL[other_chi==np.min(other_chi)]


teff_list=[]
logg_list=[]
res_list=[]
#find best fitting model
new_best=optimize.fmin(new_err.fit_func_test,(first_T,first_g,10.),args=(spectra,line_crop,model_c,0),retall=0,disp=0,xtol=1.,ftol=1.,full_output=1)
#calculate error
other_T=optimize.fmin(new_err.err_t,(first_T,first_g),args=(new_best[0][2],new_best[1],spectra,line_crop,model_c),retall=0,disp=0,xtol=1.,ftol=1.,full_output=0)

best_T=new_best[0][0]
best_g=new_best[0][1]
best_rv=new_best[0][2]
print "First solution"
print best_T, abs(best_T-other_T[0])
print best_g, abs(best_g-other_T[1])
print "rv=",best_rv
# get and save best model
lines_s,lines_m,spectra,model_n=new_err.fit_func_test((best_T,best_g,best_rv),spectra,line_crop,models=model_c,mode=1)



#repeat fit for secondary solution
second_best=optimize.fmin(new_err.fit_func_test,(other_sol[0][0],other_sol[0][1],best_rv),args=(spectra,line_crop,model_c,0),retall=0,disp=0,xtol=1.,ftol=1.,full_output=1)
other_T2=optimize.fmin(new_err.err_t,(other_sol[0][0],other_sol[0][1]),args=(best_rv,second_best[1],spectra,line_crop,model_c),retall=0,disp=0,xtol=1.,ftol=1.,full_output=0)


s_best_T=second_best[0][0]
s_best_g=second_best[0][1]
s_best_rv=second_best[0][2]

# get and save best model
lines_s_o,lines_m_o,spectra,model_n_o=new_err.fit_func_test((s_best_T,s_best_g,s_best_rv),spectra,line_crop,models=model_c,mode=1)
print" "
print "Second solution"
print s_best_T, abs(s_best_T-other_T2[0])
print s_best_g, abs(s_best_g-other_T2[1])


#=======================plotting===============================================


fig=plt.figure(figsize=(8,5))
ax = plt.subplot2grid((1,4), (0, 3))
#plt.plot(lines_s[0][:,0]-np.median(lines_s[0][:,0]),lines_s[0][:,1],color='k')#-20
#plt.plot(lines_s[0][:,0]-np.median(lines_s[0][:,0]),lines_m[0],color='r')

#plt.plot(lines_s[0][:,0]-lines_s[0][:,0][lines_s[0][:,1]==np.min(lines_s[0][:,1])][0],lines_s[0][:,1],color='k')#-20
#plt.plot(lines_s[0][:,0]-lines_s[0][:,0][lines_s[0][:,1]==np.min(lines_s[0][:,1])][0],lines_m[0],color='r')


plt.plot(lines_s[1][:,0]-lines_s[1][:,0][lines_s[1][:,1]==np.min(lines_s[1][:,1])][0],lines_s[1][:,1]+0.5,color='k')
plt.plot(lines_s[1][:,0]-lines_s[1][:,0][lines_s[1][:,1]==np.min(lines_s[1][:,1])][0],lines_m[1]+0.5,color='r')

plt.plot(lines_s[2][:,0]-lines_s[2][:,0][lines_s[2][:,1]==np.min(lines_s[2][:,1])][0],lines_s[2][:,1]+1,color='k')
plt.plot(lines_s[2][:,0]-lines_s[2][:,0][lines_s[2][:,1]==np.min(lines_s[2][:,1])][0],lines_m[2]+1,color='r')


plt.plot(lines_s[3][:,0]-lines_s[3][:,0][lines_s[3][:,1]==np.min(lines_s[3][:,1])][0],lines_s[3][:,1]+1.5,color='k')
plt.plot(lines_s[3][:,0]-lines_s[3][:,0][lines_s[3][:,1]==np.min(lines_s[3][:,1])][0],lines_m[3]+1.5,color='r')

plt.plot(lines_s[4][:,0]-lines_s[4][:,0][lines_s[4][:,1]==np.min(lines_s[4][:,1])][0],lines_s[4][:,1]+2,color='k')
plt.plot(lines_s[4][:,0]-lines_s[4][:,0][lines_s[4][:,1]==np.min(lines_s[4][:,1])][0],lines_m[4]+2,color='r')


plt.plot(lines_s[5][:,0]-lines_s[5][:,0][lines_s[5][:,1]==np.min(lines_s[5][:,1])][0],lines_s[5][:,1]+2.5,color='k')
plt.plot(lines_s[5][:,0]-lines_s[5][:,0][lines_s[5][:,1]==np.min(lines_s[5][:,1])][0],lines_m[5]+2.5,color='r')


#--------------------------------------------------------------------------------------
#plt.plot(lines_s_o[0][:,0]-lines_s_o[0][:,0][lines_s_o[0][:,1]==np.min(lines_s_o[0][:,1])][0],lines_m_o[0],color='g')
plt.plot(lines_s_o[1][:,0]-lines_s_o[1][:,0][lines_s_o[1][:,1]==np.min(lines_s_o[1][:,1])][0],lines_m_o[1]+0.5,color='g')
plt.plot(lines_s_o[2][:,0]-lines_s_o[2][:,0][lines_s_o[2][:,1]==np.min(lines_s_o[2][:,1])][0],lines_m_o[2]+1,color='g')
plt.plot(lines_s_o[3][:,0]-lines_s_o[3][:,0][lines_s_o[3][:,1]==np.min(lines_s_o[3][:,1])][0],lines_m_o[3]+1.5,color='g')
plt.plot(lines_s_o[4][:,0]-lines_s_o[4][:,0][lines_s_o[4][:,1]==np.min(lines_s_o[4][:,1])][0],lines_m_o[4]+2,color='g')
plt.plot(lines_s_o[5][:,0]-lines_s_o[5][:,0][lines_s_o[5][:,1]==np.min(lines_s_o[5][:,1])][0],lines_m_o[5]+2.5,color='g')



xticks = ax.xaxis.get_major_ticks()

ax.set_xticklabels([])
ax.set_yticklabels([])


ax2 = plt.subplot2grid((1,4), (0, 0),colspan=3)

#xdata2=[]
#ydata2=[]
#weight2= 1/(spectra[:,2]**2)
#err_data2=[]
#binsize=5
#for i in range(0,(np.size(spectra[:,0])-binsize),binsize):
#        xdata2.append(np.average(spectra[:,0][i:i+binsize]))
#        ydata2.append(np.average(spectra[:,1][i:i+binsize],weights=weight2[i:i+binsize]))
 #       err_data2.append(1/np.sqrt(np.sum(weight2[i:i+binsize])))

#spectra=np.stack((xdata2,ydata2,err_data2), axis=-1)



plt.plot(spectra[:,0],spectra[:,1],color='k')
#plt.plot(xdata2,ydata2,color='k')
model_n[np.isnan(model_n)]=0.0
model_flux=model_n[:,1]
model_wave=model_n[:,0]
spec_flux=spectra[:,1]
spec_wave=spectra[:,0]
check_f_model=model_flux[np.logical_and(model_wave>4500., model_wave<4700.)]
check_f_spec=spec_flux[np.logical_and(spec_wave>4500., spec_wave<4700.)]
adjust=np.average(check_f_model)/np.average(check_f_spec)
plt.plot(model_n[:,0]*(best_rv+300000.)/300000.,model_n[:,1]/adjust,color='r')





model_n_o[np.isnan(model_n_o)]=0.0
model_flux=model_n_o[:,1]
model_wave=model_n_o[:,0]
spec_flux=spectra[:,1]
spec_wave=spectra[:,0]
check_f_model=model_flux[np.logical_and(model_wave>4500., model_wave<4700.)]
check_f_spec=spec_flux[np.logical_and(spec_wave>4500., spec_wave<4700.)]
adjust=np.average(check_f_model)/np.average(check_f_spec)
plt.plot(model_n_o[:,0],model_n_o[:,1]/adjust,color='g')



plt.ylabel(r'F$_{\lambda}$ [erg cm$^{-2}$ s$^{-1} \AA^{-1}$]',fontsize=12)#labelpad=24
plt.xlabel(r'Wavelength $(\AA)$',fontsize=12)


#fig.tight_layout()



plt.show()
