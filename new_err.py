import numpy as np
import fitting_scripts
import sys
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import optimize

basedir='/Users/christophermanser/Storage/PhD_files/DESI'
c = 299792.458 # Speed of light in km/s
# requires:
# x: intial guess of T g and rv
# spec: spectrum to fit
# linee: list of lines to fit
# mode=0 is for finding the best fi, mode=1 is for fitting and  retriving a specific model
def fit_func_test(x,spec,linee,models='da2014',mode=0):
    T,g,rv=x
    spectra=spec
    model_test=fitting_scripts.interpolating_model_DA(T,(g/100),mod_type=models)
    try:
        model_f=model_test[:,1]
        model_w=model_test[:,0]
    except:
        print("Could not load the model")
        pass
    else:
        norm_model, m_cont_flux=fitting_scripts.norm_spectra(model_test)
       #Load Balmer line
        line_crop = linee#np.loadtxt(basedir+'/WDFitting/line_crop.dat')
        
    #Check if lines are inside spectra l
        line_crop = line_crop[(line_crop[:,0]>spectra[:,0].min()) & (line_crop[:,1]<spectra[:,0].max())]
    #Normalize the spectrum
        spectra_n, cont_flux = fitting_scripts.norm_spectra(spectra)
    #Interpolate spectrum and model onto same resolution
        m_wave_n=norm_model[:,0]*(rv+c)/c
        m_flux_n=norm_model[:,1]
        tck_l_m = interpolate.interp1d(m_wave_n,m_flux_n,kind='linear')
        spectra_n_w=spectra_n[:,0][spectra_n[:,0]>np.min(m_wave_n)]
        m_flux_n_i = tck_l_m(spectra_n[:,0])

    #Initialise: normalised models and spectra in line region, and chi2
        tmp_lines_m = []
        lines_s = []
        l_chi2 = []
        list=np.array([])
    #for each line
        for c in xrange(len(line_crop)):
        #crop model and spectra to line
        #only if line region is entirely covered by spectrum ...spectra_n[:,0] --> spectra_n_w
            if (line_crop[c,1] < spectra[:,0].max()) & (line_crop[c,0] > spectra[:,0].min()):
                l_m = m_flux_n_i.transpose()[(spectra_n[:,0]>=line_crop[c,0])&(spectra_n[:,0]<=line_crop[c,1])].transpose()
                l_s = spectra_n[(spectra_n[:,0]>=line_crop[c,0])&(spectra_n[:,0]<=line_crop[c,1])]
            #renormalise models to spectra in line region
                l_m = l_m*np.sum(l_s[:,1])/np.sum(l_m)#.reshape([len(l_m),1])
            #calculate chi2
                l_chi2.append( np.sum(((l_s[:,1]-l_m)/l_s[:,2])**2))
                list=np.concatenate((list,((l_s[:,1]-l_m)/l_s[:,2])**2),axis=0)
                tmp_lines_m.append(l_m)
                lines_s.append(l_s)
    #mean chi2 over lines
        l_chi2 = np.array(l_chi2)
        model_test=np.dstack((model_w, model_f))[0]
        if mode==1:
            return lines_s, tmp_lines_m, spectra, model_test#,l_chi2#, x_list
        elif mode==0:
           
            return np.sum(l_chi2) #this is the quantity that gets minimized 


#================================================================================================================
# This script finds errors by minimizing the function at chi+1 rather than chi
# requires:
# x: intial guess of T and g
# rv: rv value from best fit
# valore: the chi value of the best fit
# spec: spectrum to fit
# linee: list of lines to fit

def err_t(x,rv,valore,spec,linee,models='da2014'):
    T,g=x
    spectra=spec
    model_test=fitting_scripts.interpolating_model_DA(T,(g/100),mod_type=models)
    try:
        model_f=model_test[:,1]
        model_w=model_test[:,0]
    except:
        print("Could not load the model")
        pass
    else:
        norm_model, m_cont_flux=fitting_scripts.norm_spectra(model_test)
       #Load Balmer line
        line_crop = linee#np.loadtxt(basedir+'/WDFitting/line_crop.dat')
        
    #Check if lines are inside spectra l
        line_crop = line_crop[(line_crop[:,0]>spectra[:,0].min()) & (line_crop[:,1]<spectra[:,0].max())]
    #Normalize the spectrum
        spectra_n, cont_flux = fitting_scripts.norm_spectra(spectra)
    #Interpolate spectrum and model onto same resolution
        m_wave_n=norm_model[:,0]*(rv+c)/c
        m_flux_n=norm_model[:,1]
        tck_l_m = interpolate.interp1d(m_wave_n,m_flux_n,kind='linear')
        spectra_n_w=spectra_n[:,0][spectra_n[:,0]>np.min(m_wave_n)]
        m_flux_n_i = tck_l_m(spectra_n[:,0])

    #Initialise: normalised models and spectra in line region, and chi2
        tmp_lines_m = []
        lines_s = []
        l_chi2 = []
        list=np.array([])
    #for each line
        for c in xrange(len(line_crop)):
        #crop model and spectra to line
        #only if line region is entirely covered by spectrum ...spectra_n[:,0] --> spectra_n_w
            if (line_crop[c,1] < spectra[:,0].max()) & (line_crop[c,0] > spectra[:,0].min()):
                l_m = m_flux_n_i.transpose()[(spectra_n[:,0]>=line_crop[c,0])&(spectra_n[:,0]<=line_crop[c,1])].transpose()
                l_s = spectra_n[(spectra_n[:,0]>=line_crop[c,0])&(spectra_n[:,0]<=line_crop[c,1])]
            #renormalise models to spectra in line region
                l_m = l_m*np.sum(l_s[:,1])/np.sum(l_m)#.reshape([len(l_m),1])
            #calculate chi2
                l_chi2.append( np.sum(((l_s[:,1]-l_m)/l_s[:,2])**2))
                list=np.concatenate((list,((l_s[:,1]-l_m)/l_s[:,2])**2),axis=0)
                tmp_lines_m.append(l_m)
                lines_s.append(l_s)
    #mean chi2 over lines
        l_chi2 = np.array(l_chi2)
        model_test=np.dstack((model_w, model_f))[0]
        return abs(np.sum(l_chi2)-(valore+1.)) #this is the quantity that gets minimized 
