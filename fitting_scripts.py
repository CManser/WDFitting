import numpy as _np
import os
import exceptions
from pylab import *
basedir='/Users/christophermanser/Storage/PhD_files/DESI'
def _band_limits(band):
    """
    give magnitude band e.g. "sdss_r"
    return outer limits
    for use in models etc
    """
    mag = _np.loadtxt(basedir+'/sm/'+band+'.dat')
    mag = mag[mag[:,1]>0.05]
    return [mag[:,0].min(),mag[:,0].max()]


def models_normalised(quick=True, model='sdss', testing=False):
    """
    Import Normalised WD Models
    No Arguements
    Optional:
        quick=True   : Use presaved model array. Check is up to date
        model='sdss': Which model grid to use: 'sdss' (DA, fine, noIR), 'new' (DA, course, IR, new), 'old' (DA, course, IR, old), 'interp' (DA, fine++, noIR)
        testing=False      : plot testing image
    Return [model_list,model_param,orig_model_wave,orig_model_flux,tck_model,r_model]
    """
    #Use preloaded tables
    if quick:
        #
        if (model!='sdss') & (model!='new') & (model!='old') & (model!='fine++')& (model!='da2014')& (model!='pier')& (model!='pier3D')& (model!='pier3D_smooth') & (model!='pier_rad')& (model!='pier1D')& (model!='pier_smooth')& (model!='pier_rad_smooth')& (model!='pier_rad_fullres')& (model!='pier_fullres'):
            raise wdfitError('Unknown "model" in models_normalised')
        
        fn = '/WDFitting/wdfit.'+model+'.lst'
        d = '/WDFitting/WDModels_Koester.'+model+'_npy/'
        model_list = _np.loadtxt(basedir+fn, usecols=[0], dtype='string', comments='WDFitting/')
        model_param = _np.loadtxt(basedir+fn, usecols=[1,2], comments='WDFitting/')
        m_spec = _np.load(basedir+d+model_list[0])
        m_wave = m_spec[:,0]
        #
        out_m_wave = m_wave[(m_wave>=3400)&(m_wave<=13000)]
        norm_m_flux = _np.load(basedir+'/WDFitting/norm_m_flux.'+model+'.npy')
        #
        if out_m_wave.shape[0] != norm_m_flux.shape[1]:
            raise wdfitError('l and f arrays not correct shape in models_normalised')
    #Calculate
    else:
        from scipy import interpolate
        #Load models from models()
        model_list,model_param,m_wave,m_flux,tck_model,r_model = models(quick=quick, model=model)
        #Range over which DA is purely continuum
        norm_range = _np.loadtxt(basedir+'/WDFitting/wide_norm_range.dat', usecols=[0,1])
        norm_range_s = _np.loadtxt(basedir+'/WDFitting/wide_norm_range.dat', usecols=[2], dtype='string')
        #for each line, for each model, hold l and f of continuum
        cont_l = _np.empty([len(m_flux), len(norm_range)])
        cont_f = _np.empty([len(m_flux), len(norm_range)])
        #for each zone
        
        for j in xrange(len(norm_range)):
            if (norm_range[j,0] < m_wave.max()) & (norm_range[j,1] > m_wave.min()):
                #crop
                _f = m_flux.transpose()[(m_wave>=norm_range[j,0])&(m_wave<=norm_range[j,1])].transpose()
                _l = m_wave[(m_wave>=norm_range[j,0])&(m_wave<=norm_range[j,1])]
                
                #interpolate region
                print norm_range[j,0], norm_range[j,1]
                print np.size(_l)

                tck = interpolate.interp1d(_l,_f,kind='cubic')
                #interpolate onto 10* resolution
                l = _np.linspace(_l.min(),_l.max(),(len(_l)-1)*10+1)
                f = tck(l)
                #print f
                #find maxima and save
                if norm_range_s[j]=='P':
                    for k in xrange(len(f)):
                        cont_l[k,j] = l[f[k]==f[k].max()][0]
                        cont_f[k,j] = f[k].max()
                #find mean and save
                elif norm_range_s[j]=='M':
                    for k in xrange(len(f)):
                        cont_l[k,j] = _np.mean(l)
                        cont_f[k,j] = _np.mean(f[k])
                else:
                    print 'Unknown norm_range_s, ignoring'
        #Continuum
        if (norm_range.min()>3400) & (norm_range.max()<13000):
            out_m_wave = m_wave[(m_wave>=3400)&(m_wave<=13000)]
        else:
            raise wdfitError('Normalised models cropped to too small a region')
        #
        cont_m_flux = _np.empty([len(m_flux),len(out_m_wave)])
        for i in xrange(len(m_flux)):
            #not suitable for higher order fitting
            #tck = interpolate.splrep(m_wave_nr[i],m_flux_nr[i], t=[4100,4340,4900,6460], k=3)
            #print cont_f[i]
            tck = interpolate.splrep(cont_l[i],cont_f[i], t=[3885,4340,4900,6460], k=3)
            cont_m_flux[i] = interpolate.splev(out_m_wave,tck)
        #Normalised flux
        norm_m_flux = m_flux.transpose()[(m_wave>=3400)&(m_wave<=13000)].transpose()/cont_m_flux
        _np.save(basedir+'/WDFitting/norm_m_flux.'+model+'.npy', norm_m_flux)
        #
        #testing
        if testing:
            import pylab as pl
            def p():
                pl.figure(figsize=(7,9))
                ax1 = pl.subplot(211)
                pl.axvline([3885], color='g', zorder=1)
                pl.axvline([4340], color='g', zorder=1)
                pl.axvline([4900], color='g', zorder=1)
                pl.axvline([6460], color='g', zorder=1)
                pl.plot(m_wave,m_flux[i], color='grey', lw=0.8, zorder=2)
                pl.plot(out_m_wave,cont_m_flux[i], 'b-', zorder=3)
                pl.scatter(cont_l[i], cont_f[i], edgecolors='r', facecolors='none', zorder=20)
                ax2 = pl.subplot(212, sharex=ax1)
                pl.axhline([1], color='g')
                pl.plot(out_m_wave, norm_m_flux[i], 'b-')
                pl.ylim([0,2])
                pl.xlim([3400,13000])
                pl.show()
                return
            #
            for i in _np.where(model_param[:,1]==8.0)[0][::8]:
                p()
        #
    #
    return [out_m_wave,norm_m_flux,model_list,model_param]



def norm_spectra(spectra, add_infinity=True, testing=False):
    """
    Normalised spectra by DA WD continuum regions
    spectra of form array([wave,flux,error]) (err not necessary)
    only works on SDSS spectra region
    Optional:
        EDIT wide_norm_range to change whether the region is fitted for a peak or mean'd
        add_infinity=False : add a spline point at [infinity,0]
        testing=False      : plot testing image
    returns spectra, cont_flux
    """
    from scipy import interpolate
    #Range over which DA is purely continuum
    #norm_range = _np.loadtxt(basedir+'/WDFitting/wide_norm_range.dat', usecols=[0,1])
    #norm_range_s = _np.loadtxt(basedir+'/WDFitting/wide_norm_range.dat', usecols=[2], dtype='string')
    start_n=_np.array([3770.,3796.,3835.,3895.,3995.,4130.,4490.,4620.,5070.,5200.,6000.,7000.,7550.,8400.])
    end_n=_np.array([3795.,3830.,3885.,3960.,4075.,4290.,4570.,4670.,5100.,5300.,6100.,7050.,7600.,8450.])
    #norm_range=np.stack((start_n, end_n), axis=-1)
    norm_range=np.dstack((start_n, end_n))[0]
    norm_range_s=_np.array(['P','P','P','P','P','P','M','M','M','M','M','M','M','M'])
      #
    if len(spectra[0,:])>2:
        s_nr = _np.zeros([len(norm_range),3])
        spectra[:,2][spectra[:,2]==0.] = spectra[:,2].max()
    else:
        s_nr = _np.zeros([len(norm_range),2])
    #for each zone
    for j in xrange(len(norm_range)):
        if (norm_range[j,0] < spectra[:,0].max()) & (norm_range[j,1] > spectra[:,0].min()):
            #crop
            _s = spectra[(spectra[:,0]>=norm_range[j,0])&(spectra[:,0]<=norm_range[j,1])]
            #Avoids gappy spectra
            #More points than k
            k=3
            if len(_s)>k:
                #interpolate onto 10* resolution
                l = _np.linspace(_s[:,0].min(),_s[:,0].max(),(len(_s)-1)*10+1)
                if len(spectra[0,:])>2:
                    tck = interpolate.splrep(_s[:,0],_s[:,1],w=1/_s[:,2], s=1000)
                    #median errors for max/mid point
                    s_nr[j,2] = _np.median(_s[:,2]) / _np.sqrt(len(_s[:,0]))
                else:
                    tck = interpolate.splrep(_s[:,0],_s[:,1],s=0.0)
                #
                f = interpolate.splev(l,tck)
                #find maxima and save
                if norm_range_s[j]=='P':
                    s_nr[j,0] = l[f==f.max()][0]
                    s_nr[j,1] = f.max()
                #find mean and save
                elif norm_range_s[j]=='M':
                    s_nr[j,0:2] = _np.mean(l), _np.mean(f)
                else:
                    print 'Unknown norm_range_s, ignoring'
    #
    s_nr = s_nr[ s_nr[:,0] != 0 ]
    #
    #t parameter chosen by eye. Position of knots. need more near blue end
    #allow for short spectra though not as effective
    #not suitable for higher order fitting
    #
    if s_nr[:,0].max() < 4901:
        #print 'Warning: knots used for spline normalisation probably bad'
        knots=None
    elif s_nr[:,0].max() < 6460:
        #print 'Warning: knots used for spline normalisation may be bad'
        #knots = [3885,4340,4900,int(s_nr[:,0].max()-5)]
        knots = [3000,4900,4100,4340,4860,int(s_nr[:,0].max()-5)]
        linee=_np.array(knots)

#knots = [4100,4340,4900,int(s_nr[:,0].max()-5)]
    else:
        knots = [3885,4340,4900,6460]
        #knots = [4100,4340,4900,6460]
    #
    if s_nr[:,0].min() > 3885:
        print 'Warning: knots used for spline normalisation not suitable for high order fitting'
        knots=knots[1:]
    if s_nr[:,0].min() > 4340:
        print 'Warning: knots used for spline normalisation probably bad'
        knots=None
    #
    ##add a point at infinity, 0 for spline to fit too. error = mean of other errors.
    if add_infinity:
        if s_nr.shape[1] > 2:
            s_nr = _np.vstack([ s_nr, _np.array([90000.,0., _np.mean(s_nr[:,2]) ]) ])
            s_nr = _np.vstack([ s_nr, _np.array([100000.,0., _np.mean(s_nr[:,2]) ]) ])
        else:
            s_nr = _np.vstack([ s_nr, _np.array([90000.,0.]) ])
            s_nr = _np.vstack([ s_nr, _np.array([100000.,0.]) ])
    #
    #weight by errors
    try:
        if len(spectra[0,:])>2:
            tck = interpolate.splrep(s_nr[:,0],s_nr[:,1], w=1/s_nr[:,2], t=knots, k=3)
        else:
            tck = interpolate.splrep(s_nr[:,0],s_nr[:,1], t=knots, k=3)
    except ValueError:
        #print 'Warning: knots used for spline normalisation failed, now probably bad'
        knots=None
        if len(spectra[0,:])>2:
            tck = interpolate.splrep(s_nr[:,0],s_nr[:,1], w=1/s_nr[:,2], t=knots, k=3)
        else:
            tck = interpolate.splrep(s_nr[:,0],s_nr[:,1], t=knots, k=3)
    #
    cont_flux = interpolate.splev(spectra[:,0],tck)
    cont_flux = cont_flux.reshape([len(cont_flux),1])
    spectra_ret = _np.copy(spectra)
    spectra_ret[:,1:] = spectra_ret[:,1:]/cont_flux
    #
    #testing
    if testing:
        import pylab as pl
        from jg import spectra as s
        def p():
            pl.figure(figsize=(7,9))
            ax1 = pl.subplot(211)
            for bla in range(_np.size(linee)):
                pl.axvline(linee[bla], color='g', zorder=1)
            #pl.axvline([4340], color='g', zorder=1)
            #pl.axvline([4900], color='g', zorder=1)
            #pl.axvline([6460], color='g', zorder=1)
            bs =  s.bin_spc(spectra)
            pl.plot(bs[:,0], bs[:,1], color='grey', lw=0.8, zorder=2)
            pl.plot(spectra[:,0], cont_flux, 'b-', zorder=3)
            pl.scatter(s_nr[:,0], s_nr[:,1], edgecolors='r', facecolors='none', zorder=20)
            pl.ylim([0, _np.max([spectra[:,1].max(), cont_flux.max(), s_nr[:,1].max()])*1.2])
            ax2 = pl.subplot(212, sharex=ax1)
            pl.axhline([1], color='g')
            pl.plot(spectra_ret[:,0], spectra_ret[:,1], color='grey', lw=0.8)
            pl.ylim([0,2])
            pl.xlim([3400,13000])
            pl.show()
            return
        #
        p()
    #
    return spectra_ret, cont_flux


def models(quick=True, quiet=True, band='sdss_r', model='sdss'):
    """
    Import WD Models
    Optional:
        quick=True   : Use presaved model array. Check is up to date
        quiet=True   : verbose
        band='sdss_r': which mag band to calculate normalisation over (MEAN, not folded)
        model='sdss': Which model grid to use: 'sdss' (DA, fine, noIR), 'new' (DA, course, IR, new), 'old' (DA, course, IR, old), 'interp' (DA, fine++, noIR), 'db' (DB models)
    Return [model_list,model_param,orig_model_wave,orig_model_flux,tck_model,r_model]
    """
    from scipy import interpolate
    #
    if (model!='sdss') & (model!='new') & (model!='old') & (model!='fine++') & (model!='db') & (model!='da2014')& (model!='pier')& (model!='pier3D')& (model!='pier3D_smooth')&(model!='pier_rad') & (model!='pier1D')& (model!='pier_smooth')& (model!='pier_rad_smooth')& (model!='pier_rad_fullres')& (model!='pier_fullres'):
        raise wdfitError('Unknown "model" in models')
    else:
        fn = '/WDFitting/wdfit.'+model+'.lst'
        d = '/WDFitting/WDModels_Koester.'+model+'/'
    #
    #Load in table of all models
    model_list = _np.loadtxt(basedir+fn, usecols=[0], dtype='string', comments='WDFitting/')
    model_param = _np.loadtxt(basedir+fn, usecols=[1,2], comments='WDFitting/')
    orig_model_wave = _np.loadtxt(basedir+d+model_list[0], usecols=[0])
    #
    if not quiet:
        print 'Loading Models'
    if quick:
        orig_model_flux = _np.load(basedir+'/WDFitting/orig_model_flux.'+model+'.npy')
        #
        if orig_model_wave.shape[0] != orig_model_flux.shape[1]:
            raise wdfitError('l and f arrays not correct shape in models')
    else:
        orig_model_flux = _np.empty([len(model_list),len(orig_model_wave)])
        if not(model=='db'):
            for i in xrange(len(model_list)):
                if not quiet:
                    print i
                print basedir+d+model_list[i]
                orig_model_flux[i] = _np.loadtxt(basedir+d+model_list[i],usecols=[1])
            #
            _np.save(basedir+'/WDFitting/orig_model_flux.'+model+'.npy',orig_model_flux)
        else:
            from jg import spectra as _s
            for i in xrange(len(model_list)):
                if not quiet:
                    print i
                #
                tmp = _s.spectra(fn = basedir+d+model_list[i], usecols=[0,1])
                #Not uniform wavelength grid argh!
                tmp.interpolate(orig_model_wave, kind='linear', save_res=True)
                orig_model_flux[i] = tmp.f()
        #
        _np.save(basedir+'/WDFitting/orig_model_flux.'+model+'.npy',orig_model_flux)
    #Interpolate Model onto Spectra points
    #Linear
    tck_model = interpolate.interp1d(orig_model_wave,orig_model_flux,kind='linear')
    #Only calculate r model once
    band_limits = _band_limits(band)
    r_model = _np.mean(orig_model_flux.transpose()[((orig_model_wave>=band_limits[0])&(orig_model_wave<=band_limits[1]))].transpose(),axis=1)
    #
    print orig_model_wave
    print orig_model_flux
    return [model_list,model_param,orig_model_wave,orig_model_flux,tck_model,r_model]



    

# Models interpolation #
# for speed the models need to be saved specifically as numpy arrays
def interpolating_model_DA(temperature,gravity,mod_type='pier',band='none',mag=0,radius=0,distance=10,Av=0):
	"""Interpolate model atmospheres given an input Teff and logg.
	If band=='V': normalizes the model at lambda = 5500A
        If Av!=0: the Fitzpatrick 1999 law is used"""
	# PARAMETERS # 
	dir_models = basedir + '/WDFitting/WDModels_Koester.'+mod_type+'_npy/'
        if mod_type=="pier":	
            teff=np.array([1500.,1750.,2000.,2250.,2500.,2750.,3000.,3250.,3500.,3750.,4000.,4250.,4500.,4750.,5000.,5250.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,11500.,12000.,12500.,13000.,13500.,14000.,14500.,15000.,15500.,16000.,16500.,17000.,20000.,25000.,30000.,35000.,40000.,45000.,50000.,55000.,60000.,65000.,70000.,75000.,80000.,85000.,90000.])
            logg=np.array([6.50,7.00,7.50,7.75,8.00,8.25,8.50,9.00,9.50])

        elif mod_type=="pier_fullres":	
            teff=np.array([1500.,1750.,2000.,2250.,2500.,2750.,3000.,3250.,3500.,3750.,4000.,4250.,4500.,4750.,5000.,5250.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,11500.,12000.,12500.,13000.,13500.,14000.,14500.,15000.,15500.,16000.,16500.,17000.,20000.,25000.,30000.,35000.,40000.,45000.,50000.,55000.,60000.,65000.,70000.,75000.,80000.,85000.,90000.])
            logg=np.array([6.50,7.00,7.50,7.75,8.00,8.25,8.50,9.00,9.50])

        elif mod_type=="pier3D_smooth":	
            teff=np.array([1500.,1750.,2000.,2250.,2500.,2750.,3000.,3250.,3500.,3750.,4000.,4250.,4500.,4750.,5000.,5250.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,11500.,12000.,12500.,13000.,13500.,14000.,14500.,15000.,15500.,16000.,16500.,17000.,20000.,25000.,30000.,35000.,40000.,45000.,50000.,55000.,60000.,65000.,70000.,75000.,80000.,85000.,90000.])
            logg=np.array([7.00,7.50,8.00,8.50,9.00])
        elif mod_type=="pier_smooth":	
            teff=np.array([6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,11500.,12000.,12500.,13000.,13500.,14000.,14500.,15000.])
            logg=np.array([7.00,7.50,8.00,8.50,9.00])

        elif mod_type=="pier_rad":	
            teff=np.array([6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,11500.,12000.,12500.,13000.,13500.,14000.,14500.,15000.])
            logg=np.array([7.00,7.50,8.00,8.50,9.00])

        elif mod_type=="pier_rad_smooth":	
            teff=np.array([6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,11500.,12000.,12500.,13000.,13500.,14000.,14500.,15000.])
            logg=np.array([7.00,7.50,8.00,8.50,9.00])

        elif mod_type=="da2014":
            teff=np.array([6000.,6250.,6500.,6750.,7000.,7250.,7500.,7750.,8000.,8250.,8500.,8750.,9000.,9250.,9500.,9750.,10000.,10100.,10200.,10250.,10300.,10400.,10500.,10600.,10700.,10750.,10800.,10900.,11000.,11100.,11200.,11250.,11300.,11400.,11500.,11600.,11700.,11750.,11800.,11900.,12000.,12100.,12200.,12250.,12300.,12400.,12500.,12600.,12700.,12750.,12800.,12900.,13000.,13500.,14000.,14250.,14500.,14750.,15000.,15250.,15500.,15750.,16000.,16250.,16500.,16750.,17000.,17250.,17500.,17750.,18000.,18250.,18500.,18750.,19000.,19250.,19500.,19750.,20000.,21000.,22000.,23000.,24000.,25000.,26000.,27000.,28000.,29000.,30000.,35000.,40000.,45000.,50000.,55000.,60000.,65000.,70000.,75000.,80000.,90000.,100000.])
            logg=np.array([4.00,4.25,4.50,4.75,5.00,5.25,5.50,5.75,6.00,6.25,6.50,6.75,7.00,7.25,7.50,7.75,8.00,8.25,8.50,8.75,9.00,9.25,9.50])
        #print logg
        if mod_type=='pier3D':
            if temperature<6000. or temperature>90000. or gravity<6.5 or gravity>9.:
                #print 'BOOM 3D'
                return [],[]
        elif mod_type=='pier_rad':
            if temperature<6000. or temperature>15000. or gravity<7.0 or gravity>9.:
                #print 'BOOM rad'
                return [],[]
        elif mod_type=='pier_smooth':
            if temperature<6000. or temperature>15000. or gravity<7.0 or gravity>9.:
                #print 'BOOM'
                return [],[]
        elif mod_type=='da2014':
            if temperature<6000. or temperature>100000. or gravity<4.0 or gravity>9.5:
                #print 'BOOM'
                return [],[]
	
	# INTERPOLATION #
	g1,g2 = _np.max(logg[logg<=gravity]),_np.min(logg[logg>=gravity])
        if mod_type=='da2014':
            t1,t2 = str(int(max(teff[teff<=temperature]))).zfill(6),str(int(min(teff[teff>=temperature]))).zfill(6)
        else:
            t1,t2 = (int(max(teff[teff<=temperature]))),(int(min(teff[teff>=temperature])))
        #else:
         #   t1,t2 = t1.zfill(6),t2.zfill(6)
	models = []
	
        for j in [t1,t2]:
		for i in [g1,g2]:
                    #print 'WD_%.2f_%d.0.dat' % (i,j)
                    if mod_type =='da2014':
                        i=int(i*100)
                        #models.append([x for x in lista_models if x == 'da%s_%s_2.7.npy' % (j,str(int(i*100)))])
                        mod_name='da'+str(j)+'_'+str(i)+'_2.7.npy'
                        models.append(mod_name)
                    else:
                        i=('%.2f' % i)
                        mod_name='WD_'+str(i)+'_'+str(j)+'.0.npy'
                        models.append(mod_name)

	model = []
	if len(models[0])!=0:
	    m11=_np.load(dir_models+models[0])
	    model = _np.zeros((len(m11),2))
	else:
	    m11 = 0
	if len(models[1])!=0:    
	    m12=_np.load(dir_models+models[1])
	    model = _np.zeros((len(m11),2))
	else:
	    m12 = 0	
	if len(models[2])!=0:   
	    m21=_np.load(dir_models+models[2])
	    model = _np.zeros((len(m11),2))
	else:
	    m21 = 0	
	if len(models[3])!=0:    
	    m22=_np.load(dir_models+models[3])
	    model = _np.zeros((len(m11),2))
	else:
	    m22 = 0	
	
	if t1!=t2:
		t = (temperature-float(t1))/(float(t2)-float(t1))
                
	else:
		t=0
	if _np.isnan(t)==True:
		t=0	
	if g1!=g2:
		g = (gravity-g1)/(g2-g1)
	else:
		g=0
	if _np.isnan(g)==True:
		g=0	
	if len(model)!=0:
            flux_i = (1-t)*(1-g)*m11[:,1]+t*(1-g)*m21[:,1]+t*g*m22[:,1]+(1-t)*g*m12[:,1]
	    wav_i = m11[:,0]
	    model[:,0],model[:,1]=wav_i,flux_i
	else:
	    return [], []
	
		#return models,model
        return model






def corr3d(temperature,gravity,ml2a=0.8,graph=False):
	""" Determines the 3D correction (Tremblay et al. 2013, 559, A104)
	    from their atmospheric parameters
	    Optional, measure in parsecs
	    ML2/alpha = 0.8 is implemented
	    """	
	temperature, gravity = float(temperature), float(gravity)
	if temperature > 14500. or temperature < 6000. or gravity > 9.0 or gravity < 7.0:
	   #print "no correction"
	   return 0.,0.
	if ml2a==0.8:
		teff_corr = csv2rec('/storage/astro2/phsmav/data2/models/3dcorr/ml2a_08/teff_corr.csv')
		logg_corr = csv2rec('/storage/astro2/phsmav/data2/models/3dcorr/ml2a_08/logg_corr.csv')
		
	        teff = teff_corr['teff']
		logg = teff_corr['logg']
		
		# indentifying ranges Teff_corr#
		t1,t2 = _np.max(teff[teff<=temperature]),_np.min(teff[teff>=temperature])
		g1,g2 = _np.max(logg[logg<=gravity]),_np.min(logg[logg>=gravity])
		if t1!=t2:
			t1,t2 = float(t1), float(t2)
			t = (temperature-t1)/(t2-t1)
		else:
			t=0.
		if _np.isnan(t)==True:
			t=0.
		if g1!=g2:
			g1,g2 = float(g1), float(g2)
			g = (gravity-g1)/(g2-g1)
		else:
			g=0
		if _np.isnan(g)==True:
			g=0		
		m11 = teff_corr[logical_and(teff==t1,logg==g1)]['teff_corr'][0]
		m12 = teff_corr[logical_and(teff==t1,logg==g2)]['teff_corr'][0]
		m21 = teff_corr[logical_and(teff==t2,logg==g1)]['teff_corr'][0]
		m22 = teff_corr[logical_and(teff==t2,logg==g2)]['teff_corr'][0]
		teff_3dcorr = (1-t)*(1-g)*m11+t*(1-g)*m21+t*g*m22+(1-t)*g*m12
		
		m11 = logg_corr[logical_and(teff==t1,logg==g1)]['logg_corr'][0]
                m12 = logg_corr[logical_and(teff==t1,logg==g2)]['logg_corr'][0]
                m21 = logg_corr[logical_and(teff==t2,logg==g1)]['logg_corr'][0]
                m22 = logg_corr[logical_and(teff==t2,logg==g2)]['logg_corr'][0]
                logg_3dcorr = (1-t)*(1-g)*m11+t*(1-g)*m21+t*g*m22+(1-t)*g*m12
		
		if graph==True:	
			fig = figure(figsize=(7,2.5))
			fig.subplots_adjust(left=0.10,bottom=0.15,right=0.98,top=0.98,wspace=0.35)
			ax1 = fig.add_subplot(121)
			plot(teff[where(logg==7.0)],teff_corr['teff_corr'][where(logg==7.0)],'b-',label='logg = 7.0',lw=0.5)
			plot(teff[where(logg==7.5)],teff_corr['teff_corr'][where(logg==7.5)],'r-',label='logg = 7.5',lw=0.5)
			plot(teff[where(logg==8.0)],teff_corr['teff_corr'][where(logg==8.0)],'g-',label='logg = 8.0',lw=0.5)
			plot(teff[where(logg==8.5)],teff_corr['teff_corr'][where(logg==8.5)],'c-',label='logg = 8.5',lw=0.5)
			plot(teff[where(logg==9.0)],teff_corr['teff_corr'][where(logg==9.0)],'m-',label='logg = 9.0',lw=0.5)
			legend(numpoints=1,ncol=1,loc=3, fontsize=7)
			plot(temperature,teff_3dcorr,'ro',ms=2)
			xlabel('Teff')
			ylabel('Teff corr')
			xlim([6000,14500])
			ax2 = fig.add_subplot(122)
			plot(teff[where(logg==7.0)],logg_corr['logg_corr'][where(logg==7.0)],'b-',lw=0.5)
			plot(teff[where(logg==7.5)],logg_corr['logg_corr'][where(logg==7.5)],'r-',lw=0.5)
			plot(teff[where(logg==8.0)],logg_corr['logg_corr'][where(logg==8.0)],'g-',lw=0.5)
			plot(teff[where(logg==8.5)],logg_corr['logg_corr'][where(logg==8.5)],'c-',lw=0.5)
			plot(teff[where(logg==9.0)],logg_corr['logg_corr'][where(logg==9.0)],'m-',lw=0.5)
			plot(temperature,logg_3dcorr,'ro',ms=2)
			xlabel('Teff')
			ylabel('logg corr')
			xlim([6000,14500])
			show()
		return _np.round(teff_3dcorr,0),_np.round(logg_3dcorr,2)
	elif ml2a==0.8:
	   print "to be implemented"
	elif ml2a==0.7:
           print "to be implemented"




def fit_line(spectra, model_in=None, quick=True, line_lmax=None, line_lmin=None, model='sdss', diagnostic=False):
    """
    Use normalised models - can pass model_in from models_normalised() when processing many spectra
    Input Spectra
    Optional:
        model_in=None   : Given model array
        quick=True      : Use presaved model array. Check is up to date
        line_lmax=None  : Maximum wavelength to fit line profiles too e.g.=6200 ignores Ha
        line_lmin=None  : Minimum wavelength to fit line profiles too
        model='sdss': Which model grid to use: 'sdss' (DA, fine, noIR), 'new' (DA, course, IR, new), 'old' (DA, course, IR, old), 'interp' (DA, fine++, noIR)
        diagnostic=False: Return extra things
    #
    Normalise spectra
    Crop to each line
    Scale models to spectra
    Calc chi2
    return chi2,
    list of arrays of spectra at lines,
    list of arrays of best models at lines
    if diagnostic:
        return normalised spectra, continuum flux, normalised model at best chi2
    """
    from scipy import interpolate
    #min max of lines ha hb hg hd he
    line_crop = _np.loadtxt(basedir+'/WDFitting/line_crop.dat')
    if line_lmin:
        line_crop = line_crop[_np.mean(line_crop,axis=1)>line_lmin]
    if line_lmax:
        line_crop = line_crop[_np.mean(line_crop,axis=1)<line_lmax]
    #
    #Check if lines are inside spectra l
    line_crop = line_crop[(line_crop[:,0]>spectra[:,0].min()) & (line_crop[:,1]<spectra[:,0].max())]
    #load normalised models
    if model_in==None:
        m_wave_n,m_flux_n,model_list,model_param = models_normalised(quick=quick, model=model)
    else:
        m_wave_n,m_flux_n,model_list,model_param = model_in
    #
    #normalise spectra
    spectra_n, cont_flux = norm_spectra(spectra)
    #
    cont_flux = cont_flux[(spectra_n[:,0] >= m_wave_n.min()) & (spectra_n[:,0] <= m_wave_n.max())]
    spectra_n = spectra_n[(spectra_n[:,0] >= m_wave_n.min()) & (spectra_n[:,0] <= m_wave_n.max())]
    #linearly interpolate models onto spectral wavelength grid
    tck_l_m = interpolate.interp1d(m_wave_n,m_flux_n,kind='linear')
    m_flux_n_i = tck_l_m(spectra_n[:,0])
    #
    #Initialise: normalised models and spectra in line region, and chi2
    tmp_lines_m = []
    lines_s = []
    l_chi2 = []
    #for each line
    for c in xrange(len(line_crop)):
        #crop model and spectra to line
        #only if line region is entirely covered by spectrum
        if (line_crop[c,1] < spectra[:,0].max()) & (line_crop[c,0] > spectra[:,0].min()):
            l_m = m_flux_n_i.transpose()[(spectra_n[:,0]>=line_crop[c,0])&(spectra_n[:,0]<=line_crop[c,1])].transpose()
            l_s = spectra_n[(spectra_n[:,0]>=line_crop[c,0])&(spectra_n[:,0]<=line_crop[c,1])]
            #renormalise models to spectra in line region
            l_m = l_m*_np.sum(l_s[:,1])/_np.sum(l_m,axis=1).reshape([len(l_m),1])
            
            #calculate chi2
            #l_chi2[:,c] = _np.sum(((l_s[:,1]-l_m)/chi2errs)**2,axis=1)/(len(l_s)-2)
            if _np.isnan(_np.sum(((l_s[:,1]-l_m)/l_s[:,2])**2,axis=1)[0])==False:
                
                l_chi2.append( _np.sum(((l_s[:,1]-l_m)/l_s[:,2])**2,axis=1) )
            #
                tmp_lines_m.append(l_m)
                lines_s.append(l_s)
    #mean chi2 over lines
    l_chi2 = _np.array(l_chi2)
    lines_chi2 = _np.sum(l_chi2,axis=0)
    #print tmp_lines_m[c][lines_chi2==lines_chi2.min()]
    #store best model lines for output
    lines_m = []
    for c in xrange(len(line_crop)):
        if (line_crop[c,1] < spectra[:,0].max()) & (line_crop[c,0] > spectra[:,0].min()):
            lines_m.append(tmp_lines_m[c][lines_chi2==lines_chi2.min()][0])
    #
    #chi2 contour
    model_shape = [len(_np.unique(model_param[:,0])),len(_np.unique(model_param[:,1]))]
    #print _np.concatenate([model_param[:,0].reshape(model_shape), model_param[:,1].reshape(model_shape), lines_chi2.reshape(model_shape)])
    #chi2_cont = _np.concatenate([model_param[:,0].reshape(model_shape), model_param[:,1].reshape(model_shape), lines_chi2.reshape(model_shape)], axis=2)
    #
    best_TL = model_param[lines_chi2 == lines_chi2.min()][0]

    other=sorted(lines_chi2)[2]
    other_TL = model_param[lines_chi2 == other][0]

    if diagnostic:
        return  lines_s, lines_m, spectra_n, cont_flux, model_param,lines_chi2#m_flux_n_i[lines_chi2==lines_chi2.min()][0]#chi2_cont,
    else:
        return  lines_s, lines_m,best_TL,model_param,lines_chi2,other_TL
