

def simulated_survey(nsim):

    # Sorry, imports are here to make pp work
    import numpy as np
    import pandas as pd

    import astropysics.constants as C

    import cksgaia.io
    import cksgaia.fitting
    import cksgaia.completeness

    npl = nsim
    
    # radii limits to define two peaks
    r1, r2, r3, r4 = 1.0, 1.75, 1.75, 3.5
    
    # Make some fake planets
    r1frac = 0.5
    r2frac = 1-r1frac
    nstars = 5e4

    # broad distributions
    # simulated_periods = np.random.lognormal(3.5, 1.0, size=int(np.round(npl*r1frac)))
    # simulated_periods = np.append(simulated_periods, np.random.lognormal(2.5, 0.9, size=int(np.round(npl*r2frac))))
    # simulated_radii = np.random.lognormal(1.05, 0.32, size=int(np.round(npl*r1frac)))
    # simulated_radii = np.append(simulated_radii, np.random.lognormal(0.0, 0.20, size=int(np.round(npl*r2frac))))

    # narrow distributions
    #simulated_periods = np.random.lognormal(3.5, 0.7, size=int(np.round(npl*r1frac)))
    #simulated_periods = np.append(simulated_periods, np.random.lognormal(2.5, 0.9, size=int(np.round(npl*r2frac))))
    #simulated_radii = np.random.lognormal(0.867, 0.17, size=int(np.round(npl*r1frac)))
    #simulated_radii = np.append(simulated_radii, np.random.lognormal(0.26, 0.12, size=int(np.round(npl*r2frac))))

    # uniform distributions
    simulated_periods = 10**np.random.uniform(np.log(1), np.log(200), size=int(npl))
    simulated_radii = 10**(np.random.uniform(np.log(0.5), np.log(8), size=int(npl)))
    
    
    random_order = np.array(range(nsim))
    np.random.shuffle(random_order)
    simulated_periods = simulated_periods[random_order]
    simulated_radii = simulated_radii[random_order]

    Redges = np.logspace(np.log10(0.5), np.log10(20), 36)
    true_dist,_ = np.histogram(simulated_radii, bins=Redges)
    true_dist = np.array(true_dist, dtype=float) / nstars

    simdf = pd.DataFrame([])
    simdf['koi_period'] = simulated_periods
    simdf['iso_prad'] = simulated_radii
    simdf['koi_period_err1'] = simulated_periods*1e-4
    simdf['iso_prad_err1'] = simulated_radii*0.10


    # Calculate ratio of small to large
    small = simdf.query('iso_prad > %3.1f & iso_prad < %3.1f & koi_period > 1 & koi_period < 100' % (r1,r2))
    large = simdf.query('iso_prad > %3.1f & iso_prad < %3.1f & koi_period > 1 & koi_period < 100' % (r3,r4))

    true_ratio = len(small) / float(len(large))


    # Put them around stars with detected planets in my sample
    kicselect = cksgaia.io.load_table('kic')

    physmerge = cksgaia.io.load_table('fulton17')

    kois = pd.merge(physmerge, kicselect, left_on='id_kic', right_on='KICID')

    random_kic = kicselect.iloc[np.random.choice(kicselect.index, size=nsim, replace=True)]
    simdf['id_kic'] = random_kic['id_kic'].values
    simdf['iso_srad'] = random_kic['kic_srad'].values
    simdf['iso_smass'] = random_kic['kic_smass'].values


    kois = simdf.copy()

    kois['koi_sma'] = (kois['iso_smass']*(kois['koi_period']/365.)**2)**(1/3.)
    a = (C.G*kois['iso_srad']*C.Ms*((kois['koi_period']*(24*3600.))/(2*np.pi))**2)**(1/3.) * C.aupercm
    R = kois['iso_srad']*C.Rs * C.aupercm
    kois['koi_duration'] = (kois['koi_period']*24./np.pi)*np.arcsin(R/a)
    kois['koi_ror'] = (kois['iso_prad'] * (C.Re/C.Rs))/kois['iso_srad']
    kois['koi_dor'] = (kois['koi_sma'] / 0.00465047) / kois['iso_srad']

    
    # Fit CDPP values and extrapolate to transit durations
    kicselect = cksgaia.completeness.fit_cdpp(kicselect)

    # Add systematic offset to KIC stellar radii
    star_scale = np.random.normal(1.0, 0.05)  # 5% systematic for GAIA
    kicselect['kic_srad'] = kicselect['kic_srad'].values * np.clip(star_scale, 0.2, 5.0)

    # Add noise to KIC stellar radii
    star_noise = np.random.normal(0, kicselect['kic_srad'].values*0.1, size=len(kicselect))  # 10% radii unc. for GAIA
    kicselect['kic_srad'] += star_noise
    
    kois = pd.merge(kois, kicselect, on='id_kic', suffixes=['','_ks'])

    x = 1/np.sqrt(kois['koi_duration'])
    kois['koi_cdpp_dur'] = kois['kic_cdpp_fit0'] + kois['kic_cdpp_fit1'] * x + kois['kic_cdpp_fit2'] * x**2    
    kois['koi_snr'] = kois['koi_ror']**2 * (kois['koi_period']/kois['kic_tobs'])**-0.5 * 1/(kois['koi_cdpp_dur']*1e-6)


    # Calculate weights
    kois = cksgaia.completeness.get_weights(kois, kicselect)


    # Select detections
    obsdf = kois.copy()

    detected = np.ones(len(obsdf), dtype=bool)
    r = np.random.uniform(0,1, size=len(obsdf))
    detected = r <= 1/obsdf['weight']
    obsdf = obsdf[detected]
    ndet = obsdf['id_kic'].count()

    # Calculate 2D occurance grid
    Redges = np.logspace(np.log10(0.5), np.log10(20), 36)
    Pedges = np.logspace(np.log10(1.0), np.log10(100), 36)

    Rcen = 10**(np.log10(Redges[:-1]) + np.diff(np.log10(Redges)).mean()/2)
    Pcen = 10**(np.log10(Pedges[:-1]) + np.diff(np.log10(Pedges)).mean()/2)

    detections,_,_ = np.histogram2d(obsdf['koi_period'], obsdf['iso_prad'], bins=[Pedges,Redges])
    wdetections,_,_ = np.histogram2d(obsdf['koi_period'], obsdf['iso_prad'], bins=[Pedges,Redges], weights=obsdf['weight'])

    sfac = 1/obsdf['tr_prob'].mean()
    rhist = np.sum(detections, axis=0) 
    rhistn = rhist / nstars * sfac
    rerr = np.sqrt(rhist) / nstars * sfac

    whist = np.sum(wdetections, axis=0) 
    whistn = whist / nstars
    werr = rerr * (whistn/rhistn)

    
    # Calculate ratio of small to large
    small = obsdf.query('iso_prad > %3.1f & iso_prad < %3.1f & koi_period > 1 & koi_period < 100' % (r1,r2))['weight'].sum()
    large = obsdf.query('iso_prad > %3.1f & iso_prad < %3.1f & koi_period > 1 & koi_period < 100' % (r3,r4))['weight'].sum()

    final_ratio = small / large

    kx, ky = cksgaia.fitting.wkde(obsdf['iso_prad'].values, obsdf['iso_prad_err1'].values, 1+0*obsdf['weight'].values)
    wkx, wky = cksgaia.fitting.wkde(obsdf['iso_prad'].values, obsdf['iso_prad_err1'].values, obsdf['weight'].values)
    tkx, tky = cksgaia.fitting.wkde(simdf['iso_prad'].values, simdf['iso_prad_err1'].values,
                                       np.ones_like(simdf['iso_prad'].values))

    output = {}
    output['n_detections'] = ndet
    output['true_ratio'] = true_ratio
    output['final_ratio'] = final_ratio
    output['radius_bin_centers'] = Rcen
    output['observed_hist'] = rhist
    output['observed_hist_norm'] = rhistn
    output['weighted_hist_norm'] = whistn
    output['sim_data'] = obsdf
    output['kde_x'] = kx
    output['kde_y'] = ky
    output['wkde_x'] = wkx
    output['wkde_y'] = wky
    output['true_kde_x'] = tkx
    output['true_kde_y'] = tky

    
    return output

if __name__ == '__main__':
    out = simulated_survey(1000)
