#!/usr/bin/env python3
"""Generate the declarative "params" of star_psf / objects / icl / lsbg / mask / ai_services (replaces hand-written dialogs).
Rows: (name, label, type, default, group, help, extra-dict)"""
import os as _os
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.environ.get('OGF_ROOT',_os.path.dirname(_HERE))
import json,sys
R=_ROOT+'/plugins/'
def P(name,label,typ,default,group,help,**kw):
    d=dict(name=name,label=label,type=typ,default=default,group=group,help=help); d.update(kw); return d
F='float'; I='int'; S='string'; B='bool'; C='choice'
interp=['linear','cubic','nearest']
spec={}
spec['star_psf']=dict(store=dict(array='catpanel',key='psf,param,%s',save='CatalogPanelPSFParamSave'),params=[
 P('class-star-thresh','CLASS_STAR threshold',F,0.8,'Star finding','Minimum SExtractor CLASS_STAR for a PSF star candidate.',min=0,max=1),
 P('max-ellipticity','Max ellipticity',F,0.2,'Star finding','Stars rounder than this are kept.',min=0,max=1),
 P('fwhm-sigma','FWHM sigma',F,2.0,'Star finding','FWHM clipping in units of the robust scatter of the stellar locus.',min=0.1,max=20),
 P('min-flux-snr','Min flux S/N',F,10.0,'Star finding','Minimum S/N of a PSF star.',min=0,max=10000),
 P('psf-size','PSF size',I,51,'PSF','Stamp size of the empirical PSF.',min=5,max=2001,unit='pix'),
 P('ext-core-mag-min','Core mag min',F,18.0,'Extended PSF','Brightest magnitude of the core stars.',min=-10,max=40),
 P('ext-core-mag-max','Core mag max',F,22.0,'Extended PSF','Faintest magnitude of the core stars.',min=-10,max=40),
 P('ext-wing-mag-max','Wing mag max',F,16.0,'Extended PSF','Faintest magnitude of the bright stars that define the wings.',min=-10,max=40),
 P('ext-core-size','Core cutout size',I,51,'Extended PSF','Cutout size of the core PSF.',min=5,max=2001,unit='pix'),
 P('ext-wing-size','Wing cutout size',I,201,'Extended PSF','Cutout size of the wing PSF.',min=5,max=4001,unit='pix'),
 P('ext-blend-inner','Blend inner',F,20.0,'Extended PSF','Radius where the core-to-wing blend starts.',min=0,max=2000,unit='pix'),
 P('ext-blend-outer','Blend outer',F,30.0,'Extended PSF','Radius where the blend ends.',min=0,max=2000,unit='pix'),
 P('ext-saturation-limit','Saturation',F,60000.0,'Extended PSF','Pixels above this level are treated as saturated.',min=0,max=1e9,unit='ADU'),
 P('sim-telescope','Telescope',C,'auto','Simulation','Telescope for the simulated PSF (auto = from the header).',choices=['auto','jwst','hst']),
 P('sim-instrument','Instrument',S,'auto','Simulation','Instrument name, or auto.'),
 P('sim-filter','Filter',S,'auto','Simulation','Filter name, or auto.'),
 P('sim-psf-size','PSF size',I,201,'Simulation','Size of the simulated PSF.',min=5,max=4001,unit='pix'),
 P('sim-oversample','Oversample',I,1,'Simulation','Oversampling factor.',min=1,max=16,expert=1),
 P('sim-jitter-sigma','Jitter sigma',F,0.007,'Simulation','Pointing jitter (WebbPSF).',min=0,max=5,unit='arcsec',expert=1),
 P('sim-focus-offset','Focus offset',F,0.0,'Simulation','Focus offset (TinyTim).',min=-100,max=100,expert=1),
])
spec['objects']=dict(store=dict(array='catpanel',key='param,%s',save='CatalogPanelParamSave'),params=[
 P('sep-deblend-nthresh','Sub-thresholds',I,64,'Deblending','Number of deblending sub-thresholds of Separate Selected.',min=1,max=1024),
 P('sep-deblend-mincont','Min contrast',F,0.0001,'Deblending','Minimum deblending contrast.',min=0,max=1),
 P('sep-detect-thresh','Threshold',F,0.8,'Detection','Detection threshold of the cutout re-extraction.',min=0.01,max=100,unit='sigma'),
 P('sep-detect-minarea','Min area',I,3,'Detection','Minimum number of connected pixels.',min=1,max=100000,unit='pix'),
 P('sep-radius-factor','Radius factor',F,3.0,'Cutout','Cutout half-size in units of the source radius.',min=0.5,max=100),
 P('sep-back-size','Background size',I,32,'Cutout','Background mesh size of the cutout.',min=2,max=4096,unit='pix'),
])
icl_p=[
 ('expand-factor','Expand factor',F,1.5,'Masking','Source-mask expansion factor.',dict(min=0,max=50)),
 ('max-dilate-radius','Max dilate radius',I,20,'Masking','Largest dilation radius of the mask.',dict(min=0,max=10000,unit='pix')),
 ('bright-star-mag-limit','Bright star mag limit',F,18.0,'Masking','Stars brighter than this get an extended mask.',dict(min=-10,max=40)),
 ('bright-star-radius-scale','Bright star radius scale',F,10.0,'Masking','Radius scale of the bright-star masks.',dict(min=0,max=1000)),
 ('interp-method','Interpolation',C,'linear','Masking','Interpolation across masked pixels.',dict(choices=interp)),
 ('detect-thresh','Detect threshold',F,5.0,'Masking','Source detection threshold for the mask.',dict(min=0.01,max=100,unit='sigma')),
 ('bkg-order','Polynomial order',I,3,'Background','Order of the polynomial / Chebyshev background.',dict(min=0,max=20)),
 ('bkg-sigma-clip','Sigma clip',F,3.0,'Background','Sigma clipping for the background fit.',dict(min=0.5,max=20)),
 ('bkg-sep-mesh','SEP mesh size',I,256,'Background','Mesh size of the SEP background.',dict(min=4,max=8192,unit='pix')),
 ('bkg-iterative','Iterative refinement',B,0,'Background','Re-fit the background excluding sources found on the previous pass.',{}),
 ('bkg-n-iterations','N iterations',I,3,'Background','Iterations of the refinement.',dict(min=1,max=50)),
 ('bkg-convergence-tol','Convergence tol',F,0.01,'Background','Stop when the background changes less than this.',dict(min=0,max=1)),
 ('bkg-refine-thresh','Refine threshold',F,2.0,'Background','Threshold for masking sources during refinement.',dict(min=0.1,max=50,unit='sigma')),
 ('rmin','R min',F,5.0,'Profile','Innermost radius of the radial profile.',dict(min=0,max=1e6,unit='pix')),
 ('rmax','R max',F,1000.0,'Profile','Outermost radius of the radial profile.',dict(min=1,max=1e6,unit='pix')),
 ('nsteps','N steps',I,80,'Profile','Number of radial bins.',dict(min=2,max=10000)),
 ('spacing','Spacing',C,'log','Profile','Radial bin spacing.',dict(choices=['log','linear'])),
 ('ellipticity','Ellipticity',F,0.0,'Profile','Ellipticity of the elliptical annuli.',dict(min=0,max=0.99)),
 ('pa','PA',F,0.0,'Profile','Position angle of the annuli.',dict(min=-360,max=360,unit='deg')),
 ('mag-zeropoint','Mag zeropoint',F,25.0,'Calibration','Photometric zeropoint.',dict(min=-100,max=100)),
 ('pixel-scale','Pixel scale',F,0.06,'Calibration','Pixel scale.',dict(min=1e-6,max=1000,unit='arcsec/px')),
 ('mu-threshold','ICL mu threshold',F,26.5,'Calibration','Surface-brightness threshold that defines the ICL.',dict(min=0,max=40,unit='mag/arcsec2')),
 ('mu-levels','Isophotal mu levels',S,'26.0,27.0,28.0','Calibration','Comma separated surface-brightness levels of the multi-threshold measurement.',{}),
 ('measure-radius','Measure radius',F,500.0,'Calibration','Radius of the ICL measurement.',dict(min=1,max=1e6,unit='pix')),
]
spec['icl']=dict(store=dict(array='catpanel',key='icl,param,%s',save='CatalogPanelICLParamSave'),params=[P(n,l,t,d,g,h,**kw) for n,l,t,d,g,h,kw in icl_p])
lsbg_p=[
 ('mask-mag-threshold','Mask mag threshold',F,22.0,'Masking','Sources brighter than this are masked.',dict(min=-10,max=99)),
 ('mask-expand-factor','Expand factor',F,1.5,'Masking','Mask expansion factor.',dict(min=0,max=50)),
 ('max-dilate-radius','Max dilate radius',I,30,'Masking','Largest dilation radius.',dict(min=0,max=10000,unit='pix')),
 ('bright-star-mag-limit','Bright star mag limit',F,18.0,'Masking','Stars brighter than this get an extended mask.',dict(min=-10,max=40)),
 ('bright-star-radius-scale','Bright star radius scale',F,12.0,'Masking','Radius scale of the bright-star masks.',dict(min=0,max=1000)),
 ('interp-method','Interpolation',C,'linear','Masking','Interpolation across masked pixels.',dict(choices=interp)),
 ('mask-detect-thresh','Mask detect threshold',F,1.5,'Masking','Detection threshold for the mask.',dict(min=0.01,max=100,unit='sigma')),
 ('mask-detect-minarea','Mask detect min area',I,5,'Masking','Minimum area for the mask detection.',dict(min=1,max=100000,unit='pix',expert=1)),
 ('lsb-protect','LSB structure protection',B,1,'Masking','Do not mask extended low-surface-brightness structures.',{}),
 ('lsb-mu-threshold','LSB mu threshold',F,24.0,'Masking','Surface brightness above which structures are protected.',dict(min=0,max=40,unit='mag/arcsec2')),
 ('bkg-method','Method',C,'sep_large','Background','Background model.',dict(choices=['sep_large','polynomial','chebyshev'])),
 ('bkg-mesh-size','Mesh size',I,256,'Background','Mesh size of the SEP background.',dict(min=4,max=8192,unit='pix')),
 ('bkg-poly-order','Polynomial order',I,3,'Background','Order of the polynomial / Chebyshev background.',dict(min=0,max=20)),
 ('bkg-sigma-clip','Sigma clip',F,3.0,'Background','Sigma clipping.',dict(min=0.5,max=20)),
 ('bkg-n-iterations','Iterations',I,3,'Background','Refinement iterations.',dict(min=1,max=50)),
 ('bkg-refine-thresh','Refine threshold',F,2.0,'Background','Threshold for masking sources during refinement.',dict(min=0.1,max=50,unit='sigma')),
 ('bkg-rms-quantile','RMS quantile',F,0.25,'Background','Quantile used for the noise estimate.',dict(min=0.01,max=1,expert=1)),
 ('bkg-convergence-tol','Convergence tolerance',F,0.01,'Background','Stop when the background changes less than this.',dict(min=0,max=1,expert=1)),
 ('detect-thresh','Detect threshold',F,0.8,'Detection','Detection threshold.',dict(min=0.01,max=100,unit='sigma')),
 ('detect-minarea','Min area',I,50,'Detection','Minimum area of a candidate.',dict(min=1,max=1000000,unit='pix')),
 ('detect-filter-kernel','Filter kernel',C,'gauss5x5','Detection','Convolution kernel.',dict(choices=['none','gauss3x3','gauss5x5','gauss7x7','gauss9x9','tophat5','tophat7','mexhat'])),
 ('deblend-nthresh','Deblend N thresholds',I,32,'Detection','Deblending sub-thresholds.',dict(min=1,max=1024)),
 ('deblend-mincont','Deblend min contrast',F,0.005,'Detection','Deblending minimum contrast.',dict(min=0,max=1)),
 ('multiscale','Multi-scale detection',B,1,'Detection','Also detect on binned images.',{}),
 ('multiscale-factors','Scale factors',S,'1,2,4','Detection','Comma separated binning factors.',{}),
 ('sersic-fit','Sersic profile fitting',B,1,'Detection','Fit a Sersic profile to each candidate.',{}),
 ('sersic-n-min','Sersic n min (fit)',F,0.2,'Sersic fit','Lower bound of the fitted Sersic index.',dict(min=0.1,max=20,expert=1)),
 ('sersic-n-max','Sersic n max (fit)',F,10.0,'Sersic fit','Upper bound of the fitted Sersic index.',dict(min=0.1,max=20,expert=1)),
 ('sersic-re-min','Re min (fit)',F,0.5,'Sersic fit','Lower bound of the fitted effective radius.',dict(min=0,max=1e4,expert=1)),
 ('sersic-cutout-scale','Cutout scale',F,5.0,'Sersic fit','Cutout size in units of the source radius.',dict(min=1,max=100,expert=1)),
 ('sersic-max-nfev','Max function evaluations',I,500,'Sersic fit','Limit of the least-squares fit.',dict(min=10,max=1000000,expert=1)),
 ('mag-zeropoint','Mag zeropoint',F,25.0,'Calibration & filter','Photometric zeropoint.',dict(min=-100,max=100)),
 ('pixel-scale','Pixel scale',F,0.06,'Calibration & filter','Pixel scale.',dict(min=1e-6,max=1000,unit='arcsec/px')),
 ('phot-apertures','Phot apertures',S,'5,10,20,40','Calibration & filter','Comma separated aperture diameters in pixels.',{}),
 ('mu-eff-min','mu_eff min',F,24.0,'Calibration & filter','Faintest-end filter: minimum effective surface brightness.',dict(min=0,max=40,unit='mag/arcsec2')),
 ('mu-eff-max','mu_eff max',F,30.0,'Calibration & filter','Maximum effective surface brightness.',dict(min=0,max=40,unit='mag/arcsec2')),
 ('r-eff-min','R_eff min',F,2.5,'Calibration & filter','Minimum effective radius.',dict(min=0,max=1e4,unit='arcsec')),
 ('r-eff-max','R_eff max',F,60.0,'Calibration & filter','Maximum effective radius.',dict(min=0,max=1e4,unit='arcsec')),
 ('ellipticity-max','Ellipticity max',F,0.7,'Calibration & filter','Maximum ellipticity.',dict(min=0,max=1)),
 ('min-snr','Min SNR',F,2.0,'Calibration & filter','Minimum S/N.',dict(min=0,max=1e4)),
 ('sersic-n-filter-min','Sersic n min (filter)',F,0.3,'Calibration & filter','Sersic index lower limit of the filter.',dict(min=0,max=20)),
 ('sersic-n-filter-max','Sersic n max (filter)',F,6.0,'Calibration & filter','Sersic index upper limit of the filter.',dict(min=0,max=20)),
 ('sersic-chi2-max','Sersic chi2 max',F,10.0,'Calibration & filter','Maximum reduced chi2 of the Sersic fit.',dict(min=0,max=1e6)),
 ('svm-classify','SVM classification',B,0,'SVM','Classify candidates with the SVM model.',{}),
 ('svm-threshold','SVM threshold',F,0.3,'SVM','Probability threshold.',dict(min=0,max=1)),
 ('svm-checkpoint','SVM checkpoint',S,'','SVM','Model checkpoint file (empty = built-in).',{}),
]
spec['lsbg']=dict(store=dict(array='catpanel',key='lsbg,param,%s',save='CatalogPanelLSBGParamSave'),params=[P(n,l,t,d,g,h,**kw) for n,l,t,d,g,h,kw in lsbg_p])
spec['mask']=dict(store=dict(array='ogfmask',key='%s'),on_apply='OGFMaskOverlayApplied',params=[
 P('color','Colour',S,'red','Overlay','Colour name or #rrggbb of the mask overlay.'),
 P('transparency','Transparency',I,50,'Overlay','Overlay transparency in percent (0 = opaque, 100 = invisible).',min=0,max=100,unit='%'),
])
spec['ai_services']=dict(store=dict(array='ogfai',key='%s',save='OGFAIPrefSave'),on_open='OGFAIInit',on_apply='OGFAIRegRefresh',params=[
 P('sf','Service profile file',  'file','','Services','JSON profile that lists the external services (see docs/ai_services.md). Empty = built-in MOCK only.'),
 P('backend,photoz','Photo-z backend',C,'local','Backend of the local steps','local = the built-in model (unchanged behaviour); otherwise an enabled service of task photoz.',choices_proc='OGFAIBackendChoices photoz'),
 P('backend,sed_fit','SED fitting backend',C,'local','Backend of the local steps','local or an enabled sed_fit service.',choices_proc='OGFAIBackendChoices sed_fit'),
 P('backend,morphology','Galaxy morphology backend',C,'local','Backend of the local steps','local or an enabled morphology service.',choices_proc='OGFAIBackendChoices morphology'),
 P('backend,star','AI star classification backend',C,'local','Backend of the local steps','local or an enabled star_galaxy service.',choices_proc='OGFAIBackendChoices star'),
])
for pid,sp in spec.items():
    f=R+pid+'/plugin.json'; m=json.load(open(f))
    m['params']=sp['params']
    m['store']=sp['store']
    for k in('on_apply','on_open'):
        if k in sp: m[k]=sp[k]
    m['settings']='params'
    for s in m['steps']:
        if s.get('proc') in ('CatalogPanelStarPSFSettings','CatalogPanelICLSettings','CatalogPanelLSBGSettings','CatalogPanelMaskOverlaySettings'):
            s['proc']='OGFParamDialog %s'%pid
        if s.get('settings')=='CatalogPanelSeparateSettings': s['settings']='params'
        if s.get('proc')=='OGFAIRegistry': pass
    if pid=='mask':
        for s in m['steps']:
            if s['id']=='overlay_style': s['settings_only']=1
    open(f,'w').write(json.dumps(m,indent=1,ensure_ascii=False)+'\n')
    print(pid,len(sp['params']),'params')
