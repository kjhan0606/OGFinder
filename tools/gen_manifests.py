#!/usr/bin/env python3
"""One-off generator for plugins/*/plugin.json of the built-in (migrated-by-wrapper) plugins.
The generated JSON files are the source of truth afterwards (hand-edited)."""
import json, os, sys
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.environ.get('OGF_ROOT', _os.path.dirname(_HERE))
OUT = _ROOT + '/plugins'

def S(id, label, proc=None, legacy=None, session='NONE', stage=None, records=None, needs=None, **kw):
    d = {'id': id, 'label': label}
    if proc: d['proc'] = proc
    if legacy: d['legacy'] = legacy if isinstance(legacy, list) else [legacy]
    d['session'] = session
    if stage: d['stage'] = stage
    if records: d['records'] = records
    if needs: d['needs'] = needs
    d.update(kw)
    return d

def V(id, label, proc, legacy):
    return {'id': id, 'label': label, 'proc': proc, 'legacy': [legacy]}

def P(name, label, type, default, group='General', help='', **kw):
    d = {'name': name, 'label': label, 'type': type, 'default': default, 'group': group, 'help': help}
    d.update(kw)
    return d

plugins = []

# ------------------------------------------------------------------ Detect
plugins.append(dict(id='extract', name='Source extraction', tab='Detect', order=10,
  description='SExtractor-compatible detection (ds9_sextract), dual-image mode, catalog trimming.',
  requires={'binaries': ['ds9_sextract']},
  store={'array': 'catpanel', 'key': 'param,%s', 'save': 'CatalogPanelParamSave'},
  steps=[
    S('extract', 'Extract', 'CatalogPanelExtract', 'SExtractor > Extract', 'AUTO', 'detect', ['extract'], ['image'], settings='params'),
    S('dual', 'Dual-Image Extract...', 'CatalogPanelDualExtract', 'SExtractor > Dual-Image Extract...', 'AUTO', 'detect', ['analysis.dual_extract'], ['image']),
    S('trim', 'Trim catalog...', 'CatalogPanelTrimDialog', 'SExtractor > Trim...', 'AUTO', None, ['catalog.trim'], ['catalog']),
  ],
  params=[
    P('detect-thresh', 'Detect Threshold', 'float', 1.5, 'Detection', 'Detection threshold in sigma above the local background.', min=0.1, max=100, unit='sigma'),
    P('detect-minarea', 'Detect Min Area', 'int', 5, 'Detection', 'Minimum number of connected pixels above threshold.', min=1, max=100000, unit='pix'),
    P('deblend-nthresh', 'Deblend NThresh', 'int', 32, 'Detection', 'Number of deblending sub-thresholds.', min=1, max=64),
    P('deblend-mincont', 'Deblend MinCont', 'float', 0.005, 'Detection', 'Minimum contrast parameter for deblending.', min=0, max=1),
    P('conv-filter', 'Conv Filter', 'choice', 'default', 'Detection', 'Convolution filter applied before detection.', choices=['default', 'gauss5x5', 'mexhat', 'tophat']),
    P('phot-aperture', 'Phot Aperture (diam)', 'float', 5.0, 'Photometry', 'Diameter of the first fixed aperture.', min=0.1, max=1000, unit='pix'),
    P('phot-aperture-2', 'Phot Aperture 2', 'float', 4.0, 'Photometry', 'Second aperture diameter.', min=0.1, max=1000, unit='pix', expert=1),
    P('phot-aperture-3', 'Phot Aperture 3', 'float', 6.0, 'Photometry', 'Third aperture diameter.', min=0.1, max=1000, unit='pix', expert=1),
    P('phot-aperture-5', 'Phot Aperture 5', 'float', 10.0, 'Photometry', 'Fifth aperture diameter.', min=0.1, max=1000, unit='pix', expert=1),
    P('mag-zeropoint', 'Mag Zeropoint', 'float', 25.0, 'Image', 'Magnitude zeropoint of the image.'),
    P('gain', 'Gain', 'float', 0.0, 'Image', 'Detector gain (0 = ignore).', min=0),
    P('pixel-scale', 'Pixel Scale', 'float', 1.0, 'Image', 'Pixel scale in arcsec/pixel.', min=0.0001, max=100, unit='arcsec/pix'),
    P('seeing-fwhm', 'Seeing FWHM', 'float', 3.0, 'Image', 'Seeing FWHM in pixels (used for CLASS_STAR).', min=0.1, max=100, unit='pix'),
    P('back-size', 'Back Size', 'int', 64, 'Background', 'Background mesh size.', min=8, max=4096, unit='pix', expert=1),
    P('back-filtersize', 'Back Filter Size', 'int', 3, 'Background', 'Median filter size on the background mesh.', min=1, max=15, expert=1),
    P('n-workers', 'Parallel Workers (0=auto)', 'int', 0, 'Performance', 'Worker threads for the Python steps (0 = automatic).', min=0, max=256, expert=1),
  ]))

plugins.append(dict(id='objects', name='Objects / catalog editing', tab='Detect', order=20,
  description='Markers on the image, visible-only filter, interactive add/delete/separate, AI merge.',
  steps=[
    S('mark_all', 'Mark All', 'CatalogPanelMarkAll', 'Objects > Mark All', 'NONE', None, None, ['catalog']),
    S('clear_markers', 'Clear Markers', 'CatalogPanelClearMarkers', 'Objects > Clear Markers'),
    S('visible_only', 'Show Visible Only', 'CatalogPanelShowVisible', 'Objects > Show Visible Only', 'MANUAL', None, ['catalog.unrecorded'], ['catalog'], toggle='catpanel(visible_mode)'),
    S('add_objects', 'Add Objects (Click + A)', None, 'Objects > Add Objects (Click + A)', 'MANUAL', None, ['catalog.add_object'], ['catalog'], toggle='catpanel(add_objects_mode)', toggle_only=1),
    S('delete_selected', 'Delete Selected (Click + D)', 'CatalogPanelDeleteSelected', 'Objects > Delete Selected (Click + D)', 'MANUAL', None, ['catalog.delete'], ['catalog']),
    S('separate_selected', 'Separate Selected (Click + S)', 'CatalogPanelSeparateSelected', 'Objects > Separate Selected (Click + S)', 'MANUAL', None, ['catalog.separate'], ['catalog'],
      settings='CatalogPanelSeparateSettings', legacy_settings=['Objects > Separate Settings...']),
    S('save_separated', 'Save Separated Catalog...', 'CatalogPanelSeparateSave', 'Objects > Save Separated Catalog...'),
    S('load_separated', 'Load Separated Catalog...', 'CatalogPanelSeparateLoad', 'Objects > Load Separated Catalog...', 'MANUAL'),
    S('ai_merge', 'AI Merge...', 'CatalogPanelAIMerge', 'Objects > AI Merge...', 'MANUAL', None, ['catalog.merge'], ['catalog']),
  ]))

plugins.append(dict(id='bands', name='Multi-band', tab='Detect', order=30,
  description='Band registry, detection band, forced photometry in all bands, tiled view.',
  steps=[
    S('register', 'Register Current Frame as Band...', 'CatalogPanelBandsRegister', 'Bands > Register Current Frame as Band...', 'CONFIG', None, ['bands.register'], ['image']),
    S('load', 'Load Band...', 'CatalogPanelBandsLoad', 'Bands > Load Band...', 'CONFIG', None, ['bands.register']),
    S('set_detect', 'Detection Band...', 'OGFUIBandMenu detect', 'Bands > Detection Band', 'CONFIG', None, ['bands.detect']),
    S('remove', 'Remove Band...', 'OGFUIBandMenu remove', 'Bands > Remove Band', 'CONFIG', None, ['bands.remove']),
    S('list', 'List Bands...', 'CatalogPanelBandsList', 'Bands > List Bands...'),
    S('detect', 'Detect in Detection Band', 'CatalogPanelBandsDetect', 'Bands > Detect in Detection Band', 'AUTO', 'detect', ['extract'], ['image']),
    S('measure', 'Measure in All Bands...', 'CatalogPanelBandsMeasure', 'Bands > Measure in All Bands...', 'AUTO', 'detect', ['bands.measure'], ['catalog']),
    S('tile', 'Tile Bands', 'CatalogPanelBandsTile', 'Bands > Tile Bands', 'NONE', None, None, None),
    S('single', 'Single Frame View', 'CatalogPanelBandsSingleView', 'Bands > Single Frame View', 'NONE', None, None, None),
    S('copy_mask', 'Copy Mask to Other Bands', 'CatalogPanelMaskCopyToBands', 'Bands > Copy Mask to Other Bands', 'AUTO', None, ['mask.copy_to_bands']),
  ]))

# ------------------------------------------------------------------ Classify
plugins.append(dict(id='star_psf', name='Stars and PSF', tab='Classify', order=10,
  description='PSF star finding, AI star/galaxy classification, PSF building (empirical, extended, simulated), PSF I/O.',
  steps=[
    S('find_stars', 'Find Stars', None, None, 'AUTO', 'classify', ['psf.find_stars'], ['catalog'], variants=[
        V('combined', 'Combined', 'CatalogPanelFindStars combined', 'Star(PSF) > Find Stars > Combined'),
        V('class_star', 'CLASS_STAR', 'CatalogPanelFindStars class_star', 'Star(PSF) > Find Stars > CLASS_STAR'),
        V('fwhm', 'FWHM', 'CatalogPanelFindStars fwhm', 'Star(PSF) > Find Stars > FWHM')]),
    S('ai_star', 'AI Star Classification', 'CatalogPanelStarFinder', 'Star(PSF) > AI Star Classification', 'AUTO', 'classify', ['stars.classify', 'ai.run'], ['catalog']),
    S('show_stars', 'Show Stars', 'CatalogPanelShowStars', 'Star(PSF) > Show Stars'),
    S('clear_stars', 'Clear Stars', 'CatalogPanelClearStars', 'Star(PSF) > Clear Stars'),
    S('build_psf', 'Build PSF', None, None, 'AUTO', 'classify', ['psf.build'], ['catalog'], variants=[
        V('median', 'Median Stack', 'CatalogPanelBuildPSF median', 'Star(PSF) > Build PSF > Median Stack'),
        V('moffat', 'Moffat Fit', 'CatalogPanelBuildPSF moffat', 'Star(PSF) > Build PSF > Moffat Fit'),
        V('gaussian', 'Gaussian Fit', 'CatalogPanelBuildPSF gaussian', 'Star(PSF) > Build PSF > Gaussian Fit'),
        V('epsf', 'ePSF', 'CatalogPanelBuildPSF epsf', 'Star(PSF) > Build PSF > ePSF')]),
    S('ext_psf', 'Extended PSF...', 'CatalogPanelBuildExtendedPSF', 'Star(PSF) > Build PSF > Extended PSF...', 'AUTO', 'classify', ['psf.build_extended'], ['catalog']),
    S('webbpsf', 'WebbPSF (JWST)...', 'CatalogPanelSimPSFWebbPSF', 'Star(PSF) > Build PSF > WebbPSF (JWST)...', 'AUTO', 'classify', ['psf.webbpsf']),
    S('tinytim', 'TinyTim (HST)...', 'CatalogPanelSimPSFTinyTim', 'Star(PSF) > Build PSF > TinyTim (HST)...', 'AUTO', 'classify', ['psf.tinytim']),
    S('view_psf', 'View PSF', 'CatalogPanelViewPSF', 'Star(PSF) > View PSF'),
    S('save_psf', 'Save PSF...', 'CatalogPanelSavePSF', 'Star(PSF) > Save PSF...'),
    S('load_psf', 'Load PSF...', 'CatalogPanelLoadPSF', 'Star(PSF) > Load PSF...'),
    S('psf_settings', 'PSF / Star Settings...', 'CatalogPanelStarPSFSettings', 'Star(PSF) > Settings...', settings_only=1),
  ]))

plugins.append(dict(id='galaxy_model', name='Galaxy classification', tab='Classify', order=20,
  description='AI morphology classification and galaxy model fits.',
  steps=[
    S('ai_morph', 'AI Morphology Classification', 'CatalogPanelGalaxyMorphology', 'Galaxy Model > AI Morphology Classification', 'AUTO', 'classify', ['galaxy.morphology', 'ai.run'], ['catalog']),
    S('fit', 'Fit Galaxy Model', None, None, 'NONE', None, None, ['catalog'], variants=[
        V('elliptical', 'Elliptical', 'CatalogPanelGalaxyFit elliptical', 'Galaxy Model > Fit Elliptical Model'),
        V('spiral', 'Spiral', 'CatalogPanelGalaxyFit spiral', 'Galaxy Model > Fit Spiral Model')]),
    S('params', 'Extract Parameters', 'CatalogPanelGalaxyParams', 'Galaxy Model > Extract Parameters', 'NONE', None, None, ['catalog']),
  ]))

# ------------------------------------------------------------------ Measure
plugins.append(dict(id='morphology', name='Structure', tab='Measure', order=10,
  description='Sersic fitting, non-parametric morphology, bulge+disk decomposition.',
  store={'array': 'catpanel', 'key': 'bd,param,%s', 'save': 'CatalogPanelBDParamSave'},
  steps=[
    S('sersic', 'Sersic Fitting', 'CatalogPanelSersicFit', 'Analysis > Sérsic Fitting', 'AUTO', 'measure', ['analysis.sersic'], ['catalog', 'image']),
    S('morphometry', 'Non-Parametric Morphology (CAS/Gini/M20)', 'CatalogPanelMorphometry', 'Analysis > Non-Parametric Morphology (CAS/Gini/M20)', 'AUTO', 'measure', ['analysis.morphometry'], ['catalog', 'image']),
    S('bulge_disk', 'Bulge+Disk Decomposition', 'CatalogPanelBulgeDisk', 'Analysis > Bulge+Disk Decomp...', 'AUTO', 'measure', ['analysis.bulge_disk'], ['catalog', 'image'], settings='params'),
  ],
  params=[
    P('max-sources', 'Max sources', 'int', 100, 'Bulge+Disk', 'Fit at most this many sources (brightest first).', min=1, max=100000),
    P('free-bulge-n', 'Free bulge Sersic n', 'bool', 0, 'Bulge+Disk', 'Let the bulge Sersic index vary instead of fixing n=4.'),
    P('mag-zeropoint', 'Mag zeropoint', 'float', 25.0, 'Bulge+Disk', 'Magnitude zeropoint of the image.'),
    P('pixel-scale', 'Pixel scale', 'float', 0.263, 'Bulge+Disk', 'Pixel scale in arcsec/pixel.', min=0.0001, max=100, unit='arcsec/pix'),
  ]))

plugins.append(dict(id='photometry', name='Photometry', tab='Measure', order=20,
  description='PSF, multi-band and crowded-field photometry, cross-match, segmentation map, completeness.',
  steps=[
    S('psf_phot', 'PSF Photometry', 'CatalogPanelPSFPhotometry', 'Analysis > PSF Photometry', 'AUTO', 'measure', ['analysis.psf_phot'], ['catalog', 'image']),
    S('multiband', 'Multi-Band Photometry...', 'CatalogPanelMultiBand', 'Analysis > Multi-Band Photometry...', 'AUTO', 'measure', ['analysis.multiband'], ['catalog', 'image']),
    S('crowded', 'Crowded Field Photometry', 'CatalogPanelCrowdedPhot', 'Analysis > Crowded Field Photometry', 'AUTO', 'measure', ['analysis.crowded_phot'], ['catalog', 'image']),
    S('crossmatch', 'Cross-Match (VizieR)...', 'CatalogPanelCrossMatch', 'Analysis > Cross-Match (VizieR)...', 'AUTO', 'measure', ['analysis.crossmatch'], ['catalog'], network=1),
    S('segmap', 'Segmentation Map', 'CatalogPanelSegmentationMap', 'Analysis > Segmentation Map', 'AUTO', None, ['analysis.segmap'], ['image']),
    S('completeness', 'Completeness Simulation...', 'CatalogPanelCompleteness', 'Analysis > Completeness Simulation...', 'AUTO', None, ['analysis.completeness'], ['image']),
  ]))

plugins.append(dict(id='photoz_sed', name='Photo-z and SED', tab='Measure', order=30,
  description='Photometric redshift and SED fitting (local models or external AI services).',
  steps=[
    S('photoz', 'Photo-z (AI)...', 'CatalogPanelPhotoZ', 'Analysis > Photo-z (AI)...', 'AUTO', 'measure', ['analysis.photo_z', 'ai.run'], ['catalog']),
    S('sed', 'SED Fitting (AI)...', 'CatalogPanelSEDFit', 'Analysis > SED Fitting (AI)...', 'AUTO', 'measure', ['analysis.sed_fit', 'ai.run'], ['catalog']),
  ]))

plugins.append(dict(id='deconv', name='Deconvolution', tab='Measure', order=40,
  description='PSF deconvolution of the current image (Richardson-Lucy, Wiener, Tikhonov, CLEAN, MEM).',
  store={'array': 'catpanel', 'key': 'psf,param,%s', 'save': 'CatalogPanelPSFParamSave'},
  steps=[
    S('deconvolve', 'Deconvolve', None, None, 'AUTO', 'measure', ['deconv.run'], ['image'], settings='params', variants=[
        V('rl', 'Richardson-Lucy', 'CatalogPanelDeconvolve rl', 'Deconvolution > Richardson-Lucy'),
        V('rl_accelerated', 'Richardson-Lucy (Accelerated)', 'CatalogPanelDeconvolve rl_accelerated', 'Deconvolution > Richardson-Lucy (Accelerated)'),
        V('rl_tv', 'Richardson-Lucy (Regularized)', 'CatalogPanelDeconvolve rl_tv', 'Deconvolution > Richardson-Lucy (Regularized)'),
        V('wiener', 'Wiener', 'CatalogPanelDeconvolve wiener', 'Deconvolution > Wiener'),
        V('tikhonov', 'Tikhonov', 'CatalogPanelDeconvolve tikhonov', 'Deconvolution > Tikhonov'),
        V('clean', 'CLEAN', 'CatalogPanelDeconvolve clean', 'Deconvolution > CLEAN'),
        V('mem', 'Maximum Entropy (MEM)', 'CatalogPanelDeconvolve mem', 'Deconvolution > Maximum Entropy (MEM)')]),
    S('quick', 'Quick Deconvolve (RL)', 'CatalogPanelQuickDeconvolve', 'Deconvolution > Quick Deconvolve (RL)', 'AUTO', 'measure', ['deconv.run'], ['image']),
  ],
  params=[
    P('rl-iterations', 'Iterations', 'int', 30, 'Richardson-Lucy', 'Number of Richardson-Lucy iterations.', min=1, max=10000),
    P('tv-lambda', 'TV lambda', 'float', 0.001, 'Richardson-Lucy', 'Total-variation regularisation weight (regularized RL).', min=0, max=10),
    P('wiener-nsr', 'NSR', 'float', 0.01, 'Wiener', 'Noise-to-signal ratio of the Wiener filter.', min=0, max=10),
    P('tikhonov-lambda', 'Lambda', 'float', 0.001, 'Tikhonov', 'Tikhonov regularisation weight.', min=0, max=10),
    P('clean-gain', 'Gain', 'float', 0.1, 'CLEAN', 'CLEAN loop gain.', min=0.001, max=1),
    P('clean-niter', 'Iterations', 'int', 1000, 'CLEAN', 'Maximum number of CLEAN components.', min=1, max=1000000),
    P('clean-threshold', 'Threshold', 'float', 0.0, 'CLEAN', 'Stop when the residual peak is below this value (0 = off).', min=0, expert=1),
    P('mem-lambda', 'Lambda', 'float', 0.1, 'MEM', 'Entropy weight.', min=0, max=100),
    P('mem-niter', 'Iterations', 'int', 100, 'MEM', 'Number of MEM iterations.', min=1, max=100000),
  ]))

# ------------------------------------------------------------------ Low-SB
plugins.append(dict(id='mask', name='Mask', tab='Low-SB', order=10,
  description='One shared bit-flag mask per image: automatic masking, region editing, overlay, import/export.',
  steps=[
    S('auto', 'Auto Mask...', 'CatalogPanelMaskAuto', 'Mask > Auto Mask...', 'AUTO', None, ['mask.auto'], ['image']),
    S('overlay', 'Show Mask Overlay', 'CatalogPanelMaskToggleOverlay', 'Mask > Show Mask Overlay', 'NONE', None, None, None, toggle='ogfmask(overlay)'),
    S('overlay_style', 'Overlay Colour/Transparency...', 'CatalogPanelMaskOverlaySettings', 'Mask > Overlay Colour/Transparency...'),
    S('edit', 'Edit Mask', None, None, 'MANUAL', None, ['mask.add', 'mask.erase', 'mask.grow', 'mask.shrink', 'mask.invert', 'mask.clear', 'mask.undo', 'mask.redo'], ['image'], variants=[
        V('add', 'Add Regions to Mask', 'CatalogPanelMaskRegions 0', 'Mask > Add Regions to Mask'),
        V('erase', 'Erase Regions from Mask', 'CatalogPanelMaskRegions 1', 'Mask > Erase Regions from Mask'),
        V('grow', 'Grow...', 'CatalogPanelMaskGrow 0', 'Mask > Grow...'),
        V('shrink', 'Shrink...', 'CatalogPanelMaskGrow 1', 'Mask > Shrink...'),
        V('invert', 'Invert', 'CatalogPanelMaskSimple invert', 'Mask > Invert'),
        V('clear', 'Clear Mask', 'CatalogPanelMaskSimple clear', 'Mask > Clear Mask'),
        V('undo', 'Undo', 'CatalogPanelMaskSimple undo', 'Mask > Undo'),
        V('redo', 'Redo', 'CatalogPanelMaskSimple redo', 'Mask > Redo')]),
    S('save_as', 'Save Mask As...', 'CatalogPanelMaskSaveAs', 'Mask > Save As...', 'AUTO', None, ['mask.export']),
    S('import', 'Import Mask...', 'CatalogPanelMaskImport', 'Mask > Import...', 'MANUAL', None, ['mask.import']),
    S('masked', 'Show Masked Image', 'CatalogPanelMaskShowMasked', 'Mask > Show Masked Image', 'AUTO', None, ['mask.masked']),
    S('stats', 'Mask Statistics', 'CatalogPanelMaskStats', 'Mask > Statistics'),
  ]))

plugins.append(dict(id='icl', name='Intracluster light (ICL)', tab='Low-SB', order=20,
  description='ICL pipeline: masking, background model, BCG centre, profiles, measurements, decomposition.',
  steps=[
    S('mask', '1. Source Masking (shared mask)', 'CatalogPanelICLMask', 'ICL > 1. Source Masking (shared mask)', 'AUTO', None, ['mask.auto'], ['image']),
    S('view_mask', 'Show Mask Overlay', 'CatalogPanelICLViewMask', 'ICL >    Show Mask Overlay'),
    S('save_mask', 'Save Mask As...', 'CatalogPanelICLSaveMask', 'ICL >    Save Mask As...'),
    S('import_mask', 'Import Mask...', 'CatalogPanelICLImportMask', 'ICL >    Import Mask...', 'MANUAL'),
    S('background', '2. Background Model', None, None, 'AUTO', None, ['icl.background'], ['image'], variants=[
        V('polynomial', 'Polynomial Fit', 'CatalogPanelICLBackground polynomial', 'ICL > 2. Background Model > Polynomial Fit'),
        V('chebyshev', 'Chebyshev Fit', 'CatalogPanelICLBackground chebyshev', 'ICL > 2. Background Model > Chebyshev Fit'),
        V('sep_large', 'SEP Large Mesh', 'CatalogPanelICLBackground sep_large', 'ICL > 2. Background Model > SEP Large Mesh')]),
    S('view_bkg', 'View Background', 'CatalogPanelICLViewBkg', 'ICL >    View Background'),
    S('center', '3. Set BCG Center', 'CatalogPanelICLSetCenter', 'ICL > 3. Set BCG Center', 'MANUAL'),
    S('profile', 'Measure Profile', 'CatalogPanelICLProfile', 'ICL >    Measure Profile', 'AUTO', None, ['icl.profile']),
    S('sector', 'Sector Profile...', 'CatalogPanelICLSectorProfile', 'ICL >    Sector Profile...'),
    S('measure', '4. ICL Measurements', 'CatalogPanelICLMeasure', 'ICL > 4. ICL Measurements', 'AUTO', None, ['icl.measure']),
    S('measure_multi', 'Multi-Threshold ICL', 'CatalogPanelICLMeasureMulti', 'ICL >    Multi-Threshold ICL', 'AUTO', None, ['icl.measure-multi']),
    S('decompose', '5. BCG+ICL Decomposition', 'CatalogPanelICLDecompose', 'ICL > 5. BCG+ICL Decomposition', 'AUTO', None, ['icl.decompose']),
    S('color', 'Color Profile...', 'CatalogPanelICLColorProfile', 'ICL > Color Profile...', 'AUTO', None, ['icl.color']),
    S('save_profile', 'Save Profile...', 'CatalogPanelICLSaveProfile', 'ICL > Save Profile...'),
    S('load_profile', 'Load Profile...', 'CatalogPanelICLLoadProfile', 'ICL > Load Profile...'),
    S('export_script', 'Export Script...', 'CatalogPanelICLExportScript', 'ICL > Export Script...'),
    S('import_script', 'Import Script...', 'CatalogPanelICLImportScript', 'ICL > Import Script...'),
    S('icl_settings', 'ICL Settings...', 'CatalogPanelICLSettings', 'ICL > Settings...', settings_only=1),
  ]))

plugins.append(dict(id='lsbg', name='Low-surface-brightness galaxies (LSBG)', tab='Low-SB', order=30,
  description='LSBG pipeline: masking, background cleaning, detection, photometry, Sersic fit, filtering, SVM classification.',
  steps=[
    S('mask', '1. Mask Bright Sources (shared mask)', 'CatalogPanelLSBGMask', 'LSBG > 1. Mask Bright Sources (shared mask)', 'AUTO', None, ['mask.auto'], ['image']),
    S('view_mask', 'Show Masked Image', 'CatalogPanelLSBGViewMask', 'LSBG >    Show Masked Image'),
    S('save_mask', 'Save Mask As...', 'CatalogPanelLSBGSaveMask', 'LSBG >    Save Mask As...'),
    S('import_mask', 'Import Mask...', 'CatalogPanelLSBGImportMask', 'LSBG >    Import Mask...', 'MANUAL'),
    S('clean', '2. Background Model', None, None, 'AUTO', None, ['lsbg.clean'], ['image'], variants=[
        V('sep_large', 'SEP Large Mesh', 'CatalogPanelLSBGClean sep_large', 'LSBG > 2. Background Model > SEP Large Mesh'),
        V('polynomial', 'Polynomial Fit', 'CatalogPanelLSBGClean polynomial', 'LSBG > 2. Background Model > Polynomial Fit'),
        V('chebyshev', 'Chebyshev Fit', 'CatalogPanelLSBGClean chebyshev', 'LSBG > 2. Background Model > Chebyshev Fit')]),
    S('view_clean', 'View Cleaned Image', 'CatalogPanelLSBGViewClean', 'LSBG >    View Cleaned Image'),
    S('detect', '3. Detect LSBG Candidates', 'CatalogPanelLSBGDetect', 'LSBG > 3. Detect LSBG Candidates', 'AUTO', None, ['lsbg.detect']),
    S('photometry', '4. Photometry', 'CatalogPanelLSBGPhotometry', 'LSBG > 4. Photometry', 'AUTO', None, ['lsbg.photometry']),
    S('sersic', '5. Sersic Profile Fit', 'CatalogPanelLSBGSersic', 'LSBG > 5. Sérsic Profile Fit', 'AUTO', None, ['lsbg.sersic']),
    S('filter', '6. Filter + Grade', 'CatalogPanelLSBGFilter', 'LSBG > 6. Filter + Grade', 'AUTO', None, ['lsbg.filter']),
    S('svm', '7. SVM Classify', 'CatalogPanelLSBGSVMClassify', 'LSBG > 7. SVM Classify', 'AUTO', None, ['lsbg.svm-classify']),
    S('save_cat', 'Save Catalog...', 'CatalogPanelSaveCatalog', 'LSBG > Save Catalog...', 'AUTO', None, ['catalog.save'], ['catalog']),
    S('load_cat', 'Load Catalog...', 'CatalogPanelLoadCatalog', 'LSBG > Load Catalog...', 'MANUAL', None, ['catalog.load']),
    S('forced', '8. Forced Photometry (Multi-Band)', 'CatalogPanelLSBGForcedPhot', 'LSBG > 8. Forced Photometry (Multi-Band)', 'MANUAL', None, ['lsbg.forced']),
    S('run_all', 'Run Full Pipeline', 'CatalogPanelLSBGRunAll', 'LSBG > Run Full Pipeline', 'AUTO', None, ['lsbg.run']),
    S('export_script', 'Export Script...', 'CatalogPanelLSBGExportScript', 'LSBG > Export Script...'),
    S('import_script', 'Import Script...', 'CatalogPanelLSBGImportScript', 'LSBG > Import Script...'),
    S('lsbg_settings', 'LSBG Settings...', 'CatalogPanelLSBGSettings', 'LSBG > Settings...', settings_only=1),
  ]))

# ------------------------------------------------------------------ Time-domain (menu entries kept as wrappers; ogf_moving.tcl)
plugins.append(dict(id='moving', name='Moving objects and transients', tab='Time-domain', order=10,
  description='Multi-epoch asteroid / transient pipeline (reference fetch, alignment, ZOGY differencing, tracklets, orbit fit).',
  network=0, tcl=None,
  steps=[
    S('fetch', 'Fetch Reference...', 'OGFMovFetchDialog', 'Moving Objects > Fetch Reference...', 'AUTO', None, ['moving.fetch']),
    S('align', 'Align', 'OGFMovAlign', 'Moving Objects > Align', 'AUTO', None, ['moving.align']),
    S('difference', 'Difference', 'OGFMovDifference', 'Moving Objects > Difference', 'AUTO', None, ['moving.difference']),
    S('detect', 'Detect', 'OGFMovDetect', 'Moving Objects > Detect'),
    S('link', 'Link Tracklets', 'OGFMovLink', 'Moving Objects > Link Tracklets', 'AUTO', None, ['moving.link']),
    S('identify', 'Identify Known Objects', 'OGFMovIdentify', 'Moving Objects > Identify Known Objects', 'AUTO', None, ['moving.identify'], None, network=1),
    S('orbit', 'Orbit Fit...', 'OGFMovOrbitDialog', 'Moving Objects > Orbit Fit...', 'AUTO', None, ['moving.orbit']),
    S('transients', 'Transient Candidates...', 'OGFMovTransients', 'Moving Objects > Transient Candidates...', 'AUTO', None, ['moving.transients']),
    S('lightcurve', 'Light Curve', 'OGFMovLightCurve', 'Moving Objects > Light Curve', 'AUTO', None, ['moving.lightcurve']),
    S('export', 'Export...', 'OGFMovExport', 'Moving Objects > Export...', 'AUTO', None, ['moving.export']),
    S('select', 'Select Exposures...', 'OGFMovSelectFiles', 'Moving Objects > Select Exposures...', 'NONE'),
    S('results', 'Show Results Table', 'OGFMovTable', 'Moving Objects > Show Results Table'),
    S('workdir', 'Work Directory...', 'OGFMovWorkdir', 'Moving Objects > Work Directory...'),
  ]))

# ------------------------------------------------------------------ Results
plugins.append(dict(id='catalog', name='Catalog', tab='Results', order=10,
  description='Catalog input/output and interactive plots.',
  steps=[
    S('save', 'Save Catalog', 'CatalogPanelSaveCatalog', 'SExtractor > Save Catalog', 'AUTO', None, ['catalog.save'], ['catalog']),
    S('export_reg', 'Export Regions (.reg)', 'CatalogPanelExportRegions', 'SExtractor > Export Regions (.reg)', 'NONE', None, None, ['catalog']),
    S('export_fits', 'Export FITS Table', 'CatalogPanelExportFITS', 'SExtractor > Export FITS Table', 'AUTO', None, ['catalog.export_fits'], ['catalog']),
    S('load', 'Load Catalog', 'CatalogPanelLoadCatalog', 'SExtractor > Load Catalog', 'MANUAL', None, ['catalog.load']),
    S('clear', 'Clear', 'CatalogPanelClear', 'SExtractor > Clear', 'MANUAL', None, ['catalog.clear']),
    S('plot', 'Interactive Plot...', 'CatalogPanelPlotDialog', 'Analysis > Interactive Plot...', 'NONE', None, None, ['catalog']),
    S('viewer', 'Analysis Viewer...', 'CatalogPanelAnalysisViewer', 'Analysis > Analysis Viewer...'),
  ]))

plugins.append(dict(id='ai_services', name='AI services', tab='Results', order=20,
  description='Provider-agnostic bridge to external AI services (docs/ai_services.md).', network=1,
  steps=[
    S('registry', 'Service Registry...', 'OGFAIRegistry', 'Analysis > AI Services > Service Registry...', 'NONE'),
    S('run', 'Run Task on Catalog...', 'OGFAIRunDialog', 'Analysis > AI Services > Run Task on Catalog...', 'AUTO', None, ['ai.run'], ['catalog'], network=1),
    S('log', 'Show Last Run Log', 'OGFAIShowLog', 'Analysis > AI Services > Show Last Run Log'),
  ]))

plugins.append(dict(id='session', name='Session recorder', tab='Results', order=30,
  description='Export the GUI session as an input-agnostic Python pipeline script.',
  steps=[
    S('save', 'Save Session as Python Script...', 'CatalogPanelSessionSave', 'Analysis > Save Session as Python Script...'),
    S('show', 'Show Session Log', 'CatalogPanelSessionShow', 'Analysis > Show Session Log'),
    S('reset', 'Reset Session Log', 'CatalogPanelSessionReset', 'Analysis > Reset Session Log'),
  ]))


CHIP = {  # id: (short name, primary step, settings, menu)
 'extract': ('Extract', 'extract', 'params', None),
 'objects': ('Objects', 'mark_all', 'CatalogPanelSeparateSettings', None),
 'bands': ('Bands', 'measure', None, None),
 'star_psf': ('Stars/PSF', 'ai_star', 'CatalogPanelStarPSFSettings', None),
 'galaxy_model': ('Galaxies', 'ai_morph', None, None),
 'morphology': ('Structure', 'sersic', 'params', None),
 'photometry': ('Photometry', 'psf_phot', None, None),
 'photoz_sed': ('Photo-z/SED', 'photoz', None, None),
 'deconv': ('Deconvolve', 'quick', 'params', None),
 'mask': ('Mask', 'auto', 'CatalogPanelMaskOverlaySettings', None),
 'icl': ('ICL', 'measure', 'CatalogPanelICLSettings', None),
 'lsbg': ('LSBG', 'run_all', 'CatalogPanelLSBGSettings', None),
 'moving': ('Moving', 'link', None, None),
 'catalog': ('Catalog', 'save', None, None),
 'ai_services': ('AI Services', 'run', 'OGFAIRegistry', None),
 'session': ('Session', 'save', None, 'Workflow'),
}

for p in plugins:
    sh, prim, st, menu = CHIP[p['id']]
    p['short'] = sh
    if prim: p['primary'] = prim
    if st: p['settings'] = st
    if menu: p['menu'] = menu
    d = os.path.join(OUT, p['id'])
    os.makedirs(d, exist_ok=True)
    p.pop('tcl', None)
    p = {'schema': 1, **p}
    with open(os.path.join(d, 'plugin.json'), 'w') as f:
        json.dump(p, f, indent=1, ensure_ascii=False)
        f.write('\n')
print(len(plugins), 'plugins;', sum(len(p['steps']) for p in plugins), 'steps')
