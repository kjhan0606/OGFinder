"""Band registry shared by the photo-z / SED adapters: catalog band name -> pivot wavelength + the name each code uses for the same filter.

`lam` is the pivot wavelength in micron, `dl` the approximate FWHM (micron; only used by the TOY test models), `eazy` a list of regular expressions
searched (first match wins) in the names of an EAZY FILTER.RES file; `sedpy`, `cigale`, `bagpipes` are the names those codes use.  Names of the
codes were checked against the installed packages where they could be installed (docs/sed_codes.md says which); an unknown band is
skipped with a warning, never guessed.
"""
import re

B = {}


def _add(name, lam, dl, eazy, sedpy=None, cigale=None, bagpipes=None):
    B[name.upper()] = dict(name=name.upper(), lam=lam, dl=dl, eazy=eazy if isinstance(eazy, list) else [eazy], sedpy=sedpy, cigale=cigale, bagpipes=bagpipes)


for n, lam, dl in (('F435W', .4318, .0920), ('F475W', .4747, .1520), ('F555W', .5361, .1240), ('F606W', .5921, .2320), ('F625W', .6311, .1560),
                   ('F775W', .7693, .1500), ('F814W', .8057, .1540), ('F850LP', .9033, .1180)):
    _add(n, lam, dl, [r'wfc_%s_t81' % n.lower()], 'acs_wfc_%s' % n.lower(), 'hst.acs.%s' % n.upper(), 'filters/HST_ACS_WFC.%s.dat' % n.upper())
for n, lam, dl in (('F098M', .9864, .1690), ('F105W', 1.0552, .2650), ('F110W', 1.1534, .4430), ('F125W', 1.2486, .3010), ('F140W', 1.3923, .3840), ('F160W', 1.5369, .2870)):
    _add(n, lam, dl, [r'hst/wfc3/IR/%s\.dat' % n.lower()], 'wfc3_ir_%s' % n.lower(), 'hst.wfc3.%s' % n.upper(), 'filters/HST_WFC3_IR.%s.dat' % n.upper())
for n, lam, dl in (('F275W', .2710, .0405), ('F336W', .3355, .0511)):
    _add(n, lam, dl, [r'hst/wfc3/UVIS/%s\.dat' % n.lower()], 'wfc3_uvis_%s' % n.lower(), 'hst.wfc3.%s' % n.upper(), 'filters/HST_WFC3_UVIS.%s.dat' % n.upper())
for n, lam, dl in (('F090W', .9022, .1760), ('F115W', 1.1543, .2050), ('F150W', 1.5007, .3180), ('F200W', 1.9886, .4610), ('F277W', 2.7617, .6820),
                   ('F356W', 3.5682, .7240), ('F410M', 4.0822, .4280), ('F444W', 4.4036, 1.0230)):
    _add(n, lam, dl, [r'jwst_nircam_%s' % n.lower()], 'jwst_%s' % n.lower(), 'jwst.nircam.%s' % n.upper(), 'filters/JWST_NIRCam.%s.dat' % n.upper())
for n, lam, dl, ez, sp in (('SDSS_U', .3556, .0540, 'SDSS_filter_u', 'sdss_u0'), ('SDSS_G', .4702, .1380, 'SDSS_filter_g', 'sdss_g0'), ('SDSS_R', .6176, .1380, 'SDSS_filter_r', 'sdss_r0'),
                           ('SDSS_I', .7496, .1530, 'SDSS_filter_i', 'sdss_i0'), ('SDSS_Z', .8947, .1370, 'SDSS_filter_z', 'sdss_z0')):
    _add(n, lam, dl, [ez], sp, 'sdss.%s' % n[-1].lower(), 'filters/SLOAN_SDSS.%s.dat' % n[-1].lower())
_add('J', 1.25, .16, [r'VISTA/J_', r'2mass.*[_/]j'], 'twomass_J', '2mass.J', 'filters/2MASS_2MASS.J.dat')
_add('H', 1.65, .25, [r'VISTA/H_', r'2mass.*[_/]h'], 'twomass_H', '2mass.H', 'filters/2MASS_2MASS.H.dat')
_add('KS', 2.16, .26, [r'VISTA/Ks_', r'2mass.*[_/]ks?'], 'twomass_Ks', '2mass.Ks', 'filters/2MASS_2MASS.Ks.dat')
_add('IRAC1', 3.55, .75, [r'irac_tr1'], 'spitzer_irac_ch1', 'spitzer.irac.ch1', 'filters/Spitzer_IRAC.I1.dat')
_add('IRAC2', 4.49, 1.01, [r'irac_tr2'], 'spitzer_irac_ch2', 'spitzer.irac.ch2', 'filters/Spitzer_IRAC.I2.dat')
ALIASES = {'ACS_F606W': 'F606W', 'ACS_F814W': 'F814W', 'ACS_F850LP': 'F850LP', 'ACS_F435W': 'F435W', 'ACS_F775W': 'F775W', 'ACS_F475W': 'F475W', 'ACS_F555W': 'F555W',
           'WFC3_F105W': 'F105W', 'WFC3_F125W': 'F125W', 'WFC3_F160W': 'F160W', 'WFC3_F140W': 'F140W', 'WFC3_F110W': 'F110W',
           'CH1': 'IRAC1', 'CH2': 'IRAC2', '3.6': 'IRAC1', '4.5': 'IRAC2', 'KSBAND': 'KS', 'K': 'KS'}


def lookup(band):
    """Registry entry of a catalog band name (case-insensitive; 'ACS_F606W', 'F606W', 'SDSS_G'/'g'); None when unknown."""
    k = str(band).strip().upper().replace('-', '_')
    k = ALIASES.get(k, k)
    if k in B:
        return B[k]
    return None


def eazy_filter_index(entry, filter_names):
    """1-based index of the matching filter in a list of EAZY FILTER.RES names (first pattern with a hit; None if absent)."""
    for pat in entry['eazy']:
        rx = re.compile(pat, re.I)
        for i, nm in enumerate(filter_names):
            if rx.search(nm):
                return i + 1
    return None


def known_bands():
    return sorted(B, key=lambda k: B[k]['lam'])
