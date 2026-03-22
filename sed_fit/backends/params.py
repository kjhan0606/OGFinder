"""Parameter mapping between canonical names and SPS code conventions.

All backends use: Chabrier IMF, Calzetti dust law, delayed-tau SFH.
"""

import numpy as np


def canonical_to_fsps(z, log_mass, log_age, log_Z, Av, log_tau):
    """Convert canonical parameters to python-fsps format.

    Returns dict for fsps.StellarPopulation parameter setting.
    """
    # FSPS uses log(Z/Zsun) directly, age in Gyr, tau in Gyr
    age_gyr = 10 ** (log_age - 9)  # yr -> Gyr
    tau_gyr = 10 ** (log_tau - 9)
    return {
        'zred': z,
        'logzsol': log_Z,        # log(Z/Zsun)
        'dust2': Av * 0.4 / 1.0,  # FSPS dust2 ~ Av / (Rv * 1.086)
        'tage': age_gyr,          # age in Gyr
        'tau': tau_gyr,           # e-folding in Gyr
        'mass': 10 ** log_mass,
        'imf_type': 1,            # Chabrier
        'dust_type': 0,           # Calzetti
        'sfh': 4,                 # delayed-tau
    }


def canonical_to_bagpipes(z, log_mass, log_age, log_Z, Av, log_tau):
    """Convert canonical parameters to Bagpipes model_components format."""
    age_gyr = 10 ** (log_age - 9)
    tau_gyr = 10 ** (log_tau - 9)
    return {
        'redshift': z,
        'delayed': {
            'age': age_gyr,             # Gyr
            'tau': tau_gyr,             # Gyr
            'massformed': 10 ** log_mass,
            'metallicity': 10 ** log_Z * 0.02,  # absolute Z (Zsun=0.02)
        },
        'dust': {
            'type': 'Calzetti',
            'Av': Av,
        },
    }


def canonical_to_prospector(z, log_mass, log_age, log_Z, Av, log_tau):
    """Convert canonical parameters to Prospector model params."""
    age_gyr = 10 ** (log_age - 9)
    tau_gyr = 10 ** (log_tau - 9)
    return {
        'zred': z,
        'logzsol': log_Z,
        'dust2': Av * 0.4 / 1.0,
        'tage': age_gyr,
        'tau': tau_gyr,
        'mass': 10 ** log_mass,
        'imf_type': 1,
        'sfh': 4,
    }


def canonical_to_cigale(z, log_mass, log_age, log_Z, Av, log_tau):
    """Convert canonical parameters to CIGALE input format."""
    age_myr = 10 ** (log_age - 6)  # yr -> Myr
    tau_myr = 10 ** (log_tau - 6)
    return {
        'redshift': z,
        'sfhdelayed': {
            'age_main': age_myr,           # Myr
            'tau_main': tau_myr,           # Myr
        },
        'bc03': {
            'metallicity': 10 ** log_Z * 0.02,  # absolute Z
            'imf': 1,                              # Chabrier
        },
        'dustatt_calzleit': {
            'Av': Av,
        },
        'stellar_mass': 10 ** log_mass,
    }


def canonical_to_dense_basis(z, log_mass, log_age, log_Z, Av, log_tau):
    """Convert canonical parameters to Dense Basis format."""
    return {
        'z': z,
        'logM': log_mass,
        'logZsol': log_Z,
        'dust2': Av * 0.4 / 1.0,
        'tage': 10 ** (log_age - 9),
        'tau': 10 ** (log_tau - 9),
    }


# Dispatch table
_CONVERTERS = {
    'fsps': canonical_to_fsps,
    'bagpipes': canonical_to_bagpipes,
    'prospector': canonical_to_prospector,
    'cigale': canonical_to_cigale,
    'dense_basis': canonical_to_dense_basis,
}


def convert_params(backend_name, z, log_mass, log_age, log_Z, Av, log_tau):
    """Convert canonical parameters to backend-specific format.

    Parameters
    ----------
    backend_name : str
    z, log_mass, log_age, log_Z, Av, log_tau : float

    Returns
    -------
    dict : backend-specific parameter dictionary.
    """
    converter = _CONVERTERS.get(backend_name)
    if converter is None:
        raise ValueError(f"Unknown backend: {backend_name}")
    return converter(z, log_mass, log_age, log_Z, Av, log_tau)
