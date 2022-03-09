#!/usr/bin/env python
# coding: utf-8

# # Running Strong Coupling

"""Computing the running of the strong coupling within the Standard Model"""


# ## Preliminaries

import rundec


# ## Metadata

__author__ = 'Benjamin Wallisch'
__copyright__ = 'Copyright 2021'
__license__ = 'GPL'
__version__ = '1.0.0'
__maintainer__ = 'Benjamin Wallisch'
__email__ = 'bwallisch@ias.edu'
__status__ = 'Production'


# ## Running $\mathbf{\alpha}_\mathbf{s}$

def alpha_s(temp, n_loops=5):
    """Strong coupling constant at energy scale/temperature 'temp' computed from the n-loop
    renormalization group evolution within the Standard Model, including the effect of decoupling
    quarks below their masses

    Parameters
    ----------
    temp : float
        Energy scale/temperature (in GeV)

    n_loops : int
        Number of loops in the renormalization group evolution

    Returns
    -------
    float
        Strong coupling constant at 'temp'
    """
    # Parameters (all from PDG 2020)
    alphas_at_z = 0.1179  # alpha_s(m_Z)
    mass_z = 91.1876  # m_Z
    mass_top = 172.76  # (direct measurement)
    mass_bottom = 4.78  # (pole mass)
    mass_charm = 1.67  # (pole mass)

    # Initialize RunDec
    crd = rundec.CRunDec()

    if temp >= mass_top:
        # 6 active flavors above top mass
        n_active_flavors = 6
        # Run from m_Z up to top mass with decoupled top, then run to higher energies
        crd.nfMmu.nf = n_active_flavors
        crd.nfMmu.Mth = mass_top
        crd.nfMmu.muth = mass_top
        result = crd.AlL2AlH(alphas_at_z, mass_z, crd.nfMmu, temp, n_loops)

    elif temp >= mass_bottom:
        # 5 active flavors between bottom and top mass
        n_active_flavors = 5
        # Run from m_Z to larger/smaller energies
        result = crd.AlphasExact(alphas_at_z, mass_z, temp, n_active_flavors, n_loops)

    elif temp >= mass_charm:
        # 4 active flavors between charm and bottom mass
        n_active_flavors = 4
        # Run from m_Z down to bottom mass, then decouple bottom and run to smaller energies
        crd.nfMmu.nf = n_active_flavors + 1
        crd.nfMmu.Mth = mass_bottom
        crd.nfMmu.muth = mass_bottom
        result = crd.AlH2AlL(alphas_at_z, mass_z, crd.nfMmu, temp, n_loops)

    elif temp < mass_charm:
        # 3 active flavors below charm
        n_active_flavors = 3
        # Run from m_Z down to bottom mass, decouple bottom, run down to charm mass, decouple charm
        # and run to smaller energies
        crd.nfMmu.nf = n_active_flavors + 2
        crd.nfMmu.Mth = mass_bottom
        crd.nfMmu.muth = mass_bottom
        alphas_at_c =  crd.AlH2AlL(alphas_at_z, mass_z, crd.nfMmu, mass_charm, n_loops)

        crd.nfMmu.nf = n_active_flavors + 1
        crd.nfMmu.Mth = mass_charm
        crd.nfMmu.muth = mass_charm
        result = crd.AlH2AlL(alphas_at_c, mass_charm, crd.nfMmu, temp, n_loops)

    if result == 0:
        raise ValueError(f"No result; {temp} GeV likely in strong-coupling regime.")

    return result
