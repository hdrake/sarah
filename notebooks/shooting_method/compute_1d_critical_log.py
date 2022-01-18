#!/usr/bin/env python
# coding: utf-8

# ### Shooting method implementation for evaluating whether sections are hydraulically controlled

print("Importing modules")
import sys
sys.path.append("../")
from shooting import *

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

print("Test")
def f(ϕ):
    return(2. * (2*np.pi)/(60.**2 * 24.) * np.sin(np.deg2rad(ϕ)))

δρ = 0.45
ρ0 = 1028.
g = 9.81
gp = (δρ/ρ0)*g

ϕ = 62. # Latitude of Faroe Bank channel
f0 = f(ϕ)
αsill = 5.8e-6 # from Borenas and Lundberg (1988)
r = f0**2 / (gp*αsill)

# #### Define range for parameter sweep
α = -1.45
γ = 6.

nβ = np.int64(3e7)
βlim = [0.1, 2.0]
βvec = np.logspace(np.log10(βlim[0]), np.log10(βlim[1]), num=nβ)

nx = 400
x_idx = np.arange(0, nx, 1)

match = np.zeros_like(βvec)
Q = np.zeros_like(βvec)
mode = np.zeros_like(βvec)

# #### Run shooting method across parameters

for ii, β in enumerate(βvec):
    outputs = shoot_perturbations(α, β, γ, r, nx)

    if outputs is None:
        pass

    else:
        match[ii] = outputs['match']
        Q[ii] = outputs['Q']
        mode[ii] = outputs['mode']

ds = xr.Dataset()
ds.attrs['α'] = α
ds['β'] = xr.DataArray(βvec, coords=[βvec], dims=['β'])
ds.attrs['γ'] = γ
ds['match'] = xr.DataArray(match, coords=[βvec], dims=['β'], name='match')
ds['Q'] = xr.DataArray(Q, coords=[βvec], dims=['β'], name='Q')
ds['mode'] = xr.DataArray(mode, coords=[βvec], dims=['β'], name='mode')

ds.to_netcdf('../../data/critical_1d_log.nc')

