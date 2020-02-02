#!/usr/bin/env python
# coding: utf-8

# ### Shooting method implementation for evaluating whether sections are hydraulically controlled

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from utilities import *

save_perturbations = False


# #### External parameters and functions

# In[2]:


def f(ϕ): return 2. * (2*np.pi)/(60.**2 * 24.) * np.sin(np.deg2rad(ϕ))

δρ = 0.45
ρ0 = 1028.
g = 9.81
gp = (δρ/ρ0)*g

ϕ = 62. # Latitude of Faroe Bank channel
f0 = f(ϕ)
αsill = 5.8e-6 # from Borenas and Lundberg (1988)
r = f0**2 / (gp*αsill)


# #### Define range for parameter sweep

# In[15]:


nβ = np.int64(2e7)

α = -1.64
γ = 1.0

βlim = [0., 2.]
dβ = np.diff(βlim)/nβ
βvec = np.arange(βlim[0],βlim[1],dβ)

nx = 200
x_idx = np.arange(0, nx, 1)

match = np.zeros_like(βvec)
Q = np.zeros_like(βvec)
mode = np.zeros_like(βvec)

if save_perturbations:
    xarr = np.zeros((nβ, nx))
    xp = np.zeros((nβ, nx))
    zp = np.zeros((nβ, nx))


# #### Run shooting method across parameters

# In[7]:


for ii, β in enumerate(βvec):
    outputs = shoot_perturbations(α, β, γ, r, nx)

    if outputs is None:
        pass

    else:
        match[ii] = outputs['match']
        Q[ii] = outputs['Q']
        mode[ii] = outputs['mode']

        if save_perturbations:
            xp[ii,:] = outputs['xp']
            zp[ii,:] = outputs['zp']
            xarr[ii,:] = outputs['x']


# In[8]:


ds = xr.Dataset()
ds['β'] = xr.DataArray(βvec, coords=[βvec], dims=['β'])
ds['match'] = xr.DataArray(match, coords=[βvec], dims=['β'], name='match')
ds['Q'] = xr.DataArray(Q, coords=[βvec], dims=['β'], name='Q')
ds['mode'] = xr.DataArray(mode, coords=[βvec], dims=['β'], name='mode')

if save_perturbations:
    ds['x'] = xr.DataArray(xarr, coords=[βvec, x_idx], dims=['β', 'x_idx'], name='x')
    ds['d'] = calc_d(ds['x'], α, ds['β'], γ)
    ds['xp'] = xr.DataArray(xp, coords=[βvec, x_idx], dims=['β', 'x_idx'], name='xp')
    ds['zp'] = xr.DataArray(zp, coords=[βvec, x_idx], dims=['β', 'x_idx'], name='zp')

ds.to_netcdf('../../data/critical_1d.nc')

