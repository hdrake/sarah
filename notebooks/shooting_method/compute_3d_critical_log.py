#!/usr/bin/env python
# coding: utf-8

# ### Shooting method implementation for evaluating whether sections are hydraulically controlled

# In[3]:

import sys
sys.path.append("../")
from shooting import *

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

save_perturbations = False


# #### External parameters and functions

# In[4]:


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

# In[4]:


nα = 400
nβ = 400
nγ = 400

αlim = [-2.75, -1.2]
dα = np.diff(αlim)/nα
αvec = np.arange(αlim[0],αlim[1],dα)

βvec = np.logspace(-0.75, np.log10(3), num=nβ)
γvec = np.logspace(-1.5, np.log10(3**2), num=nγ)

nx = 200
x_idx = np.arange(0, nx, 1)

βarr, αarr, γarr = np.meshgrid(βvec, αvec, γvec) # for some reason this is the correct order

match = np.zeros_like(βarr)
Q = np.zeros_like(βarr)
mode = np.zeros_like(βarr)

if save_perturbations:
    xarr = np.zeros((nα, nβ, nγ, nx))
    xp = np.zeros((nα, nβ, nγ, nx))
    zp = np.zeros((nα, nβ, nγ, nx))


# #### Run shooting method across parameters

# In[5]:


for kk, γ in enumerate(γvec):
    print(kk, end=', ')
    for jj, β in enumerate(βvec):
        for ii, α in enumerate(αvec):
            
            outputs = shoot_perturbations(α, β, γ, r, nx)

            if outputs is None:
                pass

            else:
                match[ii,jj,kk] = outputs['match']
                Q[ii,jj,kk] = outputs['Q']
                mode[ii,jj,kk] = outputs['mode']

                if save_perturbations:
                    xp[ii,jj,kk,:] = outputs['xp']
                    zp[ii,jj,kk,:] = outputs['zp']
                    xarr[ii,jj,kk,:] = outputs['x']


# In[6]:


ds = xr.Dataset()
ds['α'] = xr.DataArray(αvec, coords=[αvec], dims=['α'])
ds['β'] = xr.DataArray(βvec, coords=[βvec], dims=['β'])
ds['γ'] = xr.DataArray(γvec, coords=[γvec], dims=['γ'])
ds['match'] = xr.DataArray(match, coords=[αvec, βvec, γvec], dims=['α', 'β', 'γ'], name='match')
ds['Q'] = xr.DataArray(Q, coords=[αvec, βvec, γvec], dims=['α', 'β', 'γ'], name='Q')
ds['mode'] = xr.DataArray(mode, coords=[αvec, βvec, γvec], dims=['α', 'β', 'γ'], name='mode')

if save_perturbations:
    ds['x'] = xr.DataArray(xarr, coords=[αvec, βvec, γvec, x_idx], dims=['α', 'β', 'γ', 'x_idx'], name='x')
    ds['d'] = calc_d(ds['x'], ds['α'], ds['β'], ds['γ'])
    ds['xp'] = xr.DataArray(xp, coords=[αvec, βvec, γvec, x_idx], dims=['α', 'β', 'γ', 'x_idx'], name='xp')
    ds['zp'] = xr.DataArray(zp, coords=[αvec, βvec, γvec, x_idx], dims=['α', 'β', 'γ', 'x_idx'], name='zp')

ds.to_netcdf('../../data/critical_3d_log.nc')


# In[ ]:




