#!/usr/bin/env python
# coding: utf-8

# ### Shooting method implementation for evaluating whether sections are hydraulically controlled

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from utilities import *

save_all=False

# #### External parameters and functions

# In[ ]:


def f(ϕ): return 2. * (2*np.pi)/(60.**2 * 24.) * np.sin(np.deg2rad(ϕ))

δρ = 0.45
ρ0 = 1028.
g = 9.81
gp = (δρ/ρ0)*g

ϕ = 62. # Latitude of Faroe Bank channel
f0 = f(ϕ)
αsill = 5.8e-6 # from Borenas and Lundberg (1988)
rsill = f0**2 / (gp*αsill)

r = 0.6


# #### Define range for parameter sweep

# In[ ]:


nα = 400
nβ = 400
nγ = 400

αlim = [-2.75, -1.2]
dα = np.diff(αlim)/nα
αvec = np.arange(αlim[0],αlim[1],dα)

βlim = [0., 2.1]
dβ = np.diff(βlim)/nβ
βvec = np.arange(βlim[0],βlim[1],dβ)

γlim = [0., 2.2]
dγ = np.diff(γlim)/nγ
γvec = np.arange(γlim[0],γlim[1],dγ)

nx = 200
x_idx = np.arange(0, nx, 1)

βarr, αarr, γarr = np.meshgrid(βvec, αvec, γvec) # for some reason this is the correct order

match = np.zeros_like(βarr)
Q = np.zeros_like(βarr)
mode = np.zeros_like(βarr)

if save_all:
    xarr = np.zeros((nα, nβ, nγ, nx))
    xp = np.zeros((nα, nβ, nγ, nx))
    zp = np.zeros((nα, nβ, nγ, nx))


# #### Run shooting method across parameters

# In[ ]:


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

                if save_all:
                    xp[ii,jj,kk,:] = outputs['xp']
                    zp[ii,jj,kk,:] = outputs['zp']
                    xarr[ii,jj,kk,:] = outputs['x']


ds = xr.Dataset()
ds['α'] = xr.DataArray(αvec, coords=[αvec], dims=['α'])
ds['β'] = xr.DataArray(βvec, coords=[βvec], dims=['β'])
ds['γ'] = xr.DataArray(γvec, coords=[γvec], dims=['γ'])
ds['match'] = xr.DataArray(match, coords=[αvec, βvec, γvec], dims=['α', 'β', 'γ'], name='match')
ds['Q'] = xr.DataArray(Q, coords=[αvec, βvec, γvec], dims=['α', 'β', 'γ'], name='Q')
ds['mode'] = xr.DataArray(mode, coords=[αvec, βvec, γvec], dims=['α', 'β', 'γ'], name='mode')

ds['d'] = calc_d(ds['x'], ds['α'], ds['β'], ds['γ'])

if save_all:
    ds['x'] = xr.DataArray(xarr, coords=[αvec, βvec, γvec, x_idx], dims=['α', 'β', 'γ', 'x_idx'], name='x')
    ds['xp'] = xr.DataArray(xp, coords=[αvec, βvec, γvec, x_idx], dims=['α', 'β', 'γ', 'x_idx'], name='xp')
    ds['zp'] = xr.DataArray(zp, coords=[αvec, βvec, γvec, x_idx], dims=['α', 'β', 'γ', 'x_idx'], name='zp')

ds.to_netcdf('../../data/critical_3d.nc')


# In[ ]:


tol = 0.25
c = (ds['match'].isel(α=5)>=1-tol) & (ds['match'].isel(α=5)<=1+tol)

q = c.plot()
q.set_clim([0,1])


# In[ ]:




