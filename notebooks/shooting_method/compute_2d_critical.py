#!/usr/bin/env python
# coding: utf-8

# ### Shooting method implementation for evaluating whether sections are hydraulically controlled

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from utilities import *


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


nβ = 5000
nγ = 5000

α = -2.16666666

βlim = [0., 2.]
dβ = np.diff(βlim)/nβ
βvec = np.arange(βlim[0],βlim[1],dβ)

γlim = [0., 2.]
dγ = np.diff(γlim)/nγ
γvec = np.arange(γlim[0],γlim[1],dγ)

nx = 200
x_idx = np.arange(0, nx, 1)

γarr, βarr = np.meshgrid(γvec, βvec)
match = np.zeros_like(βarr)
Q = np.zeros_like(βarr)
mode = np.zeros_like(βarr)

xarr = np.zeros((nβ, nγ, nx))
xp = np.zeros((nβ, nγ, nx))
zp = np.zeros((nβ, nγ, nx))


# #### Run shooting method across parameters

# In[ ]:


for jj, γ in enumerate(γvec):
    print(jj, end=', ')
    for ii, β in enumerate(βvec):
                
        outputs = shoot_perturbations(α, β, γ, r, nx)
        
        if outputs is None:
            pass
        
        else:
            match[ii,jj] = outputs['match']
            Q[ii,jj] = outputs['Q']
            mode[ii,jj] = outputs['mode']
            
            xp[ii,jj,:] = outputs['xp']
            zp[ii,jj,:] = outputs['zp']
            xarr[ii,jj,:] = outputs['x']


# In[ ]:


ds = xr.Dataset()
ds['β'] = xr.DataArray(βvec, coords=[βvec], dims=['β'])
ds['γ'] = xr.DataArray(γvec, coords=[γvec], dims=['γ'])
ds['match'] = xr.DataArray(match, coords=[βvec, γvec], dims=['β', 'γ'], name='match')
ds['Q'] = xr.DataArray(Q, coords=[βvec, γvec], dims=['β', 'γ'], name='Q')
ds['mode'] = xr.DataArray(mode, coords=[βvec, γvec], dims=['β', 'γ'], name='mode')

ds['x'] = xr.DataArray(xarr, coords=[βvec, γvec, x_idx], dims=['β', 'γ', 'x_idx'], name='x')
ds['d'] = calc_d(ds['x'], α, ds['β'], ds['γ'])
ds['xp'] = xr.DataArray(xp, coords=[βvec, γvec, x_idx], dims=['β', 'γ', 'x_idx'], name='xp')
ds['zp'] = xr.DataArray(zp, coords=[βvec, γvec, x_idx], dims=['β', 'γ', 'x_idx'], name='zp')

ds.to_netcdf('../../data/critical_2d.nc')


# In[ ]:


tol = 0.25
c = (ds['match']>=1-tol) & (ds['match']<=1+tol)

q = c.plot()
q.set_clim([0,1])


# In[ ]:




