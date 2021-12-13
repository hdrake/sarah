import numpy as np

def calc_h(x, r): return x**2/r;

def calc_v(x, α, β, r): return 2.*(α+1./r)*x+β
def calc_d(x, α, β, γ): return α*x**2 + β*x + γ
def calc_d_x(x, α, β): return 2.*α*x+β

def calc_xR(α, β, γ): return (- β - np.sqrt(β**2 -4.*α*γ))/(2.*α)
def calc_xL(α, β, γ): return (- β + np.sqrt(β**2 -4.*α*γ))/(2.*α)

def shoot_perturbations(α, β, γ, r, nx):
    xR = calc_xR(α, β, γ)
    xL = calc_xL(α, β, γ)

    if ((xR-xL) <= 0.) or np.isnan(xR-xL): return

    # Adapt grid spacing to width of domain!
    dx = (xR-xL)/nx
    x = np.arange(xL,xR+dx,dx)

    d = calc_d(x, α, β, γ)
    d_x = calc_d_x(x, α, β)
    v = calc_v(x, α, β, r)

    xpL = np.zeros(nx)
    xpR = np.zeros(nx)

    zpL = np.zeros(nx)
    zpR = np.zeros(nx)

    s = 2.*x/r
    sx = 2./r

    xpL[0]=1
    zpL[0]=s[0]*xpL[0]

    # from left
    for i in range(nx//2+1):
        if i==0:
            xpL[i+1] = xpL[i]+dx*(0.5/d_x[i]*(s[i]**2/v[i]**2 + sx)*xpL[i] + zpL[i]/v[i]**2)
            zpL[i+1] = zpL[i]+dx*(v[i]*0.5/d_x[i]*(s[i]**2/v[i]**2 + sx)*xpL[i])
        else:
            xpL[i+1]=xpL[i]+dx*((s[i]/d[i])*xpL[i]+1./d[i]*(d[i]/v[i]**2-1.)*zpL[i])
            zpL[i+1]=zpL[i]+dx*(s[i]*v[i]/d[i]*xpL[i]-v[i]/d[i]*zpL[i])

    xpR[-1]=1
    zpR[-1]=s[-1]*xpR[-1]

    # from right
    for j in range(nx//2+1):
        i = nx-(j+1)
        if j==0:
            xpR[i-1] = xpR[i]-dx*(0.5/d_x[i]*(s[i]**2/v[i]**2 + sx)*xpR[i] + zpR[i]/v[i]**2)
            zpR[i-1] = zpR[i]-dx*(v[i]*0.5/d_x[i]*(s[i]**2/v[i]**2 + sx)*xpR[i])
        else:
            xpR[i-1]=xpR[i]-dx*((s[i]/d[i])*xpR[i]+1/d[i]*(d[i]/v[i]**2-1)*zpR[i])
            zpR[i-1]=zpR[i]-dx*(s[i]*v[i]/d[i]*xpR[i]-v[i]/d[i]*zpR[i])

    # Redo shooting but scaled
    xpR[-1]=xpL[nx//2]/xpR[nx//2]
    zpR[-1]=s[-1]*xpR[-1]

    # from right
    for j in range(nx//2+1):
        i = nx-(j+1)
        if j==0:
            xpR[i-1] = xpR[i]-dx*(0.5/d_x[i]*(s[i]**2/v[i]**2 + sx)*xpR[i] + zpR[i]/v[i]**2)
            zpR[i-1] = zpR[i]-dx*(v[i]*0.5/d_x[i]*(s[i]**2/v[i]**2 + sx)*xpR[i])
        else:
            xpR[i-1]=xpR[i]-dx*((s[i]/d[i])*xpR[i]+1/d[i]*(d[i]/v[i]**2-1)*zpR[i])
            zpR[i-1]=zpR[i]-dx*(s[i]*v[i]/d[i]*xpR[i]-v[i]/d[i]*zpR[i])
    
    xp_tmp = np.zeros(nx)
    xp_tmp[0:nx//2] = xpL[0:nx//2]
    xp_tmp[nx//2:] = xpR[nx//2:]
    
    zp_tmp = np.zeros(nx)
    zp_tmp[0:nx//2] = zpL[0:nx//2]
    zp_tmp[nx//2:] = zpR[nx//2:]
    
    mode = np.sum((zp_tmp[1:]*zp_tmp[:-1]) < 0) # count number of zero crossings
    psi = v*d
    psi[d<0.]=0.
    
    outputs = {
        'match': np.abs(zpL[nx//2]-zpR[nx//2])/np.sqrt(zpL[nx//2]**2 + zpR[nx//2]**2),
        'Q': np.nansum(psi*dx),
        'xp': xp_tmp,
        'zp': zp_tmp,
        'x': x,
        'mode': mode
    }

    return outputs