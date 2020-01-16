
import numpy as np


# Example of critical flow

dx = 0.0015

r = 0.6
x = np.arange(-1,1+dx,dx)
h0 = 0
h = h0 + x**2/r

nalpha = 300
alpha_min = -1.75
Dalpha = 0.25
dalpha = Dalpha/nalpha
alphas = np.arange(alpha_min,alpha_min+Dalpha,dalpha)

nbeta = 300
Dbeta = 0.32
dbeta = Dbeta/nbeta
betas = np.arange(0,Dbeta,dbeta)

ngamma = 300
Dgamma = 1.
dgamma = Dgamma/ngamma
gammas = np.arange(0,Dgamma,dgamma)

Betas,Gammas,Alphas = np.meshgrid(betas,gammas,alphas)
match = np.zeros_like(Alphas)
Q = np.zeros_like(match)

for kk in range(1,nalpha):
    alpha = alphas[kk]
    print(str(kk),end=" ")
    for jj in range(1,nbeta):
        beta = betas[jj]
        for ii in range(1,ngamma):
            gamma = gammas[ii]

            d = alpha*x**2+beta*x+gamma
            d[d+h<h] = np.nan

            a = (beta - np.sqrt(beta**2-4*alpha*gamma))/(2*alpha)
            b = (- beta - np.sqrt(beta**2-4*alpha*gamma))/(2*alpha)

            x = np.arange(-a,b+dx,dx)
            nx = np.size(x)

            xpl = np.zeros(np.size(x))
            xpr = np.zeros(np.size(x))

            zpl = np.zeros(np.size(x))
            zpr = np.zeros(np.size(x))

            h = h0 + x**2/r
            s = 2*x/r
            sx = 2/r+np.zeros_like(s)
            d = alpha*x**2+beta*x+gamma
            v = 2*(alpha+1./r)*x+beta
            ddx = 2*alpha*x+beta

            xpl[0]=1
            zpl[0]=s[0]*xpl[0]

            # from left
            for i in range(np.int(np.size(xpl)/2)+1):
                if i==0:
                    xpl[i+1] = xpl[i]+dx*(0.5/ddx[i]*(s[i]**2/v[i]**2 + sx[i])*xpl[i] + zpl[i]/v[i]**2)
                    zpl[i+1] = zpl[i]+dx*(v[i]*0.5/ddx[i]*(s[i]**2/v[i]**2 + sx[i])*xpl[i])                   
                else:
                    xpl[i+1]=xpl[i]+dx*((s[i]/d[i])*xpl[i]+1/d[i]*(d[i]/v[i]**2-1)*zpl[i])
                    zpl[i+1]=zpl[i]+dx*(s[i]*v[i]/d[i]*xpl[i]-v[i]/d[i]*zpl[i])

            xpr[-1]=1
            zpr[-1]=s[-1]*xpr[-1]

            # from right
            for j in range(np.int(np.size(xpl)/2)+1):
                i = nx-(j+1)
                if j==0:
                    xpr[i-1] = xpr[i]-dx*(0.5/ddx[i]*(s[i]**2/v[i]**2 + sx[i])*xpr[i] + zpr[i]/v[i]**2)
                    zpr[i-1] = zpr[i]-dx*(v[i]*0.5/ddx[i]*(s[i]**2/v[i]**2 + sx[i])*xpr[i])                   
                else:
                    xpr[i-1]=xpr[i]-dx*((s[i]/d[i])*xpr[i]+1/d[i]*(d[i]/v[i]**2-1)*zpr[i])
                    zpr[i-1]=zpr[i]-dx*(s[i]*v[i]/d[i]*xpr[i]-v[i]/d[i]*zpr[i])


            # Redo shooting but scaled
            xpr[-1]=xpl[np.int(nx/2)]/xpr[np.int(nx/2)]
            zpr[-1]=s[-1]*xpr[-1]

            # from right
            for j in range(np.int(np.size(xpl)/2)+1):
                i = nx-(j+1)
                if j==0:
                    xpr[i-1] = xpr[i]-dx*(0.5/ddx[i]*(s[i]**2/v[i]**2 + sx[i])*xpr[i] + zpr[i]/v[i]**2)
                    zpr[i-1] = zpr[i]-dx*(v[i]*0.5/ddx[i]*(s[i]**2/v[i]**2 + sx[i])*xpr[i])                   
                else:
                    xpr[i-1]=xpr[i]-dx*((s[i]/d[i])*xpr[i]+1/d[i]*(d[i]/v[i]**2-1)*zpr[i])
                    zpr[i-1]=zpr[i]-dx*(s[i]*v[i]/d[i]*xpr[i]-v[i]/d[i]*zpr[i])

            match[ii,jj,kk] = zpl[np.int(nx/2)]/zpr[np.int(nx/2)]
            Q[ii,jj,kk] = np.nansum(v*d*dx)

np.savez("../data/critical_zoom.npz", Alphas=Alphas, Betas=Betas, Gammas=Gammas, match=match, Q=Q)
