import numpy as np
import utils

g = 9.81
ρ0 = 1028.

def add_derived_variables(d):
    d['gp'] = g*d['δρ']/ρ0
    d['f0'] = utils.f(d['ϕ'])
    d['Ld'] = np.sqrt(d['gp']*d['H'])/d['f0']

    d['rsill'] = d['f0']**2 / (d['gp']*d['αsill'])
    d['κsill'] = 2./d['rsill']

    d['λ'] = d['L']/d['Lsill'] # inverse along-channel length of sill
    return
    
# Faroe Bank Channel
fbc = {
    'δρ' : 0.45,
    'ϕ' : 62.,
    
    'H' : 500.,
    'L' : 10.e3,
    'W' : 10.e3,
    'ηinf' : 1000.,

    'αsill' : 5.8e-6,
    'Hsill' : 400.,
    'Lsill' : 100.e3,
}
add_derived_variables(fbc)

