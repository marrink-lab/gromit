#!/usr/bin/env python3

import sys

epsR, epsRF, RC, LD, LR, LC, LS = [float(i) for i in sys.argv[1:]]

Krf = (epsRF-epsR)/(2*epsRF+epsR)/(RC**3)
Crf = 3*epsRF/(2*epsRF+epsR)/RC

# Shifted power
def vf_fun(a,rc,rs):
    if rs < 0: 
        return lambda r: r**-a, lambda r: a*r**-(a+1)
    A, B = -((a+4)*rc-(a+1)*rs)/((rc-rs)**2)/(rc**(a+2)), ((a+3)*rc-(a+1)*rs)/((rc-rs)**3)/(rc**(a+2))
    C    = (rc**-a)-(a*A*(rc-rs)**3/3)-(a*B*(rc-rs)**4/4)
    V    = lambda r: r**-a - C - (r>rs and a*(A*(r-rs)**3/3+B*(r-rs)**4/4))
    F    = lambda r: a*r**-(a+1) + (r>rs and a*(A*(r-rs)**2+B*(r-rs)**3))
    return V, F

VD,FD = vf_fun(LD,LC,LS)
VR,FR = vf_fun(LR,LC,LS)

r  = [ i/1e3 or 0 for i in range(2,3001,2)]
print("""#
# Coulomb cut-off/reaction-field: epsRF = %f, epsR = %f, RC = %f
# Lennard-Jones dispersion: power=%f, cutoff=%f, %sshifted %s
# Lennard-Jones repulsion:  power=%f, cutoff=%f, %sshifted %s""" % (
    epsRF, epsR, RC,
    LD, LC, LS<0 and "not " or "", LS>=0 and "(rshift=%f)"% LS or "",
    LR, LC, LS<0 and "not " or "", LS>=0 and "(rshift=%f)"% LS or ""))

f0 = [ (1/i+Krf*i*i-Crf)/epsR for i in r ]
f1 = [ (1/(i*i)-2*Krf*i)/epsR for i in r ]

g0 = [ i<LC and -VD(i) for i in r ]
g1 = [ i<LC and -FD(i) for i in r ]

h0 = [ i<LC and  VR(i) for i in r ]
h1 = [ i<LC and  FR(i) for i in r ]

table = [(0,0,0,0,0,0,0)] + list(zip(r,f0,f1,g0,g1,h0,h1))

print("\n".join([(7*"%.10e ")%i for i in table]))


