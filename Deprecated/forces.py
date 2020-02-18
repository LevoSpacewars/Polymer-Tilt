import math


def V_chain(r,rmin, rmax, l_max, l_0, K):
    right = math.log(1-(r-l_0)**2/l_max**2,10)
    left = -K*l_max**2/2
    V = right*left
    F = K( r - l_0)
    F = F/(l_max*(1 - (r-l_0)**2/l_max**2))
    return (V,F)
def Particle_interaction(r, rmin, rmax, l_0, K):
    if r < l_0:
        k_c = K
        V = k_c *(r - l_0)**2
        F = -2*k_c*(r-l_0)
        return (V,F)
    else:
        return (0,0)
