import numpy as np
import scipy as sp
import scipy.integrate as integrate
import scipy.optimize as optimize
import matplotlib.pyplot as plt

def linearBezier(t, p0, p1):
    return [ (1.-t)*p0[0]+t*p1[0], (1.-t)*p0[1]+t*p1[1] ]

def quadraticBezier(t, p0, p1, p2):
    return linearBezier(t, linearBezier(t, p0, p1), linearBezier(t, p1, p2))

def cubicBezier(t, p0, p1, p2, p3):
    return linearBezier(t, quadraticBezier(t, p0, p1, p2), quadraticBezier(t, p1, p2, p3))

def norm(a):
    return np.sqrt(a[0]*a[0]+a[1]*a[1])

def optimizationIntegrandLow(t, p0, p1, p2, p3, a):
    b3 = cubicBezier(t, p0, p1, p2, p3)
    bhalf = cubicBezier(.5, p0, p1, p2, p3)
    b2lo = quadraticBezier(2.*t, p0, [-a[0]+bhalf[0], -a[1]+bhalf[1]], bhalf)
    return norm([b3[0]-b2lo[0],b3[1]-b2lo[1]]);

def optimizationIntegrandHigh(t, p0, p1, p2, p3, a):
    b3 = cubicBezier(t, p0, p1, p2, p3)
    bhalf = cubicBezier(t, p0, p1, p2, p3)
    b2hi = quadraticBezier(2.*(t-.5), bhalf, [bhalf[0]-a[0],bhalf[1]-a[1]], p3)
    return norm([b3[0]-b2hi[0],b3[1]-b2hi[1]]);

def optimizationIntegral(a, p0, p1, p2, p3):
    (lo, loerr) = integrate.quad(optimizationIntegrandLow, 0., .5, args=(p0,p1,p2,p3,a))
    (hi, hierr) = integrate.quad(optimizationIntegrandHigh, .5, 1., args=(p0,p1,p2,p3,a))
    return lo + hi

def cubicToQuadratics(cubic):
    a = optimize.minimize(optimizationIntegral, [0.,0.], args=(cubic[0],cubic[1],cubic[2],cubic[3])).x
    bhalf = cubicBezier(.5, cubic[0],cubic[1],cubic[2],cubic[3])
    p1prime = [-a[0]+bhalf[0], -a[1]+bhalf[1]]
    p2prime = [bhalf[0]+a[0],bhalf[1]+a[1]]
    return [ [p0, p1prime, bhalf], [bhalf, p2prime, p3] ]
    
# Example cubic spline
p0 = [ 3., -5. ]
p1 = [ 2.1, 2.9 ]
p2 = [ 3.1, -0.2 ]
p3 = [ -2.0, -1.3 ]

fig = plt.figure()
t = np.arange(0., 1., 1.e-3)
curve = [ cubicBezier(ti, p0, p1, p2, p3) for ti in t ]
x = [ curve[i][0] for i in range(len(curve)) ]
y = [ curve[i][1] for i in range(len(curve)) ]
plt.plot(x, y, 'r-')
plt.plot([p0[0]], [p0[1]], 'ro')
plt.plot([p3[0]], [p3[1]], 'ro')
plt.plot([p1[0]], [p1[1]], 'bo')
plt.plot([p2[0]], [p2[1]], 'bo')

quadratics = cubicToQuadratics([p0,p1,p2,p3])

p1prime = quadratics[0][1]
p2prime = quadratics[1][1]
bhalf = quadratics[0][2]

lcurve = [ quadraticBezier(2.*ti, p0, p1prime, bhalf) for ti in np.arange(0., .5, 1.e-3) ]
lx = [ lcurve[i][0] for i in range(len(lcurve)) ]
ly = [ lcurve[i][1] for i in range(len(lcurve)) ]

plt.plot(lx, ly, 'g-') 

hcurve = [ quadraticBezier(2.*(ti-.5), bhalf, p2prime, p3) for ti in np.arange(.5, 1., 1.e-3) ]
hx = [ hcurve[i][0] for i in range(len(hcurve)) ]
hy = [ hcurve[i][1] for i in range(len(hcurve)) ]

plt.plot(hx, hy, 'g-') 

plt.plot([p1prime[0]], [p1prime[1]], 'yo')
plt.plot([p2prime[0]], [p2prime[1]], 'yo')

fig.show()

input()
