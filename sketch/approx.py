import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def linearBezier(t, p0, p1):
    return [ (1.-t)*p0[0]+t*p1[0], (1.-t)*p0[1]+t*p1[1] ]

def quadraticBezier(t, p0, p1, p2):
    return linearBezier(t, linearBezier(t, p0, p1), linearBezier(t, p1, p2))

def cubicBezier(t, p0, p1, p2, p3):
    return linearBezier(t, quadraticBezier(t, p0, p1, p2), quadraticBezier(t, p1, p2, p3))

def optimizationIntegral(p0, p1, p2, p3):
    

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

fig.show()

input()
