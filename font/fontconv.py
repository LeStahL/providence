# Providence by Team210 - 64k intro by Team210 at Vortex 2k19
# Copyright (C) 2019  Alexander Kraus <nr4@z10.info>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import xml.etree.ElementTree as et
from svgpathtools import Path, Line, CubicBezier, parse_path
import itertools
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
    return [ [cubic[0], p1prime, bhalf], [bhalf, p2prime, cubic[3]] ]

def rescale(point, xmin, xmax, ymin, ymax):
    ret = point - complex(xmin, ymin)
    ret = complex(ret.real/abs(xmax-xmin), ret.imag/(ymax-ymin)/100.*29.)
    ret -= complex(.5,.5*29./100.)
    ret = complex(ret.real, -ret.imag)
    return ret

ast = et.parse('princess-sofia-plain.svg')
root = ast.getroot()

f = open("fontgen.py", "wt")
f.write("import numpy\n")
f.write("def glyph(char):\n")
f.write("    quads = []\n")

# Assume, that there is only one document group level present.
for g in root.findall('{http://www.w3.org/2000/svg}g'):
    for path in g:
        id = path.attrib['id'].replace('path-', '')
        print("Letter has id: ", id)
        
        d = path.attrib['d']
        
        cpath = parse_path(d)
        
        cubicBeziers = []
        quadraticBeziers = []
        lines = []
        
        xmax = -1.e9
        xmin = 1.e9
        ymax = -1.e9
        ymin = 1.e9
        
        # Find the dimensions of the canvas.
        for content in cpath:
            xmax = max(xmax, content.start.real)
            xmax = max(xmax, content.end.real)
            xmin = min(xmin, content.start.real)
            xmin = min(xmin, content.end.real)
            ymax = max(ymax, content.start.imag)
            ymax = max(ymax, content.end.imag)
            ymin = min(ymin, content.start.imag)
            ymin = min(ymin, content.end.imag)

            if(type(content).__name__ == "CubicBezier"):
                xmax = max(xmax, content.control1.real)
                xmax = max(xmax, content.control2.real)
                xmin = min(xmin, content.control1.real)
                xmin = min(xmin, content.control2.real)
                ymax = max(ymax, content.control1.imag)
                ymax = max(ymax, content.control2.imag)
                ymin = min(ymin, content.control1.imag)
                ymin = min(ymin, content.control2.imag)
                
            elif(type(content).__name__ == "QuadraticBezier"):
                xmax = max(xmax, content.control1.real)
                xmin = min(xmin, content.control1.real)
                ymax = max(ymax, content.control1.imag)
                ymin = min(ymin, content.control1.imag)
        
        print("Letter has x dimensions: ",xmin,xmax)
        print("Letter has y dimensions: ",ymin,ymax)
        
        # Assign path content to the proper arrays.
        for content in cpath:
            start = rescale(content.start, xmin, xmax, ymin, ymax)
            end = rescale(content.end, xmin, xmax, ymin, ymax)
            
            if(type(content).__name__ == "CubicBezier"):
                control1 = rescale(content.control1, xmin, xmax, ymin, ymax)
                control2 = rescale(content.control2, xmin, xmax, ymin, ymax)
                
                cubicBeziers += [ [start, control1, control2, end] ]
                print("Cubic Bezier ", cubicBeziers[-1])
                
            elif(type(content).__name__ == "QuadraticBezier"):
                control1 = rescale(content.control1, xmin, xmax, ymin, ymax)
                
                quadraticBeziers += [ [start, control1, end] ]
                print("Quadratic Bezier ", quadraticBeziers[-1])
                
            elif(type(content).__name__ == "Line"):
                lines += [ [start, end] ]
                #print("Line ", lines[-1])
                
            else:
                print("Unrecognized path control. Ignoring: ", content)
                
        
        for cubic in cubicBeziers:
            cubicarray = [ [cubic[0].real, cubic[0].imag], [cubic[1].real, cubic[1].imag], [cubic[2].real, cubic[2].imag], [cubic[3].real, cubic[3].imag] ]
            quadratics = cubicToQuadratics(cubicarray)
            quadraticBeziers += [[ complex(quadratics[0][0][0], quadratics[0][0][1]), complex(quadratics[0][1][0], quadratics[0][1][1]), complex(quadratics[0][2][0], quadratics[0][2][1]) ]]
            quadraticBeziers += [[ complex(quadratics[1][0][0], quadratics[1][0][1]), complex(quadratics[1][1][0], quadratics[1][1][1]), complex(quadratics[1][2][0], quadratics[1][2][1]) ]]
        
        # Replace magic character ids with actual characters
        if id == 'comma': id = ','
        elif id == 'dot': id = '.'
        elif id == 'xmark': id = '!'
        elif id == 'qmark': id = '?'
        
        f.write("    if char == '" + id + "':\n")
        f.write("        quads += [ ")
        for i in range(len(quadraticBeziers)-2):
            quadratic = quadraticBeziers[i]
            #print(quadratic)
            f.write("[" + str(quadratic[0].real) + "," + str(quadratic[0].imag) + "," + str(quadratic[1].real) + "," + str(quadratic[1].imag) + "," + str(quadratic[2].real) + "," + str(quadratic[2].imag) + "],");
        quadratic = quadraticBeziers[len(quadraticBeziers)-1]
        f.write("[" + str(quadratic[0].real) + "," + str(quadratic[0].imag) + "," + str(quadratic[1].real) + "," + str(quadratic[1].imag) + "," + str(quadratic[2].real) + "," + str(quadratic[2].imag) + "] ]\n");
        
        print("\n")

f.write("    return [ quads ]\n")
f.write("def pack_length(char):\n")
f.write("    gly = glyph(char)\n")
f.write("    quads = gly[0]\n")
f.write("    return 1 + len(quads)*6\n")
        
f.close()
