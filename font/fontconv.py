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

def rescale(point, xmin, xmax, ymin, ymax):
    ret = point - complex(xmin, ymin)
    ret = complex(ret.real/abs(xmax-xmin), ret.imag/(ymax-ymin)/100.*29.)
    ret -= complex(.5,.5*29./100.)
    ret = complex(ret.real, -ret.imag)
    return ret

def linearBezier(t, p0, p1):
    return (1.-t)*p0+t*p1

def quadraticBezier(t, p0, p1, p2):
    return linearBezier(t, linearBezier(t, p0, p1), linearBezier(t, p1, p2))

def cubicBezier(t, p0, p1, p2, p3):
    return linearBezier(t, quadraticBezier(t, p0, p1, p2), quadraticBezier(t, p1, p2, p3))

def approximateQuadratic(cubicBezier):
    quadraticBeziers = []
    
    
    return quadraticBeziers

ast = et.parse('princess-sofia-plain.svg')
root = ast.getroot()

f = open("font.gen.py", "wt")

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
                print("Line ", lines[-1])
                
            else:
                print("Unrecognized path control. Ignoring: ", content)
        
        # Convert the cubic array to quadratics. Do this by
        # - while sum of squared differences is too big:
        #   - add an on-curve point in the center
        #   - solve the optimization problem for the missing quad control point
        for cubic in cubicBeziers:
            #quadraticBeziers += [ 
            
        
        print("\n")
        
f.close()
