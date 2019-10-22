# Hardcyber - PC-64k-Intro by Team210 at Deadline 2k19
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

from svgpathtools import Path, Line
from svgpathtools import parse_path

path_string = 'M 30.519076,122.94061 V 87.151984 H 178.54706 v 35.788626 z m 143.690094,-4.72307 v -16.66371 h -11.93978 l 6.64809,6.88987 v 7.12767 h -6.64809 v -7.12767 h -7.22718 v 7.12767 h -7.15906 v -7.12767 h 7.15906 l -7.15906,-6.88987 v -7.149151 h 7.15906 v 7.149151 h 7.22718 v -7.149151 h 11.93978 v -2.910415 h -61.5156 v 2.910415 h 7.15911 v 7.149151 h -7.15911 v 6.88987 h 7.15911 v 7.12767 h -7.15911 l -7.14365,-7.12767 V 91.494264 H 35.170633 v 2.910415 h 10.31875 V 108.4437 h 7.2209 V 94.404679 h 7.06657 V 108.4437 l -7.06657,7.12767 h -7.2209 l -7.14378,-7.12767 v -6.88987 h -3.17497 v 16.66371 h 52.78435 V 94.404679 h 14.170427 v 7.149151 h -7.166567 v 16.66371 z m -110.992707,-16.66371 7.00411,-7.149151 h 14.16256 V 108.4437 l -7.01855,7.12767 h -14.14812 z m 14.14812,0 h -7.14401 v 6.88987 h 7.14401 z m 45.824107,14.01754 v -14.01754 l 6.81473,-7.149151 h 14.17555 v 7.149151 l -6.42779,6.88987 v 0 -6.88987 h -7.74776 v 6.88987 h 14.17555 v 7.12767 z'

path = parse_path(path_string)

# find dimensions
xmax = -1.e9
xmin = 1.e9
ymax = -1.e9
ymin = 1.e9
for line in path:
    xmax = max(xmax, line.start.real)
    xmax = max(xmax, line.end.real)
    
    xmin = min(xmin, line.start.real)
    xmin = min(xmin, line.end.real)
    
    ymax = max(ymax, line.start.imag)
    ymax = max(ymax, line.end.imag)
    
    ymin = min(ymin, line.start.imag)
    ymin = min(ymin, line.end.imag)

# rescale path
for i in range(len(path)):
    path[i].start -= complex(xmin,ymin)
    path[i].start = complex(path[i].start.real/abs(xmax-xmin), path[i].start.imag/abs(ymax-ymin)/100.*29.)
    path[i].start -= complex(.5,.5*29./100.)
    path[i].start = complex(path[i].start.real,-path[i].start.imag)
    
    path[i].end -= complex(xmin,ymin)
    path[i].end = complex(path[i].end.real/abs(xmax-xmin), path[i].end.imag/abs(ymax-ymin)/100.*29.)
    path[i].end -= complex(.5,.5*29./100.)
    path[i].end = complex(path[i].end.real,-path[i].end.imag)

# sort path
#newpath = [ path[0] ]
#del path[0]
#while len(path) > 1:
    #print(len(path))
    #for j in range(len(path)-1):
        #print(j,"/",len(path))
        #if abs(newpath[-1].end - path[j].start)<5.e-1:
            #newpath += [ path[j] ]
            #del path[j]
            #break
#path = newpath
    
with open('5711.frag', 'wt') as f:
    f.write('const int npts = ' + str(4*len(path)) + ';\n')
    f.write('const float path[npts] = float[npts](')
    
    for i in range(len(path)-1):
        line = path[i]
        f.write('{:.3f}'.format(line.start.real) + ',' + '{:.3f}'.format(line.start.imag) + ',')
        f.write('{:.3f}'.format(line.end.real) + ',' + '{:.3f}'.format(line.end.imag) + ',')
    line = path[-1]
    f.write('{:.3f}'.format(line.start.real) + ',' + '{:.3f}'.format(line.start.imag) + ',')
    f.write('{:.3f}'.format(line.end.real) + ',' + '{:.3f}'.format(line.end.imag))
    f.write(');\n')
    f.close()
