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

path_string = 'M 216.07227 508.76953 L 216.07227 563.94922 L 228.10742 563.94922 L 228.10742 544.30273 L 231.88672 544.30273 L 231.88672 563.94922 L 254.71875 563.94922 L 254.71875 508.76953 L 245.19727 508.76953 L 245.19727 527.16211 L 241.41797 527.16211 L 241.41797 508.76953 L 216.07227 508.76953 z M 258.49805 508.76953 L 258.49805 563.94922 L 268.02148 563.94922 L 268.02148 544.30273 L 271.80078 544.30273 L 271.80078 563.94922 L 307.92773 563.94922 L 307.92773 546.19336 L 294.62695 546.19336 L 294.62695 537.62305 L 277.38281 531.70312 L 294.62695 525.8125 L 294.62695 508.76953 L 258.49805 508.76953 z M 298.40625 508.76953 L 298.40625 528.51367 L 289.05273 531.71094 L 298.40625 534.92188 L 298.40625 542.41406 L 311.70703 542.41406 L 311.70703 563.94922 L 330.89258 563.94922 L 330.89258 508.76953 L 311.70703 508.76953 L 311.70703 527.16211 L 307.92773 527.16211 L 307.92773 508.76953 L 298.40625 508.76953 z M 268.02148 520.55664 L 271.80078 520.55664 L 271.80078 527.16211 L 268.02148 527.16211 L 268.02148 520.55664 z'

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
    
with open('nr4sketch.frag', 'wt') as f:
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
