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

path_string = 'm -0.50195312,277.18945 v 37.79493 H 178.25781 L 104.58789,493.5918 h 42.17773 L 219.2168,314.98438 h 42.36132 l 16.91993,-37.79493 z m 291.34570312,0 -17.67773,37.79493 H 452.5918 l -5.12305,12.92578 h -71.64844 l -66.54883,165.65039 h 113.17969 l 85.77735,-216.3711 z m 538.08398,0 -310.25195,0.25391 -15.15234,38.55273 270.12109,-1.01171 -69.57422,178.57617 h 92.82031 l 42.99219,-106.46875 -35.05078,-14.14453 -33.43945,82.81835 h -12.03906 z m 12.98438,0 -34.31641,84.79102 35.24219,13.64844 24.96289,-60.64453 h 13.78906 L 809.8125,493.5918 h 40.45117 l 87.14844,-216.40235 z m 109.39258,0 -86.98047,216.3711 h 27.96094 66.89254 l 88.709,-216.3711 z m 25.54101,37.79493 h 14.6953 l -57.7187,140.78124 h -13.57425 z m -747.36523,12.92578 -66.03711,165.65039 129,0.0293 15.71484,-37.82422 h -88.96093 l 35.90625,-90.06054 h 49.85742 l -8.14258,19.92383 -31.77734,0.70507 -13.3125,37.78125 70.24023,-0.13671 39.2793,-96.06836 z m 268.90039,0 -67.39844,165.68164 h 41.875 l 51.04297,-127.88477 h 69.23438 l -24.87891,64.4668 39.2832,-0.73047 25.86524,-63.73633 h 67.23242 L 653.92969,493.5918 h 36.99023 l 65.08203,-165.68164 z m -97.01367,37.79687 h 31.11523 l -35.70508,90.05859 h -31.58789 z'

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
    
with open('team210.frag', 'wt') as f:
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
