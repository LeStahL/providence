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

ast = et.parse('princess-sofia-plain.svg')
root = ast.getroot()

for g in root.findall('{http://www.w3.org/2000/svg}g'):
    for path in g:
        id = path.attrib['id'].replace('path-', '')
        d = path.attrib['d']
        
        cpath = parse_path(d)
        for content in cpath:
            print(content)
