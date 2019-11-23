/* Hardcyber - PC-64k-Intro by Team210 at Deadline 2k19
 * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#version 130

out vec4 gl_FragColor;

uniform vec2 iResolution;
uniform float iProgress;

const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);

void dbox(in vec2 x, in vec2 b, out float d);
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);
void stroke(in float d0, in float s, out float d);
void rand(in vec2 x, out float n);
void lfnoise(in vec2 t, out float n);
void mfnoise(in vec2 x, in float d, in float b, in float e, out float n);
void hash22(in vec2 p, out vec2 d);
void dvoronoi(in vec2 x, out float d, out vec2 p, out float control_distance);

void dleaf(in vec2 x, in float N, in float R, out float d)
{
    float phi = atan(x.y, x.x);
    d = length(x) - R*cos(N*phi);
    d *= length(x)*length(x)/R/R;
}

void leafcolor(in vec2 x, in float d, out vec3 col)
{
    float n;
    float phi = atan(x.y, x.x);
    mfnoise(phi*c.xx, 16., 400., .65, n);
    col = mix(vec3(0.97,0.97,0.97), vec3(0.79,0.80,0.78), .8*n+.5+.5*sin(1.e1*(phi+2.*6.*d)));
}

float sm(in float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void addflower(in vec2 uv, inout vec3 col)
{
    float n;
    lfnoise(1.5*uv.x*c.xx, n);
    vec3 c1;
    float d;
    
    // TODO: add center of flower
    // TODO: add lights in background or other flowers
    // TODO: make stem more realistic
    // TODO: maybe add leaves to flower
    
    // Stem
    lfnoise(1.5*uv.y*c.xx, n);
    d = mix(abs(uv.x-.1*n-.075)-.02,1.,step(0.,uv.y));
    col = mix(col, mix(vec3(0.12,0.21,0.11), mix(vec3(0.37,0.60,0.47), vec3(0.95,0.92,0.39), uv.x+.5), clamp((uv.x+.01-.1*n-.075)/.02,0.,1.)), sm(d));
    
    // Leaves
    for(float i = 0.; i < 5.; i += 1.)
    {
        float phi = pi/16.*i;
        rand(phi*c.xx, phi);
        phi *= pi;
        mat2 R = mat2(cos(phi), sin(phi), -sin(phi), cos(phi));
	    dleaf(R*uv, 6., .3, d);
    	leafcolor(R*uv, d, c1);
        col = mix(col, .8*col, sm((d-.001)/12.));
    	col = mix(col, c1, sm(d));
    }
    
    // Center of flower
    const float size = 128.;
    float phi = atan(uv.y, uv.x);
    lfnoise(5.*phi*c.xx, n);
    float da = length(uv)-.075-.01*abs(n),
        v, vc;
    vec2 vi;
    col = mix(col, vec3(0.91,0.84,0.12), sm(da));
    dvoronoi(size*uv, v, vi, vc);
    vc /= size;
    vi /= size;
    vec3 o = c.yyx,
        dir = normalize(vec3(uv-vi,0.)-o);
    c1 = 1.4*vec3(0.81,0.74,0.07);
    c1 = mix(c1, 1.5*vec3(0.55,0.36,0.00), clamp(length(uv-vi)/.01,0.,1.));
    c1 = mix(c1, .8*c1*c1, clamp(length(uv)/.075,0.,1.));
    col = mix(col, .5*col, sm(da/25.));
    col = mix(col, c1, sm(mix(length(uv-vi)-.01,1.,sm(-da))));
}

void background(in vec2 x, out vec3 col)
{
    x.x /= .5+x.y;
    col = c.yyy;
    
    col = mix(vec3(0.71,0.28,0.20), vec3(0.04,0.55,0.87),clamp(x.y,0.,1.));
    
    float n;
    mfnoise(x*vec2(1.,3.),4.,400.,.45, n);
    col = mix(col, vec3(0.65,0.22,0.95), .3+.5*n);
    float na = abs(n)-.01;
    col = mix(col, vec3(0.76,0.41,0.97), sm(na/512.));
    mfnoise(x*vec2(1.,3.),1.,400.,.55, n);
    na = abs(n+.3)-.01;
    col = mix(col, vec3(0.72,0.20,0.22), sm(na/256.));
    
    mfnoise(x*vec2(1.,3.), 8., 500., .45, n);
    col = mix(col, .4*col, .5+.5*n);
    
    col = mix(col, 2.*vec3(0.87,0.54,0.31), clamp(1.-length(x*vec2(.35,2.))/1.,0.,1.));
    
    col = sqrt(col);
}

void addwindow(in vec2 uv, inout vec3 col, in vec2 dimensions);
void addprogressbar(in vec2 uv, inout vec3 col, in vec2 dimensions, in float progress);

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = 1.5*(fragCoord-.5*iResolution.xy)/iResolution.yy;
    
    float n;
    lfnoise(1.5*uv.x*c.xx, n);
    
    vec3 col = c.yyy;
	if(uv.y > 0.) background(uv, col);
    
    for(float i=20.; i>= 3.; i -= 1.)
    {
        float na, nb;
        rand(i*c.xx, na);
        float x = mod(.3*i*(uv.x-na), 4./i)-2./i;
        rand((x-.3*i*(uv.x-na))*c.xx, nb);
    	addflower(vec2(x, .3*i*(uv.y+.3*nb)), col);
    }
    
    col *= col*col;
    col = clamp(col, 0., 1.);
    
    addwindow(uv, col, vec2(.3,.2));
    addprogressbar(uv, col, vec2(.28,.02), iProgress);
    
    // Scan lines
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.);
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
