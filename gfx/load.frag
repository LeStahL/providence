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

#version 150

out vec4 gl_FragColor;

uniform vec2 iResolution;
uniform float iProgress;

const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);

const float fsaa = 4.;

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

//distance to spline with parameter t
float dist2(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t)
{
    t = clamp(t, 0., 1.);
    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);
}

//minimum dist3ance to spline
void dspline2(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float ds)
{
    //coefficients for 0 = t^3 + a * t^2 + b * t + c
    vec2 E = x-p0, F = p2-2.*p1+p0, G = p1-p0;
    vec3 ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);

	//discriminant and helpers
    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;
    
    //triple real root
    if(dis > 0.) 
    {
        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);
        ds = dist2(p0,p1,p2,x,ui.x+ui.y-tau);
        return;
    }
    
    //three dist3inct real roots
    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;
    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;
    ds = min(
        dist2(p0,p1,p2,x, t.x),
        min(
            dist2(p0,p1,p2,x,t.y),
            dist2(p0,p1,p2,x,t.z)
        )
    );
}

// iq's smooth minimum
void smoothmin(in float a, in float b, in float k, out float dst)
{
    float h = max( k-abs(a-b), 0.0 )/k;
    dst = min( a, b ) - h*h*h*k*(1.0/6.0);
}

void dsmoothvoronoi(in vec2 x, in float sm, out float d, out vec2 z)
{
    float n;
//     lfnoise(x-iTime*c.xy, n);
    
    vec2 y = floor(x);
       float ret = 1.;
    vec2 pf=c.yy, p;
    float df=10.;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            d = length(x-p);
            
            if(d < df)
            {
                df = d;
                pf = p;
            }
        }
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            vec2 o = p - pf;
            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
            smoothmin(ret, d, sm, ret);
        }
    
    d = ret;
    z = pf;
}

void dbutterflywing(in vec2 x, out float d)
{
    x.x -= .15-.1*abs(x.y);
    
    // Butterfly main structure
    dspline2(x, vec2(.6,.4), vec2(.0,.4), c.yy, d);
  	
    float da;
    dspline2(x, c.yy, vec2(.0,-.2), vec2(.1,-.2), da);
    d = min(d, da);
    
    dspline2(x, vec2(.1,-.2), vec2(.1,-.2), vec2(.4,-.4), da);
    d = min(d, da);
    
    // Wing thickness
    da = .15*smoothstep(-.1,.2,x.y)*smoothstep(.4,.2,x.y)*cos(x.y*pi);
    da += .2*smoothstep(.1,-.2,x.y)*smoothstep(-.4,0.,x.y)*cos(x.y*pi);
    
    d = abs(d)-da-.02*smoothstep(.4,.3,abs(x.y));
    d *= 2.;
}

float sm(in float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void butterflytexture(in vec2 x, out vec3 col)
{
    float size = 4.;
    float d, v;
    vec2 vi;
    vec2 y = vec2(length(x), atan(x.y,x.x));
    dsmoothvoronoi(size*y, .5, v, vi);
    vi /= size;
    v /= size;
    rand(vi, d);
    col = mix(mix(vec3(0.01,0.02,0.59), vec3(0.49,0.00,0.87), d), c.yyy, sm(abs(v)-.03));
    col = mix(col, 2.*col, sm(1.e-3*v));
    vec3 dis;
    rand(x, dis.x);
    rand(x+1337., dis.y);
    rand(x+2337., dis.z);
    col += .2*dis;
    col = mix(col, length(col)/sqrt(3.)*c.xxx, .5);
}

void det3(in mat3 m, out float d)
{
//#ifdef determinant
    d = determinant(m);
   /*
#else
    d = m[0][0]*m[1][1]*m[2][2] + m[1][0]*m[2][1]*m[0][2] + m[2][0]*m[0][1]*m[2][2]
      - m[2][0]*m[1][1]*m[0][2] - m[1][0]*m[0][1]*m[2][2] - m[0][0]*m[2][1]*m[1][2];
#endif
*/
}

void det2(in mat2 m, out float d)
{

//#ifdef determinant
    d = determinant(m);
/*
#else
    d = m[0][0]*m[1][1] - m[1][0]*m[0][1];
#endif
*/
}

void inv(in mat3 m, out mat3 dst)
{
//#ifdef inverse
    dst = inverse(m);
/*
#else
    float d;
    det3(m, d);
    det2(mat2(m[1][1],m[2][1],m[1][2],m[2][2]), dst[0][0]);
    det2(mat2(m[0][2],m[2][2],m[0][1],m[2][1]), dst[0][1]);
    det2(mat2(m[0][1],m[1][1],m[0][2],m[1][2]), dst[0][2]);
    det2(mat2(m[
#endif
*/
}

void analytical_plane(in vec3 o, in vec3 dir, in vec3 p0, in vec3 p1, in vec3 p2, out float d, out vec2 uv)
{
    mat3 m = mat3(p1, p2, -dir),
        minv;
    det3(m, d);
    if(d == 0.) d = 1.e3;
    else
    {
        inv(m, minv);
        vec3 ltd = minv*(o-p0);
        d = ltd.z;
        uv = ltd.xy;
    }
}

void append_butterfly(in vec2 uv, inout vec3 col)
{ 
vec3 o0 = .75*c.yyx+.3*c.yzy,
        o = o0,
        r = c.xyy,
        t = c.yyy, 
        u = cross(normalize(t-o),-r),
        dir,
        n, 
        x,
        c1 = c.yyy,
        l;
    int N = 250,
        i;
    float d;
    
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);
	
    float ta = .7*sign(uv.x)*abs(pow(sin(pi* /*.3*/ 2.*iProgress),2.));
    mat2 R = mat2(cos(ta), sin(ta), -sin(ta), cos(ta));
    vec2 r1 = R*vec2(1.,0.),
        r2 = R*vec2(-1.,0.);
    
    float dp1, dp2;
    vec2 uvp1, uvp2;
    vec3 p1 = normalize(vec3(r1.x,0.,r1.y)),
        p2 = normalize(vec3(r2.x,0.,r2.y));
    analytical_plane(o, dir, c.yyy, p1, c.yxy, dp1, uvp1);
    analytical_plane(o, dir, c.yyy, p2, c.yxy, dp2, uvp2);
    float mat = -1.;
    
    if(dp1<dp2)
    {
        d = dp1;
        mat = 0.;
    }
    else
    {
        d = dp2;
        mat = 1.;
    }
    if(d > 1.e2) mat = -1.;
    
    x = o + d * dir;
    n = abs(cross(p1, c.yxy));
	l = normalize(.5*c.yyx);
    
    float db;
    
    
    if(mat == 0.)
    {
        butterflytexture(uvp2,c1);
        dbutterflywing(vec2(abs(uvp2.x),uvp2.y), db);
        c1 = mix(c1, .2*length(c1)/sqrt(3.)*c.xxx, step(-.4+.8*iProgress, uvp2.y));
        col = mix(col, c.yyy, sm(db));
    	db += .04;
        col = mix(col, c1, sm(db));
    }
    else if(mat == 1.)
    {
        dbutterflywing(vec2(abs(uvp1.x),uvp1.y), db);
        butterflytexture(uvp1,c1);
        c1 = mix(c1, .2*length(c1)/sqrt(3.)*c.xxx, step(-.4+.8*iProgress, uvp1.y));
        col = mix(col, c.yyy, sm(db));
    	db += .04;
        col = mix(col, c1, sm(db));
    }
    
    col = mix(col, 12.*vec3(0.99,0.44,0.37), sm((abs(db-.063)-.01))/12.);
    col = mix(col, 22.*c.xxy, sm((abs(db-.063)-.001))/12.);
}

// Creative Commons Attribution-ShareAlike 4.0 International Public License
// Created by David Hoskins.
// See https://www.shadertoy.com/view/4djSRW
void hash22(in vec2 p, out vec2 d)
{
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    d = fract((p3.xx+p3.yz)*p3.zy);
}

void dist(in vec2 a, in vec2 b, out float d)
{
    d = length(b-a);
}

void nearest_controlpoint(in vec2 x, out vec2 p)
{
    float dmin = 1.e5, 
        d;
    vec2 dp,
        y = floor(x);
    
    float i = 0.;
    for(float i = -1.; i <= 1.; i += 1.)
        for(float j = -1.; j <= 1.; j += 1.)
        {
            hash22(y+vec2(i,j), dp);
            dp += y+vec2(i,j);
            dist(x, dp, d);
            if(d<dmin)
            {
                dmin = d;
                p = dp;
            }
        }
}

void dvoronoi(in vec2 x, out float d, out vec2 p, out float control_distance)
{
    d = 1.e4;
    vec2 y,
        dp;
    
    nearest_controlpoint(x, p);
    y = floor(p);
    
    control_distance = 1.e4;
    
    for(float i = -2.; i <= 2.; i += 1.)
        for(float j = -2.; j <= 2.; j += 1.)
        {
            if(i==0. && j==0.) continue;
            hash22(y+vec2(i,j), dp);
            dp += y+vec2(i,j);
            vec2 o = p - dp;
            float l = length(o);
            d = min(d,abs(.5*l-dot(x-dp,o)/l));
            control_distance = min(control_distance,.5*l);
        }
}

void lfnoise(in vec2 t, out float n)
{
    vec2 i = floor(t);
    t = fract(t);
    t = smoothstep(c.yy, c.xx, t);
    vec2 v1, v2;
    rand(i, v1.x);
    rand(i+c.xy, v1.y);
    rand(i+c.yx, v2.x);
    rand(i+c.xx, v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    n = mix(v1.x, v1.y, t.x);
}

void background_geraffel(in vec2 uv, inout vec3 col)
{
	float vsize = 1.;
    float v, vc, d, r, n;
    vec2 vi, x;
    
    for(float i = 0.; i < 6.; i += 1.)
    {
        lfnoise(i*c.xx, n);
        rand(i*c.xx, vi.x);
        rand(i*c.xx+1337., vi.y);
        x = uv +(-1.+2.*vi)*vsize;
        dvoronoi(x*vsize, v, vi, vc);
        v /= vsize;
        d = length(x-vi/vsize)-.5*vc/vsize;
        rand(vi, r);
        col = mix(col, mix(col, mix(vec3(0.42,0.05,0.95),c.xxx,r*(.5+.5*n)), 2.*step(.5,r)*r), sm(d/22.));
        vsize += .4;
    }
    float phi = atan(uv.y, uv.x);
    lfnoise(150.*phi*c.xx, n);
    col = mix(col, vec3(0.69,0.01,0.37), 1.-.5*length(uv)-.1*abs(n));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord-.5*iResolution.xy)/iResolution.yy;
    uv *= 3.;
    vec3 col = vec3(0.04,0.02,0.29);
    background_geraffel(uv, col);
    append_butterfly(uv, col);
    
    col = mix(col,col*col, col);
    col = clamp(col,0.,1.);
    
    // Scan lines
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.);

    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    vec4 col = vec4(0.);
    float bound = sqrt(fsaa)-1.;
   	for(float i = -.5*bound; i<=.5*bound; i+=1.)
        for(float j=-.5*bound; j<=.5*bound; j+=1.)
        {
            vec4 c1;
            mainImage(c1, gl_FragCoord.xy+vec2(i,j)*1.5/max(bound, 1.));
     		col += c1;
        }
    col /= fsaa;
    gl_FragColor = col;
}

