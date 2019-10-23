/* 4k 
 * Copyright (C) 2019 Alexander Kraus <nr4@z10.info>
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
 
#version 130

out vec4 gl_FragColor;

uniform float iTime;
uniform vec2 iResolution;

uniform float iFader0;
uniform float iFader1;
uniform float iFader2;
uniform float iFader3;
uniform float iFader4;
uniform float iFader5;
uniform float iFader6;
uniform float iFader7;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;
const float ra = .5;
float gTime;

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
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

void hash12(in vec2 p, out float d)
{
	vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    d = fract((p3.x + p3.y) * p3.z);
}
// End of David Hoskins' code

void lfnoise(in vec2 t, out float n)
{
    vec2 i = floor(t);
    t = fract(t);
    t = smoothstep(c.yy, c.xx, t);
    vec2 v1, v2;
    rand(i, v1.x);
    rand((i+c.xy), v1.y);
    rand((i+c.yx), v2.x);
    rand((i+c.xx), v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    n = mix(v1.x, v1.y, t.x);
}
void mfnoise(in vec2 x, in float d, in float b, in float e, out float n)
{
    n = 0.;
    float a = 1., nf = 0., buf;
    for(float f = d; f<b; f *= 2.)
    {
        lfnoise(f*x, buf);
        n += a*buf;
        a *= e;
        nf += 1.;
    }
    n *= (1.-e)/(1.-pow(e, nf));
}


void dbox3(in vec3 x, in vec3 b, out float d)
{
  vec3 da = abs(x) - b;
  d = length(max(da,0.0))
         + min(max(da.x,max(da.y,da.z)),0.0);
}

void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = (sda.x<sdb.x)?sda:sdb;
}

void sub(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = (sda.x>sdb.x)?abs(sda):abs(sdb)*c.zx;
}

void dist(in vec2 a, in vec2 b, out float d)
{
    d = length(b-a);
}

void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

mat3 gR;
vec2 ind = c.yy;
void scene(in vec3 x, out vec2 sdf)
{
    x.y += .3*gTime;
    
    float nlo,
        nhi;
    
    lfnoise(1.3*x.xy, nlo);
    nlo = .5+.5*nlo;
    
    lfnoise(3.3*x.xy, nhi);
    nlo -= .1*nhi;
    
    mfnoise(x.xy, 80., 800., .3, nhi);
    
    float nd;
    lfnoise(1.e-1*x.y*c.xx, nd);
    
    sdf = vec2(x.z+(1.*nlo-.1*abs(nhi))+2.*smoothstep(.7,-.3,abs(x.x-1.2*nd)), 1.);
    sdf.x += .1;
    
    add(sdf, vec2(x.z+1., 3.), sdf);
}

void normal(in vec3 x, out vec3 n, in float dx)
{
    vec2 s, na;
    
    scene(x,s);
    scene(x+dx*c.xyy, na);
    n.x = na.x;
    scene(x+dx*c.yxy, na);
    n.y = na.x;
    scene(x+dx*c.yyx, na);
    n.z = na.x;
    n = normalize(n-s.x);
}

float sm(in float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void analytical_box(in vec3 o, in vec3 dir, in vec3 size, out vec2 d)
{
    vec3 tlo = min((size-o)/dir,(-size-o)/dir),
        thi = max((size-o)/dir,(-size-o)/dir);
    vec2 abxlo = abs(o.yz + tlo.x*dir.yz),
        abylo = abs(o.xz + tlo.y*dir.xz),
        abzlo = abs(o.xy + tlo.z*dir.xy),
        abxhi = abs(o.yz + thi.x*dir.yz),
        abyhi = abs(o.xz + thi.y*dir.xz),
        abzhi = abs(o.xy + thi.z*dir.xy);
    vec4 dn = 1.e4*c.xyyy;
    
    dn = mix(dn, vec4(tlo.x,c.xyy), float(all(lessThan(abxlo,size.yz)))*step(tlo.x,dn.x));
    dn = mix(dn, vec4(tlo.y,c.yxy), float(all(lessThan(abylo,size.xz)))*step(tlo.y,dn.x));
    dn = mix(dn, vec4(tlo.z,c.yyx), float(all(lessThan(abzlo,size.xy)))*step(tlo.z,dn.x));

    d.x = dn.r;
    
    dn = 1.e4*c.xyyy;
    dn = mix(dn, vec4(thi.x,c.xyy), float(all(lessThan(abxhi,size.yz)))*step(thi.x,dn.x));
    dn = mix(dn, vec4(thi.y,c.yxy), float(all(lessThan(abyhi,size.xz)))*step(thi.y,dn.x));
    dn = mix(dn, vec4(thi.z,c.yyx), float(all(lessThan(abzhi,size.xy)))*step(thi.z,dn.x));
    
    d.y = dn.r;
}

void analytical_sphere(in vec3 o, in vec3 dir, in float R, out vec2 d)
{
    float a = dot(dir,dir),
        b = 2.*dot(o,dir),
        cc = dot(o,o)-R*R,
        dis = b*b-4.*a*cc;
    vec2 dd = (dis<0.)?1.e4*c.xx:(c.xz*sqrt(dis)-b)/2./a;
    d = vec2(min(dd.x, dd.y), max(dd.x, dd.y));
}

void analytical_plane(in vec3 o, in vec3 dir, in vec3 plane_normal, out float d)
{
    d = -dot(o,plane_normal)/dot(dir,plane_normal);
    if(d<0.) d = 1.e4;
}

void evening_sky(in vec2 uv, out vec3 col)
{
    col = mix(vec3(0.97,0.85,0.22).rbg, vec3(1.00,0.49,0.25), uv.y);
    col = sqrt(abs(col));
    
    float n, hf_parts, na;

    lfnoise(uv, n);
    mfnoise(uv, 2.,15.e3, .45, hf_parts);
    n = mix(n, hf_parts, .5);
    
    col = mix(col, vec3(0.87,0.33,0.48), clamp(abs(.5+1.9*n)-.4,0.,1.));
//     col = mix(col, c.yyy, clamp(abs(-.5+1.9*n)-.3,0.,1.));
//     float na;
//     lfnoise(3.*uv, na);
//     mfnoise(uv, 3.,15.e3, .45, hf_parts);
//     n = mix(n, hf_parts, .5);
//     n -= .2*na;
    lfnoise(uv+.03*c.yx, n);
    mfnoise(uv+.03*c.yx, 2.,15.e3, .65, hf_parts);
    n = mix(n, hf_parts, .5);

    col = mix(col, vec3(0.45,0.14,0.04), clamp(abs(-.2+.9*n)-.3,0.,1.));

//     lfnoise(3.*uv, na);
//     n += .5*na;
    
    lfnoise(uv+.06*c.yx, n);
    mfnoise(uv+.06*c.yx, 2.,15.e3, .55, hf_parts);
    n = mix(n, hf_parts, .5);

    col = mix(col, c.xxy, clamp(abs(-.3+.9*n)-.2,0.,1.));
    col = mix(col, vec3(.01,.05,.13), clamp(abs(-.3+.9*n)-.05,0.,1.));
    
    col = mix(col, .75*vec3(0.42,0.27,0.22), sm(uv.y-.3+.3*n));
    //col = mix(col, vec3(0.11,0.05,0.11), clamp(.7+.3*n,0.,1.));
    
    col *= col;
}

void needle_tree(in vec2 uv, inout vec3 col)
{
    float d, dx;
    
    lfnoise(4.*uv.y*c.xx, dx);
    
    // Stem
    dlinesegment(uv-.02*dx*c.xy, vec2(0.,-.5), vec2(0.,.5), d);
    d = abs(d)-mix(.001,.03,.5-uv.y);
    
    col = mix(col, c.yyy, sm(d));
}

void mainImage(out vec4 fragColor, in vec2 fragCoord, in float t_)
{
    // Set up global variables
//     rot3(vec3(1.1,1.3,1.5)*iTime, gR);
    
    gTime = t_;
    
    // Set up coordinates and camera
    vec2 uv = (fragCoord.xy-.5*iResolution.xy)/iResolution.y,
        s,
        dsky_sphere;
    uv *= .5;
    vec3 col = c.yyy,
        o0 = .1*c.yyx+1.*c.yzy,
        o = o0,
        r = c.xyy,
        t = c.yyy, 
        u = cross(normalize(t-o),-r),
        dir,
        n, 
        x,
        c1 = c.yyy,
        l;
    int N = 450,
        i;
    float d,
        dlower,
        dupper;
    
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);
    
    analytical_sphere(o, dir, 20., dsky_sphere);
    analytical_plane(o, dir, c.yyx, dupper);
    analytical_plane(o+1.*c.yyx, dir, c.yyx, dlower);
    
    vec2 first_hit = vec2(dupper, 0.);
    add(first_hit, vec2(dsky_sphere.y, 1.), first_hit);
    
    d = first_hit.x;

    if(first_hit.y == 0.) // Floor
    {
        // Raymarch mountains
        for(i = 0; i<N; ++i)
        {
            x = o + d * dir;
            scene(x,s);
            if(s.x < 1.e-3) break;
//             d += s.x<1.e-2?min(s.x,1.e-2):s.x;
            d += min(s.x, 5.e-2);
        }
        
        if(i<N)
        {
            normal(x, n, 8.e-2);
            l = normalize(x + .5*n);
          
            if(s.y == 3.)
            {
                col = c.yxy;
            }
            else
            {
                col = mix(vec3(0.18,0.07,0.03), vec3(0.19,0.16,0.18), x.z);
                col = .1*col 
                    + .2*col*dot(l, n)
                    + .4*col*pow(abs(dot(reflect(l,n),dir)),2.);
            }
        }
        
        float na;
        lfnoise(128.*x.xy, na);
        
        col = mix(col, mix(vec3(0.63,0.24,0.13), c.yxy, .15+.15*na), dot(n,c.yxx));
        col = mix(col, c.yyy, dot(n,c.yzy));
        col = mix(col, vec3(0.42,0.27,0.22), tanh(.05 *(d+5.*abs(x.z))));
        col *= .5*col;
        col = mix(col, .5*vec3(0.89,0.66,0.47), abs(x.z*x.z));
        
//         col += vec3(0.42,0.27,0.22)
//         col = mix(col.grb, col, clamp(iTime/2.,0.,1.));
    }
    else if(first_hit.y == 1.) // Sky
    {
        x = o + d * dir;

        vec2 y = vec2(atan(x.y, x.x), x.z/length(x));
        evening_sky(vec2(4.,20.)*y, col);

    }
    
    // Add sun
    vec2 z = uv-.15*c.yx;
    float dsun = length(z) - .01;
    col = mix(col,5.*vec3(0.95,0.41,0.11), sm(dsun/(100.+20.*sm(1.e-3*sin(3.82*pi*atan(z.y,z.x))))));
    col = mix(col,3.*c.xxx, sm((dsun+.01)/20.));

    // Add needle tree TODO
//     needle_tree(uv, col);
    
//     z = mod(uv,.05)-.025;
//     z = abs(z)-.0001;
//     col = mix(col,mix(col, c.xxx, .2), sm(min(z.x,z.y)));
//     
//     z = mod(uv,.01)-.005;
//     z = abs(z)-.00005;
//     col = mix(col,mix(col, c.xxx, .2), sm(min(z.x,z.y)));
    
    col = mix(col, vec3(0.71,0.95,0.38), .1);
    
//     col *= 5.*col;
    
    fragColor = vec4(clamp(col,0.,1.),1.);
}

void main()
{
    vec4 cc = c.yyyy;
    for(float dt = -.05; dt <= .05; dt += .05)
    {
        mainImage(gl_FragColor, gl_FragCoord.xy, iTime + dt);
        cc += gl_FragColor;
    }
    cc /= 3.;
    gl_FragColor = cc;
}




