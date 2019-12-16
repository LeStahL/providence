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

uniform float iFSAA;
uniform vec2 iResolution;
uniform sampler2D iChannel0;
uniform float iTime;

out vec4 gl_FragColor;

const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);

const int npts = 284;
const float path[npts] = float[npts](-0.500,0.145,-0.500,0.094,-0.500,0.094,-0.329,0.094,-0.329,0.094,-0.400,-0.145,-0.400,-0.145,-0.360,-0.145,-0.360,-0.145,-0.290,0.094,-0.290,0.094,-0.250,0.094,-0.250,0.094,-0.234,0.145,-0.234,0.145,-0.500,0.145,-0.222,0.145,-0.239,0.094,-0.239,0.094,-0.068,0.094,-0.068,0.094,-0.073,0.077,-0.073,0.077,-0.141,0.077,-0.141,0.077,-0.205,-0.145,-0.205,-0.145,-0.097,-0.145,-0.097,-0.145,-0.015,0.145,-0.015,0.145,-0.222,0.145,0.291,0.145,-0.005,0.145,-0.005,0.145,-0.019,0.093,-0.019,0.093,0.238,0.094,0.238,0.094,0.172,-0.145,0.172,-0.145,0.261,-0.145,0.261,-0.145,0.302,-0.002,0.302,-0.002,0.268,0.017,0.268,0.017,0.236,-0.094,0.236,-0.094,0.225,-0.094,0.225,-0.094,0.291,0.145,0.304,0.145,0.271,0.031,0.271,0.031,0.304,0.013,0.304,0.013,0.328,0.094,0.328,0.094,0.341,0.094,0.341,0.094,0.273,-0.145,0.273,-0.145,0.311,-0.145,0.311,-0.145,0.395,0.145,0.395,0.145,0.304,0.145,0.408,0.145,0.325,-0.145,0.325,-0.145,0.352,-0.145,0.352,-0.145,0.415,-0.145,0.415,-0.145,0.500,0.145,0.500,0.145,0.408,0.145,0.432,0.094,0.446,0.094,0.446,0.094,0.391,-0.094,0.391,-0.094,0.378,-0.094,0.378,-0.094,0.432,0.094,-0.281,0.077,-0.344,-0.145,-0.344,-0.145,-0.221,-0.145,-0.221,-0.145,-0.206,-0.094,-0.206,-0.094,-0.290,-0.094,-0.290,-0.094,-0.256,0.026,-0.256,0.026,-0.209,0.026,-0.209,0.026,-0.216,-0.000,-0.216,-0.000,-0.247,-0.001,-0.247,-0.001,-0.259,-0.052,-0.259,-0.052,-0.192,-0.052,-0.192,-0.052,-0.155,0.077,-0.155,0.077,-0.281,0.077,-0.024,0.077,-0.088,-0.145,-0.088,-0.145,-0.048,-0.145,-0.048,-0.145,0.000,0.026,0.000,0.026,0.066,0.026,0.066,0.026,0.043,-0.060,0.043,-0.060,0.080,-0.059,0.080,-0.059,0.105,0.026,0.105,0.026,0.169,0.026,0.169,0.026,0.124,-0.145,0.124,-0.145,0.160,-0.145,0.160,-0.145,0.222,0.077,0.222,0.077,-0.024,0.077,-0.117,0.026,-0.087,0.026,-0.087,0.026,-0.121,-0.094,-0.121,-0.094,-0.151,-0.094,-0.151,-0.094,-0.117,0.026);

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
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

void cnoise(in float x, out float n)
{
    float x1 = floor(x), 
        x2 = ceil(x),
        h1, h2;
    rand(x1*c.xx, h1);
    rand(x2*c.xx, h2);
    x = fract(x);
    
    float a = .5*(1.-abs(h1-h2)/sqrt(3.));
    n = clamp(sqrt(3.)*mix(1.-x,x,step(h1,h2))+.5*(h1+h2-sqrt(3.)), min(h1,h2), max(h1,h2));
}

void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

void dteam210(in vec2 x, out float ret)
{
    ret = 1.;
    float da;

    float n = 0.;
    for(int i=0; i<npts/4; ++i)
    {
        vec2 ptsi = vec2(path[4*i], path[4*i+1]),
            ptsip1 = vec2(path[4*i+2], path[4*i+3]),
            k = x-ptsi, 
            d = ptsip1-ptsi;
        
        float beta = k.x/d.x,
            alpha = d.y*k.x/d.x-k.y;
        
        n += step(.0, beta)*step(beta, 1.)*step(0., alpha);
        dlinesegment(x, ptsi, ptsip1, da);
        ret = min(ret, da);
    }
    
    ret = mix(ret, -ret, mod(n, 2.));
}

float sm(in float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(d2d, abs(z)-0.5*h);
    d = min(max(w.x,w.y),0.0) + length(max(w,0.0));
}

// iq's smooth minimum
void smoothmin(in float a, in float b, in float k, out float dst)
{
    float h = max( k-abs(a-b), 0.0 )/k;
    dst = min( a, b ) - h*h*h*k*(1.0/6.0);
}

void scene(in vec3 x, out vec2 sdf)
{
    sdf = x.z*c.xy;
    const float r = .025;
    vec2 y = mod(x.xy, r)-.5*r;
    float yi = (x.y-y.y)/r;
    x.x += .5*r*mod(yi,2.);
    y = mod(x.xy, r)-.5*r;
    float d = length(y)-.25*r;
    zextrude(x.z, d, .005, d);
    smoothmin(sdf.x, d, .01, sdf.x);
    
//     float d =
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

void glanz(in vec2 uv, inout vec3 col)
{
    float phi = atan(uv.y, uv.x), 
        r = length(uv);
    col = mix(2.*c.xxx, col, smoothstep(0.,(.5+.5*pow(sin(2.2*pi*phi),3.))*.03, length(uv)));
}

void main()
{
    vec2 uv = (gl_FragCoord.xy -.5*iResolution.xy)/iResolution.y;
    vec3 col = c.yyy;
 
    vec2 s;
    vec3 o = c.yyx,
        r = c.xyy,
        t = c.yyy, 
        u = cross(normalize(t-o),-r),
        dir,
        n, 
        x,
        c1 = c.yyy,
        c2 =  texture(iChannel0, gl_FragCoord.xy/iResolution.xy).rgb,
        l;
    int N = 250,
        i;
    float d = 0.;
    
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);
    
    d = -(o.z-.01)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
        x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4) break;
//         d += s.x<1.e-2?min(s.x,1.e-2):s.x;
//         d += min(s.x, 5.e-3);
        d += s.x;
    }
    
    if(i<N)
    {
        normal(x, n, 5.e-4);
        l = normalize(x+.5*n);
        col = .7*vec3(0.16,0.02,0.40);
        col = .5*col*c2;
        float na;
        mfnoise(uv, 3.,300., .25, na);
        
        col = mix(col, 4.*col, clamp(dot(c.yxy, n),0.,1.));
        col = mix(col, 2.*col, clamp(dot(c.xxy, n),0.,1.));
        col = mix(col, .55*c.xxx, .8+.2*na);
            col = .1*col 
                + .8*col*col*dot(l, n)
                + 3.4*col*col*pow(abs(dot(reflect(l,n),dir)),1.);
        col = clamp(col, 0., 1.);
        col = 1.*col*col*col;
//         col = mix(col, col*col, );
//         col = col*col*col;
    }
    else col = c.xyy;
    float foo;
    
    // lower decoration
    cnoise(4.*uv.x, foo);
    foo = uv.y+1.-.1*foo;
    col = mix(col,.6*(col*col*col+.4*c.yyx), sm(foo-.51));

    // upper decoration
    cnoise(4.*uv.x-1337., foo);
    foo = uv.y-1.+.1*foo;
    col = mix(.6*(col*col*col+.4*c.yyx), col, sm(foo+.51));

    col = col*length(col)/sqrt(3.)*c.xxx;
    c1 = col;
    dteam210(1.1*vec2(1.3,4.)*(uv-.4375*c.yx-.38*c.xy), d);
    col = mix(col, mix(col,c.yyy,.7), sm((d-.02)/4.));
    d+=.002;
    col = mix(col, mix(col, vec3(1.,.3,.1)*4.*c1*c.xxx, .9), sm(d/4.));
    col = mix(col,2.*c.xxx, sm((abs(d)-.002))*4.);
//     col = max(c1,col);
//     d = abs(d-.01)-.005;
//     col = mix(col, c.xxx, sm(abs(d-.003)-.001));
    
    cnoise(4.*uv.x - 2337., foo);
    
    col = mix(col, c2, sm(abs(x.y+.025)-.05*foo-.375));
    col = mix(col, c.xxx, sm(abs(abs(x.y+.025)-.05*foo-.375)-.002));
    col = mix(col, 2.*vec3(0.79,0.03,0.33), sm(abs(abs(x.y+.025)-+.05*foo-.372)-.001));
 
    glanz((uv-.47*c.yx-.73*c.xy), col);
    glanz((uv-.47*c.yx-.03*c.xy), col);
//     glanz(vec2(1.3,4.)*(uv-.42*c.yx-.23*c.xy), col);
//     glanz(vec2(1.3,4.)*(uv-.43*c.yx-.53*c.xy), col);
    
    // grid lines
    vec2 y = mod(uv, .025)-.0125;
    y = abs(y);
    float rr;
    lfnoise(uv-mod(uv, .025)-iTime, rr);
    col = mix(col, mix(c.xxx,c.yyy, sm(min(y.x,y.y))),.1*rr);
    
    // Scan lines
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);
    
    gl_FragColor = vec4(clamp(col, 0.,1.), 1.);
}
