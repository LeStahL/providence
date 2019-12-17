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

const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);

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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = 1.5*(fragCoord-.5*iResolution.xy)/iResolution.yy;
    
    
    float n;
    lfnoise(1.5*uv.x*c.xx, n);
    
    vec3 col = c.yyy;
    background(uv, col);
    
    uv.x += .3*iTime;
    for(float i=20.; i>= 3.; i -= 1.)
    {
        float na, nb;
        rand(i*c.xx, na);
        float x = mod(.3*i*(uv.x-na), 4./i)-2./i;
        rand((x-.3*i*(uv.x-na))*c.xx, nb);
        vec2 dx;
        lfnoise(iTime*c.xx-i, dx.x);
        lfnoise(iTime*c.xx-1337.-i, dx.y);
    	addflower(vec2(x, .3*i*(uv.y+.3*nb))-.01*dx, col);
    }
    
    col *= col*col;
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
