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

void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = (sda.x<sdb.x)?sda:sdb;
}

// iq's smooth minimum
void smoothmin(in float a, in float b, in float k, out float dst)
{
    float h = max( k-abs(a-b), 0.0 )/k;
    dst = min( a, b ) - h*h*h*k*(1.0/6.0);
}

void scene(in vec3 x, out vec2 sdf)
{
    x.y += .3*iTime;
    
    float flo, fhi;
    mfnoise(x.xy,3.,30.,.25, flo);
    mfnoise(x.xy, 27.4e1,27.4e3, .75, fhi);
    
    sdf = vec2(x.z+.1*flo+mix(0.,.07,(pow(.5+.5*fhi,2.))), 0.);
    
    float s = 1.;
    float v, vc;
    vec2 vi;
    dvoronoi(s*x.xy, v, vi, vc);
    vi /= s;
    
    if(abs(vi.x) > .3)
    {

        vec2 y = x.xy-vi;
        float r = .05;
        smoothmin(sdf.x, length(y)-r, .2, sdf.x);

        float n;
        //mfnoise(x.z*c.xx, 1., 100., .45, n);
        lfnoise(x.z*c.xx, n);

        r = mix(r, 3.*r, clamp(x.z*x.z/(1.+.2*x.y),0.,1.));

        float phi = atan(y.y, y.x);
        float na;
        mfnoise(vec2(phi,x.z), 24.,240.,.3, na);
        s = pi/6.;
        float xi = mod(phi,s)-.5*s,
            xii = phi-xi+.52*n;
        vec2 z = y-1.5*r*vec2(cos(xii), sin(xii));
        float da = length(z)-.1*r+.01*na;
        add(sdf, vec2(da, 2.), sdf);
        smoothmin(sdf.x, da, .1, sdf.x);

        sdf.x /= 4.;
    }
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

void floor_texture(in vec2 uv, inout vec3 col)
{
    uv.y += .3*iTime;    
    
    float flo, fhi;
    lfnoise(3.*uv, flo);
    mfnoise(uv, 7.4e1,7.4e3, .75, fhi);
    float fa = .5+.5*fhi;
    
    col = mix(col, c.yyy*vec3(0.03,0.10,0.18), 4.*pow(fa,2.));
//     col = mix(col, c.xyy, abs(fa)-.25);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Set up coordinates and camera
    vec2 uv = (fragCoord.xy-.5*iResolution.xy)/iResolution.y,
        s,
        dsky_sphere;
    uv *= .5;
    vec3 col = c.yyy,
        o0 = .35*c.yyx+1.*c.yzy,
        o = o0,
        r = c.xyy,
        t = .3*c.yyx, 
        u = cross(normalize(t-o),-r),
        dir,
        n, 
        x,
        c1 = c.yyy,
        l;
    int N = 450,
        i;
    float d = 0.,
        dlower,
        dupper;
    t += uv.x * r + uv.y * u;
    dir = normalize(t-o);
    
    // Raymarch mountains
        for(i = 0; i<N; ++i)
        {
            x = o + d * dir;
            scene(x,s);
            if(s.x < 1.e-3) break;
//             d += s.x<1.e-2?min(s.x,1.e-2):s.x;
            d += min(s.x, 4.e-2);
            //d += s.x;
        }
        
        if(i<N)
        {
            normal(x, n, 8.e-2);
            l = normalize(x + .5*n);
          
            if(s.y == 3.)
            {
                col = c.yxy;
            }
//             else if(s.y == 2.)
            {
                col = mix(vec3(0.74,0.20,0.20),c.yyy,1.-smoothstep(.0,.1,x.z));
                col = mix(col, 4.*vec3(0.93,0.62,0.14), clamp(dot(n,c.yyx),0.,1.));
                col = mix(col,vec3(0.84,0.24,0.91), clamp(dot(n, c.xyy),0.,1.));
                col = mix(col, vec3(0.58,0.18,0.41), clamp(dot(n,c.yxy),0.,1.));
                
                if(s.y == 0.)floor_texture(x.xy, col);
                col = .1*col 
                    + .2*col*dot(l, n)
                    + .4*mix(col, vec3(0.39,0.68,0.55), dot(n,c.yyx))*pow(abs(dot(reflect(l,n),dir)),2.);
                //col = mix(col,vec3(0.38,1.00,0.86), clamp(dot(n,c.yzy),0.,1.));
            }
            

        }
    
    col = mix(col, mix(2.5*vec3(0.80,0.55,0.47),vec3(0.19,0.54,0.57), clamp((length(uv))/(.5+.1*cos(8.*atan(uv.y,uv.x))),0.,1.)), clamp(sqrt(d)/5.,0.,1.));
    
    col *= col;
    col = mix(col, mix(vec3(0.98,0.95,0.66),mix(c.xyy,vec3(0.98,0.95,0.66),.2), 1.-clamp(x.y/5.,0.,1.)), clamp(-x.z/.1,0.,1.));
//     if(x.z<.1)
    
    
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
