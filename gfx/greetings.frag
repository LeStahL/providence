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
        lfnoise(f*x+.1*a, buf);
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

void dpolygon(in vec2 x, in float N, out float d)
{
    d = 2.0*pi/N;
    float t = mod(acos(x.x/length(x)), d)-0.5*d;
    d = -0.5+length(x)*cos(t)/cos(0.5*d);
}

void dbeerpolygon(in vec2 x, in float N, in float s, out float d)
{
    d = 2.0*pi/N;
    float t = mod(atan(x.y,x.x), d)-0.5*d;
    d = -0.5+length(x)*cos(t)/cos(0.5*d)+.1*smoothstep(s,1.,clamp(abs(t),0.,1.));
}

//distance to spline with parameter t
float dist3(vec3 p0,vec3 p1,vec3 p2,vec3 x,float t)
{
    t = clamp(t, 0., 1.);
    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);
}

//minimum dist3ance to spline
void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds)
{
    //coefficients for 0 = t^3 + a * t^2 + b * t + c
    vec3 E = x-p0, F = p2-2.*p1+p0, G = p1-p0,
    	ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);

	//discriminant and helpers
    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;
    
    //triple real root
    if(dis > 0.) 
    {
        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);
        ds = dist3(p0,p1,p2,x,ui.x+ui.y-tau);
        return;
    }
    
    //three dist3inct real roots
    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;
    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;
    ds = min(
        dist3(p0,p1,p2,x, t.x),
        min(
            dist3(p0,p1,p2,x,t.y),
            dist3(p0,p1,p2,x,t.z)
        )
    );
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

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(d2d, abs(z)-0.5*h);
    d = min(max(w.x,w.y),0.0) + length(max(w,0.0));
}

void scene(in vec3 x, out vec2 sdf)
{
    x.y += .3*iTime;
    x.x -= .4;
    
    // Counter
    float dx, dy;
    lfnoise(x.y*c.xx, dx);
    dx *= .0005;
    mfnoise(x.y*c.xx, 1.,100., .15, dy);
    dy = .9+.1*dy;
    
    float na;
    mfnoise((x.x-dx)*c.xx*tanh(dy), 40., 4000., .75, na);
    vec2 sda = vec2(abs(x.z+.2)-.2-.001*na, 0.);
    sdf = c.xy;
    
    x.x += .5;
    
    // Beer glass body
    const float lsize = 1.;
    float y = mod(x.y, lsize)-.5*lsize,
        d,
        yi = (x.y-y)/lsize;
    rand(yi*c.xx, na);
    dbeerpolygon(vec2(x.x-.5*na, y)*8., 6.,smoothstep(-.05,.1,x.z)*smoothstep(.25,.1,x.z), d);
    float da = d, dc = d;
    zextrude(x.z-.2, d, .2, d);
    zextrude(x.z-.25, da, .25, da);
    d = max(d, -da);
    
    d = abs(d)-.01;
    d /= 12.;
    
    // Handle
    dspline3(vec3(x.x-.5*na, y, abs(x.z-.1)), vec3(.11, 0., .0), vec3(.11,0.,.06), vec3(.05, 0.,.05), da);
    da = abs(da)-.01;
    smoothmin(d, da, .02, d);
    add(sdf, vec2(d, 1.), sdf);
    
    // Beer
//     d = abs(d)-.01;
//     dc += .04;
//     dc = abs(dc
    zextrude(x.z-.075, dc, .15, d);
//     d = d-.01;
    d /= 12.;
    add(sdf, vec2(d,2.), sdf);
    
    add(sdf, sda, sdf);
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

void floor_texture(in vec2 uv, in vec3 n, inout vec3 col)
{
    uv.y += .3*iTime;
    uv.x -= .4;
    
    float dx, dy;
    lfnoise(uv.y*c.xx, dx);
    dx *= .0005;
    mfnoise(uv.y*c.xx, 1.,100., .15, dy);
    dy = .9+.1*dy;
    
    float na;
    mfnoise((uv.x-dx)*c.xx*tanh(dy), 40., 4000., .75, na);
    col = mix(vec3(0.11,0.11,0.11), vec3(0.56,0.44,0.38), .5+.5*na);
    
    col = mix(col,vec3(0.87,0.52,0.13), clamp(.1*dot(n, c.xyx),0.,1.));
}

float sm(in float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void illuminate(inout vec3 col, in vec3 dir, in vec3 l, in vec3 n, in vec3 x, in vec2 s)
{
    vec3 c1;
    if(s.y == 0.)
    {
        floor_texture(x.xy, n, c1);
        c1 = .1*c1 
            + .2*c1*dot(l, n)
            + 1.4*mix(c1, vec3(0.59,0.58,0.55),pow(abs(dot(reflect(l,n),dir)),2.));
        col = c1;
    }
    else if(s.y == 1.)
    {
        c1 = c.xxx;
        c1 = .1*c1 
            + .2*c1*dot(l, n)
            + .4*c1*pow(abs(dot(reflect(l,n),dir)),2.);
        col = mix(col, c1, .4);
    }
    else if(s.y == 2.)
    {
        c1 = mix(vec3(1.00,0.99,0.35), vec3(0.97,0.72,0.00), 1.-clamp(x.z/.2,0.,1.));
        c1 = .1*c1 
            + .2*c1*dot(l, n)
            + 2.4*c1*pow(abs(dot(reflect(l,n),dir)),2.);
        col = mix(col, c1, .4);
    }
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Set up coordinates and camera
    vec2 uv = (fragCoord.xy-.5*iResolution.xy)/iResolution.y,
        s,
        dsky_sphere;
    uv *= .5;
    vec3 col = c.xxx,
        o0 = .3*c.yyx+1.*c.yzy,
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
    float d = 0.,
        dlower,
        dupper;
    t += uv.x * r + uv.y * u;
    dir = normalize(t-o);
    
    d = -(o.z-.22)/dir.z;
    {
    // Raymarch mountains
        for(i = 0; i<N; ++i)
        {
            x = o + d * dir;
            scene(x,s);
            if(s.x < 1.e-4) break;
//             d += s.x<2.e-2?min(s.x,1.e-2):s.x;
//             d += s.x<4.e-2?min(s.x,3.e-3):min(s.x, 2.e-2);
            if(d>10.)
            {
                i = N;
                break;
            }
            d += min(s.x, 3.e-2);
//             d += s.x;
        }
        
        if(i<N)
        {
            normal(x, n, 5.e-5);
            l = normalize(.5*n);

            illuminate(col, dir, l, n, x, s);
            
            if(s.y >= 1.)
            {
                for(int j = 0; j < 4; ++j)
                {
                    o = x;
                    dir = refract(dir, n, .92);
                    d = 1.e-2;
                
                    for(i = 0; i<N; ++i)
                    {
                        x = o + d * dir;
                        scene(x,s);
                        if(s.x < 1.e-4) break;
            //             d += s.x<2.e-2?min(s.x,1.e-2):s.x;
            //             d += s.x<4.e-2?min(s.x,3.e-3):min(s.x, 2.e-2);
                        if(d>10.)
                        {
                            i = N;
                            break;
                        }
                        d += min(s.x, 5.e-2);
            //             d += s.x;
                    }
                    
                    if(i<N)
                    {
                        normal(x, n, 5.e-5);
                        l = normalize(.5*n);
                        
                        illuminate(col, dir, l ,n, x, s);
                    }
                }
            }
        }
    }

//     col = mix(col, mix(2.3*vec3(0.80,0.55,0.47),vec3(0.39,0.54,0.57), clamp((length(uv))/.3,0.,1.)), clamp(sqrt(d)/3.,0.,1.));
    
//     col *= col;
//     col = clamp(col, 0., 1.);
//     
//     // Snowflakes
//     for(float ka = 0.; ka < 3.; ka += 1.)
//     {
//         vec2 dx;
//         lfnoise(2.*uv-iTime, dx.x);
//         lfnoise(2.*uv-iTime-1337., dx.y);
//     
//         float v, vc;
//         vec2 vi;
//         float vsize = 16.+8.*ka;
//         dvoronoi(vsize*(uv+.5*iTime*c.yx-.1*dx), v, vi, vc);
//         vc /= vsize;
//         vi /= vsize;
//         col = mix(col, c.xxx, sm((length(vi-uv-.5*iTime*c.yx+.1*dx)-(.1-.2*ka)*vc)/(5.)));
//     }
    
//     col = sqrt(col);
//     col *= col; //FIXME: scale?
    col *= col*col;
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
