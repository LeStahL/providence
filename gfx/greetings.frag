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

void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
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

void rot(in float phi, out mat2 m)
{
    vec2 cs = vec2(cos(phi), sin(phi));
    m = mat2(cs.x, -cs.y, cs.y, cs.x);
}

void dcircle(in vec2 x, out float d)
{
    d = length(x)-1.0;
}

void dbox(in vec2 x, in vec2 b, out float d)
{
    vec2 da = abs(x)-b;
    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}

const int npts = 208;
const float path[npts] = float[npts](-0.500,0.145,-0.500,-0.029,-0.500,-0.029,-0.353,-0.029,-0.353,-0.029,-0.353,-0.087,-0.353,-0.087,-0.500,-0.087,-0.500,-0.087,-0.500,-0.145,-0.500,-0.145,-0.280,-0.145,-0.280,-0.145,-0.280,0.029,-0.280,0.029,-0.427,0.029,-0.427,0.029,-0.427,0.087,-0.427,0.087,-0.280,0.087,-0.280,0.087,-0.280,0.145,-0.280,0.145,-0.500,0.145,-0.500,0.145,-0.500,0.145,-0.240,0.145,-0.240,0.087,-0.240,0.087,-0.093,0.087,-0.093,0.087,-0.093,0.029,-0.093,0.029,-0.020,0.029,-0.020,0.029,-0.020,0.145,-0.020,0.145,-0.240,0.145,-0.093,0.029,-0.167,0.029,-0.167,0.029,-0.167,-0.087,-0.167,-0.087,-0.093,-0.087,-0.093,-0.087,-0.093,0.029,-0.093,0.029,-0.093,0.029,-0.167,-0.087,-0.240,-0.087,-0.240,-0.087,-0.240,-0.145,-0.240,-0.145,-0.167,-0.145,-0.167,-0.145,-0.167,-0.087,0.093,0.145,0.093,0.087,0.093,0.087,0.020,0.087,0.020,0.087,0.020,0.029,0.020,0.029,0.093,0.029,0.093,0.029,0.093,-0.087,0.093,-0.087,0.020,-0.087,0.020,-0.087,0.020,-0.145,0.020,-0.145,0.240,-0.145,0.240,-0.145,0.240,-0.087,0.240,-0.087,0.167,-0.087,0.167,-0.087,0.167,0.145,0.167,0.145,0.093,0.145,0.353,0.145,0.353,0.087,0.353,0.087,0.280,0.087,0.280,0.087,0.280,0.029,0.280,0.029,0.353,0.029,0.353,0.029,0.353,-0.087,0.353,-0.087,0.280,-0.087,0.280,-0.087,0.280,-0.145,0.280,-0.145,0.500,-0.145,0.500,-0.145,0.500,-0.087,0.500,-0.087,0.427,-0.087,0.427,-0.087,0.427,0.145,0.427,0.145,0.353,0.145);
void d5711(in vec2 x, out float ret)
{
    float d = 1.;
//     dpolygon(.5*x,6.0,d);
    
    x *= .7;
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
    
//     ret = max(d,-ret);
    ret /= .7;
}

void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

void dfarbrausch(in vec2 x, out float d)
{
//     dpolygon(.5*x,6.0,d);
    float da, d0;
    
    x += vec2(.1,0.);
    x *= 1.2;
    
    dlinesegment(x,vec2(-.65,.05),vec2(-.5,.05),d0);
    dlinesegment(x,vec2(-.5,.05),vec2(-.2,-.49),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(-.2,-.49),vec2(-.0,-.49),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(-.0,-.49),vec2(-.27,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(-.07,0.),vec2(-.27,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.2,-.49),vec2(-.07,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.4,-.49),vec2(.13,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.4,-.49),vec2(.2,-.49),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.33,0.),vec2(.13,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.33,0.),vec2(.51,-.33),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.6,-.15),vec2(.51,-.33),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.53,0.),vec2(.6,-.15),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.7,0.),vec2(.53,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.7,0.),vec2(.68,-.04),da);
    d0 = min(d0, da);
    dpolygon(5.*(x+vec2(.3,.65)),6.,da);
    d0 = min(d0, da/5.);
    dpolygon(5.*(x+vec2(-.5,.65)),6.,da);
    d0 = min(d0, da/5.);
    
    stroke(d0,.035, d0);
    d = d0;
}

void dhaujobb(in vec2 x, out float d)
{
//     dpolygon(.5*x,6.0,d);
    float da, d0;
    mat2 m;
	rot(.3,m);
    x = 1.1*x*m;
    x.x *= 1.1;
        
    x += vec2(-.05,.2);
    
    // Left leg
    dbox(x+.35*c.xx,vec2(.1,.05),d0);
    dbox(x+vec2(.3,.25),vec2(.05,.15),da);
    d0 = min(d0,da);
    dbox(x+vec2(.2,.15),vec2(.1,.05),da);
    d0 = min(d0,da);
    dbox(x+vec2(.15,.05),vec2(.05,.15),da);
    d0 = min(d0,da);
    
    // Right leg
    dbox(x-vec2(.65,.35),vec2(.05,.15),da);
    d0 = min(d0,da);

    // Torso
    rot(.2, m);
    dbox(m*(x-vec2(.25,.15)),vec2(.45,.05),da);
    d0 = min(d0,da);
    dbox(m*(x-vec2(-.15,.35)),vec2(.45,.05),da);
    d0 = min(d0,da);
    rot(pi/8.,m);
    dbox(m*(x-vec2(.0,.25)),vec2(.1,.15),da);
    d0 = min(d0,da);
    
    // Penis
    dbox(m*(x-vec2(.1,-.0)),vec2(.025,.1),da);
    d0 = min(d0,da);
    
    // Left hand
    rot(.3,m);
    dbox(m*(x-vec2(.235,.535)),vec2(.035,.15),da);
    d0 = min(d0,da);
    dbox(m*(x-vec2(.225,.7)),vec2(.075,.025),da);
    d0 = min(d0,da);
    
    // Right hand
    rot(-.2,m);
    dbox(m*(x+vec2(.585,-.2)),vec2(.0375,.1),da);
    d0 = min(d0,da);
    
    // Head
    dcircle(6.*(x-vec2(-.15,.58)),da);
    d0 = min(d0,da/6.);
    
    d0 -= .05*(abs(x.x)+abs(x.y)-.2);
    d = d0;
}

void dmercury(in vec2 x, out float d)
{
//     dpolygon(.5*x,6.0,d);
    d = 1.;
    float da;

    x += .1*c.yx;

    // Upper part
    dbox(x-.35*c.yx,vec2(.4,.35), da);
    d = min(d, da);
    dbox(x-.7*c.yx, vec2(.2,.2), da);
    d = max(d,-da);
    dbox(x-.25*c.yx,vec2(.2,.05),da);
    d = max(d,-da);
    
    // Lower part
    dbox(x+.2*c.yx,vec2(.1,.4),da);
    d = min(d, da);
    dbox(x+.2*c.yx, vec2(.4,.1),da);
    d = min(d, da);
}

void drevision(in vec2 x, in float r, out float dst)
{
    float l = length(x),
        p = atan(x.y,x.x),
	    d = abs(l-r*.07)-.02, 
        k1 = abs(l-r*.16)-.03,
        k2 = abs(l-r*.21)-.02, 
        k3 = abs(l-r*.35)-.03,
        k4 = abs(l-r*.45)-.02;
    vec4 n1;
    lfnoise(1.*c.xx-3.*iTime, n1.x);
    lfnoise(2.*c.xx-2.4*iTime, n1.y);
    lfnoise(3.*c.xx-2.9*iTime, n1.z);
    lfnoise(4.*c.xx-3.1*iTime, n1.w);
    n1 = mix(n1,c.yyyy, clamp((iTime-24.)/2.,0.,1.));
    d = min(d, mix(d, abs(l-.11)-.03, step(p, -1.71)*step(-2.73, p)));
    d = min(d, mix(d, k1, step(p+n1.x, 3.08)*step(2.82,p)));
    d = min(d, mix(d, k1, step(p+n1.x, 1.47)*step(.81,p)));
    d = min(d, mix(d, k1, step(p+n1.x, -.43)*step(-1.19,p)));
    d = min(d, mix(d, k2, step(p+n1.y, -2.88)*step(-pi,p)));
    d = min(d, mix(d, k2, step(p+n1.y, pi)*step(2.38,p)));
    d = min(d, mix(d, k2, step(p+n1.y, 2.1)*step(.51,p)));
    d = min(d, mix(d, k2, step(p+n1.y, .3)*step(-1.6,p)));
    d = min(d, abs(l-.24)-.02);
    d = min(d, mix(d, k3, step(p+n1.z, -2.18)*step(-pi, p)));
    d = min(d, mix(d, k3, step(p+n1.z, -1.23)*step(-1.7, p)));
    d = min(d, mix(d, k3, step(p+n1.z, -.58)*step(-.78, p)));
    d = min(d, mix(d, k3, step(p+n1.z, 0.)*step(-.29, p)));
    d = min(d, mix(d, k3, step(p+n1.z, 1.25)*step(1.06, p)));
    d = min(d, mix(d, k3, step(p+n1.z, 1.99)*step(.5*pi, p)));
    d = min(d, abs(l-.41)-.03);
    d = min(d, mix(d, k4, step(p+n1.w, 1.04)*step(.04, p)));
    d = min(d, mix(d, k4, step(p+n1.w, -2.2)*step(-2.34, p)));
    
    dst = d-.01;
}

// Adapted from iq, https://www.shadertoy.com/view/XsXSz4
void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst)
{
	vec2 e0 = p1 - p0;
	vec2 e1 = p2 - p1;
	vec2 e2 = p0 - p2;

	vec2 v0 = p - p0;
	vec2 v1 = p - p1;
	vec2 v2 = p - p2;

	vec2 pq0 = v0 - e0*clamp( dot(v0,e0)/dot(e0,e0), 0.0, 1.0 );
	vec2 pq1 = v1 - e1*clamp( dot(v1,e1)/dot(e1,e1), 0.0, 1.0 );
	vec2 pq2 = v2 - e2*clamp( dot(v2,e2)/dot(e2,e2), 0.0, 1.0 );
    
    float s = sign( e0.x*e2.y - e0.y*e2.x );
    vec2 d = min( min( vec2( dot( pq0, pq0 ), s*(v0.x*e0.y-v0.y*e0.x) ),
                       vec2( dot( pq1, pq1 ), s*(v1.x*e1.y-v1.y*e1.x) )),
                       vec2( dot( pq2, pq2 ), s*(v2.x*e2.y-v2.y*e2.x) ));

	dst = -sqrt(d.x)*sign(d.y);
}

void dschnappsgirls(in vec2 x, out float d)
{
//     dpolygon(.5*x,6.0,d);
    float da, d0;
    
    // Dress
    dtriangle(x, vec2(-.1,-.3), vec2(.5,-.3), vec2(.2, .6), d0);
    dlinesegment(x, vec2(-.1,.325), vec2(.5,.325), da);
    stroke(da,.06,da);
    d0 = max(d0,-da);
    
    // Head
    dcircle(7.*(x-vec2(.2,.5)), da);
    d0 = max(d0, -da+.5);
    d0 = min(d0, da/7.);
    
    // Legs
    dlinesegment(x, vec2(.125,-.3), vec2(.125,-.6), da);
    stroke(da, .06, da);
    d0 = min(d0, da);
    dlinesegment(x, vec2(.275,-.3), vec2(.275,-.6), da);
    stroke(da, .06, da);
    d0 = min(d0, da);
    
    // Shoulders
    dlinesegment(x, vec2(0.05,.25), vec2(.35,.25), da);
    stroke(da, .085, da);
    d0 = min(d0, da);
    
    // Arms
    dlinesegment(x, vec2(.385,.25), vec2(.5, -.1), da);
    stroke(da, .055, da);
    d0 = min(d0, da);
    dlinesegment(x, vec2(.017,.25), vec2(-.1, -.1), da);
    stroke(da, .055, da);
    d0 = min(d0, da);
    
    // Glass
    dtriangle(x, vec2(-.6,.3), vec2(-.4,.1), vec2(-.2,.3), da);
    stroke(da, .0125, da);
    d0 = min(d0, da);
    dlinesegment(x, vec2(-.4,.15), vec2(-.4,-.1), da);
    stroke(da, .0125, da);
    d0 = min(d0, da);
    dtriangle(x, vec2(-.5,-.15), vec2(-.3,-.15), vec2(-.4,-.1), da);
    d0 = min(d0, da);
    
    // Liquid
    dtriangle(x, vec2(-.55,.25), vec2(-.4,.1), vec2(-.25,.25), da);
    d0 = min(d0, da);
    
    // Salad
    dlinesegment(x, vec2(-.4,.1), vec2(-.2,.5), da);
    stroke(da, .01, da);
    d0 = min(d0, da);
    dcircle(24.*(x-vec2(-.3,.3)), da);
    d0 = min(d0, da/24.);
    dcircle(24.*(x-vec2(-.25,.4)), da);
    d0 = min(d0, da/24.);
    
//     d = max(d, -d0);
    d = d0;
}

void dear(in vec2 x, out float d)
{
    d = abs(2.*x.y)
        -.95+smoothstep(0.,.5,clamp(abs(x.x),0.,1.))
        -.5*min(-abs(x.x),.01);
}

void dspacepigs(in vec2 x, out float d)
{
//     dpolygon(.5*x,6.0,d);
    float da, d0;
    
    // Head
    dcircle(2.5*x,d0);
    d0 /= 2.5;
    
    // Ears
    dear(vec2(2.,5.)*x-vec2(.8,1.3), da);
    d0 = min(d0,da/10.);
    dear(vec2(2.,5.)*x+vec2(.8,-1.3), da);
    d0 = min(d0,da/10.);
    
    // Nose
    dcircle(6.*x-vec2(0.,-.5),da);
    d0 = max(d0,-da/6.);
    dcircle(24.*x-vec2(-1.5,-2.),da);
    d0 = min(d0,da/24.);
    dcircle(24.*x-vec2(1.5,-2.),da);
    d0 = min(d0,da/24.);
    
    // Eyes
    dcircle(16.*x-vec2(-3.5,2.5),da);
    d0 = max(d0,-da/16.);
    dcircle(16.*x-vec2(3.5,2.5),da);
    d0 = max(d0,-da/16.);
    dcircle(24.*x-vec2(-5.,3.5),da);
    d0 = min(d0,da/24.);
    dcircle(24.*x-vec2(5.,3.5),da);
    d0 = min(d0,da/24.);
    d = d0;
//     d = max(d, -d0);
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

void dbox3(in vec3 x, in vec3 b, out float d)
{
  vec3 da = abs(x) - b;
  d = length(max(da,0.0))
         + min(max(da.x,max(da.y,da.z)),0.0);
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
    const float lsize = .5;
    float y = mod(x.y, lsize)-.5*lsize,
        d,
        yi = (x.y-y)/lsize;
    rand(yi*c.xx, na);
//     mat2 mm = mat2(cos(pi*na), sin(pi*na), -sin(pi*na), cos(pi*na));
//     y = mm * vec2(x.x-.5*na, y).x;
    vec2 z = vec2(x.x-.5*na, y);
    
    // Bierdeckel
    d = length(z.xy+.3*c.xy)-.1;
    zextrude(x.z, d, .005, d);
    add(sdf, vec2(d, -1.), sdf);
    
    //Gläser
    float phi = 2.*pi*na;
    mat2 mm = mat2(cos(phi), sin(phi), -sin(phi), cos(phi));
    z = mm * z;
    dbeerpolygon(z*8., 6.,smoothstep(-.05,.1,x.z)*smoothstep(.25,.1,x.z), d);
    float da = d, dc = d;
    zextrude(x.z, d, .4, d);
    zextrude(x.z, da, .5, da);
    d = max(d, -da);
    
    d = abs(d)-.01;
    d /= 12.;
    
    // Handle
    dspline3(vec3(z, abs(x.z-.1)), vec3(.11, 0., .0), vec3(.11,0.,.06), vec3(.05, 0.,.05), da);
    da = abs(da)-.01;
    da = abs(da)-.001;
    smoothmin(d, da, .02, d);
    add(sdf, vec2(d, 1.), sdf);
    
    // Beer
    zextrude(x.z, dc, .3, d);
    d /= 12.;
    add(sdf, vec2(d,2.), sdf);
    
    add(sdf, sda, sdf);
    
//     dbox3(vec3(z, x.z), vec3(.1,.1,.5), da);
//     add(sdf, vec2(da, 1.), sdf);
}

void bounding_scene(in vec3 x, out float d)
{
    x.y += .3*iTime;
    x.x += .1;
    
    // Floor over Bierdeckel
    d = x.z-.005;
    
    // Bounds of beer glasses
    const float lsize = .5;
    float y = mod(x.y, lsize)-.5*lsize,
        yi = (x.y-y)/lsize,
        na;
    rand(yi*c.xx, na);
    vec2 z = vec2(x.x-.5*na, y);
    float da;
    dbox3(vec3(z, x.z), vec3(.1,.1,.5), da);
    d = min(d,da);
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

void greetings_texture(in vec2 uv, out vec3 col)
{
    uv.y += .3*iTime;
    uv.x += .1;
    
    float na;

    const float lsize = .5;
    float y = mod(uv.y, lsize)-.5*lsize,
        d,
        yi = (uv.y-y)/lsize;
    rand(yi*c.xx, na);
//     mat2 mm = mat2(cos(pi*na), sin(pi*na), -sin(pi*na), cos(pi*na));
//     y = mm * vec2(x.x-.5*na, y).x;
    vec2 z = vec2(uv.x-.5*na+.3, y);
    
    col = 2.*vec3(0.88,0.89,0.91);
    
    d = abs(length(z)-.07)-.01; // blu ring inside with gold
    col = mix(col, 3.5*vec3(0.00,0.18,0.36), sm(d));
    col = mix(col, 2.*vec3(0.64,0.61,0.42), sm(abs(d)-.002));
    
    vec3 ca;
    yi = floor(yi);
    yi = mod(yi, 7.);
    
    if(yi < .5)
    {
        ca = 1.8*vec3(0.97,0.24,0.24);
        d5711(15.*z, d);
    }
    else if(yi < 1.5)
    {
        ca = 2.8*vec3(0.03,0.39,0.35);
        dfarbrausch(15.*z, d);
    }
    else if(yi < 2.5)
    {
        ca = 2.2*sqrt(vec3(0.83,0.42,0.01));
        dhaujobb(15.*z, d);
    }
    else if(yi < 3.5)
    {
        ca = 2.2*sqrt(vec3(0.95,0.26,0.60));
        dmercury(15.*z, d);
    }
    else if(yi < 4.5)
    {
        ca = 1.*c.xxx;
        drevision(15.*z, 1., d);
    }
    else if(yi < 5.5)
    {
        ca = 2.*c.yyx;//1.8* vec3(0.72,0.62,0.05);
        dschnappsgirls(15.*z, d);
    }
    else if(yi < 6.5)
    {
        ca = c.yyy;//1.8* vec3(0.72,0.62,0.05);
        dspacepigs(15.*z, d);
    }
    else
    {
        ca = c.yxy;
        d = 1.;
    }

    
    float da = abs(length(z)-.09)-.002;
    col = mix(col, ca, sm(da)); // red outside
    
    da = length(z)-.055; // red inside
    col = mix(col, ca, sm(da));
    
    col = mix(col, 2.*vec3(0.88,0.89,0.91), sm(d/15.)); // logo
}

void illuminate(inout vec3 col, in vec3 dir, in vec3 l, in vec3 n, in vec3 x, in vec2 s)
{
    vec3 c1;
    if(s.y == -1.)
    {
        greetings_texture(x.xy, c1);
        c1 = .1*c1 
            + .2*c1*dot(l, n)
            + 1.4*c1*pow(abs(dot(reflect(l,n),dir)),2.);
        col = c1;
    }
    else if(s.y == 0.)
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
            + 1.8*c1*pow(abs(dot(reflect(l,n),dir)),2.);
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
        o0 = .4*c.yyx+1.*c.yzy,
        o = o0,
        r = c.xyy,
        t = c.yyy, 
        u = cross(normalize(t-o),-r),
        dir,
        n, 
        x,
        c1 = c.yyy,
        l;
    int N = 200,
        i;
    float d = 0.,
        dlower,
        dupper;
    t += uv.x * r + uv.y * u;
    dir = normalize(t-o);
    
    d = -(o.z-.22)/dir.z;
    
    // pretrace
//     for(i=0; i<10; ++i)
//     {
//         x = o + d * dir;
//         bounding_scene(x, s.x);
//         d += s.x;
//     }
    {
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
                for(int j = 0; j < 5; ++j)
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
//                         d += s.x;
                    }
                    
//                     if(i<N)
//                     {
                        normal(x, n, 5.e-5);
                        l = normalize(.5*n);
                        
//                     }
                        illuminate(col, dir, l ,n, x, s);
                }
            }
        }
    }

    col *= col*col*col;
    
    o = x;
    dir = normalize(-.2*c.xyz);
    d = 1.e-2;
    
    for(i = 0; i<N; ++i)
    {
        x = o + d * dir;
        scene(x,s);
//         s.x -= .01; // Soften shadow
        if(s.x < 1.e-3) break;
        if(x.z > .22) break;
        d += s.x;
    }
    
    if(d < length(-.2*c.xyz)) // We are in the dark
    {
        if(s.y <= 1.) col *= .4;
    }
    
//     greetings_texture(uv, col);
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
