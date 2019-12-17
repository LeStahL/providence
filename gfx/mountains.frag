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

void lcfnoise(in vec2 t, out float n)
{
    vec2 i = floor(t);
    t = fract(t);
//     t = mix(c.yy, c.xx, t);
    t = smoothstep(c.yy, c.xx, t);
    vec2 v1, v2;
    rand(i, v1.x);
    rand(i+c.xy, v1.y);
    rand(i+c.yx, v2.x);
    rand(i+c.xx, v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    n = mix(v1.x, v1.y, t.x);
    n = 1.-2.*abs(.5-abs(n));
}

void mfnoise(in vec2 x, in float d, in float b, in float e, out float n)
{
    n = 0.;
    float a = 1., nf = 0., buf;
    for(float f = d; f<b; f *= 2.)
    {
        lcfnoise(f*x+.1*a, buf);
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
    mfnoise(x.xy,1.,200.,.4, flo);
    flo = .5+.5*flo;
//     flo = 38.*pow(sin(.5+.5*flo),3.);
//     flo = 4.*(pow(abs(flo),1.));
    flo = 2.*smoothstep(1.,0.,flo);
    float dx;
    lfnoise(x.y*c.xx, dx);
    flo *= smoothstep(-.2,.2,abs(x.x-.3*dx));
    
    mfnoise(x.xy, 22., 2200., .75, fhi);
//     fhi = tanh(1.-pow(1.-abs(flo),2.));
//     fhi = 4.*pow(.5+.5*fhi,2.);
    
//     mfnoise(x.xy, 27.4e1,27.4e3, .75, fhi);
    
//     float d = .3*smoothstep(-1.,1.,abs(flo)-.3)+.4*flo+.01*fhi;
    sdf = vec2(x.z-.4*flo+.003*fhi-.12*abs(x.x*x.x)+.1, 0.);
//     sdf.x += .1*fhi;
    
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
    
    float flo, fhi;
    mfnoise(uv, 22., 2200., .7, fhi);
    col = c.xxx;
    col = mix(col, c.xyy, (.5+.5*fhi)*dot(n,c.yzy));

}

float sm(in float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
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
        t = .3*c.yyx, 
        u = cross(normalize(t-o),-r),
        dir,
        n, 
        x,
        c1 = c.yyy,
        l;
    int N = 850,
        i;
    float d = 0.,
        dlower,
        dupper;
    t += uv.x * r + uv.y * u;
    dir = normalize(t-o);
    
    d = -(o.z-.3)/dir.z;
    
    if(uv.y < .025+abs(uv.x*uv.x))
    {
    // Raymarch mountains
        for(i = 0; i<N; ++i)
        {
            x = o + d * dir;
            scene(x,s);
            if(s.x < 1.e-4) break;
//             d += s.x<3.e-2?min(s.x,1.e-3):s.x;
            d += s.x<4.e-2?min(s.x,3.e-3):min(s.x, 2.e-2);
            if(d>10.)
            {
                i = N;
                break;
            }
//             d += s.x;
        }
        
        if(i<N)
        {
            normal(x, n, 5.e-4);
            l = normalize(.5*n);
          
            if(s.y == 3.)
            {
                col = c.yxy;
            }
            {
                
                
                
                if(s.y == 0.)
                {
                    floor_texture(x.xy,n, col);
                    col = mix(col, c.xxx, clamp(x.z/.2,0.,1.));
                }
                col = .1*col 
                    + .2*col*dot(l, n)
                    + 1.4*mix(col, vec3(0.59,0.58,0.55), dot(n,c.yyx))*pow(abs(dot(reflect(l,n),dir)),6.);
                col = mix(col, 1.3*c.xxx, clamp(dot(n,c.xyy+.3*c.yzx),0.,1.));
                //col = mix(col,vec3(0.38,1.00,0.86), clamp(dot(n,c.yzy),0.,1.));
            }
            

        }
        else 
            col = c.xyy;
    }
    else 
    {
        col = c.xxx;
        d = 10.;
    }
    col = mix(col, mix(2.3*vec3(0.80,0.55,0.47),vec3(0.39,0.54,0.57), clamp((length(uv))/.3,0.,1.)), clamp(sqrt(d)/3.,0.,1.));
    
    col *= col;
    col = clamp(col, 0., 1.);
    
    // Snowflakes
    for(float ka = 0.; ka < 3.; ka += 1.)
    {
        vec2 dx;
        lfnoise(2.*uv-iTime, dx.x);
        lfnoise(2.*uv-iTime-1337., dx.y);
    
        float v, vc;
        vec2 vi;
        float vsize = 16.+8.*ka;
        dvoronoi(vsize*(uv+.5*iTime*c.yx-.1*dx), v, vi, vc);
        vc /= vsize;
        vi /= vsize;
        col = mix(col, c.xxx, sm((length(vi-uv-.5*iTime*c.yx+.1*dx)-(.1-.2*ka)*vc)/(5.)));
    }
    
//     col = sqrt(col);
//     col *= col; //FIXME: scale?
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
