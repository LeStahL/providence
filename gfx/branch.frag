#version 130

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

out vec4 gl_FragColor;

const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);

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

void background(in vec2 uv, out vec3 col)
{
    float n;
    lfnoise(1.4*vec2(2.,1.)*uv, n);
    col = mix(vec3(0.33,0.49,0.32), vec3(0.18,0.20,0.13), .5+.5*n);
    lfnoise(1.1*vec2(3.,2.)*uv-1337., n);
    col = mix(col, vec3(0.24,0.38,0.30), .5+.5*n);
  	lfnoise(1.1*vec2(6.,1.)*uv, n);
    col = mix(col, vec3(0.11,0.09,0.06), clamp(.2-.2*n,0.,1.));
    
    mfnoise(uv, 4.,400., .25, n);
    col = mix(col, 2.*col*col, (.5+.5*n));
}

void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = (sda.x<sdb.x)?sda:sdb;
}

void scene(in vec3 x, out vec2 sdf)
{
    x.x += .1*iTime;
    
    // Branch
    float r = .25;
    vec2 dx;
    lfnoise(x.x*c.xx, dx.x);
    lfnoise(1.3*x.x*c.xx+1337., dx.y);
    
    float n;
    mfnoise(vec2(1.,4.)*vec2(atan(x.y,x.z),x.x), 5.,500., .55, n);
    r += .02*abs(n*n);
    
    sdf = vec2(length(x.yz-vec2(.2,.1)*dx)-r, 0.);
    
    // Shrooms
    const float ssize = .4;
    vec3 y = x - vec3(0., vec2(.2,.1)*dx)-.2*c.yyx,
        z = y;
    z.x = mod(z.x, ssize)-.5*ssize;
    
    //r = smoothstep(.1, .2, clamp(y.z-.1,0.,1.));
    //r = mix(r, 0., step(abs(x.z-.1)-.4,0.));
    float theta = atan(z.x,z.y);
    r = mix(.1*smoothstep(.05,.1, clamp(abs(x.z-.19),0.,1.)),-.3,step(0.,abs(x.z-.2)-.1));
    float w = .01+.01*abs(cos(6.*theta*pi));//smoothstep(.01, .05, clamp(y.z-.1,0.,1.));
    
    //add(sdf, vec2((abs(length(z.xy)-r)-w), 1.), sdf);
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

void branch_texture(in vec2 uv, out vec3 col)
{
    uv.y += .1*iTime;
    float n;
    mfnoise(vec2(1.,3.)*uv, 5.,500., .85, n);
    col = mix(vec3(0.16,0.15,0.17), c.xxx, abs(n));
    
    n = .5+.5*n;
    col = mix(col, vec3(0.21,0.14,0.11), clamp(abs(n*n),0.,1.));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord-.5*iResolution.xy)/iResolution.yy,
        s;
    vec3 o0 = .75*c.yyx+.3*c.yzy,
        o = o0,
        r = c.xyy,
        t = c.yyy, 
        u = cross(normalize(t-o),-r),
        dir,
        n, 
        x,
        c1 = c.yyy,
        l,
        col = c.yyy;
    int N = 450,
        i;
    float d;
    
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);
    
    background(uv, col);

    // Raymarch branch
    for(i = 0; i<N; ++i)
    {
        x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4) break;
        if(x.z < -.1) 
        {
            i = N;
            break;
        }
        //d += s.x<1.e-1?min(s.x,5.e-2):s.x;
        //d += min(s.x, 5.e-3);
        d += s.x;
    }
    
    if(i<N)
    {
        normal(x, n, 5.e-4);
        l = normalize(c.yyx);
      
        if(s.y == 0.) // Branch
        {
            // Determine texture coordinates
            float r = .25;
            vec2 dx;
            lfnoise((x.x+.1*iTime)*c.xx, dx.x);
            lfnoise(1.3*(x.x+.1*iTime)*c.xx+1337., dx.y);

            vec2 texcoord = x.yz - vec2(.2,.1)*dx;
            texcoord = vec2(atan(texcoord.y, texcoord.x), x.x);
            
            branch_texture(texcoord, col);
            col = .1*col 
                + .2*col*dot(l, n)
                + 1.4*col*pow(abs(dot(reflect(l,n),dir)),2.);
            
        }
        else if(s.y == 1.) // Shroom
        {
            col = vec3(0.87,0.27,0.34);
            col = .1*col 
                + .2*col*dot(l, n)
                + 1.4*col*pow(abs(dot(reflect(l,n),dir)),2.);
        }
    }
    
    fragColor = vec4(clamp(col,0.,1.),1.);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
