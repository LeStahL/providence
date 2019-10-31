#version 130

const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);

//distance to spline with parameter t
float dist2(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t);
void dspline22(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float ds);

void dspline32(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, in vec2 p3, out float d)
{
    float da;
    vec2 pc = mix(p1,p2,.5);
    dspline2(x, p0, p1, pc, d);
    dspline2(x, pc, p2, p3, da);
    d = min(d, da);
}
