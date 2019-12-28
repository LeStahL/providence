#version 130

const vec3 c = vec3(1.,0.,-1.);

void rfloat(in float off, out float val);
void dbox(in vec2 x, in vec2 b, out float dst);
void dspline2(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float ds);
void dcircle(in vec2 x, out float d);
void dcirclesegment(in vec2 x, in float r, in float p0, in float p1, out float d);
void stroke(in float d0, in float s, out float d);
void smoothmin(in float a, in float b, in float k, out float dst);
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);

// Compute distance to glyph control points for debug purposes
void dglyphpts(in vec2 x, in float ordinal, in float size, out float dst)
{
    // Get glyph index length
    float nchars;
    rfloat(0., nchars);
    
    // Find glyph offset in glyph index
    float nglyphs, offset = 0;
    rfloat(1., nglyphs);
        
    for(float i=0.; i<nglyphs; i+=1.)
    {
        float ord;
        rfloat(2.+2.*i, ord);
//         ord = floor(ord);
        
        if(ord == ordinal)
        {
            rfloat(2.+2.*i+1., offset);
//             offset = floor(offset);
            break;
        }
    }
    
    if(offset == 0.) 
    {
        dst = 1.;
        return;
    }
    
    // Spline points
    float npoints;
    rfloat(offset, npoints);
    
    // Offsets
    float xoff = offset + 1.,
        yoff = offset + 1. + npoints;
        
    // Point distance
    dst = 1.;
    
    // Debug output of the spline control points
    for(float i=0.; i<npoints; i+=1.)
    {
        vec2 xa;
        rfloat(xoff+i,xa.x);
        rfloat(yoff+i,xa.y);
        dst = min(dst, length(x/size-xa)-5.e-2);
    }
    
    dst *= size;
}
