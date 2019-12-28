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

void dglyph(in vec2 x, in float ordinal, in float size, out float dst)
{
    // Bounding box
    float dis;
    dbox(x, 2.*size*c.xx, dis);
    if(dis > 0.)
    {
        dst = dis+.5*size;
        return;
    }

    // Find glyph offset in glyph index
    float nglyphs, offset = 0.;
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
    
    // Get distance from glyph data
    float d = 1., da = 1.;
    
    // Spline points
    float npoints;
    rfloat(offset, npoints);
//     npoints = floor(npoints);
    
    // Contours
    float ncont;
    rfloat(offset + 1. + 3. * npoints, ncont);
//     ncont = floor(ncont);
    
    // Offsets
    float xoff = offset + 1.,
        yoff = offset + 1. + npoints,
        toff = offset + 1. + 2. * npoints,
        coff = offset + 2. + 3. * npoints; 
    
    // Loop through the contours of the glyph. All of them are closed.
    for(float i=0.; i<ncont; i+=1.)
    {
        // Get the contour start and end indices from the contour array.
        float istart = 0., 
            iend;
        rfloat(coff + i, iend);
//         iend = floor(iend);
        if(i>0.)
        {
            rfloat(coff+i-1., istart);
            istart += 1.;
//             istart = floor(istart+1.);
        }
        
        // Prepare a stack
        vec2 stack[3];
        float tstack[3];
        int stacksize = 0;
        
        // Loop through the segments
        for(float j = istart; j <= iend; j += 1.)
        {
            // FIXME: round? floor?
            rfloat(toff + j, tstack[stacksize]);
//             tstack[stacksize] = round(tstack[stacksize]);
            rfloat(xoff+j, stack[stacksize].x);
            rfloat(yoff+j, stack[stacksize].y);
            ++stacksize;
            
            // Check if line segment is finished
            if(stacksize == 2)
            {
                if(tstack[0] == 1. && tstack[1] == 1.)
                {
                    dlinesegment(x/size, stack[0], stack[1], da);
                    d = min(d, da*size);
                    j -= 1.;
                    stacksize = 0;
                }
            }
            else 
            if(stacksize == 3)
            {
                if(tstack[0] == 1. && tstack[2] == 1.)
                {
                    dspline2(x/size, stack[0], stack[1], stack[2], da);
                    d = min(d, da*size);
                    j -= 1.;
                    stacksize = 0;
                }
                else
                {
                    vec2 p = mix(stack[1], stack[2], .5);
                    dspline2(x/size, stack[0], stack[1], p, da);
                    d = min(d, da*size);
                    stack[0] = p;
                    tstack[0] = 1.;
                    j -= 1.;
                    stacksize = 1;
                }
            }
        }
        rfloat(toff + istart, tstack[stacksize]);
        rfloat(xoff + istart, stack[stacksize].x);
        rfloat(yoff + istart, stack[stacksize].y);
        ++stacksize;
        if(stacksize == 2)
        {
            dlinesegment(x/size, stack[0], stack[1], da);
            d = min(d,da*size);
        }
        else if(stacksize == 3)
        {
            dspline2(x/size, stack[0], stack[1], stack[2], da);
            d = min(d, da*size);
        }
    }
    
    stroke(d,.025*size,d);

    dst = d;
}
