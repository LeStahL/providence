#version 130

const vec3 c = vec3(1.,0.,-1.);

void rfloat(in float off, out float val);
void dbox(in vec2 x, in vec2 b, out float dst);
void dspline2(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float ds);
void dcircle(in vec2 x, out float d);
void dcirclesegment(in vec2 x, in float r, in float p0, in float p1, out float d);
void stroke(in float d0, in float s, out float d);
void smoothmin(in float a, in float b, in float k, out float dst);

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
    float nglyphs, offset = 0;
    rfloat(1., nglyphs);
        
    for(float i=0.; i<nglyphs; i+=1.)
    {
        float ord;
        rfloat(2.+2.*i, ord);
        ord = floor(ord);
        
        if(ord == ordinal)
        {
            rfloat(2.+2.*i+1., offset);
            offset = floor(offset);
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
    
    // Quads
    float nquads;
    rfloat(offset, nquads);
    nquads = floor(nquads);
    offset += 1.;
    for(float i=0.; i<nquads; i+=1.)
    {
        float x1;
        rfloat(offset, x1);
        offset += 1.;
        float y1;
        rfloat(offset, y1);
        offset += 1.;
        float x2;
        rfloat(offset, x2);
        offset += 1.;
        float y2;
        rfloat(offset, y2);
        offset += 1.;
        float x3;
        rfloat(offset, x3);
        offset += 1.;
        float y3;
        rfloat(offset, y3);
        offset += 1.;
        dspline2(x/1.8, size*vec2(x1,3.*y1), size*vec2(x2,3.*y2), size*vec2(x3,3.*y3), da);
        d = min(d,da*1.8);
    }
    
    stroke(d,.05*size,d);

    dst = d;
}
