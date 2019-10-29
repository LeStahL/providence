#version 130

// layout(pixel_center_integer​) in vec4 gl_FragCoord;

uniform sampler2D iData;
uniform int iDataWidth;

out vec4 oData;

const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);

// void rshort(in float off, out float val)
// {
//     // Parity of offset determines which byte is required.
//     float hilo = mod(off, 2.);
//     // Find the pixel offset your data is in (2 unsigned shorts per pixel).
//     off = .5*off;
//     // - Determine texture coordinates.
//     //     offset = i*iDataWidth+j for (i,j) in [0,iDataWidth]^2
//     //     floor(offset/iDataWidth) = floor((i*iDatawidth+j)/iDatawidth)
//     //                              = floor(i)+floor(j/iDataWidth) = i
//     //     mod(offset, iDataWidth) = mod(i*iDataWidth + j, iDataWidth) = j
//     // - For texture coordinates (i,j) has to be rescaled to [0,1].
//     // - Also we need to add an extra small offset to the texture coordinate
//     //   in order to always "hit" the right pixel. Pixel width is
//     //     1./iDataWidth.
//     //   Half of it is in the center of the pixel.
//     vec2 ind = (vec2(mod(off, iDataWidth), floor((off)/iDataWidth))+.05)/iDataWidth;
//     // Get 4 bytes of data from the texture
//     vec4 block = texture(iData, ind);
//     // Select the appropriate word
//     vec2 data = mix(block.rg, block.ba, hilo);
//     // Convert bytes to unsigned short. The lower bytes operate on 255,
//     // the higher bytes operate on 65280, which is the maximum range 
//     // of 65535 minus the lower 255.
//     val = round(dot(vec2(255., 65280.), data));
// }

void rshort(in float off, out float val)
{
    // Parity of offset determines which byte is required.
    float hilo = mod(off, 2.);
    // Find the pixel offset your data is in (2 unsigned shorts per pixel).
    off = .5*off;
    // Find texture coordinates matching offset
    vec2 ind = vec2(mod(off, iDataWidth), floor(off/iDataWidth));
    // Determine data block
    vec4 block = texelFetch(iData, ivec2(ind), 0);
    // Select the appropriate word
    vec2 data = mix(block.rg, block.ba, hilo);
    // Convert bytes to unsigned short. The lower bytes operate on 255,
    // the higher bytes operate on 65280, which is the maximum range 
    // of 65535 minus the lower 255.
    val = round(dot(vec2(255., 65280.), data));
}

void main()
{
    int j = int(floor(gl_FragCoord.x + (gl_FragCoord.y-.5) * float(iDataWidth)));
    vec2 v;
    rshort(2.*j, v.x);
    rshort(2.*j+1, v.y);
    
//     float off = .5*float(j);
//     if(gl_FragCoord.x >= iDataWidth)oData = 1.
    
    vec2 vl = mod(v,256.0)/255.0;
    vec2 vh = floor(v/256.0)/255.0;
    oData = vec4(vl.x,vh.x,vl.y,vh.y);
}
