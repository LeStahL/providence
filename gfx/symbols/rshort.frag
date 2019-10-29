#version 130

uniform float iFontWidth;
uniform sampler2D iFont;

void rshort(in float off, out float val)
{
    // Parity of offset determines which byte is required.
    float hilo = mod(off, 2.);
    // Find the pixel offset your data is in (2 unsigned shorts per pixel).
    off = .5*off;
    // Find texture coordinates matching offset
    vec2 ind = vec2(mod(off, iFontWidth), floor(off/iFontWidth));
    // Determine data block
    vec4 block = texelFetch(iFont, ivec2(ind), 0);
    // Select the appropriate word
    vec2 data = mix(block.rg, block.ba, hilo);
    // Convert bytes to unsigned short. The lower bytes operate on 255,
    // the higher bytes operate on 65280, which is the maximum range 
    // of 65535 minus the lower 255.
    val = round(dot(vec2(255., 65280.), data));
}
