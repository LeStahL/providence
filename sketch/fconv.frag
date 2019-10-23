#version 440

precision highp float;

uniform sampler2D iData;
uniform int iDataWidth;

out vec4 oData;

const float pi = acos(-1.);

void main()
{
    int j = int(floor(gl_FragCoord.x + gl_FragCoord.y * float(iDataWidth)));
    oData.r = texture(iData, gl_FragCoord.xy/float(iDataWidth)).r;
//     oData.r = sin(2.*pi*j/12.);
}
