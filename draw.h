#ifndef DRAW_HEADER
#define DRAW_HEADER

{
    glUseProgram(shader_program_gfx_forestmountains.handle);
    glUniform1f(shader_uniform_gfx_forestmountains_iTime, t-t_forestmountains);
    glUniform2f(shader_uniform_gfx_forestmountains_iResolution, w, h);
#ifdef MIDI
    glUniform1f(shader_uniform_gfx_forestmountains_iFader0, fader0);
    glUniform1f(shader_uniform_gfx_forestmountains_iFader1, fader1);
    glUniform1f(shader_uniform_gfx_forestmountains_iFader2, fader2);
    glUniform1f(shader_uniform_gfx_forestmountains_iFader3, fader3);
    glUniform1f(shader_uniform_gfx_forestmountains_iFader4, fader4);
    glUniform1f(shader_uniform_gfx_forestmountains_iFader5, fader5);
    glUniform1f(shader_uniform_gfx_forestmountains_iFader6, fader6);
    glUniform1f(shader_uniform_gfx_forestmountains_iFader7, fader7);
#endif
}
#endif
