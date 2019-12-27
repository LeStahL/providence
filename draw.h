#ifndef DRAW_HEADER
#define DRAW_HEADER

if(t < t_forestinside)
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
else if(t < t_branch)
{
    glUseProgram(shader_program_gfx_forestinside.handle);
    glUniform1f(shader_uniform_gfx_forestinside_iTime, t-t_forestinside);
    glUniform2f(shader_uniform_gfx_forestinside_iResolution, w, h);
#ifdef MIDI
    glUniform1f(shader_uniform_gfx_forestinside_iFader0, fader0);
    glUniform1f(shader_uniform_gfx_forestinside_iFader1, fader1);
    glUniform1f(shader_uniform_gfx_forestinside_iFader2, fader2);
    glUniform1f(shader_uniform_gfx_forestinside_iFader3, fader3);
    glUniform1f(shader_uniform_gfx_forestinside_iFader4, fader4);
    glUniform1f(shader_uniform_gfx_forestinside_iFader5, fader5);
    glUniform1f(shader_uniform_gfx_forestinside_iFader6, fader6);
    glUniform1f(shader_uniform_gfx_forestinside_iFader7, fader7);
#endif
}
else if(t < t_mountains)
{
    glUseProgram(shader_program_gfx_branch.handle);
    glUniform1f(shader_uniform_gfx_branch_iTime, t-t_branch);
    glUniform2f(shader_uniform_gfx_branch_iResolution, w, h);
#ifdef MIDI
    glUniform1f(shader_uniform_gfx_branch_iFader0, fader0);
    glUniform1f(shader_uniform_gfx_branch_iFader1, fader1);
    glUniform1f(shader_uniform_gfx_branch_iFader2, fader2);
    glUniform1f(shader_uniform_gfx_branch_iFader3, fader3);
    glUniform1f(shader_uniform_gfx_branch_iFader4, fader4);
    glUniform1f(shader_uniform_gfx_branch_iFader5, fader5);
    glUniform1f(shader_uniform_gfx_branch_iFader6, fader6);
    glUniform1f(shader_uniform_gfx_branch_iFader7, fader7);
#endif
}
else if(t < t_flowers)
{
    glUseProgram(shader_program_gfx_mountains.handle);
    glUniform1f(shader_uniform_gfx_mountains_iTime, t-t_mountains);
    glUniform2f(shader_uniform_gfx_mountains_iResolution, w, h);
#ifdef MIDI
    glUniform1f(shader_uniform_gfx_mountains_iFader0, fader0);
    glUniform1f(shader_uniform_gfx_mountains_iFader1, fader1);
    glUniform1f(shader_uniform_gfx_mountains_iFader2, fader2);
    glUniform1f(shader_uniform_gfx_mountains_iFader3, fader3);
    glUniform1f(shader_uniform_gfx_mountains_iFader4, fader4);
    glUniform1f(shader_uniform_gfx_mountains_iFader5, fader5);
    glUniform1f(shader_uniform_gfx_mountains_iFader6, fader6);
    glUniform1f(shader_uniform_gfx_mountains_iFader7, fader7);
#endif
}
else if(t < t_greetings)
{
    glUseProgram(shader_program_gfx_flowers.handle);
    glUniform1f(shader_uniform_gfx_flowers_iTime, t-t_flowers);
    glUniform2f(shader_uniform_gfx_flowers_iResolution, w, h);
#ifdef MIDI
    glUniform1f(shader_uniform_gfx_flowers_iFader0, fader0);
    glUniform1f(shader_uniform_gfx_flowers_iFader1, fader1);
    glUniform1f(shader_uniform_gfx_flowers_iFader2, fader2);
    glUniform1f(shader_uniform_gfx_flowers_iFader3, fader3);
    glUniform1f(shader_uniform_gfx_flowers_iFader4, fader4);
    glUniform1f(shader_uniform_gfx_flowers_iFader5, fader5);
    glUniform1f(shader_uniform_gfx_flowers_iFader6, fader6);
    glUniform1f(shader_uniform_gfx_flowers_iFader7, fader7);
#endif
}
else {
    glUseProgram(shader_program_gfx_greetings.handle);
    glUniform1f(shader_uniform_gfx_greetings_iTime, t-t_greetings);
    glUniform2f(shader_uniform_gfx_greetings_iResolution, w, h);
#ifdef MIDI
    glUniform1f(shader_uniform_gfx_greetings_iFader0, fader0);
    glUniform1f(shader_uniform_gfx_greetings_iFader1, fader1);
    glUniform1f(shader_uniform_gfx_greetings_iFader2, fader2);
    glUniform1f(shader_uniform_gfx_greetings_iFader3, fader3);
    glUniform1f(shader_uniform_gfx_greetings_iFader4, fader4);
    glUniform1f(shader_uniform_gfx_greetings_iFader5, fader5);
    glUniform1f(shader_uniform_gfx_greetings_iFader6, fader6);
    glUniform1f(shader_uniform_gfx_greetings_iFader7, fader7);
#endif
}
#endif
