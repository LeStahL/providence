cmake_minimum_required(VERSION 3.6)
project(fconv)
add_executable(fconv fconv.c)
set_target_properties(fconv PROPERTIES LINK_FLAGS "/ENTRY:demo /SUBSYSTEM:Windows")
target_link_libraries(fconv User32 Winmm opengl32 Kernel32 ucrt Gdi32)

add_executable(loadttf loadttf.c)
target_link_libraries(loadttf User32 Winmm opengl32 Kernel32 ucrt Gdi32)
