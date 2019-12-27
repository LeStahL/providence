#ifndef SCENES_HEADER
#define SCENES_HEADER

#define t_forestmountains (0)
#define t_forestinside (32)
#define t_branch (64)
#define t_mountains (96)
#define t_flowers (128)
#define t_greetings (138)
#define duration (190)

const double start_times[] = {
    t_forestmountains,
    t_forestinside,
    t_branch,
    t_mountains,
    t_flowers,
    t_greetings,
};

const char *scene_names[] = {
    "Forest Mountains",
    "Forest Inside",
    "Branch",
    "Mountains",
    "Flowers",
    "Greetings",
};

const unsigned int nscenes = ARRAYSIZE(start_times);

// We need these two arrays to always have the same size - the following line will cause a compiler error if this is ever not the case
_STATIC_ASSERT(ARRAYSIZE(start_times) == ARRAYSIZE(scene_names));

#endif
