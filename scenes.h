#ifndef SCENES_HEADER
#define SCENES_HEADER

#define t_forestmountains (0)
#define t_forestinside (30.959999084472656)
#define t_branch (67.73999786376953)
#define t_mountains (92.9000015258789)
#define t_greetings (123.2699966430664)
#define t_flowers (144.82000732421875)
#define t_forestmountains_2 (150.3000030517578)
#define t_forestinside_2 (155.83999633789062)
#define t_branch_2 (162.8699951171875)
#define t_mountains_2 (167.44000244140625)
#define duration (188)

const double start_times[] = {
    t_forestmountains,
    t_forestinside,
    t_branch,
    t_mountains,
    t_greetings,
    t_flowers,
    t_forestmountains_2,
    t_forestinside_2,
    t_branch_2,
    t_mountains_2,
};

const char *scene_names[] = {
    "Forest Mountains",
    "Forest Inside",
    "Branch",
    "Mountains",
    "Greetings",
    "Flowers",
    "Forest Mountains II",
    "Forest Inside II",
    "Branch II",
    "Mountains II",
};

const unsigned int nscenes = ARRAYSIZE(start_times);

// We need these two arrays to always have the same size - the following line will cause a compiler error if this is ever not the case
_STATIC_ASSERT(ARRAYSIZE(start_times) == ARRAYSIZE(scene_names));

#endif
