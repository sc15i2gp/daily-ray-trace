typedef void (*bdsf_func)(spectrum, scene_point*, vec3);
typedef void (*dir_func)(vec3*, f64*, vec3, vec3);

#define BDSF(bdsf) void bdsf(spectrum, scene_point*, vec3);
#define DIRF(dirf)
#include "bdsf_list.h"
#undef DIRF
#undef BDSF

#define BDSF(bdsf) #bdsf,
#define DIRF(dirf)
const char *bdsf_name_list[] =
{
    #include "bdsf_list.h"
};
#undef DIRF
#undef BDSF

#define BDSF(bdsf) bdsf,
#define DIRF(dirf)
bdsf_func bdsf_list[] =
{
    #include "bdsf_list.h"
};
#undef DIRF
#undef BDSF
u32 num_bdsfs_defined = sizeof(bdsf_list)/sizeof(bdsf_func);

#define BDSF(bdsf)
#define DIRF(dirf) void dirf(vec3*, f64*, vec3, vec3);
#include "bdsf_list.h"
#undef DIRF
#undef BDSF

#define BDSF(bdsf)
#define DIRF(dirf) #dirf,
const char *dir_func_name_list[] =
{
    #include "bdsf_list.h"
};
#undef DIRF
#undef BDSF

#define BDSF(bdsf)
#define DIRF(dirf) dirf,
dir_func dir_func_list[] =
{
    #include "bdsf_list.h"
};
#undef DIRF
#undef BDSF
u32 num_dir_funcs_defined = sizeof(dir_func_list)/sizeof(dir_func);
