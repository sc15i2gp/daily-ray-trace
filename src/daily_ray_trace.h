#include <Windows.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <float.h>

typedef enum
{
    GEO_TYPE_NONE,
    GEO_TYPE_POINT,
    GEO_TYPE_SPHERE,
    GEO_TYPE_PLANE,
    GEO_TYPE_COUNT
}   geometry_type;

typedef struct scene_point scene_point;

#include "types.h"
#include "utils.h"
#include "spectrum.h"
#include "geometry.h"
#include "rng.h"
#include "win32_platform.h"
#include "bdsf.h"
#include "read_scene.h"

typedef struct
{
    char          name[32];
    geometry_type type;
    vec3          position;
    union
    {
        f64  radius;
        struct
        {
            vec3 normal;
            vec3 u; //Bounds vector
            vec3 v; //Bounds vector
        };
    };
} object_geometry;

typedef struct
{
    char       name[32];
    u32        is_black_body;
    u32        is_emissive;
    f64        shininess;
    spectrum   emission_spd;
    spectrum   diffuse_spd;
    spectrum   glossy_spd;
    u32        num_bdsfs;
    bdsf_func  bdsfs[16];
    dir_func   sample_direction;
} object_material;

struct scene_point
{
    vec3 position;
    vec3 normal;

    object_material *material;
    object_geometry *surface;
};

#include "utils.c"
#include "spectrum.c"
#include "geometry.c"
#include "rng.c"
#include "win32_platform.c"
#include "bdsf.c"
#include "read_scene.c"

#ifndef NAN
#error "NAN not supported, dingus!"
#endif

//Object in scene consists of surface geometry and material
//Need to iterate over geometries
//Need to iterate over emissive surfaces
//Need to be able to associate a geometry and a material
//Keep emissive materials together
//Keep geometries together

typedef struct
{
    u32 num_surfaces;
    object_geometry *surfaces;
    u32             *surface_material_indices;

    u32 num_scene_materials;
    object_material *scene_materials;
} scene_data;

typedef struct
{
    vec3 forward;
    vec3 right;
    vec3 up;
    vec3 aperture_position;
    f64  aperture_radius;
    vec3 film_bottom_left;
    f64  pixel_width;
    f64  pixel_height;
} camera_data;

void load_scene(const char *path, camera_data *camera, scene_data *scene, u32 width_px, u32 height_px);
void init_scene(scene_data *scene, scene_input_data *scene_input);
void init_camera(camera_data *, camera_input_data *);
void print_camera(camera_data *);
void print_scene(scene_data *);
void render_image(const char *output_path, const char *output_avg_path, const char *output_var_path, u32 dst_width, u32 dst_height, scene_data *scene, camera_data *camera, u32 samples);

#include "daily_ray_trace.c"
