#include <Windows.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <float.h>

#include "types.h"
#include "spectrum.h"
#include "rng.h"
#include "geometry.h"
#include "win32_platform.h"
#include "test.h"

#include "spectrum.c"
#include "rng.c"
#include "geometry.c"
#include "test.c"

#ifndef NAN
#error "NAN not supported, dingus!"
#endif

//Object in scene consists of surface geometry and material
//Need to iterate over geometries
//Need to iterate over emissive surfaces
//Need to be able to associate a geometry and a material
//Keep emissive materials together
//Keep geometries together

typedef enum
{
    GEO_TYPE_NONE,
    GEO_TYPE_SPHERE,
    GEO_TYPE_PLANE,
    GEO_TYPE_COUNT
} GEO_TYPE;

typedef struct
{
    GEO_TYPE type;
    union
    {
        struct
        {
            vec3 center;
            f64  radius;
        } sphere;
        struct
        {
            vec3 origin;
            vec3 normal;
            vec3 u; //Bounds vector
            vec3 v; //Bounds vector
        } plane;
    };
} object_geometry;

typedef struct
{
    u32 is_black_body;
    u32 is_emissive;
    spectrum emission_spd;
} object_material;

//Can use a single index to access geometries and materials if emissive ones are stored together
typedef struct
{
    u32 num_surfaces;
    object_geometry *surfaces;

    u32 num_emissive_surfaces;
    object_geometry *emissive_surfaces;

    u32 num_surface_materials;
    object_material *surface_materials;

    u32 num_emissive_materials;
    object_material *emissive_materials;
} scene_data;

typedef struct
{
    vec3 position;
    vec3 normal;

    u32  is_black_body;
    u32  is_emissive;

    spectrum *emission_spd;
} scene_point;

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

void init_camera(camera_data *, u32, u32, vec3, vec3, vec3, vec3, f64, f64, f64, f64);
void print_camera(camera_data *);
void render_image(spectrum* dst_pixels, u32 dst_width, u32 dst_height, scene_data *scene, camera_data *camera, u32 samples);

#include "daily_ray_trace.c"
