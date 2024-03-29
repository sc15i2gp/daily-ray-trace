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

typedef enum
{
    FILM_SAMPLE_NONE,
    FILM_SAMPLE_CENTER,
    FILM_SAMPLE_RANDOM,
    FILM_SAMPLE_COUNT
} film_sample_scheme;

typedef struct
{
    u32  num_pixel_samples;
    u32  max_cast_depth;
    u32  output_width;
    u32  output_height;
    f64  min_wl;
    f64  max_wl;
    f64  wl_interval;
    char input_scene[64];
    char output_spd[64];
    char average_spd[64];
    char variance_spd[64];
    char output_bmp[64];
    char average_bmp[64];
    char variance_bmp[64];
    char white_spd[64];
    char cmf_x[64];
    char cmf_y[64];
    char cmf_z[64];
    char red_spd[64];
    char green_spd[64];
    char blue_spd[64];
    char cyan_spd[64];
    char magenta_spd[64];
    char yellow_spd[64];
    film_sample_scheme pixel_scheme;
} config_arguments;

void print_config_arguments(config_arguments *config);

typedef struct
{
    u32 id;
    u32 width_in_pixels;
    u32 height_in_pixels;
    u32 number_of_wavelengths;
    u32 has_filter_values;
    f64 min_wavelength;
    f64 wavelength_interval;
} spd_file_header;

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
    f64        roughness;
    spectrum   emission_spd;
    spectrum   diffuse_spd;
    spectrum   glossy_spd;
    spectrum   mirror_spd;
    spectrum   refract_spd;
    spectrum   extinct_spd;
    u32        num_bdsfs;
    bdsf_func  bdsfs[16];
    dir_func   sample_direction;
} object_material;

struct scene_point
{
    vec3 position;
    vec3 normal;
    vec3 out; //Points back towards camera along path
    f64  on_dot; //Dot between normal and out
    f64  trans_wl;

    object_material *surface_material;
    object_material *incident_material;
    object_material *transmit_material;
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
    object_material *base_material;
    object_material *escape_material;
} scene_data;

typedef struct
{
    vec3 forward;
    vec3 right;
    vec3 up;
    vec3 aperture_position;
    f64  aperture_radius;
    f64  focal_depth;
    f64  focal_length;
    vec3 film_bottom_left;
    f64  pixel_width;
    f64  pixel_height;
} camera_data;

void load_scene(const char *path, camera_data *camera, scene_data *scene, u32 width_px, u32 height_px);
void init_scene(scene_data *scene, scene_input_data *scene_input);
void init_camera(camera_data *, camera_input_data *);
void print_camera(camera_data *);
void print_scene(scene_data *);
void render_image(config_arguments *config);

#include "daily_ray_trace.c"
