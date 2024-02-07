typedef struct
{
    vec3 position;
    vec3 up;
    vec3 right;
    vec3 forward;
    f64  fov;
    f64  fdepth;
    f64  flength;
    f64  aperture;
    u32  width_px;
    u32  height_px;
} camera_input_data;

typedef enum
{
    SPD_METHOD_NONE,
    SPD_METHOD_RGB,
    SPD_METHOD_CSV,
    SPD_METHOD_BLACKBODY,
    SPD_METHOD_COUNT
} spd_input_method;

typedef struct
{
    spd_input_method method;
    
    u32 has_scale_factor;
    f64 scale_factor;

    union
    {
        rgb_f64 rgb;
        char    csv[64];
        f64     blackbody_temp;
    };
} spd_input_data;

typedef struct
{
    char name[32];
    u32  is_black_body;
    u32  is_emissive;
    f64  shininess;

    spd_input_data emission_input;
    spd_input_data diffuse_input;
    spd_input_data glossy_input;
} material_input_data;

typedef struct
{
    char          name[32];
    geometry_type type;
    vec3          position;
    char          material_name[32];
    union
    {
        f64 radius;
        struct
        {
            vec3 normal;
            vec3 u;
            vec3 v;
        };
    };
} surface_input_data;

typedef struct
{
    u32 num_scene_materials;
    u32 num_surfaces;

    material_input_data scene_materials[16];
    surface_input_data  surfaces[16];
} scene_input_data;

void parse_scene(char *scene_contents, u32 scene_size, camera_input_data *, scene_input_data *);
