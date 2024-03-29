typedef struct
{
    vec3 target;
    vec3 position;
    f64  roll;
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
    SPD_METHOD_CONST,
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
        f64     constant;
    };
} spd_input_data;

typedef struct
{
    char name[32];
    u32  is_base_material;
    u32  is_escape_material;
    u32  is_black_body;
    u32  is_emissive;
    f64  shininess;
    f64  roughness;

    spd_input_data emission_input;
    spd_input_data diffuse_input;
    spd_input_data glossy_input;
    spd_input_data mirror_input;
    spd_input_data refract_input;
    spd_input_data extinct_input;

    u32 num_bdsfs;
    bdsf_func bdsfs[16];
    dir_func  sample_direction_function;
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
u32 load_csv_file_to_spectrum(spectrum dst, const char *path);

void parse_config(char *config_contents, u32 config_contents_size, config_arguments *config);
