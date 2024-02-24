#include "daily_ray_trace.h"

//Standards:
//  - Struct specific variants of similar functions start with lower case struct name
//      e.g minus for Vec3 is vec3_minus
//  - Words in names are separated by _
//  - Functions with a dst ptr arg have the dst ptr as the first argument

//Program description:
//  - Given a scene and sample parameters, produce an image
//  - Output bmp
//  - No windows/live preview for now

//Long term TODO:
//  - Investigate use + comparison of different RGB working spaces
//  - Investigate introduction of error in rgb/spectrum conversion
//      - e.g. which introduces more error: rgb to spectrum or rgb from spectrum
//  - SPD visualisations
//  - Spectral investigations
//      - Get a good understanding of what's going on with spectral stuff
//      - Test code
//      - Currently testing rgb -> spectrum -> rgb, could test all:
//          - rgb -> spectrum
//          - spectrum -> rgb
//          - spectrum -> xyz
//          - spectrum -> rgb
//  - Stop being ignorant about line-shape intersection methods and float precision
//  - Xorshift rng
//  - Textures
//  - Participating media
//  - Remove C std?

//TODO:
//  - Materials
//      - Mirror
//      - Gold
//      - Glass
//  - Non-pinhole camera
//  - Better memory management
//      - Have platform track allocations and open files (and free on shutdown)
//      - Memory arena(s)
//      - Remove fixed length arrays in scene structs
//  - Codify spd table entries
//      - cmfs
//      - rgb to spds
//      - white/ref white
//  - spd to bmp
//      - Separate tool spd->bmp (in case of crash or something) which can either be called in drt or
//        invoked as a standalone program
//      - Platform probably shouldn't have to know about filtering pixels or cmfs
//  - Profiling
//  - Tidy
//      - Include reference white in spd file
//      - Sort out camera
//          - Surely the camera's film dimensions shouldn't be dictated by fov?
//      - String type?
//      - Fix matrix stuff, it must not be needed
//      - Do general quality pass over code
//        - Function prototypes and struct definitions in header files
//        - Code order in files
//        - Review file include structure
//        - Reading scene files code (e.g. is_word_char vs is_letter_char || '_')
//        - Make naming consistent/good
//              - spectrum is the type, spd should be the name
//              - Add __ for global variables
//  - Test
//      - Try to aggressively test as much code as possible
//      - RNG
//          - Quality tests
//      - Camera
//  - Do general performance pass over code
//  - Output rendered scene info to files
//      - e.g. light colours/spectra, material colours etc.
//  - Better RNG
//      - Alternative methods
//      - Consider generating RNG up front
//      - Constant seed vs different seed every run
//  - Sort out units
//      - Units of spds, distances etc.
//      - Also consider attenuation and what exactly surface + emissive spds mean
//      - Actually, generally bone up on physics
//      - Maths too: Surface vs sphere integrals for light transport
//  - Visualisation/comparison methods for material parameters
//      - e.g. what effect does raising/lowering shininess within a range do for glossiness
//      - Maybe have specific output directory for a given scene
//  - Logging
//      - API
//      - How much to log?
//      - Metal
//  - Error handling
//  - Input arguments
//  - Parallelism
//  - Non-blackbody emissive sources
//  - Splitting
//      - Shadow rays/multiple direct light sampling
//      - Multiple estimate_indirect_contributions per path point
//  - Change algorithm back to newer version
//      - Some spectral calculations don't need to be done (e.g. after a ray has escaped)

//Figures to aim for:
//  - Handle images with resolutions up to HD (1920x1080)
//  - Textures up to 1024*1024
//  - Under 1GB memory allocated
//  - Handles all things I want to raytrace
//  - <speed target>
//  - <variance + statistics targets>
//  - <quality target?>

int main(int argc, char **argv)
{
    const char *config_path = "config.cfg";
    file_handle config_file = open_file(config_path, ACCESS_READ, FILE_EXISTS);
    u32 config_file_size = get_file_size(config_file);
    char *config_buffer = alloc(config_file_size);
    read_file(config_file, config_file_size, config_buffer);
    close_file(config_file);

    config_arguments args;
    memset(&args, 0, sizeof(args));
    parse_config(config_buffer, config_file_size, &args);

    printf("CONFIG ARGS:\n");
    print_config_arguments(&args);
    printf("\n");

    u32 spd_table_capacity = 32; //Should this be in config?
    init_spd_table(spd_table_capacity, args.min_wl, args.max_wl, args.wl_interval);

    camera_data camera;
    scene_data  scene;
    load_scene(args.input_scene, &camera, &scene, args.output_width, args.output_height);
    print_camera(&camera);
    print_scene(&scene);

    printf("Starting render...\n");
    render_image(args.output_spd, args.average_spd, args.variance_spd, args.output_width, args.output_height, &scene, &camera, args.num_pixel_samples);
    printf("Render complete.\n");

    printf("Converting...\n");
    spd_file_to_bmp(args.output_spd, args.output_bmp);
    spd_file_to_bmp(args.average_spd, args.average_bmp);
    spd_file_to_bmp(args.variance_spd, args.variance_bmp);
    printf("Converted.\n");

    return 0;
}
