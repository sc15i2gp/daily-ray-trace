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
//  - Config file for input args
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
//  - Materials:
//      - Mirror
//      - Gold
//      - Glass
//  - Non-pinhole camera
//  - Tidy
//      - Include reference white and cmfs in spd file?
//      - Make render_image more flexible with file writing and sample taking
//          - e.g. no reliance on sample ordering, could do the same couple of pixels
//              a few times in a row then flush
//      - Fix plane intersection method
//          - Negative normal causes problems
//          - Normals in general for planes
//          - Switching u,v in scene files causes reduced noise and vantablack shadows
//          - Should light shining on the back of a plane pass through? (Hint: Probably not)
//      - Fix points_mutually_visible
//          - Kinda complicated with float precision stuff
//          - Same with find_scene_intersection
//      - Do something better for sphere/hemisphere sampling
//          - One problem is that sphere sampling always gens positive z
//          - Actually have competing direction sampling methods (e.g uniform sphere and cos weighted)
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
    //Default args
    //Args (just set to default for now)
    const char *spectrum_output_path = "output\\output.spd";
    const char *average_output_path  = "output\\average.spd";
    const char *variance_output_path = "output\\variance.spd";
    const char *bmp_output_path = "output\\output.bmp";
    const char *avg_output_path = "output\\average.bmp";
    const char *var_output_path = "output\\variance.bmp";
    //const char *scene_input_path = "scenes\\first_scene.scn";
    //const char *scene_input_path = "scenes\\init_cornell.scn";
    const char *scene_input_path = "scenes\\cornell_plane_light.scn";
    u32 image_width_in_pixels = 800;
    u32 image_height_in_pixels = 600;
    u32 number_of_pixel_samples = 2;
    u32 number_of_image_pixels = image_width_in_pixels * image_height_in_pixels;

    init_spd_table(32, 380.0, 720.0, 5.0);

    camera_data camera;
    scene_data  scene;
    load_scene(scene_input_path, &camera, &scene, image_width_in_pixels, image_height_in_pixels);
    print_camera(&camera);
    print_scene(&scene);

    printf("Starting render...\n");
    render_image(spectrum_output_path, average_output_path, variance_output_path, image_width_in_pixels, image_height_in_pixels, &scene, &camera, number_of_pixel_samples);
    printf("Render complete.\n");

    printf("Converting...\n");
    spd_file_to_bmp(spectrum_output_path, bmp_output_path);
    spd_file_to_bmp(average_output_path, avg_output_path);
    spd_file_to_bmp(variance_output_path, var_output_path);
    printf("Converted.\n");

    return 0;
}
