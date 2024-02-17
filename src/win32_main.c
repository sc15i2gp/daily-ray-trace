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
//  - Variance measuring
//      - Move writing to output spd file to render_image
//      - Remove file pointer assumption from spd to bmp
//      - Move spd file to bmp out of main
//      - Write filter values to output spd file
//      - Reduce memory used for file output
//      - Compute variances and write them to their own spd file
//  - Input arguments
//  - Russian roulette
//  - Tidy
//      - Function prototypes and struct definitions in header files
//      - Include reference white and cmfs in spd file?
//      - Remove fixed length arrays in scene structs
//      - Track open files and allocations in platform, properly free stuff upon program completion
//      - Make naming consistent/good
//          - spectrum is the type, spd should be the name
//          - Add __ for global variables
//      - Fix plane intersection method
//          - Negative normal causes problems
//          - Normals in general for planes
//          - Switching u,v in scene files causes reduced noise and vantablack shadows
//          - Should light shining on the back of a plane pass through? (Hint: Probably not)
//      - Fix points_mutually_visible
//          - Kinda complicated with float precision stuff
//          - Same with find_scene_intersection
//      - I don't think is_blackbody and is_emissive are necessary in scene_point
//      - Do something better for sphere/hemisphere sampling
//          - One problem is that sphere sampling always gens positive z
//          - Actually have competing direction sampling methods (e.g uniform sphere and cos weighted)
//      - Review csv loading
//          - Maybe move it out of spectrum
//      - Good memory management (or at least better)
//      - Reading scene files code (e.g. is_word_char vs is_letter_char || '_')
//      - String type?
//      - Sort out camera
//          - Surely the camera's film dimensions shouldn't be dictated by fov?
//          - Non-pinhole
//      - Review file include structure
//      - Fix matrix stuff, it must not be needed
//      - Do general quality pass over code
//      - Move some functions to utils file
//          - Maths functions mostly such as clamp or max
//      - Better error handling
//  - Test
//      - Try to aggressively test as much code as possible
//      - RNG
//          - Quality tests
//      - Camera
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
//  - Parallelism
//  - Separate tool to do spd->bmp
//  - Non-blackbody emissive sources
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
    const char *bmp_output_path = "output\\output.bmp";
    //const char *scene_input_path = "scenes\\first_scene.scn";
    //const char *scene_input_path = "scenes\\init_cornell.scn";
    const char *scene_input_path = "scenes\\cornell_plane_light.scn";
    u32 image_width_in_pixels = 800;
    u32 image_height_in_pixels = 600;
    u32 number_of_pixel_samples = 32;
    u32 number_of_image_pixels = image_width_in_pixels * image_height_in_pixels;

    init_spd_table(32, 69, 380.0, 720.0, 5.0);

    camera_data camera;
    scene_data  scene;
    load_scene(scene_input_path, &camera, &scene, image_width_in_pixels, image_height_in_pixels);
    print_camera(&camera);
    print_scene(&scene);

    printf("Starting render...\n");
    render_image(spectrum_output_path, image_width_in_pixels, image_height_in_pixels, &scene, &camera, number_of_pixel_samples);
    printf("Render complete.\n");

    spd_file_to_bmp(spectrum_output_path, bmp_output_path);

    return 0;
}
