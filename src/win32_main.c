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
//  - Input arguments
//  - Russian roulette
//  - Tidy
//      - Remove fixed length arrays in scene structs
//      - Use less memory for output spd
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

f64 clamp(f64 f, f64 min, f64 max)
{
    f = (f < min) ? min : f;
    f = (f > max) ? max : f;
    return f;
}

rgb_u8 rgb_f64_to_rgb_u8(rgb_f64 in_rgb)
{
    rgb_u8 out_rgb;
    in_rgb.r = clamp(in_rgb.r, 0.0, 1.0);
    in_rgb.g = clamp(in_rgb.g, 0.0, 1.0);
    in_rgb.b = clamp(in_rgb.b, 0.0, 1.0);
    out_rgb.r = (u8)(in_rgb.r * 255.0);
    out_rgb.g = (u8)(in_rgb.g * 255.0);
    out_rgb.b = (u8)(in_rgb.b * 255.0);

    return out_rgb;
}

//Assume handle starts at pixel data in spd_file
void spd_file_to_bmp(file_handle spd_file, spd_file_header *header, const char *bmp_path, spectrum cmf_x, spectrum cmf_y, spectrum cmf_z, spectrum ref_white)
{
    u32 number_of_pixels = header->width_in_pixels * header->height_in_pixels;
    u32 pixels_size_rgb_u8 = number_of_pixels * sizeof(rgb_u8);
    rgb_u8 *pixels_rgb_u8 = alloc(pixels_size_rgb_u8);
    spectrum pixel_spd = alloc_spd();
    for(u32 pixel = 0; pixel < number_of_pixels; pixel += 1)
    {
        read_file(spd_file, spectrum_size, pixel_spd.samples);
        rgb_f64 pixel_f64 = spectrum_to_rgb_f64(pixel_spd, cmf_x, cmf_y, cmf_z, ref_white);
        pixels_rgb_u8[pixel] = rgb_f64_to_rgb_u8(pixel_f64);
    }

    write_pixels_to_bmp(pixels_rgb_u8, header->width_in_pixels, header->height_in_pixels, bmp_path);
    free_spd(pixel_spd);
    unalloc(pixels_rgb_u8, 0);
}

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
    u32 number_of_pixel_samples = 2;
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

    //spd file -> bmp
    spectrum ref_white = alloc_spd();
    spectrum cmf_x     = alloc_spd();
    spectrum cmf_y     = alloc_spd();
    spectrum cmf_z     = alloc_spd();

    const_spectrum(ref_white, 1.0);
    load_csv_file_to_spectrum(cmf_x, "spectra\\cmf_x.csv");
    load_csv_file_to_spectrum(cmf_y, "spectra\\cmf_y.csv");
    load_csv_file_to_spectrum(cmf_z, "spectra\\cmf_z.csv");

    spd_file_header header;
    file_handle spectrum_output_file = open_file(spectrum_output_path, ACCESS_READ, FILE_EXISTS);
    read_file(spectrum_output_file, sizeof(header), &header);
    spd_file_to_bmp(spectrum_output_file, &header, bmp_output_path, cmf_x, cmf_y, cmf_z, ref_white);
    close_file(spectrum_output_file);

    return 0;
}
