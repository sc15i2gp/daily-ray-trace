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

//TODO:
//  - Render blinn-phong sphere lit by point light
//  - Shape and direction sampling
//  - Geometry intersection
//  - Add func decls to geometry.h
//  - Figure out better strategy for spectra which doesn't leave massive gaps
//      - (without recompiling)
//  - Variance
//  - Sort out camera
//      - Surely the camera's film dimensions shouldn't be dictated by fov?
//      - Non-pinhole

//Figures to aim for:
//  - Handle images with resolutions up to HD (1920x1080)
//  - Textures up to 1024*1024
//  - Under 1GB memory allocated
//  - Handles all things I want to raytrace
//  - <speed target>
//  - <variance + statistics targets>
//  - <quality target?>

rgb_u8 rgb_f64_to_rgb_u8(rgb_f64 in_rgb)
{
    rgb_u8 out_rgb;
    out_rgb.r = (u8)(in_rgb.r * 255.0);
    out_rgb.g = (u8)(in_rgb.g * 255.0);
    out_rgb.b = (u8)(in_rgb.b * 255.0);

    return out_rgb;
}

u32 write_bmp_to_file(windows_bmp *bmp, const char *file_path)
{
    DWORD  bytes_written = 0;
    HANDLE output_file   = CreateFile(file_path, GENERIC_WRITE, FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);

    if(output_file == INVALID_HANDLE_VALUE) return 0;

    BOOL success = WriteFile(output_file, (void*)bmp->file_header, bmp->file_header->bfSize, &bytes_written, NULL);
    CloseHandle(output_file);

    return success;
}

u32 write_pixels_to_bmp(rgb_u8 *pixels, u32 width, u32 height, const char *path)
{
    u32 pixels_size = width * height * sizeof(u32);
    u32 bmp_size = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + pixels_size;
    windows_bmp bmp = {0};
    u8 *raw_bmp = (u8*)VirtualAlloc(0, bmp_size, MEM_COMMIT, PAGE_READWRITE);
    bmp.file_header = (BITMAPFILEHEADER*)raw_bmp;
    bmp.info_header = (BITMAPINFOHEADER*)(raw_bmp + sizeof(BITMAPFILEHEADER));
    bmp.pixels = (rgb_u8*)(raw_bmp + sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER));

    memcpy(bmp.pixels, pixels, pixels_size);
    
    bmp.file_header->bfType = 0x4d42;
    bmp.file_header->bfSize = bmp_size;
    bmp.file_header->bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
    bmp.info_header->biSize = sizeof(BITMAPINFOHEADER);
    bmp.info_header->biWidth = width;
    bmp.info_header->biHeight = height;
    bmp.info_header->biBitCount = 32;
    bmp.info_header->biCompression = BI_RGB;
    bmp.info_header->biSizeImage = 0;
    bmp.info_header->biXPelsPerMeter = 3780;
    bmp.info_header->biYPelsPerMeter = 3780;
    bmp.info_header->biClrUsed = 0;
    bmp.info_header->biClrImportant = 0;

    u32 success = write_bmp_to_file(&bmp, path);
    VirtualFree(raw_bmp, bmp_size, MEM_RELEASE);

    return success;
}

typedef struct
{
    u32 id;
    u32 width_in_pixels;
    u32 height_in_pixels;
    u32 number_of_wavelengths;
    f64 min_wavelength;
    f64 wavelength_interval;
} spd_file_header;

//Assume handle starts at pixel data in spd_file
void spd_file_to_bmp(HANDLE spd_file, spd_file_header *header, const char *bmp_path, spectrum *cmf_x, spectrum *cmf_y, spectrum *cmf_z, spectrum *ref_white)
{
    u32 number_of_pixels = header->width_in_pixels * header->height_in_pixels;
    u32 pixels_size_rgb_u8 = number_of_pixels * sizeof(rgb_u8);
    rgb_u8 *pixels_rgb_u8 = (rgb_u8*)VirtualAlloc(NULL, pixels_size_rgb_u8, MEM_COMMIT, PAGE_READWRITE);
    DWORD bytes_read;

    spectrum pixel_spd;
    for(u32 pixel = 0; pixel < number_of_pixels; ++pixel)
    {
        ReadFile(spd_file, &pixel_spd, sizeof(pixel_spd), &bytes_read, NULL);
        rgb_f64 pixel_f64 = spectrum_to_rgb_f64(&pixel_spd, cmf_x, cmf_y, cmf_z, ref_white);
        pixels_rgb_u8[pixel] = rgb_f64_to_rgb_u8(pixel_f64);
    }

    write_pixels_to_bmp(pixels_rgb_u8, header->width_in_pixels, header->height_in_pixels, bmp_path);
    VirtualFree(pixels_rgb_u8, pixels_size_rgb_u8, MEM_RELEASE);
}

int main(int argc, char **argv)
{
    //call_test_funcs();

    //Default args
    const char *default_spectrum_output_path = "output\\output.spd";
    const char *default_bmp_output_path = "output\\output.bmp";
    u32 default_number_of_spectrum_samples = 69;
    f64 default_smallest_wavelength = 380.0;
    f64 default_largest_wavelength = 720.0;
    f64 default_sample_interval = 5.0;
    u32 default_image_width = 800;
    u32 default_image_height = 600;
    u32 default_spd_pixel_data_size = default_image_width * default_image_height * sizeof(spectrum);
    u32 default_number_of_pixel_samples = 1;

    //Args (just set to default for now)
    number_of_spectrum_samples = default_number_of_spectrum_samples;
    smallest_wavelength = default_smallest_wavelength;
    largest_wavelength = default_largest_wavelength;
    sample_interval = default_sample_interval;
    const char *spectrum_output_path = default_spectrum_output_path;
    const char *bmp_output_path = default_bmp_output_path;
    u32 image_width_in_pixels = default_image_width;
    u32 image_height_in_pixels = default_image_height;
    u32 number_of_image_pixels = image_width_in_pixels * image_height_in_pixels;
    u32 number_of_pixel_samples = default_number_of_pixel_samples;
    u32 spd_pixel_data_size = default_spd_pixel_data_size;

    //Spectrum file contents:
    //- Number of spectra/pixels (dims)
    //- Number of samples per spectrum
    //- Wavelength low, high and interval
    //- Spectral data
    HANDLE spectrum_output_file = CreateFile(spectrum_output_path, GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);

    spd_file_header header;
    memset(&header, 0, sizeof(header));
    header.id = 0xedfeefbe;
    header.width_in_pixels = image_width_in_pixels;
    header.height_in_pixels = image_height_in_pixels;
    header.number_of_wavelengths = number_of_spectrum_samples;
    header.min_wavelength = smallest_wavelength;
    header.wavelength_interval = sample_interval;

    //Write spd file header
    DWORD bytes_written;
    WriteFile(spectrum_output_file, &header, sizeof(header), &bytes_written, NULL);

    //Alloc a certain amount of memory for pixel data
    u32 allocated_pixels = spd_pixel_data_size / sizeof(spectrum);
    spectrum *spd_pixels = (spectrum*)VirtualAlloc(NULL, spd_pixel_data_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    printf("Size = %u\n", spd_pixel_data_size);

#if 1
    camera_data camera;
    vec3 camera_pos     = {0.0, 0.0, 8.0};
    vec3 camera_up      = {0.0, 1.0, 0.0};
    vec3 camera_right   = {1.0, 0.0, 0.0};
    vec3 camera_forward = {0.0, 0.0, -1.0};
    f64 fov          = 90.0;
    f64 focal_depth  = 8.0;
    f64 focal_length = 0.5;
    f64 aperture_rad = 0.0;
    init_camera(&camera, image_width_in_pixels, image_height_in_pixels, camera_pos, camera_up, camera_right, camera_forward, fov, focal_depth, focal_length, aperture_rad);
    print_camera(&camera);
    render_image(spd_pixels, image_width_in_pixels, image_height_in_pixels, NULL, &camera, 1);
#else
    spectrum d;
    load_csv_file_to_spectrum(&d, "spectra\\red_rgb_to_spd.csv");
    for(u32 sample = 0; sample < number_of_pixel_samples; ++sample)
    {
        for(u32 pixel = 0; pixel < number_of_image_pixels; ++pixel)
        {
            u32 spd_pixel = pixel % allocated_pixels;
            copy_spectrum(&spd_pixels[spd_pixel], &d);
            if(spd_pixel == allocated_pixels - 1)
            {
                WriteFile(spectrum_output_file, spd_pixels, spd_pixel_data_size, &bytes_written, NULL);
            }
        }
        SetFilePointer(spectrum_output_file, sizeof(header), NULL, FILE_BEGIN);
    }
#endif
    WriteFile(spectrum_output_file, spd_pixels, spd_pixel_data_size, &bytes_written, NULL);
    CloseHandle(spectrum_output_file);
    VirtualFree(spd_pixels, spd_pixel_data_size, MEM_RELEASE);

    //spd file -> bmp
    spectrum ref_white;
    spectrum cmf_x, cmf_y, cmf_z;
    const_spectrum(&ref_white, 1.0);
    load_csv_file_to_spectrum(&cmf_x, "spectra\\cmf_x.csv");
    load_csv_file_to_spectrum(&cmf_y, "spectra\\cmf_y.csv");
    load_csv_file_to_spectrum(&cmf_z, "spectra\\cmf_z.csv");

    spectrum_output_file = CreateFile(spectrum_output_path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    DWORD bytes_read;
    ReadFile(spectrum_output_file, &header, sizeof(header), &bytes_read, NULL);
    spd_file_to_bmp(spectrum_output_file, &header, bmp_output_path, &cmf_x, &cmf_y, &cmf_z, &ref_white);
    CloseHandle(spectrum_output_file);

    return 0;
}
