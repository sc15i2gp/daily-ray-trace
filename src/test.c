#include "daily_ray_trace.h"

void write_test_bmp(u32 colour, const char *test_bmp_path)
{
    printf("Writing test bmp...");
    windows_bmp test_bmp              = {0};
    u32         test_bmp_width        = 300;
    u32         test_bmp_height       = 200;
    u32         number_of_test_pixels = test_bmp_width * test_bmp_height;
    u32         test_pixels_size      = number_of_test_pixels * sizeof(u32);
    u32         test_headers_size     = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
    u32         test_bmp_size         = test_headers_size + test_pixels_size;
    u8          *raw_test_bmp         = (u8*)VirtualAlloc(0, test_bmp_size, MEM_COMMIT, PAGE_READWRITE);

    test_bmp.file_header = (BITMAPFILEHEADER*)raw_test_bmp;
    test_bmp.info_header = (BITMAPINFOHEADER*)(raw_test_bmp + sizeof(BITMAPFILEHEADER));
    test_bmp.pixels      = (rgb_u8*)(raw_test_bmp + test_headers_size);

    for(u32 pixel = 0; pixel < number_of_test_pixels; ++pixel)
    {
        test_bmp.pixels[pixel].bgra[colour] = 255;
    }

    test_bmp.file_header->bfType          = 0x4d42;
    test_bmp.file_header->bfSize          = test_bmp_size;
    test_bmp.file_header->bfOffBits       = test_headers_size;
    test_bmp.info_header->biSize          = sizeof(BITMAPINFOHEADER);
    test_bmp.info_header->biWidth         = test_bmp_width;
    test_bmp.info_header->biHeight        = test_bmp_height;
    test_bmp.info_header->biBitCount      = 32;
    test_bmp.info_header->biCompression   = BI_RGB;
    test_bmp.info_header->biSizeImage     = 0;
    test_bmp.info_header->biXPelsPerMeter = 3780;
    test_bmp.info_header->biYPelsPerMeter = 3780;
    test_bmp.info_header->biClrUsed       = 0;
    test_bmp.info_header->biClrImportant  = 0;

    u32 test_success = write_bmp_to_file(&test_bmp, test_bmp_path);
    if(!test_success)  printf("Writing test bmp failed.\n");
    else               printf("Writing test bmp succeeded.\n");

    VirtualFree(raw_test_bmp, 0, MEM_RELEASE);
}

u32 test_rgb_spectrum_conversion(f64 low_rgb, f64 high_rgb, f64 rgb_iterand)
{
    //Steps:
    //  - Load spectra
    //  - Malloc for rgb inputs/outputs
    //  - Do the conversion for each rgb and store results
    //  - Compare results and get a figure for error

    printf("Testing spectral conversion...\n");
    
    spectrum white   = alloc_spd();
    spectrum red     = alloc_spd();
    spectrum green   = alloc_spd();
    spectrum blue    = alloc_spd();
    spectrum cyan    = alloc_spd();
    spectrum magenta = alloc_spd();
    spectrum yellow  = alloc_spd();
    spectrum cmf_x   = alloc_spd();
    spectrum cmf_y   = alloc_spd();
    spectrum cmf_z   = alloc_spd();
    spectrum s       = alloc_spd();
    const_spectrum(white, 1.0);
    zero_spectrum(s);
    load_csv_file_to_spectrum(red, "spectra\\red_rgb_to_spd.csv");
    load_csv_file_to_spectrum(green, "spectra\\green_rgb_to_spd.csv");
    load_csv_file_to_spectrum(blue, "spectra\\blue_rgb_to_spd.csv");
    load_csv_file_to_spectrum(cyan, "spectra\\cyan_rgb_to_spd.csv");
    load_csv_file_to_spectrum(magenta, "spectra\\magenta_rgb_to_spd.csv");
    load_csv_file_to_spectrum(yellow, "spectra\\yellow_rgb_to_spd.csv");
    load_csv_file_to_spectrum(cmf_x, "spectra\\cmf_x.csv");
    load_csv_file_to_spectrum(cmf_y, "spectra\\cmf_y.csv");
    load_csv_file_to_spectrum(cmf_z, "spectra\\cmf_z.csv");

    u32 rgb_it = ((u32)(high_rgb - low_rgb)/rgb_iterand) + 1U;
    u32 number_of_rgb_values = rgb_it * rgb_it * rgb_it;
    rgb_f64 *inputs = (rgb_f64*)VirtualAlloc(0, number_of_rgb_values * sizeof(rgb_f64), MEM_COMMIT, PAGE_READWRITE);
    rgb_f64 *results = (rgb_f64*)VirtualAlloc(0, number_of_rgb_values * sizeof(rgb_f64), MEM_COMMIT, PAGE_READWRITE);

    rgb_f64 rgb = {0.0, 0.0, 0.0};
    rgb_f64 rgb_result = {0.0, 0.0, 0.0};
    u32 rgb_value_index = 0;
    for(f64 r = low_rgb; r <= high_rgb; r += rgb_iterand)
    {
        rgb.r = r;
        for(f64 g = low_rgb; g <= high_rgb; g += rgb_iterand)
        {
            rgb.g = g;
            for(f64 b = low_rgb; b <= high_rgb; b += rgb_iterand)
            {
                rgb.b = b;
                rgb_f64_to_spectrum(rgb, s, white, red, green, blue, cyan, magenta, yellow);
                results[rgb_value_index] = spectrum_to_rgb_f64(s, cmf_x, cmf_y, cmf_z, white);
                inputs[rgb_value_index] = rgb;
                ++rgb_value_index;
            }
        }
    }


    //Info of interest:
    //  - Smallest, largest, average errors for both overall rgb and individual rgb channels
    //  - Locations + actual values corresponding to smallest, largest average errors
    f64 smallest_error = DBL_MAX;
    f64 largest_error = 0.0;
    f64 average_error = 0.0;
    rgb_f64 smallest_rgb_errors = {DBL_MAX, DBL_MAX, DBL_MAX};
    rgb_f64 largest_rgb_errors = {0.0, 0.0, 0.0};
    rgb_f64 average_rgb_errors = {0.0, 0.0, 0.0};
    for(u32 i = 0; i < number_of_rgb_values; ++i)
    {
        rgb_f64 input = inputs[i];
        rgb_f64 result = results[i];

        rgb_f64 rgb_error;
        rgb_error.r = fabs(result.r - input.r);
        rgb_error.g = fabs(result.g - input.g);
        rgb_error.b = fabs(result.b - input.b);
        
        f64 error = (rgb_error.r*rgb_error.r) + (rgb_error.g*rgb_error.g) + (rgb_error.b*rgb_error.b);
        error = sqrt(error);
        if(error < smallest_error)
        {
            smallest_error = error;
        }
        if(error > largest_error)
        {
            largest_error = error;
        }
        average_error += error;
        for(u32 j = 0; j < 3; ++j)
        {
            if(rgb_error.rgb[j] < smallest_rgb_errors.rgb[j])
            {
                smallest_rgb_errors.rgb[j] = rgb_error.rgb[j];
            }
            if(rgb_error.rgb[j] > largest_rgb_errors.rgb[j])
            {
                largest_rgb_errors.rgb[j] = rgb_error.rgb[j];
            }
            average_rgb_errors.rgb[j] += rgb_error.rgb[j];
        }
    }
    average_error /= ((f64)(number_of_rgb_values));
    average_rgb_errors.r /= ((f64)(number_of_rgb_values));
    average_rgb_errors.g /= ((f64)(number_of_rgb_values));
    average_rgb_errors.b /= ((f64)(number_of_rgb_values));
    printf("Smallest overall error: %f\n", smallest_error);
    printf("Largest overall error:  %f\n", largest_error);
    printf("Average overall error:  %f\n", average_error);
    printf("Smallest rgb errors: %f %f %f\n", smallest_rgb_errors.r, smallest_rgb_errors.g, smallest_rgb_errors.b);
    printf("Largest rgb errors:  %f %f %f\n", largest_rgb_errors.r, largest_rgb_errors.g, largest_rgb_errors.b);
    printf("Average rgb errors:  %f %f %f\n", average_rgb_errors.r, average_rgb_errors.g, average_rgb_errors.b);

    printf("Done testing spectral conversion.\n");

    u32 success = 1;
    return success;
}

u32 test_spectral_operations()
{
    u32 success = 1;
    
    //Test basic ops (sum, mul)
    
    //Test csv loading
    //  - Exact
    //  - Lerp
    //  - Spectral range greater than file
    //  - Spectral range shorter than file
    //  - Overlapping ranges both left and right

    return success;
}

void test_shape_intersection(u32 film_width, u32 film_height)
{
    printf("Testing shape intersection...\n");
    u32 film_pixel_count = film_width * film_height;
    u32 film_size = film_pixel_count * sizeof(f64);

    vec3 film_bottom_left = {-2.0, -2.0, 1.0};
    vec3 film_top_right = {2.0, 2.0, 1.0};
    f64 *film = (f64*)VirtualAlloc(0, film_size, MEM_COMMIT, PAGE_READWRITE);

    printf("Testing line-sphere intersection...\n");
    vec3 sphere_c = {0.0, 0.0, 0.0};
    f64 sphere_r = 1.0;

    f64 film_p_width = (film_top_right.x - film_bottom_left.x)/((f64)film_width);
    f64 film_p_height = (film_top_right.y - film_bottom_left.y)/((f64)film_height);
    vec3 line_d = {0.0, 0.0, -1.0};
    for(u32 y = 0; y < film_height; ++y)
    {
        for(u32 x = 0; x < film_width; ++x)
        {
            f64 line_o_x = film_bottom_left.x + ((f64)x + 0.5) * film_p_width;
            f64 line_o_y = film_bottom_left.y + ((f64)y + 0.5) * film_p_height;
            vec3 line_o = {line_o_x, line_o_y, film_bottom_left.z};
            f64 l = line_sphere_intersection(line_o, line_d, sphere_c, sphere_r);
            if(l != INFINITY) film[film_width*y + x] = 1.0 - l;
            else              film[film_width*y + x] = INFINITY;
        }
    }

    u32 film_pixels_size = film_pixel_count * sizeof(rgb_u8);
    rgb_u8 *film_pixels = (rgb_u8*)VirtualAlloc(0, film_pixels_size, MEM_COMMIT, PAGE_READWRITE);
    for(u32 i = 0; i < film_pixel_count; ++i)
    {
        if(film[i] != INFINITY)
        {
            film_pixels[i].r = (u8)(255.0 * film[i]);
            film_pixels[i].g = (u8)(255.0 * film[i]);
            film_pixels[i].b = (u8)(255.0 * film[i]);
            film_pixels[i].a = 255;

        }
        else
        {
            film_pixels[i].r = 255;
            film_pixels[i].g = 105;
            film_pixels[i].b = 184;
            film_pixels[i].a = 255;
        }
    }

    write_pixels_to_bmp(film_pixels, film_width, film_height, "output\\sphere_test.bmp");

    printf("Testing line-plane interaction...\n");
    memset(film, 0, film_size);
    memset(film_pixels, 0, film_pixels_size);

    vec3 plane_p, plane_u, plane_v, plane_n;
    vec3 p = {-0.5, -0.5, 0.0};
    vec3 u = {0.5, -0.5, 0.0};
    vec3 v = {-0.5, 0.5, 0.0};
    create_plane_from_points(p, u, v, &plane_p, &plane_u, &plane_v, &plane_n);

    for(u32 y = 0; y < film_height; ++y)
    {
        for(u32 x = 0; x < film_width; ++x)
        {
            f64 line_o_x = film_bottom_left.x + ((f64)x + 0.5) * film_p_width;
            f64 line_o_y = film_bottom_left.y + ((f64)y + 0.5) * film_p_height;
            vec3 line_o = {line_o_x, line_o_y, film_bottom_left.z};
            f64 l = line_plane_intersection(line_o, line_d, plane_p, plane_n, plane_u, plane_v);
            film[film_width*y + x] = l;
        }
    }
    for(u32 i = 0; i < film_pixel_count; ++i)
    {
        if(film[i] != INFINITY)
        {
            film_pixels[i].r = (u8)(255.0 * film[i]);
            film_pixels[i].g = (u8)(255.0 * film[i]);
            film_pixels[i].b = (u8)(255.0 * film[i]);
            film_pixels[i].a = 255;

        }
        else
        {
            film_pixels[i].r = 255;
            film_pixels[i].g = 105;
            film_pixels[i].b = 184;
            film_pixels[i].a = 255;
        }
    }
    write_pixels_to_bmp(film_pixels, film_width, film_height, "output\\plane_test.bmp");

    VirtualFree(film, 0, MEM_RELEASE);
    VirtualFree(film_pixels, 0, MEM_RELEASE);
}

void call_test_funcs()
{
    write_test_bmp(0, "output\\test_output_blue.bmp");
    write_test_bmp(1, "output\\test_output_green.bmp");
    write_test_bmp(2, "output\\test_output_red.bmp");

    test_rgb_spectrum_conversion(0.0, 1.0, 0.1);
    test_shape_intersection(1600, 1600);
}

int main()
{
    init_spd_table(16, 380.0, 720.0, 5.0);
    call_test_funcs();
}
