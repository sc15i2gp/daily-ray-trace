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

//TODO:
//  - Xorshift rng
//  - Raytrace algorithm
//      - Sample pixel
//      - Sample lens(?)
//      - Cast ray
//      - Visibility between two points
//      - Pixel filtering

// dst += (src1 * d)
// Acc: dst += src
// Mul: dst *= d

//8 bit rgb values 0-255
typedef struct
{
    union
    {
        u32 value;
        struct
        {
            u8 b;
            u8 g;
            u8 r;
            u8 a;
        };
        u8  bgra[4];
    };
}   rgb_u8;

rgb_u8 rgb_f64_to_rgb_u8(rgb_f64 in_rgb)
{
    rgb_u8 out_rgb;
    out_rgb.r = (u8)(in_rgb.r * 255.0);
    out_rgb.g = (u8)(in_rgb.g * 255.0);
    out_rgb.b = (u8)(in_rgb.b * 255.0);

    return out_rgb;
}

typedef struct
{
    BITMAPFILEHEADER *file_header;
    BITMAPINFOHEADER *info_header;
    rgb_u8             *pixels; //bgra
}   windows_bmp;

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

    VirtualFree(raw_test_bmp, test_bmp_size, MEM_RELEASE);
}

void zero_f64_array(f64 *dst_array, u32 dst_number_of_elements)
{
    for(u32 i = 0; i < dst_number_of_elements; ++i) dst_array[i] = 0.0;
}

u32 test_rgb_spectrum_conversion(f64 low_rgb, f64 high_rgb, f64 rgb_iterand)
{
    //Steps:
    //  - Load spectra
    //  - Malloc for rgb inputs/outputs
    //  - Do the conversion for each rgb and store results
    //  - Compare results and get a figure for error

    printf("Testing spectral conversion...\n");
    
    spectrum white;
    spectrum red, green, blue;
    spectrum cyan, magenta, yellow;
    spectrum cmf_x, cmf_y, cmf_z;
    const_spectrum(&white, 1.0);
    load_csv_file_to_spectrum(&red, "spectra\\red_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&green, "spectra\\green_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&blue, "spectra\\blue_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&cyan, "spectra\\cyan_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&magenta, "spectra\\magenta_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&yellow, "spectra\\yellow_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&cmf_x, "spectra\\cmf_x.csv");
    load_csv_file_to_spectrum(&cmf_y, "spectra\\cmf_y.csv");
    load_csv_file_to_spectrum(&cmf_z, "spectra\\cmf_z.csv");

    u32 rgb_it = ((u32)(high_rgb - low_rgb)/rgb_iterand) + 1U;
    u32 number_of_rgb_values = rgb_it * rgb_it * rgb_it;
    rgb_f64 *inputs = (rgb_f64*)VirtualAlloc(0, number_of_rgb_values * sizeof(rgb_f64), MEM_COMMIT, PAGE_READWRITE);
    rgb_f64 *results = (rgb_f64*)VirtualAlloc(0, number_of_rgb_values * sizeof(rgb_f64), MEM_COMMIT, PAGE_READWRITE);

    rgb_f64 rgb = {0.0, 0.0, 0.0};
    rgb_f64 rgb_result = {0.0, 0.0, 0.0};
    spectrum s;
    zero_spectrum(&s);
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
                rgb_f64_to_spectrum(rgb, &s, &white, &red, &green, &blue, &cyan, &magenta, &yellow);
                results[rgb_value_index] = spectrum_to_rgb_f64(&s, &cmf_x, &cmf_y, &cmf_z, &white);
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

void print_spectrum(spectrum *spd)
{
    printf("Wavelength | Value\n");
    for(u32 i = 0; i < number_of_spectrum_samples; ++i)
    {
        f64 wl = smallest_wavelength + (((f64)i) * sample_interval);
        printf("%f | %f\n", wl, spd->samples[i]);
    }
}

typedef struct
{
    f64 x;
    f64 y;
    f64 z;
} vec3;

vec3 vec3_sum(vec3 a, vec3 b)
{
    vec3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

vec3 vec3_sub(vec3 a, vec3 b)
{
    vec3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}

f64 vec3_dot(vec3 a, vec3 b)
{
    f64 result = a.x*b.x + a.y*b.y + a.z*b.z;
    return result;
}

vec3 vec3_mul_by_f64(vec3 v, f64 f)
{
    vec3 result;
    result.x = f * v.x;
    result.y = f * v.y;
    result.z = f * v.z;
    return result;
}

vec3 vec3_div_by_f64(vec3 v, f64 f)
{
    vec3 result;
    result.x = v.x / f;
    result.y = v.y / f;
    result.z = v.z / f;
    return result;
}

f64 vec3_length(vec3 v)
{
    f64 result = vec3_dot(v, v);
    result = sqrt(result);
    return result;
}

vec3 vec3_normalise(vec3 v)
{
    vec3 result = v;
    result = vec3_div_by_f64(v, vec3_length(v));
    return result;
}

f64 point_to_line_distance(vec3 p, vec3 l_0, vec3 l_1)
{
    vec3 normalised_l = vec3_normalise(vec3_sub(l_0, l_1));
    vec3 q = vec3_sum(l_1, vec3_mul_by_f64(normalised_l, vec3_dot(vec3_sub(p, l_1), normalised_l)));
    vec3 q_to_p = vec3_sub(p, q);
    f64 distance = vec3_dot(q_to_p, q_to_p);
    return distance;
}

//line_o: line origin
//line_d: normalised line direction
//sphere_c: sphere center
//sphere_r: sphere radius
//Returns the smallest positive intersection point from line_o along line_d through the sphere
//Returns NaN if there is no intersection point, or if there are only negative intersection points
f64 line_sphere_intersection(vec3 line_o, vec3 line_d, vec3 sphere_c, f64 sphere_r)
{
    vec3 c_to_o = vec3_sub(line_o, sphere_c);

    f64 a               = 1.0;
    f64 b               = -2.0 * vec3_dot(c_to_o, line_d);
    f64 c               = vec3_dot(c_to_o, c_to_o) - sphere_r*sphere_r;
    f64 discriminant    = b*b - 4.0 * a * c;

    if(discriminant < 0.0) return NAN;

    f64 sqrt_discriminant = sqrt(discriminant);

    f64 a_2 = 2.0 * a;

    f64 solution_0 = (b + sqrt_discriminant)/a_2;
    f64 solution_1 = (b - sqrt_discriminant)/a_2;

    if      (solution_0 < 0.0 && solution_1 < 0.0)  return NAN;         //No positive solutions
    else if (solution_0 >= 0.0 && solution_1 < 0.0) return solution_0;  //Only solution_0 positive
    else if (solution_1 >= 0.0 && solution_0 < 0.0) return solution_1;  //Only solution_1 positive
    else if (solution_0 <= solution_1)              return solution_0;  //Smallest solution_0
    else                                            return solution_1;  //Smallest solution_1
}

//plane consists of a reference point as origin, a normal and two bounding vectors extending from reference
//line_o: line origin
//line_d: normalised line direction
//plane_p: plane reference point
//plane_n: normalised plane normal
//plane_u: first plane bounds vector
//plane_v: second plane bounds vector
//Returns positive intersection point from line_o along line_d through the plane
//Returns NaN if there is no positive intersection point
f64 line_plane_intersection(vec3 line_o, vec3 line_d, vec3 plane_p, vec3 plane_n, vec3 plane_u, vec3 plane_v)
{
    if(vec3_dot(line_d, plane_n) == 0.0) return NAN; //Line parallel to plane
    
    vec3    o_to_p  = vec3_sub(plane_p, line_o);
    f64     l       = vec3_dot(o_to_p, plane_n) / vec3_dot(line_d, plane_n);
    vec3    i       = vec3_sum(line_o, vec3_mul_by_f64(line_d, l));
    vec3    j       = vec3_sub(i, plane_p);

    f64 plane_u_length          = vec3_length(plane_u);
    f64 plane_v_length          = vec3_length(plane_v);
    f64 j_dot_plane_u           = vec3_dot(j, plane_u);
    f64 j_dot_plane_v           = vec3_dot(j, plane_v);

    if( 0.0 <= j_dot_plane_u && j_dot_plane_u <= plane_u_length &&
        0.0 <= j_dot_plane_v && j_dot_plane_v <= plane_v_length)
    {
        return l;
    }

    return NAN;
}

//triangle consists of three points in CCW order when viewed down from the normal
//line_o: line origin
//line_d: normalised line direction
f64 line_triangle_intersection(vec3 line_o, vec3 line_d, vec3 triangle_a, vec3 triangle_b, vec3 triangle_c, vec3 triangle_n)
{
    vec3 o_to_a = vec3_sub(triangle_a, line_o);
    f64 l = vec3_dot(o_to_a, triangle_n)/vec3_dot(line_d, triangle_n);
    vec3 i = vec3_sum(line_o, vec3_mul_by_f64(line_d, l));

    f64 a = point_to_line_distance(i, triangle_b, triangle_c)/point_to_line_distance(triangle_a, triangle_b, triangle_c);
    f64 b = point_to_line_distance(i, triangle_c, triangle_a)/point_to_line_distance(triangle_b, triangle_c, triangle_a);
    f64 c = point_to_line_distance(i, triangle_a, triangle_b)/point_to_line_distance(triangle_c, triangle_a, triangle_b);
    f64 d = a + b + c;

    if(d > 1.0) return NAN;

    return l;
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
            film[film_width*y + x] = 1.0 - l;
        }
    }

    u32 film_pixels_size = film_pixel_count * sizeof(rgb_u8);
    rgb_u8 *film_pixels = (rgb_u8*)VirtualAlloc(0, film_pixels_size, MEM_COMMIT, PAGE_READWRITE);
    for(u32 i = 0; i < film_pixel_count; ++i)
    {
        if(!isnan(film[i]))
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

    vec3 plane_p = {-0.5, -0.5, 0.0};
    vec3 plane_u = {1.0, 0.0, 0.0};
    vec3 plane_v = {0.0, 1.0, 0.0};
    vec3 plane_n = {0.0, 0.0, 1.0};

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
        if(!isnan(film[i]))
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

    VirtualFree(film, film_size, MEM_RELEASE);
    VirtualFree(film_pixels, film_pixels_size, MEM_RELEASE);
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

void call_test_funcs()
{
    write_test_bmp(0, "output\\test_output_blue.bmp");
    write_test_bmp(1, "output\\test_output_green.bmp");
    write_test_bmp(2, "output\\test_output_red.bmp");

    test_rgb_spectrum_conversion(0.0, 1.0, 0.1);
    test_shape_intersection(1600, 1600);
}

int main(int argc, char **argv)
{
    call_test_funcs();

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
    u32 spd_pixel_data_size = 400 * 300 * sizeof(spectrum);

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

    spectrum d;
    //for(u32 i = 0; i < number_of_spectrum_samples; ++i) d.samples[i] = 0.8;
    load_csv_file_to_spectrum(&d, "spectra\\red_rgb_to_spd.csv");
    for(u32 pixel = 0; pixel < number_of_image_pixels; ++pixel)
    {
        u32 spd_pixel = pixel % allocated_pixels;
        copy_spectrum(&spd_pixels[spd_pixel], &d);
        if(spd_pixel == allocated_pixels - 1)
        {
            WriteFile(spectrum_output_file, spd_pixels, spd_pixel_data_size, &bytes_written, NULL);
        }
    }
    CloseHandle(spectrum_output_file);
    VirtualFree(spd_pixels, spd_pixel_data_size, MEM_RELEASE);

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
