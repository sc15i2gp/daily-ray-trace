#include "daily_ray_trace.h"

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

//TODO:
//  - PRNG
//  - Line/shape intersection
//  - Raytrace algorithm

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
    //  - 
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

int main(int argc, char **argv)
{
    write_test_bmp(0, "output\\test_output_blue.bmp");
    write_test_bmp(1, "output\\test_output_green.bmp");
    write_test_bmp(2, "output\\test_output_red.bmp");

    test_rgb_spectrum_conversion(0.0, 1.0, 0.1);
    return 0;
}
