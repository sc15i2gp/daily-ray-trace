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

//TODO:
//  - Separate spectral code to another file
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
    rgb_u8 p8 = {.a = 0, .r = 1, .g = 2, .b = 3};
    rgb_f64 p64 = {.r = 0.1, .g = 0.2, .b = 0.3};
    printf("Hello World\n");
    printf("%u, %u, %u, %u\n", p8.a, p8.r, p8.g, p8.b);
    printf("%f, %f, %f\n", p64.r, p64.g, p64.b);
    printf("lerp(5.0, 3.0, 6.0, 5.0, 2.0) = %f, should be %f\n", lerp(5.0, 3.0, 6.0, 5.0, 2.0), 3.0);
    printf("Loading d65.csv...\n");

    spectrum white_spd, white_rgb_spd;
    spectrum red_spd, green_spd, blue_spd;
    spectrum cyan_spd, magenta_spd, yellow_spd;
    spectrum cmf_x, cmf_y, cmf_z;
    const_spectrum(&white_spd, 1.0);
    load_csv_file_to_spectrum(&white_rgb_spd, "spectra\\white_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&red_spd, "spectra\\red_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&green_spd, "spectra\\green_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&blue_spd, "spectra\\blue_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&cyan_spd, "spectra\\cyan_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&magenta_spd, "spectra\\magenta_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&yellow_spd, "spectra\\yellow_rgb_to_spd.csv");
    load_csv_file_to_spectrum(&cmf_x, "spectra\\cmf_x.csv");
    load_csv_file_to_spectrum(&cmf_y, "spectra\\cmf_y.csv");
    load_csv_file_to_spectrum(&cmf_z, "spectra\\cmf_z.csv");

    spectrum s;
    zero_spectrum(&s);
    print_spectrum(&s);
    rgb_f64_to_spectrum(p64, &s, &white_rgb_spd, &red_spd, &green_spd, &blue_spd, &cyan_spd, &magenta_spd, &yellow_spd);
    print_spectrum(&s);
    rgb_f64 f = spectrum_to_rgb_f64(&s, &cmf_x, &cmf_y, &cmf_z, &white_spd);

    printf("First rgb = (%f, %f, %f)\n", p64.r, p64.g, p64.b);
    printf("Final rgb = (%f, %f, %f)\n", f.r, f.g, f.b);
    return 0;
}
