#include <stdio.h>
#include <stdint.h>
#include <Windows.h>
#include <math.h>

//Program description:
//  - Given a scene and sample parameters, produce an image
//  - Output bmp
//  - No windows/live preview for now

//Long term TODO:
//  - Investigate use + comparison of different RGB working spaces

//TODO:
//  - Simple spectral operations (sum, multiply etc.)
//      - + tests
//  - RGB to/from spectrum
//      - rgb to spectrum
//      - spectrum to xyz
//      - xyz to rgb
//  - Blackbody spd
//  - Separate spectral code to another file
//  - PRNG
//  - Line/shape intersection
//  - Raytrace algorithm

// dst += (src1 * d)
// Acc: dst += src
// Mul: dst *= d

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float    f32;
typedef double   f64;

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

//64 bit rgb values 0.0-1.0
typedef struct
{
    union
    {
        struct
        {
            f64 r;
            f64 g;
            f64 b;
        };
        struct
        {
            f64 x;
            f64 y;
            f64 z;
        };
        f64 rgb[3];
        f64 xyz[3];
    };
}   rgb_f64;

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

#define MAX_NUM_SPECTRUM_VALUES 128
//TODO: Is this info nececssary outside of i/o?
u32 number_of_spectrum_samples = 69;
f64 smallest_wavelength       = 380.0; //In nm
f64 largest_wavelength        = 720.0; //In nm
f64 sample_interval           = 5.0;   //In nm

//Spectrum as point samples over a range of wavelengths
//All spectra have the same number of samples
//All samples "sample[i]" are at the same wavelength (smallest_wavelength + i*sample_interval)
typedef struct
{
    f64 samples[MAX_NUM_SPECTRUM_VALUES];
}   spectrum;

void zero_spectrum_full(spectrum *dst)
{
    for(u32 i = 0; i < MAX_NUM_SPECTRUM_VALUES; ++i) dst->samples[i] = 0.0;
}

void zero_spectrum(spectrum *dst)
{
    for(u32 i = 0; i < number_of_spectrum_samples; ++i)
    {
        dst->samples[i] = 0.0;
    }
}

void zero_f64_array(f64 *dst_array, u32 dst_number_of_elements)
{
    for(u32 i = 0; i < dst_number_of_elements; ++i) dst_array[i] = 0.0;
}

void const_spectrum(spectrum *dst, f64 value)
{
    for(u32 i = 0; i < number_of_spectrum_samples; ++i) dst->samples[i] = value;
}

void spectral_mul_by_spectrum(spectrum *dst, spectrum *src0, spectrum *src1)
{
    for(u32 i = 0; i < number_of_spectrum_samples; ++i)
    {
        dst->samples[i] = src0->samples[i] * src1->samples[i];
    }
}

void spectral_mul_by_scalar(spectrum *dst, spectrum *src, f64 d)
{
    for(u32 i = 0; i < number_of_spectrum_samples; ++i)
    {
        dst->samples[i] = src->samples[i] * d;
    }
}

void spectral_acc(spectrum *dst, spectrum *src)
{
    for(u32 i = 0; i < number_of_spectrum_samples; ++i)
    {
        dst->samples[i] += src->samples[i];
    }
}

void spectral_acc_mul_by_scalar(spectrum *dst, spectrum *src, f64 d)
{
    for(u32 i = 0; i < number_of_spectrum_samples; ++i)
    {
        dst->samples[i] += src->samples[i] * d;
    }
}

void spectral_sum(spectrum *dst, spectrum *src0, spectrum *src1)
{
    for(u32 i = 0; i < number_of_spectrum_samples; ++i)
    {
        dst->samples[i] = src0->samples[i] + src1->samples[i];
    }
}

rgb_f64 spectrum_to_xyz(spectrum *spd, spectrum *cmf_x, spectrum *cmf_y, spectrum *cmf_z, spectrum *ref_white)
{
    rgb_f64 xyz = {0.0, 0.0, 0.0};
    f64 n = 0.0;
    for(u32 i = 0; i < number_of_spectrum_samples; ++i)
    {
        n += (cmf_y->samples[i] * ref_white->samples[i]);
    }
    n *= sample_interval;

    for(u32 i = 0; i < number_of_spectrum_samples; ++i)
    {
        xyz.x += (cmf_x->samples[i] * spd->samples[i] * ref_white->samples[i]);
        xyz.y += (cmf_y->samples[i] * spd->samples[i] * ref_white->samples[i]);
        xyz.z += (cmf_z->samples[i] * spd->samples[i] * ref_white->samples[i]);
    }
    xyz.x *= (sample_interval/n);
    xyz.y *= (sample_interval/n);
    xyz.z *= (sample_interval/n);

    return xyz;
}

rgb_f64 spectrum_to_rgb_f64(spectrum *spd, spectrum *ref_white, spectrum *cmf_x, spectrum *cmf_y, spectrum *cmf_z)
{
    rgb_f64 spd_xyz = spectrum_to_xyz(spd, cmf_x, cmf_y, cmf_z, ref_white);
    
    rgb_f64 linear_rgb;
    linear_rgb.r = (2.3706743 * spd_xyz.x) - (0.9000405 * spd_xyz.y) - (0.4706338 * spd_xyz.z);
    linear_rgb.g = (-0.5138850 * spd_xyz.x) + (1.4253036 * spd_xyz.y) + (0.0885814 * spd_xyz.z);
    linear_rgb.b = (0.0052982 * spd_xyz.x) - (0.0146949 * spd_xyz.y) + (1.0093968 * spd_xyz.z);

    return linear_rgb;
}

//RGB to spectrum
//  - smallest_rgb
//  - mid_rgb
//  - largest_rgb
//  - rgb_to_spd is spectrum of largest_rgb to spd
//  - cym_to_spd is:
//      - If R is smallest_rgb then cyan
//      - If G is smallest_rgb then magenta
//      - If B is smallest_rgb then yellow
//  - spectrum = smallest_rgb * white_rgb_spd + cym_spd * (mid_rgb - smallest_rgb) + rgb_spd * (largest_rgb - mid_rgb)

void rgb_f64_to_spectrum(   rgb_f64 rgb,         spectrum *dst, 
                            spectrum *white_spd, spectrum *red_spd, 
                            spectrum *green_spd, spectrum *blue_spd,
                            spectrum *cyan_spd,  spectrum *magenta_spd,
                            spectrum *yellow_spd)
{
    spectrum *rgb_spectra[] = {red_spd,  green_spd,   blue_spd};
    spectrum *cmy_spectra[] = {cyan_spd, magenta_spd, yellow_spd};
    f64 rgb_values[] =        {rgb.r,    rgb.g,       rgb.b};
    u32 indices[] = {0, 1, 2};
    u32 tmp_val;

    if(rgb_values[indices[0]] > rgb_values[indices[1]])
    {
        tmp_val = indices[1];
        indices[1] = indices[0];
        indices[0] = tmp_val;
    }
    if(rgb_values[indices[1]] > rgb_values[indices[2]])
    {
        tmp_val = indices[2];
        indices[2] = indices[1];
        indices[1] = tmp_val;
    }
    if(rgb_values[indices[0]] > rgb_values[indices[1]])
    {
        tmp_val = indices[1];
        indices[1] = indices[0];
        indices[0] = tmp_val;
    }

    u32 small_index = indices[0];
    u32 mid_index = indices[1];
    u32 large_index = indices[2];
    f64 rgb_diff_mid_small = rgb_values[mid_index] - rgb_values[small_index];
    f64 rgb_diff_large_mid = rgb_values[large_index] - rgb_values[mid_index];
    spectral_mul_by_scalar(dst, white_spd, rgb_values[small_index]); //dst = white_rgb_spd * smallest_rgb
    spectral_acc_mul_by_scalar(dst, cmy_spectra[small_index], rgb_diff_mid_small); //dst += cym_spd * (mid - small)
    spectral_acc_mul_by_scalar(dst, rgb_spectra[large_index], rgb_diff_large_mid); //dst += rgb_spd * (large - mid)
}

char *find_next_newline(char *c)
{
    while(1)
    {
        if(*c == '\n') return c;
        if(*c == EOF || *c == 0) return NULL;
        ++c;
    }
}

char *find_next_number(char *c)
{
    while(1)
    {
        if(*c >= '0' && *c <= '9') return c;
        if(*c == EOF || *c == 0) return NULL;
        ++c;
    }
}

char *find_next_char(char *c, char desired)
{
    while(1)
    {
        if(*c == desired) return c;
        if(*c == EOF || *c == 0) return NULL;
        ++c;
    }
}

f64 lerp(f64 x, f64 x0, f64 x1, f64 y0, f64 y1)
{
    return y0 + ((x - x0) * ((y1 - y0) / (x1 - x0)));
}

//CSV file structure:
//  - First line ignored (for headings)
//  - Each line is a single sample
//  - A sample consists of two csvs: "wavelength, value"
u32 load_csv_file_to_spectrum(spectrum *dst, const char *csv_path)
{
    HANDLE csv_file_handle = CreateFile(csv_path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if(csv_file_handle == INVALID_HANDLE_VALUE) return 0;

    printf("Opening csv file %s...\n", csv_path);
    DWORD csv_bytes_read = 0;
    DWORD csv_file_size = GetFileSize(csv_file_handle, NULL) + 1;
    char *csv_file_buffer = (char*)VirtualAlloc(0, csv_file_size, MEM_COMMIT, PAGE_READWRITE);
    printf("Reading csv file %s...\n", csv_path);
    u32 success = ReadFile(csv_file_handle, csv_file_buffer, csv_file_size, &csv_bytes_read, NULL);
    csv_file_buffer[csv_file_size] = 0;
    CloseHandle(csv_file_handle);
    printf("Closed csv file %s.\n", csv_path);

    f64 file_wavelengths[MAX_NUM_SPECTRUM_VALUES];
    f64 file_wavelength_values[MAX_NUM_SPECTRUM_VALUES];
    zero_f64_array(file_wavelengths, MAX_NUM_SPECTRUM_VALUES);
    zero_f64_array(file_wavelength_values, MAX_NUM_SPECTRUM_VALUES);
    u32 file_number_of_samples = 0;

    //TODO: Add support for files with wavelengths in um
    //NOTE: Currently only supporting wavelengths in nm
    
    //Count number of samples in csv file
    printf("Counting number of samples in char buffer...\n");
    for(char *c = find_next_newline(csv_file_buffer); c != NULL; c = find_next_newline(c)) 
    {
        ++c;
        if(find_next_number(c)) ++file_number_of_samples;
    }
    //Read wavelengths and sample values
    printf("Reading sample data in char buffer...\n");
    char *c = find_next_newline(csv_file_buffer) + 1;
    for(u32 i = 0; i < file_number_of_samples; ++i)
    {
        c = find_next_number(c);
        file_wavelengths[i] = atof(c);
        c = find_next_char(c, ',');
        c = find_next_number(c);
        file_wavelength_values[i] = atof(c);
        c = find_next_newline(c);
    }
    printf("Freeing char buffer...\n");
    VirtualFree(csv_file_buffer, csv_file_size, MEM_RELEASE);

    printf("Number of samples read = %u\n", file_number_of_samples);

    //Lerp values and write them to spectrum
    u32 file_wl_index = 0;
    for(u32 sample = 0; sample < number_of_spectrum_samples; ++sample)
    {
        f64 sample_wavelength = smallest_wavelength + ((f64)sample) * sample_interval;
        for(; file_wavelengths[file_wl_index+1] < sample_wavelength; ++file_wl_index);
        dst->samples[sample] = lerp(sample_wavelength, file_wavelengths[file_wl_index], file_wavelengths[file_wl_index + 1], file_wavelength_values[file_wl_index], file_wavelength_values[file_wl_index+1]);
    }

    /*
    printf("Wavelength | Value\n");
    for(u32 i = 0; i < file_number_of_samples; ++i)
    {
        f64 sample_wl = smallest_wavelength + ((f64)i) * sample_interval;
        printf("%f | %f\n", sample_wl, dst->samples[i]);
    }
*/
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
    rgb_f64_to_spectrum(p64, &s, &white_spd, &red_spd, &green_spd, &blue_spd, &cyan_spd, &magenta_spd, &yellow_spd);
    print_spectrum(&s);
    rgb_f64 f = spectrum_to_rgb_f64(&s, &white_spd, &cmf_x, &cmf_y, &cmf_z);

    printf("First rgb = (%f, %f, %f)\n", p64.r, p64.g, p64.b);
    printf("Final rgb = (%f, %f, %f)\n", f.r, f.g, f.b);
    return 0;
}
