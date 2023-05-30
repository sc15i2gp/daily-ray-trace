#include <stdio.h>
#include <stdint.h>
#include <Windows.h>

//Program description:
//  - Given a scene and sample parameters, produce an image
//  - Output bmp
//  - No windows/live preview for now

//TODO:
//  - Simple spectral operations (sum, multiply etc.)
//      - + tests
//  - RGB to/from spectrum
//      - Functions
//      - Analyse error/correctness
//      - Alternative methods (e.g. table rep of matching funcs vs gaussian approx)
//      - Tests
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
    f64 r;
    f64 g;
    f64 b;
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

void zero_f64_array(f64 *dst_array, u32 dst_number_of_elements)
{
    for(u32 i = 0; i < dst_number_of_elements; ++i) dst_array[i] = 0.0;
}

void const_spectrum(spectrum *dst, f64 value)
{
    for(u32 i = 0; i < number_of_spectrum_samples; ++i) dst->samples[i] = value;
}

void spectral_sum(spectrum *dst, spectrum *src)
{
    for(u32 sample = 0; sample < number_of_spectrum_samples; ++sample)
    {
        dst->samples[sample] += src->samples[sample];
    }
}

void spectral_mul_by_scalar(spectrum *dst, f64 scalar)
{
    for(u32 sample = 0; sample < number_of_spectrum_samples; ++sample)
    {
        dst->samples[sample] *= src->samples[sample];
    }
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

    printf("Wavelength | Value\n");
    for(u32 i = 0; i < file_number_of_samples; ++i)
    {
        f64 sample_wl = smallest_wavelength + ((f64)i) * sample_interval;
        printf("%f | %f\n", sample_wl, dst->samples[i]);
    }

    return success;
}

u32 test_spectral_operations()
{
    u32 success = 1;

    return success;
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
    spectrum test_spectrum_0;
    zero_spectrum_full(&test_spectrum_0);
    u32 d65_success = load_csv_file_to_spectrum(&test_spectrum_0, "spectra\\d65.csv");
    if(d65_success) printf("Loaded d65.csv\n");
    else printf("Failed to load d65.csv\n");
    spectrum test_spectrum_1;
    zero_spectrum_full(&test_spectrum_1);

    return 0;
}
