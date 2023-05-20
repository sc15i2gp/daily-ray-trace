#include <stdio.h>
#include <stdint.h>
#include <Windows.h>

//Program description:
//  - Given a scene and sample parameters, produce an image
//  - Output bmp
//  - No windows/live preview for now

//TODO:
//  - Write to bmp
//  - Spectral stuff
//  - Line/shape intersection
//  - Raytrace algorithm

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float  f32;
typedef double f64;

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
        u8 bgra[4];
    };
} rgb8;

//64 bit rgb values 0.0-1.0
typedef struct
{
    f64 r;
    f64 g;
    f64 b;
} rgb64;

typedef struct
{
    BITMAPFILEHEADER *file_header;
    BITMAPINFOHEADER *info_header;
    rgb8             *pixels; //bgra
} windows_bmp;

u32 write_bmp_to_file(windows_bmp *bmp, const char *file_path)
{
    DWORD bytes_written = 0;
    HANDLE output_file = CreateFile(file_path, GENERIC_WRITE, FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if(output_file == INVALID_HANDLE_VALUE) return 0;

    BOOL success = WriteFile(output_file, (void*)bmp->file_header, bmp->file_header->bfSize, &bytes_written, NULL);
    CloseHandle(output_file);

    return success;
}

void write_test_bmp(u32 colour, const char *test_bmp_path)
{
    printf("Writing test bmp...");
    windows_bmp test_bmp = {0};
    u32 test_bmp_width = 300;
    u32 test_bmp_height = 200;
    u32 number_of_test_pixels = test_bmp_width * test_bmp_height;
    u32 test_pixels_size = number_of_test_pixels * sizeof(u32);
    u32 test_headers_size = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
    u32 test_bmp_size = test_headers_size + test_pixels_size;
    u8 *raw_test_bmp = (u8*)VirtualAlloc(0, test_bmp_size, MEM_COMMIT, PAGE_READWRITE);

    test_bmp.file_header = (BITMAPFILEHEADER*)raw_test_bmp;
    test_bmp.info_header = (BITMAPINFOHEADER*)(raw_test_bmp + sizeof(BITMAPFILEHEADER));
    test_bmp.pixels = (rgb8*)(raw_test_bmp + test_headers_size);

    for(u32 pixel = 0; pixel < number_of_test_pixels; ++pixel)
    {
        test_bmp.pixels[pixel].bgra[colour] = 255;
    }

    test_bmp.file_header->bfType = 0x4d42;
    test_bmp.file_header->bfSize = test_bmp_size;
    test_bmp.file_header->bfOffBits = test_headers_size;
    test_bmp.info_header->biSize = sizeof(BITMAPINFOHEADER);
    test_bmp.info_header->biWidth = test_bmp_width;
    test_bmp.info_header->biHeight = test_bmp_height;
    test_bmp.info_header->biBitCount = 32;
    test_bmp.info_header->biCompression = BI_RGB;
    test_bmp.info_header->biSizeImage = 0;
    test_bmp.info_header->biXPelsPerMeter = 3780;
    test_bmp.info_header->biYPelsPerMeter = 3780;
    test_bmp.info_header->biClrUsed = 0;
    test_bmp.info_header->biClrImportant = 0;

    u32 test_success = write_bmp_to_file(&test_bmp, test_bmp_path);
    if(!test_success) printf("Writing test bmp failed.\n");
    else printf("Writing test bmp succeeded.\n");

    VirtualFree(raw_test_bmp, test_bmp_size, MEM_RELEASE);
}

int main(int argc, char **argv)
{
    write_test_bmp(0, "output\\test_output_blue.bmp");
    write_test_bmp(1, "output\\test_output_green.bmp");
    write_test_bmp(2, "output\\test_output_red.bmp");
    rgb8 p8 = {.a = 0, .r = 1, .g = 2, .b = 3};
    rgb64 p64 = {.r = 0.1, .g = 0.2, .b = 0.3};
    printf("Hello World\n");
    printf("%u, %u, %u, %u\n", p8.a, p8.r, p8.g, p8.b);
    printf("%f, %f, %f\n", p64.r, p64.g, p64.b);
    return 0;
}
