typedef struct
{
    BITMAPFILEHEADER *file_header;
    BITMAPINFOHEADER *info_header;
    rgb_u8             *pixels; //bgra
}   windows_bmp;

u32 write_bmp_to_file(windows_bmp*, const char*);
u32 write_pixels_to_bmp(rgb_u8*, u32, u32, const char*);
