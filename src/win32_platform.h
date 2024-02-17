typedef struct
{
    BITMAPFILEHEADER *file_header;
    BITMAPINFOHEADER *info_header;
    rgb_u8             *pixels; //bgra
}   windows_bmp;

typedef HANDLE file_handle;
//TODO: Use bit flags
typedef enum
{
    ACCESS_NONE,
    ACCESS_READ,
    ACCESS_WRITE,
    ACCESS_READWRITE
} file_access_type;

typedef enum
{
    OPEN_NONE,
    FILE_NEW,
    FILE_EXISTS
} file_open_type;

typedef struct
{
    u32 id;
    u32 width_in_pixels;
    u32 height_in_pixels;
    u32 number_of_wavelengths;
    f64 min_wavelength;
    f64 wavelength_interval;
} spd_file_header;

void *alloc(u32);
void unalloc(void *, u32 size);
void *alloc_pages(u32);

file_handle open_file(const char *path, file_access_type access, file_open_type open_type);
void close_file(file_handle file);
u32 get_file_size(file_handle file);
void read_file(file_handle file, u32 read_size, void *buffer);
void write_file(file_handle file, u32 write_size, void *buffer);
void set_file_pointer(file_handle file, u32 loc);

u32 write_bmp_to_file(windows_bmp*, const char*);
u32 write_pixels_to_bmp(rgb_u8*, u32, u32, const char*);
void spd_file_to_bmp(const char *spd_path, const char *bmp_path);
