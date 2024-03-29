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
    LARGE_INTEGER start;
    LARGE_INTEGER stop;
} timer;

f64 pc_frequency;

void start_timer(timer*);
void stop_timer(timer*);
f64  time_elapsed_in_ms(timer*);

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
void write_rgb_f64_pixels_to_bmp(rgb_f64*, u32, u32, const char *);
