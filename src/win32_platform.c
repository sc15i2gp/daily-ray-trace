u32 write_bmp_to_file(windows_bmp *bmp, const char *file_path)
{
    file_handle output_file = open_file(file_path, ACCESS_WRITE, FILE_NEW);

    write_file(output_file, bmp->file_header->bfSize, bmp->file_header);
    close_file(output_file);
}

u32 write_pixels_to_bmp(rgb_u8 *pixels, u32 width, u32 height, const char *path)
{
    u32 pixels_size = width * height * sizeof(u32);
    u32 bmp_size = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + pixels_size;
    windows_bmp bmp = {0};
    u8 *raw_bmp = alloc(bmp_size);
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
    unalloc(raw_bmp, 0);

    return success;
}

void *alloc_pages(u32 num_pages)
{
    SYSTEM_INFO sys_inf;
    GetSystemInfo(&sys_inf);
    u32    page_size = sys_inf.dwPageSize;
    u32    size = num_pages * page_size;
    void   *p   = alloc(size);
    return p;
}

void *alloc(u32 size)
{
    void   *p = VirtualAlloc(NULL, size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    return p;
}

void unalloc(void *p, u32 p_size)
{
    VirtualFree(p, p_size, MEM_RELEASE);
}

file_handle open_file(const char *path, file_access_type access_type, file_open_type open_type)
{
    DWORD win_access_type = 0;
    DWORD win_share_mode  = 0;
    DWORD win_open_type   = 0;
    switch(access_type)
    {
        case ACCESS_READ:
        {
            win_access_type = GENERIC_READ;
            win_share_mode  = FILE_SHARE_READ;
            break;
        }
        case ACCESS_WRITE:
        {
            win_access_type = GENERIC_WRITE;
            win_share_mode  = FILE_SHARE_WRITE;
            break;
        }
        case ACCESS_READWRITE:
        {
            win_access_type = GENERIC_READ    | GENERIC_WRITE;
            win_share_mode  = FILE_SHARE_READ | FILE_SHARE_WRITE;
            break;
        }
    }
    switch(open_type)
    {
        case FILE_NEW:
        {
            win_open_type = CREATE_ALWAYS;
            break;
        }
        case FILE_EXISTS:
        {
            win_open_type = OPEN_EXISTING;
            break;
        }
    }
    file_handle file = CreateFile(path, win_access_type, win_share_mode, NULL, win_open_type, FILE_ATTRIBUTE_NORMAL, NULL);

    return file;
}

void close_file(file_handle file)
{
    CloseHandle(file);
}

u32 get_file_size(file_handle file)
{
    DWORD file_size = GetFileSize(file, NULL);
    return (u32)file_size;
}

void read_file(file_handle file, u32 read_size, void *dst)
{
    DWORD bytes_read;
    ReadFile(file, dst, read_size, &bytes_read, NULL);
}

void write_file(file_handle file, u32 write_size, void *src)
{
    DWORD bytes_written;
    WriteFile(file, src, write_size, &bytes_written, NULL);
}

void set_file_pointer(file_handle file, u32 loc)
{
    SetFilePointer(file, loc, NULL, FILE_BEGIN);
}

