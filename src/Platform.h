#include <Windows.h>

#define BYTES(n) n
#define KILOBYTES(n) 1024 * BYTES(n)
#define MEGABYTES(n) 1024 * KILOBYTES(n)
#define GIGABYTES(n) 1024 * MEGABYTES(n)

void* alloc(int);
void dealloc(void*);

char* read_file_contents(const char* path);
void write_file_contents(const char* path, char* contents, int contents_size);
