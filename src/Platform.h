#include <Windows.h>

#define BYTES(n) n
#define KILOBYTES(n) 1024 * BYTES(n)
#define MEGABYTES(n) 1024 * KILOBYTES(n)
#define GIGABYTES(n) 1024 * MEGABYTES(n)

int window_width(HWND);
int window_height(HWND);

void* alloc(int);
void dealloc(void*);
void zero_mem(void*, long unsigned int);

char* read_file_contents(const char* path);
void write_file_contents(const char* path, char* contents, int contents_size);

struct Timer
{
	LARGE_INTEGER start_time;
	LARGE_INTEGER stop_time;
};

unsigned int get_pc_time();
void query_pc_frequency();
void start_timer(Timer*);
void stop_timer(Timer*);
double cycles_to_s(long unsigned int);
double cycles_to_ms(long unsigned int);
double elapsed_time_in_s(Timer*);
double elapsed_time_in_ms(Timer*);
long unsigned int elapsed_time_in_cycles(Timer*);
