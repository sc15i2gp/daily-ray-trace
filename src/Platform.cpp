#include "Platform.h"

RECT window_rect(HWND window)
{
	RECT rect = {};
	GetClientRect(window, &rect);
	return rect;
}

int window_width(HWND window)
{
	RECT rect = window_rect(window);
	return rect.right - rect.left;
}

int window_height(HWND window)
{
	RECT rect = window_rect(window);
	return rect.bottom - rect.top;
}

void* alloc(int size)
{
	return VirtualAlloc(0, size, MEM_COMMIT, PAGE_READWRITE);
}

void dealloc(void* ptr)
{
	VirtualFree(ptr, 0, MEM_RELEASE);
}

char* read_file_contents(const char* path)
{
	HANDLE file = CreateFile(path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	DWORD size = GetFileSize(file, NULL);
	char* contents = (char*)alloc(size);
	ZeroMemory(contents, size);
	DWORD bytes_read = 0;
	ReadFile(file, contents, size, &bytes_read, NULL);

	CloseHandle(file);
	return contents;
}

void write_file_contents(const char* path, char* contents, int contents_size)
{
	DWORD bytes_written = 0;
	HANDLE file = CreateFile(path, GENERIC_WRITE, FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	WriteFile(file, contents, contents_size, &bytes_written, NULL);
	CloseHandle(file);
}

double pc_frequency = 0.0; //In counts/s

void query_pc_frequency()
{
	LARGE_INTEGER pc = {};
	QueryPerformanceFrequency(&pc); //NOTE: FAILURE POINT if function fails (returns false)
	pc_frequency = (double)pc.QuadPart;
}

void start_timer(Timer* t)
{
	QueryPerformanceCounter(&t->start_time);
}

void stop_timer(Timer* t)
{
	QueryPerformanceCounter(&t->stop_time);
}

double elapsed_time_in_s(Timer* t)
{
	double elapsed = (double)(t->stop_time.QuadPart - t->start_time.QuadPart);
	return elapsed/pc_frequency;
}

double elapsed_time_in_ms(Timer* t)
{
	return elapsed_time_in_s(t) * 1000.0;
}
