#include <Windows.h>
#include <stdio.h>
#include <stdint.h>
#include "Maths.h"

//TODO: LONGTERM
//	- Timing/performance
//	- Implement printf/sprintf/etc.
//	- Implement maths functions (trig, pow, ln etc.)
//	- Remove CRT

//TODO: Solve rendering equation
//	- Recursive raytrace
//	- Physically(ish) based
//	- Output result to window

//TODO: Intersection tests
//	- If eye ray intersects sphere/plane, colour pixel black

//TODO: Image plane

//DOING: Vector maths
//	- Test RemedyDBG


bool running = true;
void* back_buffer_mem = NULL;
int back_buffer_width = 0;
int back_buffer_height = 0;
BITMAPINFO back_buffer_info = {};

inline
void set_pixel_channel(uint32_t* pixel, uint8_t channel, uint8_t value)
{
	*(((uint8_t*)pixel) + channel) = value;
}

void render_that_good_shit_right_there()
{
	uint32_t* pixels = (uint32_t*)back_buffer_mem;
	uint32_t stride = back_buffer_width;
	for(int y = 0; y < back_buffer_height; ++y)
	{
		for(int x = 0; x < back_buffer_width; ++x)
		{
			uint32_t* pixel = pixels + y*stride + x;
			set_pixel_channel(pixel, 0, (uint8_t)(x));
			set_pixel_channel(pixel, 1, (uint8_t)(y));
			set_pixel_channel(pixel, 2, 0);
		}
	}
}

void size_window_back_buffer(int new_back_buffer_width, int new_back_buffer_height)
{
	if(back_buffer_mem) VirtualFree(back_buffer_mem, 0, MEM_RELEASE);

	back_buffer_info.bmiHeader.biSize = sizeof(back_buffer_info.bmiHeader);
	back_buffer_info.bmiHeader.biWidth = new_back_buffer_width;
	back_buffer_info.bmiHeader.biHeight = -new_back_buffer_height;
	back_buffer_info.bmiHeader.biPlanes = 1;
	back_buffer_info.bmiHeader.biBitCount = 32;
	back_buffer_info.bmiHeader.biCompression = BI_RGB;

	back_buffer_width = new_back_buffer_width;
	back_buffer_height = new_back_buffer_height;

	int bytes_per_pixel = 4;
	int back_buffer_size = bytes_per_pixel * back_buffer_width * back_buffer_height;
	back_buffer_mem = VirtualAlloc(0, back_buffer_size, MEM_COMMIT, PAGE_READWRITE);
}

void update_window_front_buffer(HWND window)
{
	HDC window_device_context = GetDC(window);
	RECT front_buffer_rect;
	GetClientRect(window, &front_buffer_rect);
	int front_buffer_width = front_buffer_rect.right - front_buffer_rect.left;
	int front_buffer_height = front_buffer_rect.bottom - front_buffer_rect.top;
	StretchDIBits
	(
		window_device_context, 
		0, 0, front_buffer_width, front_buffer_height, 
		0, 0, back_buffer_width, back_buffer_height, 
		back_buffer_mem, &back_buffer_info,
		DIB_RGB_COLORS, SRCCOPY
	);
	ReleaseDC(window, window_device_context);
}

LRESULT CALLBACK window_event_callback(HWND window, UINT message, WPARAM wparam, LPARAM lparam)
{
	LRESULT result = 0;
	switch(message)
	{
		case WM_SIZE:
		{
			RECT window_rect;
			GetClientRect(window, &window_rect);
			int window_width = window_rect.right - window_rect.left;
			int window_height = window_rect.bottom - window_rect.top;
			size_window_back_buffer(window_width, window_height);
			break;
		}
		case WM_CLOSE:
			running = false;
			break;
		default:
			result = DefWindowProc(window, message, wparam, lparam);
			break;
	}
	return result;
}

int WINAPI WinMain(HINSTANCE instance, HINSTANCE prev_instance, LPSTR cmd_line, int show_cmd_line)
{
	Vec2 v = {1.0, 2.0};
	printf("%f\n", dot(v, v));
	WNDCLASS window_class = {};
	window_class.style = CS_HREDRAW | CS_VREDRAW;
	window_class.lpfnWndProc = window_event_callback;
	window_class.lpszClassName = "RaytraceClass";

	if(!RegisterClass(&window_class))
	{
		printf("PROGRAM FAILED: Could not register window class\n");
		return 1;
	}

	HWND window = 	CreateWindowEx
			(
				0, window_class.lpszClassName, "Raytrace", WS_OVERLAPPEDWINDOW | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT,
				800, 600, 0, 0, prev_instance, NULL
			);
	if(!window)
	{
		printf("PROGRAM FAILED: Could not create window\n");
		return 2;
	}

	while(running)
	{
		MSG message;
		while(PeekMessage(&message, 0, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&message);
			DispatchMessage(&message);
		}

		//Render stuff
		render_that_good_shit_right_there();

		update_window_front_buffer(window);
	}

	return 0;
}
