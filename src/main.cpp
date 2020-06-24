#include <Windows.h>
#include <stdio.h>
#include <stdint.h>
#include "Maths.h"

//TODO: LONGTERM
//	- Timing/performance
//	- Implement printf/sprintf/etc.
//	- Implement maths functions (trig, pow, ln etc.)
//	- Remove CRT
//	- Direct MVSC build output files (that aren't exe) to build folder
//	- Gamma correction

//TODO: Solve rendering equation
//	- Recursive raytrace
//	- Physically(ish) based
//	- Output result to window

//TODO: Intersection tests
//	- If eye ray intersects sphere/plane, colour pixel black

//DOING:

struct RGB8
{
	uint8_t B;
	uint8_t G;
	uint8_t R;
};
typedef Vec3 RGB64;

RGB8 rgb64_to_rgb8(RGB64 rgb64)
{
	RGB8 rgb8 = {};

	if(rgb64.R < 0.0) rgb64.R = 0.0;
	if(rgb64.R > 1.0) rgb64.R = 1.0;
	if(rgb64.G < 0.0) rgb64.G = 0.0;
	if(rgb64.G > 1.0) rgb64.G = 1.0;
	if(rgb64.B < 0.0) rgb64.B = 0.0;
	if(rgb64.B > 1.0) rgb64.B = 1.0;

	rgb8.R = (uint8_t)(rgb64.R * 255.0);
	rgb8.G = (uint8_t)(rgb64.G * 255.0);
	rgb8.B = (uint8_t)(rgb64.B * 255.0);

	return rgb8;
}

struct Render_Buffer
{
	uint32_t* pixels;
	int width;
	int height;
};

enum Pixel_Channel
{
	P_B = 0,
	P_G,
	P_R
};

inline
void set_pixel_channel(uint32_t* pixel, Pixel_Channel channel, uint8_t value)
{
	*(((uint8_t*)pixel) + channel) = value;
}

inline
uint32_t* get_pixel(Render_Buffer* r_buffer, int x, int y)
{
	return r_buffer->pixels + y*r_buffer->width + x;
}


void set_render_buffer_pixel_colour(Render_Buffer* r_buffer, int x, int y, RGB8 colour)
{
	uint32_t* pixel = get_pixel(r_buffer, x, y);
	uint32_t* c = (uint32_t*)(&colour);
	*pixel = *c;
}

void clear_render_buffer(Render_Buffer* r_buffer, RGB8 colour)
{
	for(int y = 0; y < r_buffer->height; ++y)
	{
		for(int x = 0; x < r_buffer->width; ++x) set_render_buffer_pixel_colour(r_buffer, x, y, colour);
	}
}







BITMAPINFO __back_buffer_info__ = {};
Render_Buffer __window_back_buffer__ = {};

bool running = true;

/*
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
*/

void size_window_back_buffer(Render_Buffer* back_buffer, int new_back_buffer_width, int new_back_buffer_height)
{
	if(back_buffer->pixels) VirtualFree(back_buffer->pixels, 0, MEM_RELEASE);

	__back_buffer_info__.bmiHeader.biSize = sizeof(__back_buffer_info__.bmiHeader);
	__back_buffer_info__.bmiHeader.biWidth = new_back_buffer_width;
	__back_buffer_info__.bmiHeader.biHeight = -new_back_buffer_height;
	__back_buffer_info__.bmiHeader.biPlanes = 1;
	__back_buffer_info__.bmiHeader.biBitCount = 32;
	__back_buffer_info__.bmiHeader.biCompression = BI_RGB;

	back_buffer->width = new_back_buffer_width;
	back_buffer->height = new_back_buffer_height;

	int bytes_per_pixel = 4;
	int back_buffer_size = bytes_per_pixel * back_buffer->width * back_buffer->height;
	back_buffer->pixels = (uint32_t*)VirtualAlloc(0, back_buffer_size, MEM_COMMIT, PAGE_READWRITE);
}

void update_window_front_buffer(HWND window, Render_Buffer* back_buffer)
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
		0, 0, back_buffer->width, back_buffer->height, 
		back_buffer->pixels, &__back_buffer_info__,
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
			size_window_back_buffer(&__window_back_buffer__, window_width, window_height);
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

RGB64 cast_ray(Ray ray)
{
	return ray.direction;
}

void raytrace_scene(Render_Buffer* render_target, double fov, double near_plane)
{
	Vec3 eye = {0.0, 0.0, 1.0};
	Vec3 forward = {0.0, 0.0, -1.0};
	Vec3 right = {1.0, 0.0, 0.0};
	Vec3 up = {0.0, 1.0, 0.0};
	Vec3 pixel_center = {};

	int image_plane_width_px = render_target->width;
	int image_plane_height_px = render_target->height;
	double aspect_ratio = (double)(image_plane_width_px)/(double)(image_plane_height_px);
	double image_plane_width = 2.0 * near_plane * tan_deg(fov/2.0);
	double image_plane_height = image_plane_width / aspect_ratio;
	double pixel_width = image_plane_width/(double)(image_plane_width_px);
	double pixel_height = image_plane_height/(double)(image_plane_height_px);

	Vec3 image_plane_top_left = eye + near_plane * forward - 0.5 * image_plane_width * right + 0.5 * image_plane_height * up;
	for(int y = 0; y < image_plane_height_px; ++y)
	{
		for(int x = 0; x < image_plane_width_px; ++x)
		{
			Vec3 pixel_center = image_plane_top_left + ((double)x + 0.5)*pixel_width*right - ((double)y + 0.5)*pixel_height*up;
			Ray eye_ray = {};
			eye_ray.origin = eye;
			eye_ray.direction = normalise(pixel_center - eye_ray.origin);
			RGB64 raycast_result = cast_ray(eye_ray);
			set_render_buffer_pixel_colour(render_target, x, y, rgb64_to_rgb8(raycast_result));
		}
	}
}

int WINAPI WinMain(HINSTANCE instance, HINSTANCE prev_instance, LPSTR cmd_line, int show_cmd_line)
{
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

		RGB8 clear_colour = {105, 192, 255};
		clear_render_buffer(&__window_back_buffer__, clear_colour);

		raytrace_scene(&__window_back_buffer__, 90.0, 1.0);

		update_window_front_buffer(window, &__window_back_buffer__);
	}

	return 0;
}
