#include <Windows.h>
#include <stdio.h>
#include <stdint.h>
#include "Maths.h"

//TODO: LONGTERM
//	- Recursive raytrace
//	- Timing/performance
//	- Implement printf/sprintf/etc.
//	- Implement maths functions (trig, pow, ln etc.)
//	- Remove CRT
//	- Direct MVSC build output files (that aren't exe) to build folder
//	- Gamma correction
//	- Volumetric effects
//	- Different camera models/lenses (fisheye etc.)
//	- UI stuff

//TODO: NOW
//	- Spectral power distribution (SPD) representation of light
//	- SPD of D65 illuminant
//	- SPD -> CIE XYZ -> RGB
//	- RGB -> SPD


void* alloc(int size)
{
	return VirtualAlloc(0, size, MEM_COMMIT, PAGE_READWRITE);
}

void dealloc(void* ptr)
{
	VirtualFree(ptr, 0, MEM_RELEASE);
}

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

void size_window_back_buffer(Render_Buffer* back_buffer, int new_back_buffer_width, int new_back_buffer_height)
{
	if(back_buffer->pixels) dealloc(back_buffer->pixels);

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
	back_buffer->pixels = (uint32_t*)alloc(back_buffer_size);
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
	RGB64 red = {1.0, 0.0, 0.0};
	RGB64 black = {};
	double t = 0.0;
	Vec3 p = Vec3{-0.2, 0.5, -0.2};
	Vec3 u = 0.5*normalise(Vec3{1.0, 0.0, 0.0});
	Vec3 v = normalise(Vec3{0.0, -0.6, 0.5});
	Vec3 n = normalise(cross(u, v));
	if(ray_intersects_plane(ray, p, n, u, v, &t))
	{
		return black;
	}
	else
	{
		return red;
	}
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

long double c = 2.99792458e8L; //Speed of light
long double h = 6.626176e-34L; //Planck constant
long double k = 1.380662e-23L; //Boltzmann constant
//Uses Planck's formula to compute the power of a black body radiator's emission at a given temperature and wavelength
//Temperature in kelvin and wavelength in meters
long double compute_black_body_power(long double temperature, long double wavelength)
{
	long double numerator = 2.0L * PI * h * c * c;

	long double lambda_5 = powl(wavelength, 5.0L);
	long double e_power_numerator = (h * c) / k;
	long double e_power_denominator = temperature * wavelength;
	long double e_power = e_power_numerator / e_power_denominator;
	long double e_term = expl(e_power);

	long double denominator = lambda_5 * (e_term - 1.0L);

	return numerator/denominator;
}

#define BYTES(n) n
#define KILOBYTES(n) 1024 * BYTES(n)
#define MEGABYTES(n) 1024 * KILOBYTES(n)
#define GIGABYTES(n) 1024 * MEGABYTES(n)

//start_wavelength is the wavelength for the first value in spd
//wavelength_sample_size is the range of wavelengths that a single sample in spd covers
//To find the wavelength of spd[i], start_wavelength + i * wavelength_sample_size
//Prints values to csv file
void output_black_body_curve_to_csv(long double* spds, long double* temperatures, int number_of_spds, long double start_wavelength, long double wavelength_sample_size, int number_of_samples)
{
	//Allocate resources to open file
	HANDLE csv_file = CreateFile("black_body.csv", GENERIC_WRITE, FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	char* csv_contents = (char*)alloc(MEGABYTES(64));

	//Convert spd to csv_contents
	//First, write the header row to contents
	int contents_size = sprintf(csv_contents, "%s, %dK, %dK, %dK, %dK\n", "Wavelength", (int)temperatures[0], (int)temperatures[1], (int)temperatures[2], (int)temperatures[3]);
	//Next write spectrum entries
	for(int i = 0; i < number_of_samples; ++i)
	{
		long double wavelength = start_wavelength + (long double)(i) * wavelength_sample_size;
		long double spd_samples[4] = {};
		for(int j = 0; j < 4; ++j) spd_samples[j] = (spds + j*number_of_samples)[i];
		contents_size += sprintf(csv_contents + contents_size, "%Lg, %Lg, %Lg, %Lg, %Lg\n", wavelength, spd_samples[0], spd_samples[1], spd_samples[2], spd_samples[3]);
	}
	
	//Write contents to csv_file
	DWORD bytes_written = 0;
	WriteFile(csv_file, csv_contents, contents_size+1, &bytes_written, NULL);
	
	//Deallocate resources
	CloseHandle(csv_file);
	dealloc(csv_contents);
}

#define __USE_MINGW_ANSI_STDIO 1
int WINAPI WinMain(HINSTANCE instance, HINSTANCE prev_instance, LPSTR cmd_line, int show_cmd_line)
{
	WNDCLASS window_class = {};
	window_class.style = CS_HREDRAW | CS_VREDRAW;
	window_class.lpfnWndProc = window_event_callback;
	window_class.lpszClassName = "RaytraceClass";


	//long double temperatures[4] = {4500.0L, 4000.0L, 3500.0L, 3000.0L};
	long double temperatures[4] = {4000.0L, 4000.0L, 4000.0L, 2856.0L};
	int number_of_spds = 4;
	long double start_wavelength = 300.0e-9L;
	long double end_wavelength = 780.0e-9L;
	long double wavelength_sample_size = 5.0e-9L;
	int number_of_spd_samples = (int)((end_wavelength - start_wavelength)/wavelength_sample_size) + 1;
	int spd_mem_size = number_of_spd_samples * sizeof(long double);
	//Generate spd for black body at 6500K from 380nm to 700nm at 5nm intervals
	long double* spds = (long double*)alloc(spd_mem_size * number_of_spds);
	for(int i = 0; i < number_of_spds; ++i)
	{
		long double* spd = spds + i * number_of_spd_samples;
		for(int j = 0; j < number_of_spd_samples; ++j)
		{
			long double wavelength = start_wavelength + (long double)(j) * wavelength_sample_size;
			spd[j] = compute_black_body_power(temperatures[i], wavelength);
		}
	}
	output_black_body_curve_to_csv(spds, temperatures, number_of_spds, start_wavelength, wavelength_sample_size, number_of_spd_samples);
	dealloc(spds);

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

		raytrace_scene(&__window_back_buffer__, 90.0, 0.1);

		update_window_front_buffer(window, &__window_back_buffer__);
	}

	return 0;
}
