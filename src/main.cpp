#include <stdio.h>
#include <stdint.h>
#include "Platform.h"
#include "Maths.h"
#include "Colour.h"

//TODO: LONGTERM
//	- Output to jpg/png
//	- Investigate standardising spectra wavelength ranges in the program
//	- Investigate different SPD representations
//	- Kirschoff's law for non-black body sources
//	- Different RGB -> SPD method: "Physically Meaningful Rendering using Tristimulus Colours"
//	- Texturing (images + surface properties)
//	- Different camera models/lenses (fisheye etc.)
//	- Direct MVSC build output files (that aren't exe) to build folder
//	- Invesitgate use of unity build
//	- Remove CRT
//		- Implement printf/sprintf/etc.
//		- Implement maths functions (trig, pow, ln etc.)
//	- Add error handling to platform functions

//TODO: NOW
//	- Recursive raytrace
//		- Define scene (Cornell box with sphere and point light)
//		- Direct lighting
//		- Indirect lighting
//		- Monte carlo integration
//			- Importance sampling?
//		- Area lighting
//	- More geometry models
//	- Output raytraced scene to file
//		- Output raytraced image as bmp file
//		- OpenGL for rendering scene preview and UI
//		- UI for camera control, raytrace beginning
//		- Maybe screen shows progress of raytrace
//	- Better Sampling
//		- Integration importance sampling
//		- Image plane super sampling
//	- Profiling
//	- Other reflectance models
//		- Fresnel equation specular reflection
//	- Volumetric transport

//SPD -> RGB:
//	- First, SPD -> XYZ
//	- XYZ is a colour space which describes how the eye reacts to spectra.
//	- Each of its coordinates is the integral of a spectrum transformed by colour matching functions (aka standard observer functions).
//	- Colour matching functions describe how the eye reacts to different wavelengths within ~300 - ~800 nm
//	- The 3 (x, y, z) were found experimentally at different points in time through different methods (as well as deriving new curves by averaging old ones)
//	- Multiplying a spectrum by one of these 3 results in a function describing the eye's response to light.
//	- The integrals of these curves are represented by XYZ.
//	- The data set this project will initially use for the colour functions is CIE 1931 2 degree observer.
//	- The XYZ colour space has no meaning without a reference white point
//	- The white point is the whitest point in the colour space. The colour space can be normalised by division by the white point

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

char* find_next_line(char* c)
{
	for(; *c != '\n' && *c != 0; ++c);
	++c;
	return c;
}

char* find_next_character(char* c, char character)
{
	for(; *c != character && *c != 0; ++c);
	return c;
}

bool is_number_character(char c)
{
	return (c >= '0' && c <= '9') || c == '-';
}

char* find_next_number(char* c)
{
	for(; !is_number_character(*c) && *c != 0; ++c);
	return c;
}

Spectrum reference_white = {};
Spectrum colour_matching_functions[3] = {};

Spectrum white_rgb_to_spd = {};
Spectrum cyan_rgb_to_spd = {};
Spectrum magenta_rgb_to_spd = {};
Spectrum yellow_rgb_to_spd = {};
Spectrum red_rgb_to_spd = {};
Spectrum green_rgb_to_spd = {};
Spectrum blue_rgb_to_spd = {};

Spectrum load_spd(const char* spd_path)
{
	Spectrum spd = {};
	char* spd_contents = read_file_contents(spd_path);
	
	//Count number of samples
	int number_of_samples = 0;
	for(char* c = spd_contents; *c != 0; c = find_next_line(c))
	{
		//Increment number of samples
		++number_of_samples;
	}

	//printf("Number of samples = %d\n", number_of_samples);
	spd.number_of_samples = number_of_samples;

	char* c = spd_contents;
	for(int i = 0; i < number_of_samples; ++i)
	{
		if(i == 0)
		{
			spd.start_wavelength = (long double) atof(c);
		}
		else if(i == number_of_samples - 1)
		{
			spd.end_wavelength = (long double) atof(c);
		}

		c = find_next_number(find_next_character(c, ','));
		spd.samples[i] = (long double) atof(c);

		c = find_next_line(c);
	}

	dealloc(spd_contents);

	return spd;
}

void load_colour_matching_functions()
{
	colour_matching_functions[0] = load_spd("cmf_x.csv");
	colour_matching_functions[1] = load_spd("cmf_y.csv");
	colour_matching_functions[2] = load_spd("cmf_z.csv");
}

void load_d65_illuminant()
{
	reference_white = load_spd("d65.csv");
	//normalise(reference_white);
	for(int i = 0; i < reference_white.number_of_samples; ++i) reference_white.samples[i] /= 100.0L;
}

void load_e_illuminant()
{
	long double wavelength_interval = 5.0L;
	reference_white.start_wavelength = 380.0L;
	reference_white.end_wavelength = 720.0L;
	reference_white.number_of_samples = ((reference_white.end_wavelength-reference_white.start_wavelength)/ wavelength_interval) + 1;
	for(int i = 0; i < reference_white.number_of_samples; ++i) reference_white.samples[i] = 1.0L;
}

void load_rgb_to_spd_functions()
{
	white_rgb_to_spd = load_spd("white_rgb_to_spd.csv");
	cyan_rgb_to_spd = load_spd("cyan_rgb_to_spd.csv");
	magenta_rgb_to_spd = load_spd("magenta_rgb_to_spd.csv");
	yellow_rgb_to_spd = load_spd("yellow_rgb_to_spd.csv");
	red_rgb_to_spd = load_spd("red_rgb_to_spd.csv");
	green_rgb_to_spd = load_spd("green_rgb_to_spd.csv");
	blue_rgb_to_spd = load_spd("blue_rgb_to_spd.csv");
}

Vec3 spectrum_to_xyz(Spectrum spd)
{
	//SPD -> XYZ
	Spectrum spd_xyz[3] = {};
	spd_xyz[0] = spd * colour_matching_functions[0];
	spd_xyz[1] = spd * colour_matching_functions[1];
	spd_xyz[2] = spd * colour_matching_functions[2];

	Vec3 xyz = {};
	long double l = (spd.end_wavelength - spd.start_wavelength)/(long double)(spd.number_of_samples);
	for(int i = 0; i < 3; ++i)
	{
		for(int j = 0; j < spd.number_of_samples; ++j)
		{
			xyz[i] += (double)spd_xyz[i].samples[j];
		}
	}
	xyz *= l;
	
	return xyz;
}

RGB64 gamma_linear_rgb(RGB64 linear_rgb)
{
	RGB64 rgb = {};
	for(int i = 0; i < 3; ++i)
	{
		rgb[i] = (linear_rgb[i] <= 0.0031308) ? 12.92 * linear_rgb[i] : (1.055 * pow(linear_rgb[i], 1.0/2.4)) - 0.055;
	}
	return rgb;
}

RGB64 inverse_gamma_rgb(RGB64 rgb)
{
	RGB64 linear_rgb = {};
	for(int i = 0; i < 3; ++i)
	{
		linear_rgb[i] = (rgb[i] <= 0.04045) ? rgb[i]/12.92 : pow(((rgb[i] + 0.055)/1.055), 2.4);
	}
	return linear_rgb;
}

RGB64 spectrum_to_RGB64(Spectrum spd)
{
	Vec3 xyz = spectrum_to_xyz(spd);
	Vec3 white_xyz = spectrum_to_xyz(reference_white);

	xyz /= white_xyz.y;

	//XYZ -> RGB
	RGB64 linear_rgb = {};
	linear_rgb.R = 3.24097 * xyz.x - 1.53738 * xyz.y - 0.49831 * xyz.z;
	linear_rgb.G = -0.96934 * xyz.x + 1.87596 * xyz.y + 0.04156 * xyz.z;
	linear_rgb.B = 0.05563 * xyz.x - 0.20398 * xyz.y + 1.05697 * xyz.z;

	RGB64 rgb = gamma_linear_rgb(linear_rgb);
	return rgb;
}

Spectrum RGB64_to_spectrum(RGB64 rgb)
{
	Spectrum spd = {};
	spd.start_wavelength = white_rgb_to_spd.start_wavelength;
	spd.end_wavelength = white_rgb_to_spd.end_wavelength;
	spd.number_of_samples = white_rgb_to_spd.number_of_samples;
	if(rgb.R <= rgb.G && rgb.R <= rgb.B)
	{
		spd += rgb.R * white_rgb_to_spd;
		if(rgb.G <= rgb.B)
		{
			spd += (rgb.G - rgb.R)*cyan_rgb_to_spd;
			spd += (rgb.B - rgb.G)*blue_rgb_to_spd;
		}
		else
		{
			spd += (rgb.B - rgb.R)*cyan_rgb_to_spd;
			spd += (rgb.G - rgb.B)*green_rgb_to_spd;
		}
	}
	else if(rgb.G <= rgb.R && rgb.G <= rgb.B)
	{
		spd += rgb.G * white_rgb_to_spd;
		if(rgb.R <= rgb.B)
		{
			spd += (rgb.R - rgb.G)*magenta_rgb_to_spd;
			spd += (rgb.B - rgb.R)*blue_rgb_to_spd;
		}
		else
		{
			spd += (rgb.B - rgb.G)*magenta_rgb_to_spd;
			spd += (rgb.R - rgb.B)*red_rgb_to_spd;
		}
	}
	else
	{
		spd += rgb.B * white_rgb_to_spd;
		if(rgb.R <= rgb.G)
		{
			spd += (rgb.R - rgb.B)*yellow_rgb_to_spd;
			spd += (rgb.G - rgb.R)*green_rgb_to_spd;
		}
		else
		{
			spd += (rgb.G - rgb.B)*yellow_rgb_to_spd;
			spd += (rgb.R - rgb.G)*red_rgb_to_spd;
		}
	}

	return spd;
}

struct Surface_Point
{
	Spectrum lambertian_spd;
	Spectrum specular_spd;
	Vec3 normal;
	Vec3 position;
};


Spectrum lambertian_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	return p.lambertian_spd / PI;
}

double d_max(double a, double b)
{
	return (a < b) ? b : a;
}

Spectrum specular_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	Vec3 bisector = (incoming + outgoing)/2.0;
	double specular_coefficient = pow(d_max(0.0, dot(p.normal, bisector)), 80.0);
	return specular_coefficient * p.specular_spd;
}

RGB64 cast_ray(Ray ray, Spectrum light_spd, Spectrum sphere_spd, Spectrum specular_spd)
{
	RGB64 black = {0.0, 0.0, 0.0};
	double t = 0.0;
	Vec3 p = Vec3{-0.2, 0.5, -0.2};
	Vec3 u = 0.5*normalise(Vec3{1.0, 0.0, 0.0});
	Vec3 v = normalise(Vec3{0.0, -0.6, 0.5});
	Vec3 n = normalise(cross(u, v));

	Vec3 light_pos = {0.0, 1.0, 1.0};
	Vec3 sphere_pos = {};
	if(ray_intersects_sphere(ray, sphere_pos, 0.3f, &t))
	{
		Vec3 intersection = ray.origin + t * ray.direction;
		Vec3 incoming = light_pos - intersection;
		Vec3 outgoing = -ray.direction;
		Surface_Point point = {};
		point.lambertian_spd = sphere_spd;
		point.specular_spd = specular_spd;
		point.position = intersection;
		point.normal = normalise(intersection - sphere_pos);
		Spectrum result = dot(incoming, point.normal) * light_spd * (lambertian_bsdf(point, incoming, outgoing) + specular_bsdf(point, incoming, outgoing));
		return spectrum_to_RGB64(result);
	}
	else
	{
		return black;
	}
}

void raytrace_scene(Render_Buffer* render_target, double fov, double near_plane, Spectrum light_spd, Spectrum sphere_spd, Spectrum specular_spd)
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
			RGB64 raycast_result = cast_ray(eye_ray, light_spd, sphere_spd, specular_spd);
			set_render_buffer_pixel_colour(render_target, x, y, rgb64_to_rgb8(raycast_result));
		}
	}
}

void output_spd_to_csv(Spectrum spd, const char* path = "black_body.csv")
{
	char* csv_contents = (char*)alloc(MEGABYTES(64));

	int contents_size = sprintf(csv_contents, "%s, %s\n", "Wavelength", "Value");
	long double wavelength_sample_size = (spd.end_wavelength - spd.start_wavelength)/(long double)(spd.number_of_samples - 1);
	for(int i = 0; i < spd.number_of_samples; ++i)
	{
		long double wavelength = spd.start_wavelength + (long double)(i) * wavelength_sample_size;
		contents_size += sprintf(csv_contents + contents_size, "%Lg, %Lg\n", wavelength, spd.samples[i]);
	}

	//Write contents to csv_file
	write_file_contents(path, csv_contents, contents_size+1);
	dealloc(csv_contents);
}

#define __USE_MINGW_ANSI_STDIO 1
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

	long double l = 4000.0L;
	//long double l = 10000.0L;
	Spectrum light_spd = generate_black_body_spd(l, 380.0L, 720.0L);
	normalise(light_spd);
	load_colour_matching_functions();
	//load_d65_illuminant();
	load_e_illuminant();
	load_rgb_to_spd_functions();

	Spectrum mat_spd = RGB64_to_spectrum(RGB64{0.25, 0.8, 0.4});
	Spectrum specular_spd = load_spd("d65.csv");
	while(running)
	{
		MSG message;
		while(PeekMessage(&message, 0, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&message);
			DispatchMessage(&message);
		}

		RGB8 clear_colour = {};
		clear_render_buffer(&__window_back_buffer__, clear_colour);

		raytrace_scene(&__window_back_buffer__, 90.0, 0.1, light_spd, mat_spd, specular_spd);

		update_window_front_buffer(window, &__window_back_buffer__);
	}

	return 0;
}
