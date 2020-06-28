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
//	- Invesitgate use of unity build

//TODO: NOW
//	- SPD of D65 illuminant
//	- Single sphere and single point light interaction
//		- Choose emitted colour of point light
//			- Generate light spd
//			- SPD -> RGB
//			- Show emitted colour of light at eye ray intersection
//		- Choose material colour
//			- Generate/find material spd
//				- Find spd source on internet
//				- RGB -> SPD
//			- Show material spectrum at eye ray intersection
//		- Multiply spds of light and material at eye ray intersection
//	- Kirschoff's law for non-black body sources

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

#define SPECTRUM_SAMPLE_MAX 1024

struct Spectrum
{
	long double samples[SPECTRUM_SAMPLE_MAX];
	long double start_wavelength; //In nm
	long double end_wavelength; //In nm
	int number_of_samples;
};

//NOTE: Currently assumes both spectra have same wavelength ranges and intervals
Spectrum operator*(Spectrum spd_0, Spectrum spd_1)
{
	Spectrum spd = {};
	spd.start_wavelength = spd_0.start_wavelength;
	spd.end_wavelength = spd_0.end_wavelength;
	spd.number_of_samples = spd_0.number_of_samples;

	for(int i = 0; i < spd.number_of_samples; ++i)
	{
		spd.samples[i] = spd_0.samples[i] * spd_1.samples[i];
	}

	return spd;
}

void normalise(Spectrum& spd)
{
	long double highest_value = 0.0L;
	for(int i = 0; i < spd.number_of_samples; ++i) if(spd.samples[i] > highest_value) highest_value = spd.samples[i];

	for(int i = 0; i < spd.number_of_samples; ++i) spd.samples[i] /= highest_value;
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

inline
long double m_to_nm(long double m)
{
	long double nm = m * 1e9;
	return nm;
}

inline
long double nm_to_m(long double nm)
{
	long double m = nm * 1e-9;
	return m;
}

long double sample_interval(Spectrum spd)
{
	return (spd.end_wavelength - spd.start_wavelength)/(long double)(spd.number_of_samples - 1);
}

Spectrum generate_black_body_spd(long double temperature, long double start_wavelength = 300.0, long double end_wavelength = 780.0, long double sample_interval = 5.0)
{
	Spectrum spd = {};
	spd.number_of_samples = (int)((end_wavelength - start_wavelength)/sample_interval) + 1;
	for(int i = 0; i < spd.number_of_samples; ++i)
	{
		long double wavelength = nm_to_m(start_wavelength + (long double)(i) * sample_interval);
		spd.samples[i] = compute_black_body_power(temperature, wavelength);
	}
	spd.start_wavelength = start_wavelength;
	spd.end_wavelength = end_wavelength;
	return spd;
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

Spectrum load_spd(const char* spd_path)
{
	Spectrum spd = {};
	HANDLE spd_file = CreateFile(spd_path, GENERIC_READ, FILE_SHARE_WRITE, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	DWORD spd_size = GetFileSize(spd_file, NULL);
	char* spd_contents = (char*)alloc(spd_size);
	DWORD bytes_read = 0;
	ReadFile(spd_file, spd_contents, spd_size, &bytes_read, NULL);
	
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
	CloseHandle(spd_file);

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

RGB64 spectrum_to_RGB64(Spectrum spd)
{
	Vec3 xyz = spectrum_to_xyz(spd);
	Vec3 white_xyz = spectrum_to_xyz(reference_white);

	xyz /= white_xyz.y;

	printf("%f %f %f\n", xyz.R, xyz.G, xyz.B);
	//XYZ -> RGB
	RGB64 linear_rgb = {};
	linear_rgb.R = 3.24097 * xyz.x - 1.53738 * xyz.y - 0.49831 * xyz.z;
	linear_rgb.G = -0.96934 * xyz.x + 1.87596 * xyz.y + 0.04156 * xyz.z;
	linear_rgb.B = 0.05563 * xyz.x - 0.20398 * xyz.y + 1.05697 * xyz.z;

	RGB64 rgb = gamma_linear_rgb(linear_rgb);
	return rgb;
}

RGB64 cast_ray(Ray ray, Spectrum light_spd)
{
	RGB64 black = {0.0, 0.0, 0.0};
	double t = 0.0;
	Vec3 p = Vec3{-0.2, 0.5, -0.2};
	Vec3 u = 0.5*normalise(Vec3{1.0, 0.0, 0.0});
	Vec3 v = normalise(Vec3{0.0, -0.6, 0.5});
	Vec3 n = normalise(cross(u, v));
	if(ray_intersects_sphere(ray, Vec3{}, 0.3f, &t))
	{
		return spectrum_to_RGB64(light_spd);
	}
	else
	{
		return black;
	}
}

void raytrace_scene(Render_Buffer* render_target, double fov, double near_plane, Spectrum light_spd)
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
			RGB64 raycast_result = cast_ray(eye_ray, light_spd);
			set_render_buffer_pixel_colour(render_target, x, y, rgb64_to_rgb8(raycast_result));
		}
	}
}

#define BYTES(n) n
#define KILOBYTES(n) 1024 * BYTES(n)
#define MEGABYTES(n) 1024 * KILOBYTES(n)
#define GIGABYTES(n) 1024 * MEGABYTES(n)

//start_wavelength is the wavelength for the first value in spd
//wavelength_sample_size is the range of wavelengths that a single sample in spd covers
//To find the wavelength of spd[i], start_wavelength + i * wavelength_sample_size
//Prints values to csv file
void output_black_body_curves_to_csv(Spectrum* spds, long double* temperatures)
{
	long double wavelength_sample_size = (spds[0].end_wavelength - spds[0].start_wavelength)/(long double)(spds[0].number_of_samples - 1);
	//Allocate resources to open file
	HANDLE csv_file = CreateFile("black_body.csv", GENERIC_WRITE, FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	char* csv_contents = (char*)alloc(MEGABYTES(64));

	//Convert spd to csv_contents
	//First, write the header row to contents
	int contents_size = sprintf(csv_contents, "%s, %dK, %dK, %dK, %dK\n", "Wavelength", (int)temperatures[0], (int)temperatures[1], (int)temperatures[2], (int)temperatures[3]);
	//Next write spectrum entries
	for(int i = 0; i < spds[0].number_of_samples; ++i)
	{
		long double wavelength = spds[0].start_wavelength + (long double)(i) * wavelength_sample_size;
		long double spd_samples[4] = {};
		for(int j = 0; j < 4; ++j) spd_samples[j] = spds[j].samples[i];
		contents_size += sprintf(csv_contents + contents_size, "%Lg, %Lg, %Lg, %Lg, %Lg\n", wavelength, spd_samples[0], spd_samples[1], spd_samples[2], spd_samples[3]);
	}
	
	//Write contents to csv_file
	DWORD bytes_written = 0;
	WriteFile(csv_file, csv_contents, contents_size+1, &bytes_written, NULL);
	
	//Deallocate resources
	CloseHandle(csv_file);
	dealloc(csv_contents);
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
	HANDLE csv_file = CreateFile(path, GENERIC_WRITE, FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	DWORD bytes_written = 0;
	WriteFile(csv_file, csv_contents, contents_size+1, &bytes_written, NULL);
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

	//long double l = 4000.0L;
	long double l = 10000.0L;
	Spectrum light_spd = generate_black_body_spd(l, 380.0L, 780.0L);
	normalise(light_spd);
	output_spd_to_csv(light_spd);
	load_colour_matching_functions();
	load_d65_illuminant();
	output_spd_to_csv(colour_matching_functions[0], "cmf_x_1.csv");
	output_spd_to_csv(colour_matching_functions[1], "cmf_y_1.csv");
	output_spd_to_csv(colour_matching_functions[2], "cmf_z_1.csv");
	output_spd_to_csv(reference_white, "white.csv");
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

		raytrace_scene(&__window_back_buffer__, 90.0, 0.1, light_spd);

		update_window_front_buffer(window, &__window_back_buffer__);
	}

	return 0;
}
