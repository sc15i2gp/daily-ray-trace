#include "daily_ray_trace.h"

//CODE STANDARDS:
//	- No templates (they kinda suck and are hard to use)
//	- Keep to C features as much as possible, except for:
//		- Operator overloading
//		- Constructors (for profiling)

//NOTES on bugs and problems:
//	- Incorrect colour in indirect lighting tint near cornell box walls
//		- I'm fairly sure it's a problem with the bsdf function
//		- It uses blinn-phong glossy, which means an indirect ray's contribution
//		- will only be significant when the incident and incoming rays are close
//		- to the surface normal
//		- (SOLVED) It was due to not reversing the outgoing direction in the indirect
//		- lighting contribution function. Fixed by reversing the ougoing direction in cast_ray()
//	- Program currently halts after a while due to windows (not responding)
//		- Can fix that with separate thread showing image progress
//		- (PARTIALLY SOLVED) Raytrace happens on separate thread, however it still halts a bit
//		- when copying raytrace output buffer to window back buffer
//	- Black pixels after adding cosine weighted hemisphere sampling
//		- (SOLVED) It's due to the cos_weighted_sample_hemisphere returning v such that dot(v, surface_normal) = 0
//		- which is used to compute the direction pdf value
//	- Using e illuminant on sphere at origin in cornell box, eye rays intersecting with light return black radiance (but other lighting is correct)
//		- (SOLVED) Algorithm wasn't taken into account emissive contributions
//	- Artifacts near specular direction of eye ray to light source
//		- Specifically sphere source
//		- Not sure yet what's causing it
//		- The artifacts are symmetrical about the center of the sphere along the plane
//		- (SOLVED) Was caused by floating point precision error
//	- Had difficulty debugging when looking for problems causing black pixels
//		- Turns out I'd forgotten that the lens model changed which flips the image output buffer

//NOTES: Parameters I've wanted to change/view but had to rebuild for:
//	- Emission/reflectance spectra
//	- Number of render samples to take
//	- Different sampling strategies
//		- Uniform vs cosine weighted hemisphere sampling
//	- Material parameters
//		- Gloss shininess
//	- Scene objects
//		- Generally: positions, spectra, light types etc.
//	- Image plane
//		- Position

//TODO: LONGTERM
//	- Dispersion
//	- Change reference white + find more precise spectra (maybe parameterise)
//	- Animation
//	- CSG

//TODO: Alternative methods
//	- Trowbridge-Reitz distribution for microfacets
//	- Oren-Nayar for diffuse materials
//	- Ward BSDF for microfacets
//	- Fourier BSDFs
//	- Investigate different SPD representations
//	- Different camera models/lenses (fisheye etc.)
//	- Different RGB -> SPD method: "Physically Meaningful Rendering using Tristimulus Colours"
//	- Consider using explicitly sized types
//	- Maybe use vec4 for everything (to make matrix operations slightly easier and 4 numbers easier to optimise than 3)
//	- Kirschoff's law for non-black body sources
//	- Better RNG
//	- Remove CSTDLIB
//		- Implement printf/sprintf/etc.
//		- Implement maths functions (trig, pow, ln etc.)
//		- Numeric limits (e.g. DBL_MAX)
//	- Bidirectional methods

//TODO: Necessary
//	- Textures
//		- Colour/image
//		- Normal map/bump map/displacement map
//		- Frosted glass roughness
//	- Move platform code to platform file
//	- Optimise
//		- Reduce amount of code
//		- Profiling
//		- Write faster code
//	- Quality of life
//		- Output time taken in seconds, minutes, hours etc.
//		- Resizing screen and how rendering should be handled/updated etc.
//		- Functionality for producing comparison images (e.g. sampling strategies)
//		- Output to jpg/bmp/png
//		- Try and get rid of reliance on plane normal direction (hard to use else)
//		- Invesitgate use of unity build
//		- Note failure points/error cases and handle
//		- Direct build output files (that aren't exe) to a particular location
//		- Add error handling to platform functions
//		- Add citations in code for origin of algorithm/technique/bsdf etc.
//	- Scene editing
//		- OpenGL for rendering scene preview and UI
//		- UI for camera control, raytrace control, parameterisation etc.
//		- Maybe screen shows progress of raytrace
//	- Subsurface scattering
//	- Allow objects within transmission media (will not work correctly as of yet)
//	- Generic triangulation rendering
//	- Skybox/infinite light/infinite geometry (such as infinite ground plane)
//	- Volumetric transport

//TODO: NOW
//	- Optimise memory stuff
//		- Reduce number of calls to memcpy and memset
//		- Improve cache performance
//	- Optimisation/Cleaning + ease of debugging
//		- Sort out floating point precision issues
//		- Remove superfluous code 
//		- Profiling
//		- Optimise slow methods
//		- Maybe multithreading
//	- Fun things to render
//		- Water in a box
//		- Frosted glass with earth texture for frostiness
//		- Torus

//Reduce time spent copying:
//	- Majority of time now spent copying material data
//	- Don't store any spectra on the stack
//	- Have pool of spectra which can be accessed

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
	void* ptr = VirtualAlloc(0, size, MEM_COMMIT, PAGE_READWRITE);
	ZeroMemory(ptr, size);
	return ptr;
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

double cycles_to_s(long unsigned int cycles)
{
	return (double)(cycles)/pc_frequency;
}

double cycles_to_ms(long unsigned int cycles)
{
	return cycles_to_s(cycles) * 1000.0;
}

long unsigned int elapsed_time_in_cycles(Timer* t)
{
	long unsigned int elapsed = (long unsigned int)(t->stop_time.QuadPart - t->start_time.QuadPart);
	return elapsed;
}

double elapsed_time_in_s(Timer* t)
{
	return cycles_to_s(elapsed_time_in_cycles(t));
}

double elapsed_time_in_ms(Timer* t)
{
	return cycles_to_ms(elapsed_time_in_cycles(t));
}

void output_to_ppm(const char* path, Texture r_buffer)
{
	char* contents = (char*)alloc(MEGABYTES(10));

	int contents_size = sprintf(contents, "P%d\n", 3);
	contents_size += sprintf(contents + contents_size, "%d %d\n", r_buffer.width, r_buffer.height);
	contents_size += sprintf(contents + contents_size, "255\n");
	for(int y = 0; y < r_buffer.height; ++y)
	{
		for(int x = 0; x < r_buffer.width; ++x)
		{
			RGB8 p = *TEXTURE_READ(RGB8, r_buffer, x, y);
			contents_size += sprintf(contents + contents_size, "%d %d %d\n", p.R, p.G, p.B);
		}
	}
	
	write_file_contents(path, contents, contents_size);
	dealloc(contents);
}

void invert_render_buffer(Texture r_buffer)
{
	for(int y = 0; y < r_buffer.height/2; ++y)
	{
		for(int x = 0; x < r_buffer.width; ++x)
		{
			uint32_t* p_0 = TEXTURE_READ(uint32_t, r_buffer, x,  y);
			uint32_t* p_1 = TEXTURE_READ(uint32_t, r_buffer, x, r_buffer.height - 1 - y);
			uint32_t p_temp = *p_0;
			*p_0 = *p_1;
			*p_1 = p_temp;
		}
	}

	for(int y = 0; y < r_buffer.height; ++y)
	{
		for(int x = 0; x < r_buffer.width/2; ++x)
		{
			uint32_t* p_0 = TEXTURE_READ(uint32_t, r_buffer, x, y);
			uint32_t* p_1 = TEXTURE_READ(uint32_t, r_buffer, r_buffer.width - 1 - x, y);
			uint32_t p_temp = *p_0;
			*p_0 = *p_1;
			*p_1 = p_temp;
		}
	}
}

Texture __window_back_buffer__ = {};

bool running = true;

LRESULT CALLBACK window_event_callback(HWND window, UINT message, WPARAM wparam, LPARAM lparam)
{
	LRESULT result = 0;
	switch(message)
	{
		case WM_CLOSE:
			running = false;
			break;
		default:
			result = DefWindowProc(window, message, wparam, lparam);
			break;
	}
	return result;
}

DWORD WINAPI render_image_task(LPVOID param)
{
	bool* completed_raytrace = (bool*)param;
	render_image(&__window_back_buffer__);
	*completed_raytrace = true;
	return 0;
}

#define RENDER_TARGET_WIDTH 800
#define RENDER_TARGET_HEIGHT 600

#define __USE_MINGW_ANSI_STDIO 1
int WINAPI WinMain(HINSTANCE instance, HINSTANCE prev_instance, LPSTR cmd_line, int show_cmd_line)
{
	query_pc_frequency();

	HWND window = {};
	
	WNDCLASS window_class = {};
	window_class.style = CS_HREDRAW | CS_VREDRAW;
	window_class.lpfnWndProc = window_event_callback;
	window_class.lpszClassName = "RaytraceClass";

	if(!RegisterClass(&window_class))
	{
		printf("PROGRAM FAILED: Could not register window class\n");
		return 1;
	}

	window = 	CreateWindowEx
			(
				0, window_class.lpszClassName, "Raytrace", WS_OVERLAPPEDWINDOW | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT,
				RENDER_TARGET_WIDTH, RENDER_TARGET_HEIGHT, 0, 0, prev_instance, NULL
			);
	if(!window)
	{
		printf("PROGRAM FAILED: Could not create window\n");
		return 2;
	}

	RGB8 clear_colour = {};
	TEXTURE_CLEAR(__window_back_buffer__, clear_colour);

	BITMAPINFO __back_buffer_info__ = {};
	__back_buffer_info__.bmiHeader.biSize = sizeof(__back_buffer_info__.bmiHeader);
	__back_buffer_info__.bmiHeader.biWidth = window_width(window);
	__back_buffer_info__.bmiHeader.biHeight = -window_height(window);
	__back_buffer_info__.bmiHeader.biPlanes = 1;
	__back_buffer_info__.bmiHeader.biBitCount = 32;
	__back_buffer_info__.bmiHeader.biCompression = BI_RGB;

	__window_back_buffer__.width = window_width(window);
	__window_back_buffer__.height = window_height(window);

	int bytes_per_pixel = 4;
	int back_buffer_size = bytes_per_pixel * __window_back_buffer__.width * __window_back_buffer__.height;
	__window_back_buffer__.pixels = (uint8_t*)alloc(back_buffer_size);
	__window_back_buffer__.pixel_size = sizeof(uint32_t);

	init_profiling();

	bool completed_raytrace = false;
	HANDLE raytrace_thread = CreateThread(NULL, 0, render_image_task, &completed_raytrace, 0, NULL);

	while(running && !completed_raytrace)
	{
		MSG message;
		while(PeekMessage(&message, 0, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&message);
			DispatchMessage(&message);
		}

		HDC window_device_context = GetDC(window);
		StretchDIBits
		(
			window_device_context, 
			0, 0, window_width(window), window_height(window), 
			__window_back_buffer__.width, __window_back_buffer__.height, -__window_back_buffer__.width, -__window_back_buffer__.height, 
			__window_back_buffer__.pixels, &__back_buffer_info__,
			DIB_RGB_COLORS, SRCCOPY
		);
		ReleaseDC(window, window_device_context);
	}
	
	TerminateThread(raytrace_thread, 0);
	
	print_render_profile();
	print_profile();

	printf("Writing back buffer to file\n");
	invert_render_buffer(__window_back_buffer__);
	output_to_ppm("output.ppm", __window_back_buffer__);

	CloseHandle(raytrace_thread);
	return 0;
}
