#include "daily_ray_trace.h"

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
//	- Going multithreaded caused incorrect image
//		- (POSSIBLE CAUSE) RNG sequences always the same and rand is not thread safe

//TODO:
//	- Features
//		- Objects inside other objects (transmission media)
//		- Volumetric transport
//		- Not fixed depth, biased statistical method
//		- Scene file, not hard coded
//			- Geometry
//			- Material BDSFs
//			- Spectra
//		- Bump maps
//		- Subsurface scattering
//		- Skybox
//		- More (maybe improved) geometry
//			- Torus
//			- CSG
//		- Support for 3d models
//		- Alternate spd representations
//			- RGB vs spectrum
//			- Basis (or other) functions
//		- Multiple platforms
//			- Windows
//			- Linux
// 	- Testing and quality
//		- Regressions
//		- Prevent bugs
//		- Edge cases
//			- Make sure camera works from all orientations
//		- Comparison and analysis (such as variance, convergence etc.)
//		- Performance
//			- Code performance
//			- Algorithm performance (reduce variance)
//		- Floating point stuff
//			- Precision
//			- Double vs float
//		- Program robustness
//		- Clean up project directory (maybe consider dev + release sort of structure)
//		- Improve texture access interface
//		- Improve RGB -> SPD (and vice versa) methods
//			- "Physically Meaningful Rendering using Tristimulus Colours"
//		- Robust RNG
//		- Betture texture sampling
//		- Better sampling of image plane (subdivision strategies, filter function etc.)
//		- Parameterise reference spectra
//		- Try and get rid of reliance on plane normal direction (hard to use else)
//		- Alt methods:
//			- Trowbridge-Reitz distribution for microfacets
//			- Oren-Nayar for diffuse materials
//			- Ward BSDF for microfacets
//			- Fourier BSDFs
//			- Investigate different SPD representations
//			- Different camera models/lenses (fisheye etc.)
//			- Bidirectional methods
//		- Improve program output
//		- Parameterise number of render samples to take
//		- Paper and code citations
//		- Multiple platforms
//		- Remove cstdlib
//	- Optimise
//		- Better and more disciplined memory layout
//		- Reduce code size
//		- Reduce memory footprint
//		- Do similar things together (e.g. spectral calculations, geometric calculations)
//		- Improve integration algorithm (reduce variance, biased statistical method to reduce depth of path in scene)
//	- Fun things to render
//		- Water in a box
//		- Frosted glass with earth texture for frostiness
//		- Torus
//		- Cylinder with two shadows: One showing a circle, one a square
//		- Light bulb with filament
//		- Spotlight

//DOING:
// - Get program working as baseline for further development
//		- Get profiling working
//		- Get debugging working
//		- Remove multithreading for now
//		- Figure out why debug build sample numbers aren't increasing
// - Clean up project directory structure
// - Remove/alter hard coded parameters
//		- Scene
//		- Number of samples to take
//		- Number of samples in spectra
// - Bone up on maths
//		- Light transport
//		- Monte Carlo integration
//		- Signal processing
// - Remove fixed depth casting
// - Transmission media
//		- Nested solids
//		- Volumetric transport
//		- Subsurface scattering
// - Reduce memory footprint
//		- Use file streaming and small amount of mem in cache
// - Sort out float precision issue

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

unsigned int get_pc_time()
{
	LARGE_INTEGER time;
	QueryPerformanceCounter(&time);
	unsigned int pc_time = (unsigned int)time.QuadPart;
	return pc_time;
}

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

void output_to_bmp(const char* path, Texture r_buffer)
{
	BITMAPFILEHEADER bmp_file_header = {};
	BITMAPINFOHEADER bmp_info_header = {};

	int texture_size = r_buffer.width * r_buffer.height * r_buffer.pixel_size;
	bmp_file_header.bfType = 0x4d42;
	bmp_file_header.bfSize = sizeof(bmp_file_header) + sizeof(bmp_info_header) + texture_size;
	bmp_file_header.bfOffBits = bmp_file_header.bfSize - texture_size;

	bmp_info_header.biSize = sizeof(bmp_info_header);
	bmp_info_header.biWidth = r_buffer.width;
	bmp_info_header.biHeight = r_buffer.height; //NOTE: THIS MAY NEED TO CHANGE IF PROBLEMS PRINTING FILE
	bmp_info_header.biPlanes = 1;
	bmp_info_header.biBitCount = 32;
	bmp_info_header.biCompression = BI_RGB;

	char* contents = (char*)alloc(bmp_file_header.bfSize);
	char* current_contents = contents;
	memcpy(current_contents, &bmp_file_header, sizeof(bmp_file_header));
	current_contents += sizeof(bmp_file_header);
	memcpy(current_contents, &bmp_info_header, sizeof(bmp_info_header));
	current_contents += sizeof(bmp_info_header);
	memcpy(current_contents, r_buffer.pixels, texture_size);

	write_file_contents(path, contents, bmp_file_header.bfSize);
	dealloc(contents);
}

/*
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
*/


#define RENDER_TARGET_WIDTH 800
#define RENDER_TARGET_HEIGHT 600

#define __USE_MINGW_ANSI_STDIO 1
int WINAPI WinMain(HINSTANCE instance, HINSTANCE prev_instance, LPSTR cmd_line, int show_cmd_line)
{
	query_pc_frequency();

	RGB8 clear_colour = {};
	Texture __window_back_buffer__ = {};
	TEXTURE_CLEAR(__window_back_buffer__, clear_colour);

	__window_back_buffer__.width = RENDER_TARGET_WIDTH;
	__window_back_buffer__.height = RENDER_TARGET_HEIGHT;

	int bytes_per_pixel = 4;
	int back_buffer_size = bytes_per_pixel * __window_back_buffer__.width * __window_back_buffer__.height;
	__window_back_buffer__.pixels = (uint8_t*)alloc(back_buffer_size);
	__window_back_buffer__.pixel_size = sizeof(uint32_t);

	init_profiling();

	render_image(&__window_back_buffer__);

	print_render_profile();
	print_profile();

	printf("Writing back buffer to file\n");
	output_to_bmp("output.bmp", __window_back_buffer__);

	return 0;
}
