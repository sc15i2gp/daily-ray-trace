#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include "Platform.h"
#include "Maths.h"
#include "Colour.h"

//CODE STANDARDS:
//	- No templates (they kinda suck and are hard to use)
//	- Keep to C features as much as possible, except for:
//		- Operator overloading

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
//	- Add debug info to surface points
//	- Texturing
//		- Colour/image
//		- Normal map
//		- Displacement map
//	- Skybox/infinite light/infinite geometry
//		- Skybox
//		- Sun
//		- Infinite ground plane
//	- Allow objects within transmission media
//		- Distinguish between air and vacuum
//		- Stack-like structure for tracking current medium
//	- Work out how to do sampling with time/animation stuff
//		- Go through monte carlo integration
//		- Go through light calculations and units
//		- Decide whether it is worth doing at this time
//	- Optimisation/Cleaning
//		- Remove superfluous code 
//		- Profiling
//		- Optimise slow methods

struct Pixel_Render_Buffer
{
	uint32_t* pixels;
	int width;
	int height;
};

inline
uint32_t* get_pixel(Pixel_Render_Buffer* r_buffer, int x, int y)
{
	return r_buffer->pixels + y*r_buffer->width + x;
}


void set_render_buffer_pixel_colour(Pixel_Render_Buffer* r_buffer, int x, int y, RGB8 colour)
{
	uint32_t* pixel = get_pixel(r_buffer, x, y);
	uint32_t* c = (uint32_t*)(&colour);
	*pixel = *c;
}

void clear_render_buffer(Pixel_Render_Buffer* r_buffer, RGB8 colour)
{
	for(int y = 0; y < r_buffer->height; ++y)
	{
		for(int x = 0; x < r_buffer->width; ++x) set_render_buffer_pixel_colour(r_buffer, x, y, colour);
	}
}

void output_to_ppm(const char* path, Pixel_Render_Buffer* r_buffer)
{
	char* contents = (char*)alloc(MEGABYTES(10));

	int contents_size = sprintf(contents, "P%d\n", 3);
	contents_size += sprintf(contents + contents_size, "%d %d\n", r_buffer->width, r_buffer->height);
	contents_size += sprintf(contents + contents_size, "255\n");
	for(int y = 0; y < r_buffer->height; ++y)
	{
		for(int x = 0; x < r_buffer->width; ++x)
		{
			RGB8 p = *((RGB8*)get_pixel(r_buffer, x, y));
			contents_size += sprintf(contents + contents_size, "%d %d %d\n", p.R, p.G, p.B);
		}
	}
	
	write_file_contents(path, contents, contents_size);
	dealloc(contents);
}

void invert_render_buffer(Pixel_Render_Buffer* r_buffer)
{
	for(int y = 0; y < r_buffer->height/2; ++y)
	{
		for(int x = 0; x < r_buffer->width; ++x)
		{
			uint32_t* p_0 = get_pixel(r_buffer, x,  y);
			uint32_t* p_1 = get_pixel(r_buffer, x, r_buffer->height - 1 - y);
			uint32_t p_temp = *p_0;
			*p_0 = *p_1;
			*p_1 = p_temp;
		}
	}

	for(int y = 0; y < r_buffer->height; ++y)
	{
		for(int x = 0; x < r_buffer->width/2; ++x)
		{
			uint32_t* p_0 = get_pixel(r_buffer, x, y);
			uint32_t* p_1 = get_pixel(r_buffer, r_buffer->width - 1 - x, y);
			uint32_t p_temp = *p_0;
			*p_0 = *p_1;
			*p_1 = p_temp;
		}
	}
}

struct Spectrum_Render_Buffer
{
	Spectrum* pixels;
	int width;
	int height;
};

inline
Spectrum* get_pixel(Spectrum_Render_Buffer* r_buffer, int x, int y)
{
	return r_buffer->pixels + y*r_buffer->width + x;
}

void set_render_buffer_pixel_spectrum(Spectrum_Render_Buffer* r_buffer, int x, int y, Spectrum s)
{
	Spectrum* pixel = get_pixel(r_buffer, x, y);
	*pixel = s;
}

void clear_render_buffer(Spectrum_Render_Buffer* r_buffer, Spectrum s)
{
	for(int y = 0; y < r_buffer->height; ++y)
	{
		for(int x = 0; x < r_buffer->width; ++x) set_render_buffer_pixel_spectrum(r_buffer, x, y, s);
	}
}

void write_spectrum_render_buffer_to_pixel_render_buffer(Spectrum_Render_Buffer* spectrum_buffer, Pixel_Render_Buffer* pixel_buffer)
{
	for(int y = 0; y < spectrum_buffer->height; ++y)
	{
		for(int x = 0; x < spectrum_buffer->width; ++x)
		{
			Spectrum* pixel_spectrum = get_pixel(spectrum_buffer, x, y);
			RGB64 pixel_colour_64 = spectrum_to_RGB64(*pixel_spectrum);
			RGB8 pixel_colour = rgb64_to_rgb8(pixel_colour_64);
			set_render_buffer_pixel_colour(pixel_buffer, x, y, pixel_colour);
		}
	}
}

BITMAPINFO __back_buffer_info__ = {};
Pixel_Render_Buffer __window_back_buffer__ = {};

bool running = true;

void size_window_back_buffer(Pixel_Render_Buffer* back_buffer, int new_back_buffer_width, int new_back_buffer_height)
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

void update_window_front_buffer(HWND window, Pixel_Render_Buffer* back_buffer)
{
	HDC window_device_context = GetDC(window);
	StretchDIBits
	(
		window_device_context, 
		0, 0, window_width(window), window_height(window), 
		back_buffer->width, back_buffer->height, -back_buffer->width, -back_buffer->height, 
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
			size_window_back_buffer(&__window_back_buffer__, window_width(window), window_height(window));
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

struct Surface_Point; //Forward decl
typedef Spectrum (*REFLECTION_MODEL_FUNCTION)(Surface_Point, Vec3, Vec3);
typedef double (*DISTRIBUTION_FUNCTION)(Surface_Point, Vec3, Vec3);
typedef Vec3 (*INDIRECT_SAMPLE_FUNCTION)(Surface_Point, Vec3, double*);

enum BSDF_TYPE
{
	BSDF_TYPE_DIFFUSE = 1 << 0,
	BSDF_TYPE_SPECULAR = 1 << 1
};

struct BSDF
{
	BSDF_TYPE type;
	//Computes proportion of light reflected at given point with given incident and reflection directions
	REFLECTION_MODEL_FUNCTION bsdf;
	//Returns probability of sampled direction
	DISTRIBUTION_FUNCTION pdf;
	//Samples a random direction at a given point with a given reflection direction
	//This is included here so that, for a reflection model, importance sampling can take place
	INDIRECT_SAMPLE_FUNCTION sample_direction;
};

enum Material_Type
{
	MAT_TYPE_NONE = 0,
	MAT_TYPE_DIELECTRIC,
	MAT_TYPE_CONDUCTOR
};

#define MAT_BSDF_MAX 8
struct Material
{
	Material_Type type;
	int number_of_bsdfs;
	BSDF bsdfs[MAT_BSDF_MAX];

	//Plastic data
	Spectrum diffuse_spd;
	Spectrum glossy_spd;
	double shininess;

	//Conductor data
	Spectrum refract_index;
	Spectrum extinct_index;
	bool is_dielectric;

	double roughness;
};

struct Surface_Point
{
	char* name;
	Spectrum diffuse_spd;
	Spectrum glossy_spd;
	Spectrum emission_spd;
	Spectrum incident_refract_index;
	Spectrum transmit_refract_index;
	Vec3 normal;
	Vec3 position;
	Material material;
	bool exists;
	bool is_emissive;
};

//REFLECTION MODELS
Spectrum diffuse_phong_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	return p.material.diffuse_spd / PI;
}

Spectrum glossy_phong_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	Vec3 bisector = (incoming + outgoing)/2.0;
	double specular_coefficient = pow(d_max(0.0, dot(p.normal, bisector)), p.material.shininess);
	return specular_coefficient * p.material.glossy_spd;
}

//TODO: Have specular reflectance and transmittance in material
Spectrum perfect_specular_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	if(incoming == reflect_vector(-outgoing, p.normal)) 
	{
		return generate_constant_spd(1.0);
	}
	else return Spectrum{};
}

double fresnel_reflectance_dielectric(double incident_refraction_index, double transmit_refraction_index, double incident_cos, double transmit_cos)
{
	double parallel_reflectance = 
		(transmit_refraction_index * incident_cos - incident_refraction_index * transmit_cos) / 
		(transmit_refraction_index * incident_cos + incident_refraction_index * transmit_cos);
	double perpendicular_reflectance = 
		(incident_refraction_index * incident_cos - transmit_refraction_index * transmit_cos) / 
		(incident_refraction_index * incident_cos + transmit_refraction_index * transmit_cos);

	double reflectance = 0.5 * (parallel_reflectance*parallel_reflectance + perpendicular_reflectance*perpendicular_reflectance);
	return reflectance;
}

double fresnel_reflectance_dielectric(double incident_refraction_index, double transmit_refraction_index, double incident_cos)
{
	double relative_refract_index = incident_refraction_index / transmit_refraction_index;
	double incident_sin_sq = d_max(0.0, 1.0 - incident_cos * incident_cos);
	double transmit_sin_sq = relative_refract_index * relative_refract_index * incident_sin_sq;
	
	//Total internal reflection
	if(transmit_sin_sq >= 1.0) 
	{
		return 1.0;
	}

	double transmit_cos = sqrt(d_max(0.0, 1.0 - transmit_sin_sq * transmit_sin_sq));

	return fresnel_reflectance_dielectric(incident_refraction_index, transmit_refraction_index, incident_cos, transmit_cos);
}


Spectrum fresnel_reflectance_dielectric(Spectrum incident_refract_index, Spectrum transmit_refract_index, double incident_cos)
{
	Spectrum reflectance = {};
	for(int i = 0; i < number_of_samples; ++i)
	{
		reflectance.samples[i] = fresnel_reflectance_dielectric(incident_refract_index.samples[i], transmit_refract_index.samples[i], incident_cos);
	}
	return reflectance;
}

double fresnel_transmittance_dielectric(double incident_refraction_index, double transmit_refraction_index, double incident_cos, double transmit_cos)
{
	double transmittance = 1.0 - fresnel_reflectance_dielectric(incident_refraction_index, transmit_refraction_index, incident_cos, transmit_cos);
	return transmittance;
}

double fresnel_reflectance_conductor(double incident_refract_index, double transmit_refract_index, double transmit_extinct_index, double incident_cos)
{
	double relative_refract_index = transmit_refract_index / incident_refract_index;
	double relative_extinct_index = transmit_extinct_index / incident_refract_index;

	double incident_cos_sq = incident_cos * incident_cos;
	double incident_sin_sq = 1.0 - incident_cos_sq;
	double relative_refract_index_sq = relative_refract_index * relative_refract_index;
	double relative_extinct_index_sq = relative_extinct_index * relative_extinct_index;
	double r = relative_refract_index_sq - relative_extinct_index_sq - incident_sin_sq;
	double a_sq_plus_b_sq = sqrt(r * r + 4.0 * relative_refract_index_sq * relative_extinct_index_sq);
	double a = sqrt(0.5 * (a_sq_plus_b_sq + r));
	double s = a_sq_plus_b_sq + incident_cos_sq;
	double t = 2.0 * a * incident_cos;
	double u = incident_cos_sq * a_sq_plus_b_sq + incident_sin_sq * incident_sin_sq;
	double v = t * incident_sin_sq;
	double parallel_reflectance = (s - t) / (s + t);
	double perpendicular_reflectance = parallel_reflectance * (u - v) / (u + v);

	return 0.5 * (parallel_reflectance + perpendicular_reflectance);
}

Spectrum fresnel_reflectance_conductor(Spectrum incident_refract_index, Spectrum transmit_refract_index, Spectrum transmit_extinct_index, double incident_cos)
{
	Spectrum reflectance = {};
	for(int i = 0; i < number_of_samples; ++i)
	{
		reflectance.samples[i] = fresnel_reflectance_conductor(incident_refract_index.samples[i], transmit_refract_index.samples[i], transmit_extinct_index.samples[i], incident_cos);
	}
	return reflectance;
}

Spectrum fresnel_transmittance_dielectric(Spectrum incident_refract_index, Spectrum transmit_refract_index, double incident_cos, double transmit_cos)
{
	Spectrum transmittance = {};
	for(int i = 0; i < number_of_samples; ++i)
	{
		transmittance.samples[i] = fresnel_transmittance_dielectric(incident_refract_index.samples[i], transmit_refract_index.samples[i], incident_cos, transmit_cos);
	}
	return transmittance;
}

//NOTE: The index of refraction chosen is the spectrum's value at 630 nm
//	Most materials' given refraction indices are the single values at 633 nm
//	630 is the closest value which can be used to this which is likely to work for any interval length
//	This is temporary and a better solution will be implemented as necessary.
#define TRANSMISSION_WAVELENGTH 630

Vec3 transmit_vector(double incident_refract_index, double transmit_refract_index, Vec3 normal, Vec3 outgoing)
{
	//Handles entering and leaving medium
	if(dot(outgoing, normal) < 0.0) normal = -normal;

	//Parallel and perpendicular to p.normal
	Vec3 transmit_perpendicular_component = (incident_refract_index/transmit_refract_index) * (dot(outgoing, normal)*normal -outgoing);
	Vec3 transmit_parallel_component = -sqrt(1.0 - dot(transmit_perpendicular_component, transmit_perpendicular_component))*normal;
	
	Vec3 transmit_direction = transmit_perpendicular_component + transmit_parallel_component;
	return transmit_direction;
}

Vec3 transmit_vector(Spectrum& incident_refract_index_spd, Spectrum& transmit_refract_index_spd, Vec3 normal, Vec3 outgoing)
{

	//Handles entering and leaving medium
	if(dot(outgoing, normal) < 0.0) normal = -normal;

	double incident_refract_index = spd_value_at_wavelength(incident_refract_index_spd, TRANSMISSION_WAVELENGTH);
	double transmit_refract_index = spd_value_at_wavelength(transmit_refract_index_spd, TRANSMISSION_WAVELENGTH);
	return transmit_vector(incident_refract_index, transmit_refract_index, normal, outgoing);
}

//TODO: Make it so specular bsdfs don't have to check if incoming is specular reflection vector
Spectrum fresnel_specular_reflection_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	if(incoming == reflect_vector(-outgoing, p.normal))
	{
		double incident_cos = abs(dot(outgoing, p.normal));
		double reflectance_cos = abs(dot(incoming, p.normal));

		if(p.material.type == MAT_TYPE_CONDUCTOR) return fresnel_reflectance_conductor(p.incident_refract_index, p.material.refract_index, p.material.extinct_index, incident_cos) / reflectance_cos;
		else if(p.material.type == MAT_TYPE_DIELECTRIC) return fresnel_reflectance_dielectric(p.incident_refract_index, p.transmit_refract_index, incident_cos) / reflectance_cos;
	}
	return Spectrum{};
}


Spectrum fresnel_transmission_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	Spectrum transmittance = {};
	if(incoming == transmit_vector(p.incident_refract_index, p.transmit_refract_index, p.normal, outgoing))
	{
		double incident_cos = abs(dot(outgoing, p.normal));
		double transmit_cos = abs(dot(incoming, p.normal));
		return fresnel_transmittance_dielectric(p.incident_refract_index, p.transmit_refract_index, incident_cos, transmit_cos) / transmit_cos;
	}
	return transmittance;
}

Spectrum fresnel_reflection_transmission_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	//TODO: Make bsdf fresnel specular names more consistent
	Spectrum f = fresnel_specular_reflection_bsdf(p, incoming, outgoing) + fresnel_transmission_bsdf(p, incoming, outgoing);
	return f;
}

//MATERIAL BSDFs
Spectrum plastic_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	if(dot(p.normal, incoming) > 0.0)
	{
		Spectrum diffuse = diffuse_phong_bsdf(p, incoming, outgoing);
		Spectrum glossy = glossy_phong_bsdf(p, incoming, outgoing);
		return diffuse + glossy;
	}
	return Spectrum{};
}

Spectrum mirror_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	if(dot(p.normal, incoming) > 0.0)
	{
		return perfect_specular_bsdf(p, incoming, outgoing);
	}
	return Spectrum {};
}

//NOTE: Isotropic roughness
double microfacet_distribution(Surface_Point p, Vec3 v)
{
	double cos_th_sq = dot(p.normal, v) * dot(p.normal, v);
	double sin_th_sq = 1.0 - cos_th_sq;
	double sin_th = sqrt(sin_th_sq);
	double tan_th_sq = sin_th_sq / cos_th_sq;
	double cos_ph = (sin_th == 0.0) ? 1.0 : clamp(v.x / sin_th, -1.0, 1.0);
	double sin_ph = (sin_th == 0.0) ? 0.0 : clamp(v.y / sin_th, -1.0, 1.0);
	double cos_ph_sq = cos_ph * cos_ph;
	double sin_ph_sq = sin_ph * sin_ph;

	double r_sq = p.material.roughness * p.material.roughness;
	double d = exp(-tan_th_sq / r_sq) / (PI * r_sq * cos_th_sq * cos_th_sq);

	return d;
}

double lambda(Surface_Point p, Vec3 v)
{
	double cos_th = dot(p.normal, v);
	double sin_th = sqrt(1.0 - cos_th * cos_th);
	double abs_tan_th = abs(sin_th/cos_th);
	double cos_ph = (sin_th == 0.0) ? 1.0 : clamp(v.x / sin_th, -1.0, 1.0);
	double sin_ph = (sin_th == 0.0) ? 0.0 : clamp(v.y / sin_th, -1.0, 1.0);
	double cos_ph_sq = cos_ph * cos_ph;
	double sin_ph_sq = sin_ph * sin_ph;

	double r_sq = p.material.roughness * p.material.roughness;
	double alpha = sqrt(cos_ph_sq * r_sq + sin_ph_sq * r_sq); //alpha = sqrt(2) * sigma
	double a = 1.0 / (alpha * abs_tan_th);
	
	if(a >= 1.6) return 1.0;

	double d = (1.0 - 1.259*a + 0.396*a*a) / (3.535*a + 2.181*a*a);
	return d;
}

double geometric_attenuation(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	return 1.0 / (1.0 + lambda(p, outgoing) + lambda(p, incoming));
}

double geometric_attenuation(Surface_Point p, Vec3 v)
{
	return 1.0 / (1.0 + lambda(p, v));
}

//Cite: Microfacet models for refraction through rough surfaces
//TODO: Compare performance of working out trig values from identities vs inverse functions
//	e.g. sin(acos(n_h_dot) = sqrt(1 - n_h_dot*n_h_dot) 
double beckmann_distribution(Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double n_h_dot = dot(surface_normal, microfacet_normal);
	if(n_h_dot <= 0.0) return 0.0;
	else
	{
		double a_sq = lobe_width * lobe_width;
		double cos_h_sq = n_h_dot * n_h_dot;
		double cos_h_quart = cos_h_sq * cos_h_sq;
		double tan_h_sq = (1.0/cos_h_sq) - 1.0;
		double d = 1.0/(PI * a_sq * cos_h_quart);
		double e = exp(-tan_h_sq/a_sq);

		return d * e;
		/*
		double nh_dot_sq = n_h_dot * n_h_dot;
		double e_quot = (nh_dot_sq - 1.0) / (lobe_width*lobe_width*nh_dot_sq);
		double denominator = PI * lobe_width * lobe_width * nh_dot_sq * nh_dot_sq;
		double d = exp(e_quot) / denominator;
		return d;
		*/
		
	}

}

double ggx_distribution(Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double nh_dot = dot(surface_normal, microfacet_normal);
	if(nh_dot <= 0.0) return 0.0;
	else
	{
		double a_sq = lobe_width * lobe_width;
		double cos_h_sq = nh_dot * nh_dot;
		double cos_h_quart = cos_h_sq * cos_h_sq;
		double tan_h_sq = (1.0/cos_h_sq) - 1.0;

		double d = a_sq / (PI * cos_h_quart * (a_sq + tan_h_sq) * (a_sq + tan_h_sq));
		return d;
	}
}

double ggx_g_1(Vec3 v, Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double vh_dot = dot(v, microfacet_normal);
	double vn_dot = dot(v, surface_normal);
	double vnh_dot_quot = abs(vh_dot / vn_dot);

	if(vnh_dot_quot <= 0.0) return 0.0;
	else
	{
		double a_sq = lobe_width * lobe_width;
		double tan_vn_sq = (1.0/(vn_dot*vn_dot)) - 1.0;
		double d = 2.0 / (1.0 + sqrt(1.0 + a_sq*tan_vn_sq));
		return d;
	}
}

//Cite: Microfacet models for refraction through rough surfaces
//TODO: Make these maths variable names more consistent/better thought out
double g_1(Vec3 v, Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double v_h_dot = dot(v, microfacet_normal);
	double v_n_dot = dot(v, surface_normal);
	double v_h_n_dot_quot = abs(v_h_dot / v_n_dot);
	double tan_v_n = sqrt((1.0/(v_n_dot*v_n_dot)) - 1.0);
	double a = 1.0/(lobe_width * tan_v_n);
	
	if(v_h_n_dot_quot <= 0.0) return 0.0;
	else if(a < 1.6)
	{
		double d = v_h_n_dot_quot * (3.535*a + 2.181*a*a)/(1.0 + 2.276*a + 2.577*a*a);
		return d;
	}
	else return 1.0;
}

double geometric_attenuation(Vec3 incoming, Vec3 outgoing, Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double g = ggx_g_1(outgoing, surface_normal, microfacet_normal, lobe_width) * ggx_g_1(incoming, surface_normal, microfacet_normal, lobe_width);
	return g;
}

//Cite: Microfacet models for refraction through rough surfaces
Spectrum cook_torrance_reflectance_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	Vec3 microfacet_normal = normalise(outgoing + incoming);
	Vec3 surface_normal = p.normal;
	double cos_sn_out = abs(dot(outgoing, surface_normal));
	double cos_sn_in = abs(dot(incoming, surface_normal));

	double cos_mn_out = abs(dot(outgoing, microfacet_normal));
	double cos_mn_in = abs(dot(incoming, microfacet_normal));

	Spectrum fr = {};

	if(p.material.type == MAT_TYPE_CONDUCTOR) fr = fresnel_reflectance_conductor(p.incident_refract_index, p.material.refract_index, p.material.extinct_index, cos_mn_out) / cos_mn_in;
	else if(p.material.type == MAT_TYPE_DIELECTRIC) fr = fresnel_reflectance_dielectric(p.incident_refract_index, p.transmit_refract_index, cos_mn_out) / cos_mn_in;

	Spectrum reflectance = ggx_distribution(surface_normal, microfacet_normal, p.material.roughness) * geometric_attenuation(incoming, outgoing, surface_normal, microfacet_normal, p.material.roughness) * fr;
	reflectance /= (4.0 * cos_sn_out * cos_sn_in);

	return reflectance;
}

//TODO: Function to get transmit_cos without actually computing transmission vector if possible
Spectrum torrance_sparrow_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	double cos_th_out = dot(outgoing, p.normal);
	double cos_th_in = dot(incoming, p.normal);

	Spectrum fr = {};
	double incident_cos = abs(cos_th_out);
	double reflectance_cos = abs(cos_th_in);

	if(p.material.type == MAT_TYPE_CONDUCTOR) fr = fresnel_reflectance_conductor(p.incident_refract_index, p.material.refract_index, p.material.extinct_index, incident_cos) / reflectance_cos;
	else if(p.material.type == MAT_TYPE_DIELECTRIC) fr = fresnel_reflectance_dielectric(p.incident_refract_index, p.transmit_refract_index, incident_cos) / reflectance_cos;

	Vec3 halfway = (outgoing + incoming)/2.0;
	Spectrum reflectance = microfacet_distribution(p, halfway) * geometric_attenuation(p, incoming, outgoing) * fr;
	reflectance /= (4.0 * cos_th_in * cos_th_out);

	return reflectance;
}

//GENERAL BSDF METHOD
Spectrum bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	Spectrum reflectance = {};
	BSDF* material_bsdfs = p.material.bsdfs;
	for(int i = 0; i < p.material.number_of_bsdfs; ++i)
	{
		reflectance += material_bsdfs[i].bsdf(p, incoming, outgoing);
	}
	return reflectance;
}

//DIRECTION SAMPLING METHODS
double diffuse_pdf(Surface_Point p, Vec3 outgoing, Vec3 sampled)
{
	return cos_weighted_sample_hemisphere_pdf(p.normal, sampled);
}

Vec3 sample_diffuse_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	return cos_weighted_sample_hemisphere(p.normal, pdf_value);
}

//NOTE: Although these two sample functions are the same now, the glossy one may be changed soon
Vec3 sample_glossy_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	/*
	*pdf_value = 1.0;
	Vec3 dir = reflect_vector(-outgoing, p.normal);
	return dir;
	*/
	return sample_diffuse_direction(p, outgoing, pdf_value);
}

Vec3 sample_cook_torrance_reflection_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	double f = uniform_sample();
	double g = uniform_sample();

	double a = p.material.roughness;
	double a_sq = a * a;
#if 0
	double l = log(1.0 - f);
	if(1.0 - f <= 0.0) l = 0.0;

	double tan_th_mn_sq = -a_sq * l;
#else
	double tan_th_mn = (a * sqrt(f)) / sqrt(1.0 - f);
	double tan_th_mn_sq = tan_th_mn * tan_th_mn;
#endif
	double cos_th_mn = 1.0 / sqrt(1.0 + tan_th_mn_sq);
	double sin_th_mn = sqrt(1.0 - cos_th_mn * cos_th_mn);
	double phi_mn = 2.0 * PI * g;

	Vec3 microfacet_normal = {sin_th_mn * cos(phi_mn), sin_th_mn * sin(phi_mn), cos_th_mn};
	Mat3x3 r = find_rotation_between_vectors(Vec3{0.0, 0.0, 1.0}, p.normal);
	microfacet_normal = r * microfacet_normal;
	double mn_sn_dot = dot(microfacet_normal, p.normal);
	if(mn_sn_dot < 0.0) 
	{
		microfacet_normal = -microfacet_normal;
		mn_sn_dot = -mn_sn_dot;
	}
	
	Vec3 reflection_direction = reflect_vector(-outgoing, microfacet_normal);
	//*pdf_value = dot(outgoing, microfacet_normal) * geometric_attenuation(reflection_direction, outgoing, p.normal, microfacet_normal, p.material.roughness) / (dot(outgoing, p.normal) * dot(microfacet_normal, p.normal));
	double d = ggx_distribution(p.normal, microfacet_normal, p.material.roughness) * mn_sn_dot;

	*pdf_value = d * (1.0/(4.0 * dot(outgoing, microfacet_normal)));

	return reflection_direction;
}

Vec3 sample_torrance_sparrow_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	double f = uniform_sample();
	double g = uniform_sample();

	double l = log(1.0 - f);
	if(1.0 - f <= 0.0) l = 0.0;
	
	double tan_th_sq = - p.material.roughness * p.material.roughness * l;
	double phi = 2.0 * PI * g;

	double cos_th = 1.0 / sqrt(1.0 + tan_th_sq);
	double sin_th = sqrt(1.0 - cos_th * cos_th);
	
	Vec3 halfway = {sin_th * cos(phi), sin_th * sin(phi), cos_th};
	if(dot(halfway, p.normal) < 0.0) halfway = -halfway;
	*pdf_value = microfacet_distribution(p, halfway) * geometric_attenuation(p, outgoing) * abs(dot(outgoing, halfway)) / abs(dot(outgoing, p.normal));

	return reflect_vector(-outgoing, halfway);
}

Vec3 sample_specular_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	*pdf_value = 1.0;
	return reflect_vector(-outgoing, p.normal);
}

Vec3 sample_specular_transmission_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	*pdf_value = 1.0;
	return transmit_vector(p.incident_refract_index, p.transmit_refract_index, p.normal, outgoing);
}


//TODO: Total internal reflection
//TODO: Fix this function
//	- transmit_cos doesn't need the transmission vector to work it out
//	- fresnel_reflectance should return 1.0 with total internal reflection
Vec3 sample_specular_reflection_or_transmission_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	double incident_refract_index = spd_value_at_wavelength(p.incident_refract_index, TRANSMISSION_WAVELENGTH);
	double transmit_refract_index = spd_value_at_wavelength(p.transmit_refract_index, TRANSMISSION_WAVELENGTH);
	double incident_cos = abs(dot(outgoing, p.normal));

	double reflectance = fresnel_reflectance_dielectric(incident_refract_index, transmit_refract_index, incident_cos);

	double f = uniform_sample();

	Vec3 sampled_direction = {};
	if(f < reflectance)
	{
		sampled_direction = reflect_vector(-outgoing, p.normal);
		*pdf_value = reflectance;
	}
	else
	{
		sampled_direction = transmit_vector(incident_refract_index, transmit_refract_index, p.normal, outgoing);
		*pdf_value = 1.0 - reflectance;
	}
	
	return sampled_direction;
}

Material create_plastic(Spectrum diffuse_spd, Spectrum glossy_spd, double shininess)
{
	Material plastic = {};
	//Reflectances
	plastic.type = MAT_TYPE_DIELECTRIC;
	plastic.diffuse_spd = diffuse_spd;
	plastic.glossy_spd = glossy_spd;
	plastic.shininess = shininess;
	
	//BSDFs
	plastic.number_of_bsdfs = 2;
	BSDF diffuse_reflection = {};
	diffuse_reflection.type = BSDF_TYPE_DIFFUSE;
	diffuse_reflection.pdf = diffuse_pdf;
	diffuse_reflection.bsdf = diffuse_phong_bsdf;
	diffuse_reflection.sample_direction = sample_diffuse_direction;

	BSDF glossy_reflection = {};
	glossy_reflection.type = BSDF_TYPE_DIFFUSE;
	glossy_reflection.pdf = diffuse_pdf;
	glossy_reflection.bsdf = glossy_phong_bsdf;
	glossy_reflection.sample_direction = sample_glossy_direction;
	plastic.bsdfs[0] = diffuse_reflection;
	plastic.bsdfs[1] = glossy_reflection;

	return plastic;
}

Material create_mirror()
{
	Material mirror = {};
	mirror.type = MAT_TYPE_CONDUCTOR;
	mirror.number_of_bsdfs = 1;
	mirror.bsdfs[0].type = BSDF_TYPE_SPECULAR;
	mirror.bsdfs[0].bsdf = perfect_specular_bsdf;
	mirror.bsdfs[0].sample_direction = sample_specular_direction;

	return mirror;
}

Material create_conductor(Spectrum refract_index, Spectrum extinct_index, double roughness)
{
	Material conductor = {};
	conductor.type = MAT_TYPE_CONDUCTOR;
	conductor.number_of_bsdfs = 1;

	if(roughness >= 0.001)
	{
		conductor.bsdfs[0].type = BSDF_TYPE_SPECULAR;
		conductor.bsdfs[0].bsdf = cook_torrance_reflectance_bsdf;
		conductor.bsdfs[0].sample_direction = sample_cook_torrance_reflection_direction;
	}
	else
	{
		conductor.bsdfs[0].type = BSDF_TYPE_SPECULAR;
		conductor.bsdfs[0].bsdf = fresnel_specular_reflection_bsdf;
		conductor.bsdfs[0].sample_direction = sample_specular_direction;
	}

	conductor.refract_index = refract_index;
	conductor.extinct_index = extinct_index;

	conductor.roughness = roughness;

	return conductor;
}

Material create_dielectric(Spectrum refract_index)
{
	Material dielectric = {};
	dielectric.type = MAT_TYPE_DIELECTRIC;
	dielectric.number_of_bsdfs = 1;
	dielectric.bsdfs[0].type = BSDF_TYPE_SPECULAR;
	dielectric.bsdfs[0].bsdf = fresnel_reflection_transmission_bsdf;
	dielectric.bsdfs[0].sample_direction = sample_specular_reflection_or_transmission_direction;

	dielectric.bsdfs[1].type = BSDF_TYPE_DIFFUSE;
	dielectric.bsdfs[1].pdf = diffuse_pdf;
	dielectric.bsdfs[1].bsdf = glossy_phong_bsdf;
	dielectric.bsdfs[1].sample_direction = sample_glossy_direction;

	dielectric.glossy_spd = generate_constant_spd(1.0);
	dielectric.shininess = 16.0;

	dielectric.refract_index = refract_index;

	return dielectric;
}

typedef Spectrum Radiance;

enum Light_Type
{
	LIGHT_TYPE_NONE = 0,
	LIGHT_TYPE_POINT,
	LIGHT_TYPE_AREA
};

struct Scene_Light
{
	Light_Type type;
	Spectrum emission_spd;
	int geometry;
	double area;
};

enum Geometry_Type
{
	GEO_TYPE_NONE = 0,
	GEO_TYPE_POINT,
	GEO_TYPE_SPHERE,
	GEO_TYPE_PLANE,
	GEO_TYPE_MODEL
};

struct Point
{
	Vec3 position;
};

struct Model_Vertex
{
	Vec3 position;
	Vec3 normal;
	Vec2 texture_coords;
};

struct Model
{
	int number_of_vertices;
	Model_Vertex* vertices;
};

struct Scene_Geometry
{
	Geometry_Type type;
	union
	{
		Point point;
		Sphere sphere;
		Plane plane;
		Model model;
	};
};

#define SCENE_OBJECT_MAX 32

struct Scene_Object
{
	char* name;
	Scene_Geometry geometry;
	Material material;
	Light_Type light_type;
	Spectrum emission_spd;
	bool is_emissive;
};

struct Scene
{
	int number_of_objects;
	Scene_Object objects[SCENE_OBJECT_MAX];
};

void add_object_to_scene(Scene* scene, Scene_Object object)
{
	scene->objects[scene->number_of_objects++] = object;
}

void add_point_light_to_scene(Scene* scene, Point p, Spectrum emission_spd)
{
	Scene_Object point_light = {};
	Scene_Geometry point_geometry = {};
	point_geometry.type = GEO_TYPE_POINT;
	point_geometry.point = p;
	
	point_light.geometry = point_geometry;
	point_light.emission_spd = emission_spd;
	point_light.light_type = LIGHT_TYPE_POINT;
	point_light.is_emissive = true;

	add_object_to_scene(scene, point_light);
}

void add_sphere_light_to_scene(Scene* scene, Sphere s, Spectrum emission_spd)
{
	Scene_Object sphere_light = {};
	Scene_Geometry sphere_geometry = {};
	sphere_geometry.type = GEO_TYPE_SPHERE;
	sphere_geometry.sphere = s;

	sphere_light.geometry = sphere_geometry;
	sphere_light.emission_spd = emission_spd;
	sphere_light.light_type = LIGHT_TYPE_AREA;
	sphere_light.is_emissive = true;

	add_object_to_scene(scene, sphere_light);
}

void add_plane_light_to_scene(Scene* scene, Plane p, Spectrum emission_spd, char* name)
{
	Scene_Object plane_light = {};
	Scene_Geometry plane_geometry = {};
	plane_geometry.type = GEO_TYPE_PLANE;
	plane_geometry.plane = p;

	plane_light.name = name;
	plane_light.geometry = plane_geometry;
	plane_light.emission_spd = emission_spd;
	plane_light.light_type = LIGHT_TYPE_AREA;
	plane_light.is_emissive = true;

	add_object_to_scene(scene, plane_light);
}

void add_sphere_to_scene(Scene* scene, Sphere s, Material material, char* name)
{
	Scene_Object sphere = {};
	Scene_Geometry sphere_geometry = {};
	sphere_geometry.type = GEO_TYPE_SPHERE;
	sphere_geometry.sphere = s;

	sphere.name = name;
	sphere.geometry = sphere_geometry;
	sphere.material = material;

	add_object_to_scene(scene, sphere);
}

void add_plane_to_scene(Scene* scene, Plane p, Material material, char* name)
{
	Scene_Object plane = {};
	Scene_Geometry plane_geometry = {};
	plane_geometry.type = GEO_TYPE_PLANE;
	plane_geometry.plane = p;

	plane.name = name;
	plane.geometry = plane_geometry;
	plane.material = material;

	add_object_to_scene(scene, plane);
}

void add_model_to_scene(Scene* scene, Model m, Material material, char* name)
{
	Scene_Object model = {};
	Scene_Geometry model_geometry = {};
	model_geometry.type = GEO_TYPE_MODEL;
	model_geometry.model = m;

	model.name = name;
	model.geometry = model_geometry;
	model.material = material;

	add_object_to_scene(scene, model);
}

bool ray_intersects_model(Ray ray, Model m, double* t, int* intersecting_triangle_index, double* a, double* b, double* c)
{
	for(int i = 0; i < m.number_of_vertices; i += 3)
	{
		Vec3 triangle_normal = (m.vertices[i].normal + m.vertices[i+1].normal + m.vertices[i+2].normal)/3.0;
		if(dot(ray.direction, triangle_normal) != 0.0)
		{
			double distance_to_triangle = dot(m.vertices[i].position - ray.origin, triangle_normal)/dot(ray.direction, triangle_normal);
			Vec3 p = ray.origin + distance_to_triangle * ray.direction;

			*a = sqrt(point_to_line_distance_sq(p, m.vertices[i+1].position, m.vertices[i+2].position)/point_to_line_distance_sq(m.vertices[i].position, m.vertices[i+1].position,m.vertices[i+2].position));
			*b = sqrt(point_to_line_distance_sq(p, m.vertices[i+2].position, m.vertices[i].position)/point_to_line_distance_sq(m.vertices[i+1].position, m.vertices[i+2].position, m.vertices[i].position));
			*c = sqrt(point_to_line_distance_sq(p, m.vertices[i].position, m.vertices[i+1].position)/point_to_line_distance_sq(m.vertices[i+2].position, m.vertices[i].position, m.vertices[i+1].position));

			double d = *a + *b + *c;
			if(d <= 1.00001)
			{
				*intersecting_triangle_index = i;
				*t = distance_to_triangle;
				return distance_to_triangle > 0.0;
			}
		}
	}
	return false;
}

Surface_Point find_ray_scene_intersection(Scene* scene, Ray ray)
{
	int intersecting_object = -1;
	double length_along_ray = DBL_MAX;
	Vec3 surface_normal = {};
	int intersect_triangle_index = -1;

	//Barycentric coordinates
	double bc_a = 0.0;
	double bc_b = 0.0;
	double bc_c = 0.0;
	for(int i = 0; i < scene->number_of_objects; ++i)
	{
		Scene_Geometry* ith_object_geometry = &(scene->objects[i].geometry);
		if(ith_object_geometry->type != GEO_TYPE_POINT)
		{
			bool ray_intersects_ith_object = false;
			double ith_length_along_ray = DBL_MAX;

			switch(ith_object_geometry->type)
			{
				case GEO_TYPE_SPHERE: 
				{
					ray_intersects_ith_object = ray_intersects_sphere(ray, ith_object_geometry->sphere, &ith_length_along_ray);
					break;
				}
				case GEO_TYPE_PLANE: 
				{
					ray_intersects_ith_object = ray_intersects_plane(ray, ith_object_geometry->plane, &ith_length_along_ray);
					break;
				}
				case GEO_TYPE_MODEL:
				{
					ray_intersects_ith_object = ray_intersects_model(ray, ith_object_geometry->model, &ith_length_along_ray, &intersect_triangle_index, &bc_a, &bc_b, &bc_c);
					break;
				}
			}

			if(ray_intersects_ith_object && ith_length_along_ray < length_along_ray)
			{
				intersecting_object = i;
				length_along_ray = ith_length_along_ray;
			}
		}
	}

	Surface_Point p = {};
	if(intersecting_object >= 0)
	{
		Vec3 intersection = ray.origin + length_along_ray * ray.direction;
		Scene_Object* object = scene->objects + intersecting_object;
		Scene_Geometry* geometry = &(object->geometry);
		switch(geometry->type)
		{
			case GEO_TYPE_SPHERE: 
			{
				surface_normal = normalise(intersection - geometry->sphere.center);
				break;
			}
			case GEO_TYPE_PLANE: 
			{
				surface_normal = geometry->plane.n;
				break;
			}
			case GEO_TYPE_MODEL:
			{
				Model m = geometry->model;
				int i = intersect_triangle_index;
				surface_normal = bc_a * m.vertices[i].normal + bc_b * m.vertices[i+1].normal + bc_c * m.vertices[i+2].normal;
				break;
			}
		}

		//If object transmits and the light is incident to the object
		//	p.incident_refract_index = vacuum_refract_index
		//	p.transmit_refract_index = object_refract_index
		//Else if object transmits and the light is transmitting through the object
		//	p.incident_refract_index = object_refract_index
		//	p.transmit_refract_index = vacuum_refract_index
		if(dot(ray.direction, surface_normal) < 0.0)
		{
			p.incident_refract_index = generate_constant_spd(1.0);
			p.transmit_refract_index = object->material.refract_index;
		}
		else
		{
			p.incident_refract_index = object->material.refract_index;
			p.transmit_refract_index = generate_constant_spd(1.0);
		}
		p.name = object->name;
		p.material = object->material;
		p.emission_spd = object->emission_spd;
		p.normal = surface_normal;
		p.position = intersection;
		p.exists = true;
		p.is_emissive = object->is_emissive;
	}

	return p;
}

double compute_area(Scene_Geometry geometry)
{
	switch(geometry.type)
	{
		case GEO_TYPE_POINT: return 0.0;
		case GEO_TYPE_SPHERE: return area(geometry.sphere);
		case GEO_TYPE_PLANE: return area(geometry.plane);
	}
	return 0.0;
}

Vec3 uniform_sample_geometry(Scene_Geometry geometry, Vec3 p, double* pdf)
{
	switch(geometry.type)
	{
		case GEO_TYPE_POINT: 
		{
			*pdf = 1.0;
			return geometry.point.position;
		}
		case GEO_TYPE_SPHERE: 
		{
			return uniform_sample_sphere_subtended(geometry.sphere, p, pdf);
		}
		case GEO_TYPE_PLANE: 
		{
			return uniform_sample_plane(geometry.plane, pdf);
		}
	}
	return Vec3{};
}

//TODO: Move this into direct_light_contribution
Vec3 sample_light_point(Scene_Object l, Vec3 p, double* pdf)
{
	Scene_Geometry light_geometry = l.geometry;
	double light_area = compute_area(light_geometry);
	Vec3 light_point = uniform_sample_geometry(light_geometry, p, pdf);
	return light_point;
}

Radiance direct_light_contribution(Scene* scene, Surface_Point p, Ray outgoing)
{
	Radiance contribution = {};
	if(!p.is_emissive)
	{
		for(int i = 0; i < scene->number_of_objects; ++i)
		{
			if(scene->objects[i].is_emissive)
			{
				double light_pdf = 1.0;
				Scene_Object light = scene->objects[i];
				Vec3 light_point = sample_light_point(light, p.position, &light_pdf);
				
				Ray shadow_ray = {};
				shadow_ray.direction = normalise(light_point - p.position);
				shadow_ray.origin = p.position + 0.001*shadow_ray.direction;

				Surface_Point shadow_test_point = find_ray_scene_intersection(scene, shadow_ray);
				double shadow_ray_length = length(shadow_test_point.position - p.position);
				double ray_length = length(light_point - p.position);
				double difference_between_ray_lengths = shadow_ray_length - ray_length;
				bool light_point_visible = !(shadow_test_point.exists && difference_between_ray_lengths < -5e-12);
				//bool light_point_visible = shadow_test_point.position == light_point;

				if(light_point_visible && light_pdf > 0.0)
				{
					Vec3 incoming = light_point - p.position;
					double dist = length(incoming);
					incoming = normalise(incoming);

					double attenuation_factor = (light.light_type == LIGHT_TYPE_POINT) ? 1.0/(4.0 * PI * dist * dist) : 1.0;
					double d = abs(dot(incoming, p.normal));
					double f = attenuation_factor * d / light_pdf;

					contribution += f * bsdf(p, incoming, outgoing.direction) * light.emission_spd;
				}
			}
		}
	}

	return contribution;
}

//MIRROR CHANGE: Probability of sampling specular direction depends on material
Vec3 choose_incoming_direction(Surface_Point p, Vec3 outgoing, double* pdf_value, bool* consider_emissive)
{
	double r = uniform_sample();
	double n = (double)p.material.number_of_bsdfs;
	int chosen_bsdf = (int)floor(r * n);
	BSDF_TYPE chosen_bsdf_type = p.material.bsdfs[chosen_bsdf].type;
	*consider_emissive = chosen_bsdf_type & BSDF_TYPE_SPECULAR;
	Vec3 direction = p.material.bsdfs[chosen_bsdf].sample_direction(p, outgoing, pdf_value);

	//Sum probabilities for BSDFs with same distributions (that aren't specular)
	for(int i = 0; i < p.material.number_of_bsdfs; ++i)
	{
		if(i != chosen_bsdf && chosen_bsdf_type != BSDF_TYPE_SPECULAR && p.material.bsdfs[i].type == chosen_bsdf_type)
		{
			*pdf_value += p.material.bsdfs[i].pdf(p, outgoing, direction);
		}
	}
	*pdf_value /= n;
	return direction;
}

int max_depth = 8; //NOTE: Arbitrarily chosen

Radiance cast_ray(Scene* scene, Ray eye_ray)
{
	Radiance eye_ray_radiance = {};
	Radiance direct_contribution = {};
	Ray outgoing = eye_ray;
	Ray incoming = {};
	Spectrum f = generate_constant_spd(1.0);
	double dir_pdf = 0.0;
	bool consider_emissive = true;
	bool collided_once = false;
	for(int depth = 0; depth < max_depth; ++depth)
	{
		Surface_Point p = find_ray_scene_intersection(scene, outgoing);
		if(p.exists && p.is_emissive && consider_emissive)
		{
			eye_ray_radiance += f * p.emission_spd;
			break;
		}
		else if(p.exists && !p.is_emissive)
		{
			collided_once = true;
			outgoing.direction = -outgoing.direction; //Reverse for bsdf computation, needs to start other way round for intersection test
			eye_ray_radiance += f * direct_light_contribution(scene, p, outgoing);
			
			//Choose new incoming direction
			incoming.direction = choose_incoming_direction(p, outgoing.direction, &dir_pdf, &consider_emissive);
			incoming.origin = p.position + 0.001*incoming.direction;
			
			if(dir_pdf <= 0.0) break;

			//Compute new direction pdf value
			f *= (abs(dot(p.normal, incoming.direction)/dir_pdf)) * bsdf(p, incoming.direction, outgoing.direction);
			
			outgoing = incoming;
			outgoing.origin += 0.001*outgoing.direction;
		}
		else break;
	}

	if(eye_ray_radiance.samples[0] < 0.0) 
	{
		printf("RADIANCE IS NEGATIVE\n");
	}
	else if(isnan(eye_ray_radiance.samples[0])) 
	{
		printf("RADIANCE IS NAN\n");
	}
	return eye_ray_radiance;
}

void raytrace_scene(Spectrum_Render_Buffer* render_target, Scene* scene, double fov, double focal_length, double focal_depth, double aperture_radius, int i)
{
	Vec3 image_plane_position = {0.0, 0.0, 8.0};
	Vec3 forward = {0.0, 0.0, -1.0};
	Vec3 right = {1.0, 0.0, 0.0};
	Vec3 up = {0.0, 1.0, 0.0};

	int image_plane_width_px = render_target->width;
	int image_plane_height_px = render_target->height;

	double aspect_ratio = (double)(image_plane_width_px)/(double)(image_plane_height_px);
	
	double image_plane_aperture_distance = (focal_length * focal_depth) / (focal_length + focal_depth);
	double image_plane_width = 2.0 * image_plane_aperture_distance * tan_deg(fov/2.0);
	double image_plane_height = image_plane_width / aspect_ratio;
	
	double pixel_width = image_plane_width/(double)(image_plane_width_px);
	double pixel_height = image_plane_height/(double)(image_plane_height_px);

	Vec3 image_plane_top_left = image_plane_position - 0.5 * image_plane_width * right + 0.5 * image_plane_height * up;

	Vec3 pinhole_position = image_plane_position + image_plane_aperture_distance * forward;

	Spectrum pixel_spectrum = {};
	Vec3 pixel_top_left = {};
	Vec3 sampled_pixel_point = {};
	Ray eye_ray = {};

	for(int y = 0; y < image_plane_height_px; ++y)
	{
		for(int x = 0; x < image_plane_width_px; ++x)
		{
			pixel_top_left = image_plane_top_left + ((double)x)*pixel_width*right - ((double)y)*pixel_height*up;
			Vec3 x_pixel_sample = (0.5 * pixel_width + (double)(i%2) * pixel_width) * Vec3{1.0, 0.0, 0.0};
			Vec3 y_pixel_sample = (0.5 * pixel_height + (double)(i / 2) * pixel_height) * Vec3{0.0, 1.0, 0.0};
			sampled_pixel_point = pixel_top_left + x_pixel_sample + y_pixel_sample;
			if(aperture_radius > 0)
			{
				Vec3 plane_of_focus_point = sampled_pixel_point + focal_depth * normalise(pinhole_position - sampled_pixel_point);
				Mat3x3 r = find_rotation_between_vectors(forward, Vec3{0.0, 0.0, 1.0});
				Vec3 sampled_lens_point = pinhole_position + r * ((focal_length / aperture_radius) * uniform_sample_disc());
				eye_ray.origin = sampled_lens_point;
				eye_ray.direction = normalise(plane_of_focus_point - eye_ray.origin);
			}
			else
			{
				eye_ray.origin = sampled_pixel_point;
				eye_ray.direction = normalise(pinhole_position - eye_ray.origin);
			}
			pixel_spectrum = cast_ray(scene, eye_ray);			
			set_render_buffer_pixel_spectrum(render_target, x, y, pixel_spectrum);
		}
	}
}

#define RENDER_TARGET_WIDTH 800
#define RENDER_TARGET_HEIGHT 600

Model create_triangle_model()
{
	Model triangle = {};
	triangle.number_of_vertices = 3;
	triangle.vertices = (Model_Vertex*)malloc(triangle.number_of_vertices * sizeof(Model_Vertex));

	double h = 2.0;
	double z = 2.0;

	Vec3 position_0 = {0.0f, h, z-4.0};
	Vec3 position_1 = {-h, -h, z};
	Vec3 position_2 = {h, -h, z};
	Vec3 normal = normalise(cross(position_2 - position_0, position_1 - position_0));

	//Top vertex
	triangle.vertices[0] =
	{
		position_0,
		normal,
		{0.5f, 1.0f}
	};

	//Bottom left vertex
	triangle.vertices[1] = 
	{
		position_1,
		normal,
		{0.0f, 0.0f}
	};

	//Bottom right vertex
	triangle.vertices[2] = 
	{
		position_2,
		normal,
		{1.0f, 0.0f}
	};

	return triangle;

}

void load_scene(Scene* scene)
{
	//Spectrum light_spd = generate_black_body_spd(4000.0);
	//normalise(light_spd);
	Spectrum light_spd = generate_constant_spd(1.0);
	Spectrum white_diffuse_spd = RGB64_to_spectrum(RGB64{0.55, 0.55, 0.55});
	Spectrum white_glossy_spd = RGB64_to_spectrum(RGB64{0.7, 0.7, 0.7});
	Spectrum red_diffuse_spd = RGB64_to_spectrum(RGB64{0.5, 0.0, 0.0});
	Spectrum red_glossy_spd = RGB64_to_spectrum(RGB64{0.7, 0.6, 0.6});
	Spectrum green_diffuse_spd = RGB64_to_spectrum(RGB64{0.1, 0.35, 0.1});
	Spectrum green_glossy_spd = RGB64_to_spectrum(RGB64{0.45, 0.55, 0.45});
	Spectrum blue_diffuse_spd = RGB64_to_spectrum(RGB64{0.2, 0.2, 0.8});
	Spectrum blue_glossy_spd = RGB64_to_spectrum(RGB64{0.8, 0.8, 0.9});
	Spectrum orange_diffuse_spd = RGB64_to_spectrum(RGB64{1.0, 0.64, 0.0});
	Spectrum orange_glossy_spd = RGB64_to_spectrum(RGB64{1.0, 0.8, 0.0});

	Spectrum sphere_diffuse_spd = RGB64_to_spectrum(RGB64{0.2, 0.8, 0.8});
	Spectrum sphere_glossy_spd = RGB64_to_spectrum(RGB64{0.8, 0.9, 0.9});

	Material white_material = create_plastic(white_diffuse_spd, white_glossy_spd, 32.0);
	Material red_material = create_plastic(red_diffuse_spd, red_glossy_spd, 32.0);
	Material green_material = create_plastic(green_diffuse_spd, green_glossy_spd, 32.0);
	Material blue_material = create_plastic(blue_diffuse_spd, blue_glossy_spd, 32.0);
	Material orange_material = create_plastic(orange_diffuse_spd, orange_glossy_spd, 32.0);
	Material sphere_material = create_plastic(sphere_diffuse_spd, sphere_glossy_spd, 50.0);
	Material mirror = create_mirror();
	Spectrum gold_refract_index = load_spd("au_spec_n.csv");
	Spectrum gold_extinct_index = load_spd("au_spec_k.csv");
	Material gold = create_conductor(gold_refract_index, gold_extinct_index, 0.344);
	Spectrum glass_refract_index = load_spd("glass.csv");
	Material glass = create_dielectric(glass_refract_index);

	Sphere glass_sphere =
	{
		{1.5, -1.8, 2.0}, 1.0
	};

	Sphere gold_sphere = 
	{
		{-2.0, -2.2, -0.5}, 0.75
	};

	Sphere mirror_sphere =
	{
		{0.0, 1.0, -1.0}, 1.0
	};

	Sphere plastic_sphere =
	{
		{2.0, -1.0, -2.0}, 1.0
	};

	Plane mirror_plane = create_plane_from_points(Vec3{-1.0, 1.0, -2.4}, Vec3{1.0, 1.0, -2.4}, Vec3{-1.0, -1.0, -2.9});

	Point light_p = {{0.0, 2.0, 2.0}};
	Sphere light_s = 
	{
		{0.0, 2.5, 2.0}, 0.5
	};

	double h = 3.0;
	Plane back_wall = create_plane_from_points(Vec3{-h, h, -h}, Vec3{h, h, -h}, Vec3{-h, -h, -h});
	Plane left_wall = create_plane_from_points(Vec3{-h, h, h}, Vec3{-h, h, -h}, Vec3{-h, -h, h});
	Plane right_wall = create_plane_from_points(Vec3{h, h, -h}, Vec3{h, h, h}, Vec3{h, -h, -h});
	Plane front_wall = create_plane_from_points(Vec3{h, h, h}, Vec3{-h, h, h}, Vec3{h, -h, h});
	Plane floor = create_plane_from_points(Vec3{-h, -h, -h}, Vec3{h, -h, -h}, Vec3{-h, -h, h});
	Plane ceiling = create_plane_from_points(Vec3{-h, h, h}, Vec3{h, h, h}, Vec3{-h, h, -h});

	Plane light_plane = create_plane_from_points(Vec3{-0.5, 2.9, 0.5}, Vec3{0.5, 2.9, 0.5}, Vec3{-0.5, 2.9, -0.5});

	Model triangle = create_triangle_model();
	//add_sphere_light_to_scene(scene, light_s, 10.0 * light_spd);
	add_plane_light_to_scene(scene, light_plane, light_spd, "Scene light");
	//add_point_light_to_scene(scene, light_p, 64.0 * light_spd);
	//add_sphere_to_scene(scene, mirror_sphere, mirror);
	//add_sphere_to_scene(scene, glass_sphere, glass);
	//add_sphere_to_scene(scene, gold_sphere, gold);
	//add_sphere_to_scene(scene, plastic_sphere, sphere_material);
	//add_plane_to_scene(scene, mirror_plane, mirror);

	add_plane_to_scene(scene, back_wall, blue_material, "Back wall");
	add_plane_to_scene(scene, left_wall, red_material, "Left wall");
	add_plane_to_scene(scene, right_wall, green_material, "Right wall");
	add_plane_to_scene(scene, floor, white_material, "Floor");
	add_plane_to_scene(scene, ceiling, white_material, "Ceiling");

	add_model_to_scene(scene, triangle, orange_material, "Model");
}

bool ready_to_display_spectrum_buffer = false;
bool completed_raytrace = false;
Spectrum_Render_Buffer final_image_buffer = {};
HWND window = {};

int number_of_render_samples = 25256;
double total_render_time = 0.0;
double average_sample_render_time = 0.0;
double max_sample_render_time = 0.0;
double min_sample_render_time = DBL_MAX;

void print_render_profile()
{
	printf("Time to render image: %fms\n", total_render_time);
	printf("Average sample render time: %fms\n", average_sample_render_time);
	printf("Max sample render time: %fms\n", max_sample_render_time);
	printf("Min sample render time: %fms\n", min_sample_render_time);
}

DWORD WINAPI render_image(LPVOID param)
{
	printf("Starting render...\n");
	Timer timer = {};

	Spectrum_Render_Buffer spectrum_buffer = {};
	spectrum_buffer.pixels = (Spectrum*)alloc(RENDER_TARGET_WIDTH * RENDER_TARGET_HEIGHT * sizeof(Spectrum));
	spectrum_buffer.width = window_width(window);
	spectrum_buffer.height = window_height(window);

	Scene scene = {};
	load_scene(&scene);

	for(int pass = 0; pass < number_of_render_samples; ++pass)
	{
		printf("Pass %d/%d\r", pass+1, number_of_render_samples);
		start_timer(&timer);
		
		//raytrace_scene(&spectrum_buffer, &scene, 90.0, 0.5, 80.0 1.4, pass % 4);
		raytrace_scene(&spectrum_buffer, &scene, 90.0, 0.5, 8.0, 0.0, pass % 4);

		for(int j = 0; j < final_image_buffer.width * final_image_buffer.height; ++j)
		{
			final_image_buffer.pixels[j] = ((double)pass * final_image_buffer.pixels[j] + spectrum_buffer.pixels[j])/(double)(pass+1);
		}

		stop_timer(&timer);

		double elapsed = elapsed_time_in_ms(&timer);
		total_render_time += elapsed;
		average_sample_render_time = ((double)pass * average_sample_render_time + elapsed)/(double)(pass + 1);
		max_sample_render_time = d_max(max_sample_render_time, elapsed);
		min_sample_render_time = d_min(min_sample_render_time, elapsed);

		ready_to_display_spectrum_buffer = true;
	}
	printf("Render completed\n");
	print_render_profile();

	completed_raytrace = true;
	return 0;
}

#define __USE_MINGW_ANSI_STDIO 1
int WINAPI WinMain(HINSTANCE instance, HINSTANCE prev_instance, LPSTR cmd_line, int show_cmd_line)
{
	srand(100000);

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

	query_pc_frequency();

	load_colour_data();

	RGB8 clear_colour = {};
	Spectrum clear_spectrum = {};

	final_image_buffer.pixels = (Spectrum*)alloc(RENDER_TARGET_WIDTH * RENDER_TARGET_HEIGHT * sizeof(Spectrum));
	final_image_buffer.width = window_width(window);
	final_image_buffer.height = window_height(window);
	clear_render_buffer(&__window_back_buffer__, clear_colour);
	clear_render_buffer(&final_image_buffer, clear_spectrum);

	HANDLE raytrace_thread = CreateThread(NULL, 0, render_image, NULL, 0, NULL);

	while(running && !completed_raytrace)
	{
		MSG message;
		while(PeekMessage(&message, 0, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&message);
			DispatchMessage(&message);
		}

		if(ready_to_display_spectrum_buffer)
		{
			write_spectrum_render_buffer_to_pixel_render_buffer(&final_image_buffer, &__window_back_buffer__);
			ready_to_display_spectrum_buffer = false;
		}

		update_window_front_buffer(window, &__window_back_buffer__);
	}
	
	TerminateThread(raytrace_thread, 0);
	//WaitForSingleObject(raytrace_thread, INFINITE);
	
	if(!completed_raytrace) print_render_profile();

	printf("Writing back buffer to file\n");
	invert_render_buffer(&__window_back_buffer__);
	output_to_ppm("output.ppm", &__window_back_buffer__);

	CloseHandle(raytrace_thread);
	return 0;
}
