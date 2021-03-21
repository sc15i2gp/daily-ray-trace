#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include "Platform.h"
#include "Maths.h"
#include "Colour.h"
#include "bsdf.h"

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
//	- Texturing
//		- Spectrum any material spectrum
//		- Double shininess and roughness
//		- Vector normal map
//	- Tetrahedron
//	- Frosted glass
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
//		- Sort out floating point precision issues
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

struct Geometry_Intersection_Point
{
	//General intersection data
	int scene_object;
	Vec3 position;
	
	//Model geometry intersection data
	int model_triangle_index;
	double bc_a;
	double bc_b;
	double bc_c;
};

Geometry_Intersection_Point find_ray_scene_intersection(Scene* scene, Ray ray)
{
	Geometry_Intersection_Point intersection_point;
	intersection_point.scene_object = -1;

	double length_along_ray = DBL_MAX;
	double ith_length_along_ray = DBL_MAX;
	for(int i = 0; i < scene->number_of_objects; ++i)
	{
		Scene_Geometry* ith_object_geometry = &(scene->objects[i].geometry);
		if(ith_object_geometry->type)
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
					ray_intersects_ith_object = ray_intersects_model(ray, ith_object_geometry->model, &ith_length_along_ray, &intersection_point.model_triangle_index, &intersection_point.bc_a, &intersection_point.bc_b, &intersection_point.bc_c);
					break;
				}
			}

			if(ray_intersects_ith_object && ith_length_along_ray < length_along_ray)
			{
				intersection_point.scene_object = i;
				length_along_ray = ith_length_along_ray;
			}
		}

	}
	intersection_point.position = ray.origin + length_along_ray * ray.direction;

	return intersection_point;
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
				Vec3 light_point = {};
				switch(light.geometry.type)
				{
					case GEO_TYPE_POINT: 
					{
						light_pdf = 1.0;
						light_point = light.geometry.point.position;
					}
					case GEO_TYPE_SPHERE: 
					{
						light_point = uniform_sample_sphere_subtended(light.geometry.sphere, p.position, &light_pdf);
					}
					case GEO_TYPE_PLANE: 
					{
						light_point = uniform_sample_plane(light.geometry.plane, &light_pdf);
					}
				}
				
				Ray shadow_ray = {};
				shadow_ray.direction = normalise(light_point - p.position);
				shadow_ray.origin = p.position + 0.001*shadow_ray.direction;

				Geometry_Intersection_Point shadow_test_point = find_ray_scene_intersection(scene, shadow_ray);

				double shadow_ray_length = length(shadow_test_point.position - p.position);
				double ray_length = length(light_point - p.position);
				double difference_between_ray_lengths = shadow_ray_length - ray_length;
				bool light_point_visible = !(difference_between_ray_lengths < -5e-12);

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
	for(int depth = 0; depth < max_depth; ++depth)
	{
		Geometry_Intersection_Point gp = find_ray_scene_intersection(scene, outgoing);
		Surface_Point p = {};
		if(gp.scene_object >= 0)
		{
			Vec3 surface_normal = {};
			Vec2 texture_coordinates = {};
			Scene_Object* object = scene->objects + gp.scene_object;
			Scene_Geometry* geometry = &(object->geometry);
			switch(geometry->type)
			{
				case GEO_TYPE_SPHERE: 
				{
					surface_normal = normalise(gp.position - geometry->sphere.center);
					//TODO: texture_coords = ?
					break;
				}
				case GEO_TYPE_PLANE: 
				{
					surface_normal = geometry->plane.n;
					//TODO: texture_coords = ?
					break;
				}
				case GEO_TYPE_MODEL:
				{
					Model m = geometry->model;
					int i = gp.model_triangle_index;
					surface_normal = gp.bc_a * m.vertices[i].normal + gp.bc_b * m.vertices[i+1].normal + gp.bc_c * m.vertices[i+2].normal;
					texture_coordinates = gp.bc_a * m.vertices[i].texture_coords + gp.bc_b * m.vertices[i+1].texture_coords + gp.bc_c * m.vertices[i+2].texture_coords;
					break;
				}
			}

			//If object transmits and the light is incident to the object
			if(dot(outgoing.direction, surface_normal) < 0.0)
			{
				p.incident_refract_index = generate_constant_spd(1.0);
				p.transmit_refract_index = object->material.refract_index;
			}
			//Else if object transmits and the light is transmitting through the object
			else
			{
				p.incident_refract_index = object->material.refract_index;
				p.transmit_refract_index = generate_constant_spd(1.0);
			}
			p.name = object->name;
			p.material = object->material;
			p.emission_spd = object->emission_spd;
			p.normal = surface_normal;
			p.position = gp.position;
			p.exists = true;
			p.is_emissive = object->is_emissive;
		}

		if(p.exists && p.is_emissive && consider_emissive)
		{
			eye_ray_radiance += f * p.emission_spd;
			break;
		}
		else if(p.exists && !p.is_emissive)
		{
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
