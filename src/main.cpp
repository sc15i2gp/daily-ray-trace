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
//	- Glossy lobe sampling
//	- Easily compare sampling strategies
//	- Maybe use vec4 for everything (to make matrix operations slightly easier and 4 numbers easier to optimise than 3)
//	- Note failure points/error cases and handle
//	- Change reference white + find more precise spectra (maybe parameterise)
//	- Consider using explicitly sized types
//	- Scene editing
//		- OpenGL for rendering scene preview and UI
//		- UI for camera control, raytrace beginning
//		- Maybe screen shows progress of raytrace
//	- Volumetric transport
//	- Skybox/infinite light/infinite geometry (such as infinite ground plane)
//	- Better RNG
//	- Output to bmp/jpg/png
//	- Investigate different SPD representations
//	- Kirschoff's law for non-black body sources
//	- Different RGB -> SPD method: "Physically Meaningful Rendering using Tristimulus Colours"
//	- Texturing (images + surface properties)
//	- Different camera models/lenses (fisheye etc.)
//	- Try and get rid of reliance on plane normal direction (hard to use else)
//	- Direct MVSC build output files (that aren't exe) to build folder
//	- Invesitgate use of unity build
//	- Remove CSTDLIB
//		- Implement printf/sprintf/etc.
//		- Implement maths functions (trig, pow, ln etc.)
//		- Numeric limits (e.g. DBL_MAX)
//	- Add error handling to platform functions

//TODO: NOW
//	- More materials/effects
//		- Mirror
//		- Metal
//		- Glass
//	- Reduce variance
//		- Image plane super sampling
//		- Russian roulette with/without max depth
//	- Move platform code to platform file
//	- Optimise
//	- Quality of life

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
typedef Vec3 (*INDIRECT_SAMPLE_FUNCTION)(Surface_Point, Vec3, double*);

struct BSDF
{
	//Needed to check if the emissive component of an object in the reflected direction should be considered
	bool is_mirror_reflection; 
	//Computes proportion of light reflected at given point with given incident and reflection directions
	REFLECTION_MODEL_FUNCTION bsdf; 
	//Samples a random direction at a given point with a given reflection direction
	//This is included here so that, for a reflection model, importance sampling can take place
	INDIRECT_SAMPLE_FUNCTION sample_direction;
};

#define MAT_BSDF_MAX 8
struct Material
{
	int number_of_bsdfs;
	BSDF bsdfs[MAT_BSDF_MAX];
	Spectrum diffuse_spd;
	Spectrum glossy_spd;
	double shininess;
};

struct Surface_Point
{
	Spectrum diffuse_spd;
	Spectrum glossy_spd;
	Spectrum emission_spd;
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
	if(incoming == reflect_vector(outgoing, p.normal)) return generate_constant_spd(1.0);
	else return Spectrum{};
}

//MATERIAL BSDFs
Spectrum plastic_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	Spectrum diffuse = diffuse_phong_bsdf(p, incoming, outgoing);
	Spectrum glossy = glossy_phong_bsdf(p, incoming, outgoing);
	return diffuse + glossy;
}

Spectrum mirror_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	return perfect_specular_bsdf(p, incoming, outgoing);
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
Vec3 sample_diffuse_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	return cos_weighted_sample_hemisphere(p.normal, pdf_value);
}

//NOTE: Although these two sample functions are the same now, the glossy one may be changed soon
Vec3 sample_glossy_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	*pdf_value = 1.0;
	return reflect_vector(outgoing, p.normal);
}

Vec3 sample_specular_direction(Surface_Point p, Vec3 outgoing, double* pdf_value)
{
	*pdf_value = 1.0;
	return reflect_vector(outgoing, p.normal);
}

Material create_plastic(Spectrum diffuse_spd, Spectrum glossy_spd, double shininess)
{
	Material plastic = {};
	//Reflectances
	plastic.diffuse_spd = diffuse_spd;
	plastic.glossy_spd = glossy_spd;
	plastic.shininess = shininess;
	
	//BSDFs
	plastic.number_of_bsdfs = 2;
	BSDF diffuse_reflection = {};
	diffuse_reflection.bsdf = diffuse_phong_bsdf;
	diffuse_reflection.sample_direction = sample_diffuse_direction;
	BSDF glossy_reflection = {};
	glossy_reflection.bsdf = glossy_phong_bsdf;
	glossy_reflection.sample_direction = sample_glossy_direction;
	plastic.bsdfs[0] = diffuse_reflection;
	plastic.bsdfs[1] = glossy_reflection;

	return plastic;
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
	GEO_TYPE_TRIANGULATION
};

struct Point
{
	Vec3 position;
};

struct Scene_Geometry
{
	Geometry_Type type;
	union
	{
		Point point;
		Sphere sphere;
		Plane plane;
	};
};

#define SCENE_OBJECT_MAX 32

struct Scene_Object
{
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

void add_plane_light_to_scene(Scene* scene, Sphere s, Spectrum emission_spd)
{
	//TODO
}

void add_sphere_to_scene(Scene* scene, Sphere s, Material material)
{
	Scene_Object sphere = {};
	Scene_Geometry sphere_geometry = {};
	sphere_geometry.type = GEO_TYPE_SPHERE;
	sphere_geometry.sphere = s;

	sphere.geometry = sphere_geometry;
	sphere.material = material;

	add_object_to_scene(scene, sphere);
}

void add_plane_to_scene(Scene* scene, Plane p, Material material)
{
	Scene_Object plane = {};
	Scene_Geometry plane_geometry = {};
	plane_geometry.type = GEO_TYPE_PLANE;
	plane_geometry.plane = p;

	plane.geometry = plane_geometry;
	plane.material = material;

	add_object_to_scene(scene, plane);
}

Surface_Point find_ray_scene_intersection(Scene* scene, Ray ray)
{
	int intersecting_object = -1;
	double length_along_ray = DBL_MAX;
	Vec3 surface_normal = {};
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
		}

		p.material = object->material;
		p.emission_spd = object->emission_spd;
		p.normal = surface_normal;
		p.position = intersection;
		p.exists = true;
		p.is_emissive = object->is_emissive;
	}

	return p;
}

bool points_mutually_visible(Scene* scene, Vec3 p_0, Vec3 p_1)
{
	Ray shadow_ray = {};
	shadow_ray.direction = normalise(p_1 - p_0);
	shadow_ray.origin = p_0 + 0.001*shadow_ray.direction;

	Surface_Point shadow_test_point = find_ray_scene_intersection(scene, shadow_ray);
	double shadow_ray_length = length(shadow_test_point.position - p_0);
	double ray_length = length(p_1 - p_0);
	double difference_between_ray_lengths = shadow_ray_length - ray_length;
	return !(shadow_test_point.exists && difference_between_ray_lengths < -1e-12);
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
				
				if(points_mutually_visible(scene, p.position, light_point) && light_pdf > 0.0)
				{
					Vec3 incoming = light_point - p.position;
					double dist = length(incoming);
					incoming = normalise(incoming);

					double attenuation_factor = (light.light_type == LIGHT_TYPE_POINT) ? 1.0/(4.0 * PI * dist * dist) : 1.0;
					double d = dot(incoming, p.normal);
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
	*consider_emissive = p.material.bsdfs[chosen_bsdf].is_mirror_reflection;
	return p.material.bsdfs[chosen_bsdf].sample_direction(p, outgoing, pdf_value);
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
		Surface_Point p = find_ray_scene_intersection(scene, outgoing);
		if(p.exists && p.is_emissive && consider_emissive)
		{
			eye_ray_radiance += p.emission_spd;
		}
		else if(p.exists && !p.is_emissive)
		{
			outgoing.direction = -outgoing.direction; //Reverse for bsdf computation, needs to start other way round for intersection test
			eye_ray_radiance += f * direct_light_contribution(scene, p, outgoing);
			
			//Choose new incoming direction
			incoming.direction = choose_incoming_direction(p, outgoing.direction, &dir_pdf, &consider_emissive);
			incoming.origin = p.position + 0.001*incoming.direction;

			//Compute new direction pdf value
			f *= (dot(p.normal, incoming.direction)/dir_pdf) * bsdf(p, outgoing.direction, incoming.direction);
			
			outgoing = incoming;
			outgoing.origin += 0.001*outgoing.direction;
		}
		else break;
	}

	return eye_ray_radiance;
}

void raytrace_scene(Spectrum_Render_Buffer* render_target, double fov, double near_plane, Scene* scene)
{
	Vec3 eye = {0.0, 0.0, 8.0};
	Vec3 forward = {0.0, 0.0, -1.0};
	Vec3 right = {1.0, 0.0, 0.0};
	Vec3 up = {0.0, 1.0, 0.0};

	int image_plane_width_px = render_target->width;
	int image_plane_height_px = render_target->height;
	double aspect_ratio = (double)(image_plane_width_px)/(double)(image_plane_height_px);
	double image_plane_width = 2.0 * near_plane * tan_deg(fov/2.0);
	double image_plane_height = image_plane_width / aspect_ratio;
	double pixel_width = image_plane_width/(double)(image_plane_width_px);
	double pixel_height = image_plane_height/(double)(image_plane_height_px);

	Vec3 image_plane_top_left = eye + near_plane * forward - 0.5 * image_plane_width * right + 0.5 * image_plane_height * up;
	Spectrum pixel_spectrum = {};
	Vec3 pixel_center = {};
	Ray eye_ray = {};
	for(int y = 0; y < image_plane_height_px; ++y)
	{
		for(int x = 0; x < image_plane_width_px; ++x)
		{
			pixel_center = image_plane_top_left + ((double)x + 0.5)*pixel_width*right - ((double)y + 0.5)*pixel_height*up;
			eye_ray.origin = eye;
			eye_ray.direction = normalise(pixel_center - eye_ray.origin);
			pixel_spectrum = cast_ray(scene, eye_ray);			
			set_render_buffer_pixel_spectrum(render_target, x, y, pixel_spectrum);
		}
	}
}

#define RENDER_TARGET_WIDTH 800
#define RENDER_TARGET_HEIGHT 600

void load_scene(Scene* scene)
{
	//Spectrum light_spd = generate_black_body_spd(4000.0);
	//normalise(light_spd);
	Spectrum light_spd = generate_constant_spd(1.0);
	Spectrum white_diffuse_spd = RGB64_to_spectrum(RGB64{0.8, 0.8, 0.8});
	Spectrum white_glossy_spd = RGB64_to_spectrum(RGB64{0.9, 0.9, 0.9});
	Spectrum red_diffuse_spd = RGB64_to_spectrum(RGB64{0.8, 0.2, 0.2});
	Spectrum red_glossy_spd = RGB64_to_spectrum(RGB64{0.9, 0.8, 0.8});
	Spectrum green_diffuse_spd = RGB64_to_spectrum(RGB64{0.2, 0.8, 0.2});
	Spectrum green_glossy_spd = RGB64_to_spectrum(RGB64{0.8, 0.9, 0.8});

	Spectrum sphere_diffuse_spd = RGB64_to_spectrum(RGB64{0.2, 0.2, 0.8});
	Spectrum sphere_glossy_spd = RGB64_to_spectrum(RGB64{0.8, 0.8, 0.9});

	Material white_material = create_plastic(white_diffuse_spd, white_glossy_spd, 28.0);
	Material red_material = create_plastic(red_diffuse_spd, red_glossy_spd, 100.0);
	Material green_material = create_plastic(green_diffuse_spd, green_glossy_spd, 100.0);
	Material sphere_material = create_plastic(sphere_diffuse_spd, sphere_glossy_spd, 100.0);

	Sphere sphere = 
	{
		{-0.5, -1.0, 1.0}, 0.5
	};
	Point light_p = {};
	Sphere light_s = 
	{
		{0.0, 2.0, 2.0}, 0.5
	};
	double h = 3.0;
	Plane back_wall = create_plane_from_points(Vec3{-h, h, -h}, Vec3{h, h, -h}, Vec3{-h, -h, -h});
	Plane left_wall = create_plane_from_points(Vec3{-h, h, h}, Vec3{-h, h, -h}, Vec3{-h, -h, h});
	Plane right_wall = create_plane_from_points(Vec3{h, h, -h}, Vec3{h, h, h}, Vec3{h, -h, -h});
	Plane floor = create_plane_from_points(Vec3{-h, -h, -h}, Vec3{h, -h, -h}, Vec3{-h, -h, h});
	Plane ceiling = create_plane_from_points(Vec3{-h, h, h}, Vec3{h, h, h}, Vec3{-h, h, -h});
	add_sphere_light_to_scene(scene, light_s, 2.0 * light_spd);
	//add_point_light_to_scene(scene, light_p, 64.0 * light_spd);
	add_sphere_to_scene(scene, sphere, sphere_material);
	add_plane_to_scene(scene, back_wall, white_material);
	add_plane_to_scene(scene, left_wall, red_material);
	add_plane_to_scene(scene, right_wall, green_material);
	add_plane_to_scene(scene, floor, white_material);
	add_plane_to_scene(scene, ceiling, white_material);

}

bool ready_to_display_spectrum_buffer = false;
bool completed_raytrace = false;
Spectrum_Render_Buffer final_image_buffer = {};
HWND window = {};

int number_of_render_samples = 16;
double total_render_time = 0.0;
double average_sample_render_time = 0.0;
double max_sample_render_time = 0.0;
double min_sample_render_time = DBL_MAX;

void print_render_profile()
{
	printf("Time to render %d samples: %fms\n", number_of_render_samples, total_render_time);
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
		printf("Starting pass %d\n", pass);
		start_timer(&timer);
		
		raytrace_scene(&spectrum_buffer, 90.0, 0.1, &scene);

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
		printf("Completed pass %d\n", pass);
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
	output_to_ppm("output.ppm", &__window_back_buffer__);

	CloseHandle(raytrace_thread);
	return 0;
}
