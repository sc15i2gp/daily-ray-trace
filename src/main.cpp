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

//TODO: LONGTERM
//	- Output to jpg/png
//	- Investigate different SPD representations
//	- Kirschoff's law for non-black body sources
//	- Different RGB -> SPD method: "Physically Meaningful Rendering using Tristimulus Colours"
//	- Texturing (images + surface properties)
//	- Different camera models/lenses (fisheye etc.)
//	- Direct MVSC build output files (that aren't exe) to build folder
//	- Invesitgate use of unity build
//	- Remove CSTDLIB
//		- Implement printf/sprintf/etc.
//		- Implement maths functions (trig, pow, ln etc.)
//		- Numeric limits (e.g. DBL_MAX)
//	- Add error handling to platform functions

//TODO: NOW
//	- Recursive raytrace
//		- Area lighting
//		- Indirect lighting
//	- Profiling
//	- Output raytraced scene to file
//		- Output raytraced image as bmp file
//		- OpenGL for rendering scene preview and UI
//		- UI for camera control, raytrace beginning
//		- Maybe screen shows progress of raytrace
//	- Reduce variance
//		- Integration importance sampling
//		- Image plane super sampling
//		- Russian roulette with/without max depth
//	- More geometry models + materials
//		- Metal
//		- Glass
//		- Mirror
//	- Volumetric transport
//	- Skybox/infinite light/infinite geometry (such as infinite ground plane)
//	- Investigate standardising spectra wavelength ranges in the program
//		- Make sure all spectra computed are consistent (have same wavelength range + number of samples)
//		- Maybe parameterise these

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

struct Surface_Point
{
	Spectrum diffuse_spd;
	Spectrum glossy_spd;
	Spectrum emission_spd;
	Vec3 normal;
	Vec3 position;
	bool exists;
	bool is_emissive;
};


Spectrum diffuse_phong_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	return p.diffuse_spd / PI;
}

double d_max(double a, double b)
{
	return (a < b) ? b : a;
}

Spectrum glossy_phong_bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	Vec3 bisector = (incoming + outgoing)/2.0;
	double specular_coefficient = pow(d_max(0.0, dot(p.normal, bisector)), 32.0);
	return specular_coefficient * p.glossy_spd;
}

Spectrum bsdf(Surface_Point p, Vec3 incoming, Vec3 outgoing)
{
	Spectrum diffuse = diffuse_phong_bsdf(p, incoming, outgoing);
	Spectrum glossy = glossy_phong_bsdf(p, incoming, outgoing);
	return diffuse + glossy;
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
	Spectrum diffuse_spd;
	Spectrum glossy_spd;
	Spectrum emission_spd;
	Light_Type light_type;
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

void add_sphere_to_scene(Scene* scene, Sphere s, Spectrum diffuse_spd, Spectrum glossy_spd)
{
	Scene_Object sphere = {};
	Scene_Geometry sphere_geometry = {};
	sphere_geometry.type = GEO_TYPE_SPHERE;
	sphere_geometry.sphere = s;

	sphere.geometry = sphere_geometry;
	sphere.diffuse_spd = diffuse_spd;
	sphere.glossy_spd = glossy_spd;

	add_object_to_scene(scene, sphere);
}

void add_plane_to_scene(Scene* scene, Plane p, Spectrum diffuse_spd, Spectrum glossy_spd)
{
	Scene_Object plane = {};
	Scene_Geometry plane_geometry = {};
	plane_geometry.type = GEO_TYPE_PLANE;
	plane_geometry.plane = p;

	plane.geometry = plane_geometry;
	plane.diffuse_spd = diffuse_spd;
	plane.glossy_spd = glossy_spd;

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

		p.diffuse_spd = object->diffuse_spd;
		p.glossy_spd = object->glossy_spd;
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
	return !(shadow_test_point.exists && length(shadow_test_point.position - p_0) < length(p_1 - p_0));
}

Radiance cast_ray(Scene*, Ray, bool consider_emissive, int depth);

Radiance indirect_light_contribution(Scene* scene, Surface_Point p, Ray outgoing, int depth)
{
	Radiance contribution = {};
	contribution.start_wavelength = 380.0;
	contribution.end_wavelength = 720.0;
	contribution.number_of_samples = 69;
	if(!p.is_emissive)
	{
		//NOTE: Does not account for specular reflection
		Vec3 new_direction = {};
		double direction_pdf = 0.0;
			
		new_direction = uniform_sample_hemisphere(p.normal);
		direction_pdf = 1.0/(2.0 * PI);
		Ray new_ray = {};
		new_ray.direction = normalise(new_direction);
		new_ray.origin = p.position + 0.001*new_ray.direction;

		Spectrum indirect_contribution = cast_ray(scene, new_ray, false, depth+1);
		contribution += dot(new_ray.direction, outgoing.direction) * (bsdf(p, new_ray.direction, outgoing.direction) * indirect_contribution)/direction_pdf;
	}

	return contribution;
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

Vec3 uniform_sample_geometry(Scene_Geometry geometry)
{
	switch(geometry.type)
	{
		case GEO_TYPE_POINT: return geometry.point.position;
		case GEO_TYPE_SPHERE: return uniform_sample_sphere(geometry.sphere);
		case GEO_TYPE_PLANE: return uniform_sample_plane(geometry.plane);
	}
	return Vec3{};
}

Vec3 sample_light_point(Scene_Object l, double* pdf)
{
	Scene_Geometry light_geometry = l.geometry;
	double light_area = compute_area(light_geometry);
	Vec3 light_point = uniform_sample_geometry(light_geometry);
	*pdf = (l.light_type == LIGHT_TYPE_POINT) ? 1.0 : (1.0 / light_area);
	return light_point;
}

Radiance direct_light_contribution(Scene* scene, Surface_Point p, Ray outgoing)
{
	Radiance contribution = {};
	contribution.start_wavelength = 380.0;
	contribution.end_wavelength = 720.0;
	contribution.number_of_samples = 69;
	if(!p.is_emissive)
	{
		outgoing.direction = -outgoing.direction;
		for(int i = 0; i < scene->number_of_objects; ++i)
		{
			if(scene->objects[i].is_emissive)
			{
				double light_pdf = 1.0;
				Scene_Object light = scene->objects[i];
				Vec3 light_point = sample_light_point(light, &light_pdf);
				
				if(points_mutually_visible(scene, p.position, light_point))
				{
					Vec3 incoming = light_point - p.position;
					double dist = length(incoming);
					incoming = normalise(incoming);

					double attenuation_factor = (light.light_type == LIGHT_TYPE_POINT) ? 1.0/(4.0 * PI * dist * dist) : 1.0;
					double f = attenuation_factor * dot(incoming, p.normal) / light_pdf;

					contribution += f * bsdf(p, incoming, outgoing.direction) * light.emission_spd;
				}
			}
		}
	}

	return contribution;
}

int max_depth = 3; //NOTE: Arbitrarily chosen
Radiance cast_ray(Scene* scene, Ray ray, bool consider_emissive, int depth)
{
	Radiance ray_radiance = {};
	ray_radiance.start_wavelength = 380.0;
	ray_radiance.end_wavelength = 720.0;
	ray_radiance.number_of_samples = 69;
	
	if(depth < max_depth)
	{
		Surface_Point p = find_ray_scene_intersection(scene, ray);
		if(p.exists)
		{
			if(consider_emissive)
			{//If eye ray or specular reflection ray
				ray_radiance += p.emission_spd;
			}
			ray_radiance += direct_light_contribution(scene, p, ray) + indirect_light_contribution(scene, p, ray, depth);
		}
	}

	return ray_radiance;
}

void raytrace_scene(Render_Buffer* render_target, double fov, double near_plane, Scene* scene)
{
	Vec3 eye = {0.0, 0.0, 5.0};
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
			RGB64 raycast_result = spectrum_to_RGB64(cast_ray(scene, eye_ray, true, 0));
			set_render_buffer_pixel_colour(render_target, x, y, rgb64_to_rgb8(raycast_result));
		}
	}
}


RGB64 cast_ray(Ray ray, Spectrum light_spd, Spectrum sphere_spd, Spectrum specular_spd)
{
	RGB64 black = {0.0, 0.0, 0.0};
	double t = 0.0;
	Plane p = {};
	p.p = Vec3{-0.2, 0.5, -0.2};
	p.u = 0.5*normalise(Vec3{1.0, 0.0, 0.0});
	p.v = normalise(Vec3{0.0, -0.6, 0.5});
	p.n = normalise(cross(p.u, p.v));

	Sphere s = {};
	s.radius = 0.3;

	Vec3 light_pos = {0.0, 1.0, 1.0};
	if(ray_intersects_sphere(ray, s, &t))
	{
		Vec3 intersection = ray.origin + t * ray.direction;
		Vec3 incoming = light_pos - intersection;
		Vec3 outgoing = -ray.direction;
		Surface_Point point = {};
		point.diffuse_spd = sphere_spd;
		point.glossy_spd = sphere_spd;
		point.position = intersection;
		point.normal = normalise(intersection - s.center);
		//point.normal = -p.n;
		Spectrum result = dot(incoming, point.normal) * light_spd * bsdf(point, incoming, outgoing); 
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
	srand(NULL);
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


	load_colour_data();

	long double l = 4000.0L;
	//long double l = 10000.0L;
	Spectrum light_spd = generate_black_body_spd(l, 380.0L, 720.0L);
	normalise(light_spd);
	Spectrum mat_spd = RGB64_to_spectrum(RGB64{0.25, 0.8, 0.4});
	Spectrum glossy_spd = RGB64_to_spectrum(RGB64{0.6, 0.9, 0.8});
	Spectrum white_diffuse_spd = RGB64_to_spectrum(RGB64{0.8, 0.8, 0.8});
	Spectrum white_glossy_spd = RGB64_to_spectrum(RGB64{0.9, 0.9, 0.9});
	Spectrum red_diffuse_spd = RGB64_to_spectrum(RGB64{0.8, 0.2, 0.2});
	Spectrum red_glossy_spd = RGB64_to_spectrum(RGB64{0.9, 0.8, 0.8});
	Spectrum green_diffuse_spd = RGB64_to_spectrum(RGB64{0.2, 0.8, 0.2});
	Spectrum green_glossy_spd = RGB64_to_spectrum(RGB64{0.8, 0.9, 0.8});

	Sphere sphere = {};
	sphere.radius = 0.5;
	double h = 2.0;
	Plane back_wall = create_plane_from_points(Vec3{-h, h, -h}, Vec3{h, h, -h}, Vec3{-h, -h, -h});
	Plane left_wall = create_plane_from_points(Vec3{-h, h, h}, Vec3{-h, h, -h}, Vec3{-h, -h, h});
	Plane right_wall = create_plane_from_points(Vec3{h, h, -h}, Vec3{h, h, h}, Vec3{h, -h, -h});
	Plane floor = create_plane_from_points(Vec3{-h, -h, -h}, Vec3{h, -h, -h}, Vec3{-h, -h, h});
	Plane ceiling = create_plane_from_points(Vec3{-h, h, h}, Vec3{h, h, h}, Vec3{-h, h, -h});
	Scene scene = {};
	add_sphere_light_to_scene(&scene, sphere, light_spd);
	add_plane_to_scene(&scene, back_wall, white_diffuse_spd, white_glossy_spd);
	add_plane_to_scene(&scene, left_wall, red_diffuse_spd, red_glossy_spd);
	add_plane_to_scene(&scene, right_wall, green_diffuse_spd, green_glossy_spd);
	add_plane_to_scene(&scene, floor, white_diffuse_spd, white_glossy_spd);
	add_plane_to_scene(&scene, ceiling, white_diffuse_spd, white_glossy_spd);

	RGB8 clear_colour = {};
	clear_render_buffer(&__window_back_buffer__, clear_colour);
	//raytrace_scene(&__window_back_buffer__, 90.0, 0.1, light_spd, mat_spd, specular_spd);
	raytrace_scene(&__window_back_buffer__, 90.0, 0.1, &scene);
	while(running)
	{
		MSG message;
		while(PeekMessage(&message, 0, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&message);
			DispatchMessage(&message);
		}



		update_window_front_buffer(window, &__window_back_buffer__);
	}

	return 0;
}
