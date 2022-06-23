#include <float.h>

typedef Spectrum Radiance;

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
		Sphere sphere;
		Plane plane;
		Model model;
	};
};

#define SCENE_OBJECT_MAX 32

struct Scene_Object
{
	Scene_Geometry geometry;
	Material material;
	char* name;
};

struct Scene
{
	Scene_Object objects[SCENE_OBJECT_MAX];
	Material air_material;
	int number_of_objects;
};

struct Geometry_Intersection_Point
{
	//General intersection data
	int scene_object;
	
	//Model geometry intersection data
	int model_triangle_index;
	Vec3 position;
	Vec3 barycentric_coordinates;
};

struct Camera
{
	//Program input
	//Needed by:
	//	- Render image
	//Orientation
	Vec3 forward;
	Vec3 right;
	Vec3 up;
	
	//Lens
	//Needed by:
	//	- Render image
	//	- Sample direction
	Vec3 pinhole_position;
	double focal_depth;
	double focal_length;
	double aperture_radius;
	
	//Computed from program input
	//Film
	//Needed by:
	//	- Render image
	Vec3 film_top_left;
	double film_pixel_width; //Width of a pixel on the image plane
	double film_pixel_height; //Height of a pixel on the image plane
	int film_width_px; //Width of film in pixels
	int film_height_px; //Height of film in pixels
};


void render_image(Texture* render_target, int number_of_samples);

struct NEW_Spectrum_Buffer
{
	//NOTE: 32 chosen because the max number of spectra needed rn is 22
	Spectrum spectra[32]; 
};

struct NEW_Camera
{
	Vec3 pinhole_position;
	Vec3 film_normal;
	double focal_depth;
	double focal_length;
	double aperture_radius;
};

struct NEW_Scene_Geometry
{
	Geometry_Type type;
	union
	{
		Sphere* sphere;
		Plane* plane;
		Model* model;
	};
};

struct NEW_Scene_Material
{
	bool is_emissive;
	double roughness;
	double shininess;
	Spectrum refract_index_spd;
	Spectrum extinct_index_spd;
	Spectrum glossy_spd;
	Spectrum diffuse_spd;
	Spectrum emission_spd;
	int number_of_bdsfs;
	NEW_BDSF bdsfs[16];
};

struct NEW_Scene
{
	NEW_Camera camera;
	Vec3 camera_position;
	Vec3 camera_up;
	Vec3 camera_right;
	double camera_fov;

	//Geometry buffers
	Sphere spheres[16];
	Plane planes[16];
	//Model models[16];
	
	int number_of_geometries;
	NEW_Scene_Geometry geometries[16];
	int number_of_materials;
	NEW_Scene_Material film_material;
	NEW_Scene_Material air_material;
	NEW_Scene_Material null_material;
	NEW_Scene_Material materials[16];

	int number_of_objects;
	NEW_Scene_Geometry* object_geometries[16];
	NEW_Scene_Material* object_materials[16];
	
	int number_of_light_sources;
	NEW_Scene_Geometry* light_source_geometries[16];
	NEW_Scene_Material* light_source_materials[16];
};

struct NEW_Geometry_Sample
{
	Vec3 position;
	//Vec3 barycentric_coordinates;
	double pdf;
	int object;
};

struct NEW_Scene_Point
{
	//Material info
	double roughness;
	NEW_Scene_Material* incident_material;
	NEW_Scene_Material* transmit_material;
	NEW_Scene_Material* material;

	//Spatial info
	int number_of_light_contributions;
	Vec3 normal;
	Vec3 position;
	//Vec3 barycentric_coordinates;
	double pdf;
	int object;
	Ray in_ray; //Both rays point away from the scene point
	Ray out_ray;
	NEW_Geometry_Sample light_contributions[16];
};

void NEW_render_image(RGB64* render_target, int target_width, int target_height, int render_samples);
/*



void NEW_render_image(RGB64*, int, int, int);
*/
