#include <float.h>

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
	bool is_emissive;
};

struct Scene
{
	int number_of_objects;
	Scene_Object objects[SCENE_OBJECT_MAX];
};

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

void render_image(Texture* render_target);
void print_render_profile();
