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

void render_image(Texture* render_target, int number_of_samples);
