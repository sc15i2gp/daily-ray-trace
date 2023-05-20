enum Geometry_Type
{
	GEO_TYPE_NONE = 0,
	GEO_TYPE_POINT,
	GEO_TYPE_SPHERE,
	GEO_TYPE_PLANE,
	GEO_TYPE_MODEL
};

struct Camera
{
	Vec3 pinhole_position;
	Vec3 film_normal;
	double focal_depth;
	double focal_length;
	double aperture_radius;
};

struct Scene_Geometry
{
	Geometry_Type type;
	union
	{
		Sphere* sphere;
		Plane* plane;
		//Model* model;
	};
};

//Properties of a surface at all points
struct Scene_Material
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
	BDSF bdsfs[16];
};

struct Scene
{
	Camera camera;
	Vec3 camera_position;
	Vec3 camera_up;
	Vec3 camera_right;
	double camera_fov;

	//Geometry buffers
	Sphere spheres[16];
	Plane planes[16];
	//Model models[16];
	
	int number_of_geometries;
	Scene_Geometry geometries[16];
	int number_of_materials;
	Scene_Material film_material;
	Scene_Material air_material;
	Scene_Material null_material;
	Scene_Material materials[16];

	int number_of_objects;
	Scene_Geometry* object_geometries[16];
	Scene_Material* object_materials[16];
	
	int number_of_light_sources;
	Scene_Geometry* light_source_geometries[16];
	Scene_Material* light_source_materials[16];
};

struct Light_Contribution
{
	Spectrum* emission_spd;
	Scene_Geometry* geometry;
	Vec3 position; //On the light's surface
};

struct Scene_Path_Point
{
	//Surface point info
	Surface_Point properties; //Properties of a surface at a particular point

	//Scene point info
	Vec3 in_direction; //Both rays point away from the scene point
	Vec3 out_direction;
	int number_of_light_contributions;
	bool is_emissive;
	Light_Contribution light_contributions[16];
	int number_of_bdsfs;
	BDSF* bdsfs;
};

void render_image(RGB64* render_target, int target_width, int target_height, int render_samples, bool invert_final_image = false);
