struct Surface_Point; //Forward decl
typedef void (*REFLECTION_MODEL_FUNCTION)(Surface_Point&, Vec3, Vec3, Spectrum&, Spectrum&);
typedef double (*DISTRIBUTION_FUNCTION)(Surface_Point&, Vec3, Vec3);
typedef Vec3 (*INDIRECT_SAMPLE_FUNCTION)(Surface_Point&, Vec3, double*);

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

	//Textures
	Texture emission_spd_texture;
	Texture diffuse_spd_texture;
	Texture glossy_spd_texture;
	Texture shininess_texture;
	Texture refract_index_texture;
	Texture extinct_index_texture;
	Texture roughness_texture;
};

struct Surface_Point
{
	char* name;
	Vec3 normal;
	Vec3 position;
	Vec2 texture_coordinates;
	Material* surface_material;
	Material* incident_material;
	Material* transmit_material;
	bool exists;
	bool is_emissive;
};

void bsdf(Surface_Point&, Vec3 incoming, Vec3 outgoing, Spectrum&);
Material create_plastic(Spectrum diffuse_spd, Spectrum glossy_spd, double shininess);
Material create_textured_plastic(Texture diffuse_spd_texture, Spectrum glossy_spd, double shininess);
Material create_mirror();
Material create_conductor(Spectrum refract_index, Spectrum extinct_index, double roughness);
Material create_dielectric(Spectrum refract_index);
