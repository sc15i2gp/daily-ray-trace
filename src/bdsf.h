struct Surface_Point; //Forward decl
typedef void (*REFLECTION_MODEL_FUNCTION)(Surface_Point&, Vec3, Vec3, Spectrum&, Spectrum&);
typedef double (*DISTRIBUTION_FUNCTION)(Surface_Point&, Vec3, Vec3);
typedef Vec3 (*INDIRECT_SAMPLE_FUNCTION)(Surface_Point&, Vec3, double*);

enum BDSF_TYPE
{
	BDSF_TYPE_DIFFUSE = 1 << 0,
	BDSF_TYPE_SPECULAR = 1 << 1
};

struct BDSF
{
	BDSF_TYPE type;
	//Computes proportion of light reflected at given point with given incident and reflection directions
	REFLECTION_MODEL_FUNCTION bdsf;
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

#define MAT_BDSF_MAX 8
struct Material
{
	Material_Type type;
	int number_of_bdsfs;
	BDSF bdsfs[MAT_BDSF_MAX];

	bool is_emissive;

	//Textures
	Texture emission_spd_texture; //Spectral, can vary over surface
	Texture diffuse_spd_texture; //Spectral, can vary over surface
	Texture glossy_spd_texture; //Spectral, can vary over surface
	Texture refract_index_texture; //v
	Texture extinct_index_texture; //Is spectral but isn't treated as such
	Texture shininess_texture; //Same as roughness
	Texture roughness_texture; //Never spectral, value is a statistical measure of roughness
};

struct Surface_Point
{
	char* name;
	bool exists;
	Vec3 normal;
	Vec3 position;
	Vec2 texture_coordinates;
	Material* surface_material;
	Material* incident_material;
	Material* transmit_material;
};

void bdsf(Surface_Point&, Vec3 incoming, Vec3 outgoing, Spectrum&);
Material create_plastic(Spectrum diffuse_spd, Spectrum glossy_spd, double shininess);
Material create_textured_plastic(Texture diffuse_spd_texture, Spectrum glossy_spd, double shininess);
Material create_mirror();
Material create_conductor(Spectrum refract_index, Spectrum extinct_index, double roughness);
Material create_dielectric(Spectrum refract_index);
