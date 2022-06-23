
enum Material_Type
{
	MAT_TYPE_NONE = 0,
	MAT_TYPE_DIELECTRIC,
	MAT_TYPE_CONDUCTOR
};

struct Surface_Point; //Forward decl
typedef void (*REFLECTION_MODEL_FUNCTION)(Surface_Point&, Vec3, Vec3, Spectrum&, Spectrum&);
typedef Vec3 (*INDIRECT_SAMPLE_FUNCTION)(Surface_Point&, Vec3, double*);
typedef double (*DISTRIBUTION_FUNCTION)(Surface_Point&, Vec3, Vec3);

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
	Texture refract_index_texture; //v, doesn't vary over surface
	Texture extinct_index_texture; //Is spectral but isn't treated as such, doesn't vary over surface
	Texture shininess_texture; //Same as roughness
	Texture roughness_texture; //Never spectral, value is a statistical measure of roughness, can vary over surface
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

struct NEW_Surface_Point;
typedef void (*NEW_REFLECTANCE_FUNCTION)(NEW_Surface_Point*, Vec3, Vec3, Spectrum*);
typedef Vec3 (*NEW_DIRECTION_SAMPLE_FUNCTION)(NEW_Surface_Point*, Vec3, double* pdf);

struct NEW_BDSF
{
	NEW_REFLECTANCE_FUNCTION reflectance;
	NEW_DIRECTION_SAMPLE_FUNCTION sample_direction;
};

struct NEW_Camera; //Forward decl
struct NEW_Surface_Point
{
	Vec3 position;
	Vec3 normal;
	double shininess;
	double roughness;
	NEW_Camera* camera;
	Spectrum* incident_refract_index_spd;
	Spectrum* incident_extinct_index_spd;
	Spectrum* transmit_refract_index_spd;
	Spectrum* transmit_extinct_index_spd;
	Spectrum* glossy_spd;
	Spectrum* diffuse_spd;
};

Vec3 NEW_sample_camera_lens_direction(NEW_Surface_Point*, Vec3, double*);
Vec3 NEW_sample_camera_pinhole_direction(NEW_Surface_Point*, Vec3, double*);
Vec3 NEW_sample_diffuse_direction(NEW_Surface_Point*, Vec3, double*);
Vec3 NEW_sample_glossy_direction(NEW_Surface_Point*, Vec3, double*);

void NEW_const_1_reflectance(NEW_Surface_Point*, Vec3, Vec3, Spectrum* fr);
void NEW_diffuse_phong_reflectance(NEW_Surface_Point*, Vec3, Vec3, Spectrum* fr);
void NEW_glossy_phong_reflectance(NEW_Surface_Point*, Vec3, Vec3, Spectrum* fr);
