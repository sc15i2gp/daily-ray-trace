
struct Surface_Point; //Forward decl
typedef void (*REFLECTANCE_FUNCTION)(Surface_Point*, Vec3, Vec3, Spectrum*);
typedef Vec3 (*DIRECTION_SAMPLE_FUNCTION)(Surface_Point*, Vec3);
typedef double (*PDF_VALUE_FUNCTION)(Surface_Point*, Vec3, Vec3);

struct BDSF
{
	REFLECTANCE_FUNCTION reflectance;
	DIRECTION_SAMPLE_FUNCTION sample_direction;
	PDF_VALUE_FUNCTION pdf;
};

struct Camera; //Forward decl
//Properties of a surface at a point
struct Surface_Point
{
	Vec3 position;
	Vec3 normal;
	double shininess;
	double roughness;
	Camera* camera;
	Spectrum* incident_refract_index_spd;
	Spectrum* incident_extinct_index_spd;
	Spectrum* transmit_refract_index_spd;
	Spectrum* transmit_extinct_index_spd;
	Spectrum* glossy_spd;
	Spectrum* diffuse_spd;
};

Vec3 sample_camera_lens_direction(Surface_Point*, Vec3);
Vec3 sample_camera_pinhole_direction(Surface_Point*, Vec3);
Vec3 sample_diffuse_direction(Surface_Point*, Vec3);
Vec3 sample_glossy_direction(Surface_Point*, Vec3);
Vec3 sample_specular_reflection_direction(Surface_Point*, Vec3);
Vec3 sample_cook_torrance_reflection_direction(Surface_Point*, Vec3);

void const_1_reflectance(Surface_Point*, Vec3, Vec3, Spectrum* fr);
void diffuse_phong_reflectance(Surface_Point*, Vec3, Vec3, Spectrum* fr);
void glossy_phong_reflectance(Surface_Point*, Vec3, Vec3, Spectrum* fr);
void mirror_reflectance(Surface_Point*, Vec3, Vec3, Spectrum* fr);

double const_1_pdf(Surface_Point*, Vec3, Vec3);
double cos_weighted_hemisphere_pdf(Surface_Point*, Vec3, Vec3);
double cook_torrance_reflection_pdf(Surface_Point*, Vec3, Vec3);
