
struct NEW_Surface_Point; //Forward decl
typedef void (*NEW_REFLECTANCE_FUNCTION)(NEW_Surface_Point*, Vec3, Vec3, Spectrum*);
typedef Vec3 (*NEW_DIRECTION_SAMPLE_FUNCTION)(NEW_Surface_Point*, Vec3);
typedef double (*NEW_PDF_VALUE_FUNCTION)(NEW_Surface_Point*, Vec3, Vec3);

struct NEW_BDSF
{
	NEW_REFLECTANCE_FUNCTION reflectance;
	NEW_DIRECTION_SAMPLE_FUNCTION sample_direction;
	NEW_PDF_VALUE_FUNCTION pdf;
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

Vec3 NEW_sample_camera_lens_direction(NEW_Surface_Point*, Vec3);
Vec3 NEW_sample_camera_pinhole_direction(NEW_Surface_Point*, Vec3);
Vec3 NEW_sample_diffuse_direction(NEW_Surface_Point*, Vec3);
Vec3 NEW_sample_glossy_direction(NEW_Surface_Point*, Vec3);

void NEW_const_1_reflectance(NEW_Surface_Point*, Vec3, Vec3, Spectrum* fr);
void NEW_diffuse_phong_reflectance(NEW_Surface_Point*, Vec3, Vec3, Spectrum* fr);
void NEW_glossy_phong_reflectance(NEW_Surface_Point*, Vec3, Vec3, Spectrum* fr);

double NEW_const_1_pdf(NEW_Surface_Point*, Vec3, Vec3);
double NEW_cos_weighted_hemisphere_pdf(NEW_Surface_Point*, Vec3, Vec3);
