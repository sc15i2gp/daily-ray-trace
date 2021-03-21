#pragma once
#include "Colour.h"
#include "Maths.h"

struct Surface_Point; //Forward decl
typedef Spectrum (*REFLECTION_MODEL_FUNCTION)(Surface_Point, Vec3, Vec3);
typedef double (*DISTRIBUTION_FUNCTION)(Surface_Point, Vec3, Vec3);
typedef Vec3 (*INDIRECT_SAMPLE_FUNCTION)(Surface_Point, Vec3, double*);

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

	//Plastic data
	Spectrum diffuse_spd;
	Spectrum glossy_spd;
	double shininess;

	//Fresnel data
	Spectrum refract_index;
	Spectrum extinct_index;
	double roughness;
};

struct Surface_Point
{
	char* name;
	Spectrum diffuse_spd;
	Spectrum glossy_spd;
	Spectrum emission_spd;
	Spectrum incident_refract_index;
	Spectrum transmit_refract_index;
	Vec3 normal;
	Vec3 position;
	Material material;
	bool exists;
	bool is_emissive;
};

Spectrum bsdf(Surface_Point, Vec3 incoming, Vec3 outgoing);
Material create_plastic(Spectrum diffuse_spd, Spectrum glossy_spd, double shininess);
Material create_mirror();
Material create_conductor(Spectrum refract_index, Spectrum extinct_index, double roughness);
Material create_dielectric(Spectrum refract_index);
