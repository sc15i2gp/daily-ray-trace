#pragma once
#include <immintrin.h>
#include <stdint.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

double uniform_sample();

#define PI 3.141592653L

double sin_deg(double t);
double cos_deg(double t);
double tan_deg(double t);

double d_max(double, double);
double d_min(double, double);

double clamp(double, double low, double high);

double lerp(double x, double x_0, double x_1, double y_0, double y_1);

/*	Vectors	*/

struct Vec2
{
	union
	{
		double xy[2];
		struct
		{
			double x;
			double y;
		};
	};
	double& operator[](int);
};

struct Vec3
{
	union
	{
		double xyz[3];
		struct
		{
			double x;
			double y;
			double z;
		};
		struct
		{
			double R;
			double G;
			double B;
		};
		struct
		{
			Vec2 xy;
			double _z;
		};
		struct
		{
			double _x;
			Vec2 yz;
		};
	};
	double& operator[](int);
};

struct Vec4
{
	union
	{
		double xyzw[4];
		struct
		{
			double x;
			double y;
			double z;
			double w;
		};
		struct
		{
			Vec3 xyz;
			double _w;
		};
		struct
		{
			double _x;
			Vec3 yzw;
		};
		struct
		{
			Vec2 xy;
			Vec2 zw;
		};
	};
	double& operator[](int);
};

//Operators:
//	- Standard: +, -, *, /, +=, -=, *=, /=
//	- Magnitude
//	- Dot, Cross


/*	Vec2	*/

Vec2 operator+(Vec2, Vec2);
Vec2 operator-(Vec2, Vec2);
Vec2 operator-(Vec2);
Vec2 operator*(Vec2, Vec2);
Vec2 operator*(double, Vec2);
Vec2 operator/(Vec2, Vec2);
Vec2 operator/(Vec2, double);
void operator+=(Vec2&, Vec2);
void operator-=(Vec2&, Vec2);
void operator*=(Vec2&, Vec2);
void operator*=(Vec2&, double);
void operator/=(Vec2&, Vec2);
void operator/=(Vec2&, double);

bool operator==(Vec2, Vec2);
bool operator!=(Vec2, Vec2);

double dot(Vec2, Vec2);
double length(Vec2);
Vec2 normalise(Vec2);

/*	Vec3	*/

Vec3 operator+(Vec3, Vec3);
Vec3 operator-(Vec3, Vec3);
Vec3 operator-(Vec3);
Vec3 operator*(Vec3, Vec3);
Vec3 operator*(double, Vec3);
Vec3 operator/(Vec3, Vec3);
Vec3 operator/(Vec3, double);
void operator+=(Vec3&, Vec3);
void operator-=(Vec3&, Vec3);
void operator*=(Vec3&, Vec3);
void operator*=(Vec3&, double);
void operator/=(Vec3&, Vec3);
void operator/=(Vec3&, double);

bool operator==(Vec3, Vec3);
bool operator!=(Vec3, Vec3);

Vec3 cross(Vec3, Vec3);
Vec3 normalise(Vec3);
double dot(Vec3, Vec3);
double length(Vec3);

/*	Vec4	*/

Vec4 operator+(Vec4, Vec4);
Vec4 operator-(Vec4, Vec4);
Vec4 operator-(Vec4);
Vec4 operator*(Vec4, Vec4);
Vec4 operator*(double, Vec4);
Vec4 operator/(Vec4, Vec4);
Vec4 operator/(Vec4, double);
void operator+=(Vec4&, Vec4);
void operator-=(Vec4&, Vec4);
void operator*=(Vec4&, Vec4);
void operator*=(Vec4&, double);
void operator/=(Vec4&, Vec4);
void operator/=(Vec4&, double);

bool operator==(Vec4, Vec4);
bool operator!=(Vec4, Vec4);

Vec4 cross(Vec4, Vec4);
Vec4 normalise(Vec4);
double dot(Vec4, Vec4);
double length(Vec4);

/*	Matrices	*/

struct Mat3x3
{
	Vec3 columns[3];

	Vec3 row(int);
	Vec3 column(int);
	Vec3& operator[](int);
};

Mat3x3 operator+(Mat3x3, Mat3x3);
Vec3 operator*(Mat3x3, Vec3);
Mat3x3 operator*(Mat3x3, Mat3x3);
Mat3x3 operator*(double, Mat3x3);
Mat3x3 identity3x3();

//TODO: Rotation from v to w?
Mat3x3 find_rotation_between_vectors(Vec3 v, Vec3 w);

Vec3 reflect_vector(Vec3 v, Vec3 n); //Reflects v through n

/*	Geometry	*/

struct Ray
{
	Vec3 origin;
	Vec3 direction;
};

struct Sphere
{
	Vec3 center;
	double radius;
};

double area(Sphere);
Vec3 uniform_sample_sphere(Sphere, Vec3 p, double* pdf_value);
Vec3 uniform_sample_sphere_subtended(Sphere,Vec3 p, double* pdf_value); //Samples from points on sphere visible from p
Vec3 uniform_sample_hemisphere(Vec3 normal, double* pdf_value);
Vec3 uniform_sample_disc();
double cos_weighted_sample_hemisphere_pdf(Vec3 normal, Vec3 v);
Vec3 cos_weighted_sample_hemisphere(Vec3 normal, double* pdf_value);

struct Plane
{
	Vec3 p; //Corner point of plane
	Vec3 n; //Normal
	//Boundary vectors are not normalised so that their lengths can be used as boundary lengths
	Vec3 u; //First boundary vector
	Vec3 v; //Second boundary vector
};

Plane create_plane_from_bounds(Vec3 p, Vec3 u, Vec3 v);
Plane create_plane_from_points(Vec3 p, Vec3 u, Vec3 v);
double area(Plane);
Vec3 uniform_sample_plane(Plane, double* pdf_value);

bool ray_intersects_sphere(Ray, Sphere, double* t);
bool ray_intersects_plane(Ray ray, Plane, double* t);
