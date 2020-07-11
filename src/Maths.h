#pragma once
#include <immintrin.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

double uniform_sample();

#define PI 3.141592653L

double sin_deg(double t);
double cos_deg(double t);
double tan_deg(double t);

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

Vec4 cross(Vec4, Vec4);
Vec4 normalise(Vec4);
double dot(Vec4, Vec4);
double length(Vec4);

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
Vec3 uniform_sample_sphere(Sphere);
Vec3 uniform_sample_hemisphere(Vec3 normal);

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
Vec3 uniform_sample_plane(Plane);

bool ray_intersects_sphere(Ray, Sphere, double* t);
bool ray_intersects_plane(Ray ray, Plane, double* t);
