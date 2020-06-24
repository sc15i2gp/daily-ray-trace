#pragma once
#include <immintrin.h>
#include <stdint.h>
#include <math.h>

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
double dot(Vec4, Vec4);
double length(Vec4);

struct Ray
{
	Vec3 origin;
	Vec3 direction;
};
