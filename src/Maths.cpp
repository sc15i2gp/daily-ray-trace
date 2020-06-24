#include "Maths.h"

double sin_deg(double t)
{
	return sin(t * (PI/180.0));
}

double cos_deg(double t)
{
	return cos(t * (PI/180.0));
}

double tan_deg(double t)
{
	return tan(t * (PI/180.0));
}

/*	Vectors	*/

double& Vec2::operator[](int index)
{
	return xy[index];
}

double& Vec3::operator[](int index)
{
	return xyz[index];
}

double& Vec4::operator[](int index)
{
	return xyzw[index];
}

/*	Vec2	*/

Vec2 operator+(Vec2 v, Vec2 w)
{
	v.x += w.x;
	v.y += w.y;
	return v;
}

Vec2 operator-(Vec2 v, Vec2 w)
{
	v.x -= w.x;
	v.y -= w.y;
	return v;
}

Vec2 operator-(Vec2 v)
{
	v.x = -v.x;
	v.y = -v.y;
	return v;
}

Vec2 operator*(Vec2 v, Vec2 w)
{
	v.x *= w.x;
	v.y *= w.y;
	return v;
}

Vec2 operator*(double d, Vec2 v)
{
	v.x *= d;
	v.y *= d;
	return v;
}

Vec2 operator/(Vec2 v, Vec2 w)
{
	v.x /= w.x;
	v.y /= w.y;
	return v;
}

Vec2 operator/(Vec2 v, double d)
{
	v.x /= d;
	v.y /= d;
	return v;
}

void operator+=(Vec2& v, Vec2 w)
{
	v = v + w;
}

void operator-=(Vec2& v, Vec2 w)
{
	v = v - w;
}

void operator*=(Vec2& v, Vec2 w)
{
	v = v * w;
}

void operator*=(Vec2& v, double d)
{
	v = d * v;
}

void operator/=(Vec2& v, Vec2 w)
{
	v = v / w;
}

void operator/=(Vec2& v, double d)
{
	v = v / d;
}

double dot(Vec2 v, Vec2 w)
{
	return v.x * w.x + v.y * w.y;
}

double length(Vec2 v)
{
	return sqrt(dot(v,v));
}

Vec2 normalise(Vec2 v)
{
	return v/length(v);
}

/*	Vec3	*/

Vec3 operator+(Vec3 v, Vec3 w)
{
	v.x += w.x;
	v.y += w.y;
	v.z += w.z;
	return v;
}

Vec3 operator-(Vec3 v, Vec3 w)
{
	v.x -= w.x;
	v.y -= w.y;
	v.z -= w.z;
	return v;
}

Vec3 operator-(Vec3 v)
{
	v.x = -v.x;
	v.y = -v.y;
	v.z = -v.z;
	return v;
}

Vec3 operator*(Vec3 v, Vec3 w)
{
	v.x *= w.x;
	v.y *= w.y;
	v.z *= w.z;
	return v;
}

Vec3 operator*(double d, Vec3 v)
{
	v.x *= d;
	v.y *= d;
	v.z *= d;
	return v;
}

Vec3 operator/(Vec3 v, Vec3 w)
{
	v.x /= w.x;
	v.y /= w.y;
	v.z /= w.z;
	return v;
}

Vec3 operator/(Vec3 v, double d)
{
	v.x /= d;
	v.y /= d;
	v.z /= d;
	return v;
}

void operator+=(Vec3& v, Vec3 w)
{
	v = v + w;
}

void operator-=(Vec3& v, Vec3 w)
{
	v = v - w;
}

void operator*=(Vec3& v, Vec3 w)
{
	v = v * w;
}

void operator*=(Vec3& v, double d)
{
	v = d * v;
}

void operator/=(Vec3& v, Vec3 w)
{
	v = v / w;
}

void operator/=(Vec3& v, double d)
{
	v = v / d;
}

Vec3 cross(Vec3 v, Vec3 w)
{
	Vec3 u = {};
	u.x = v.y*w.z - v.z*w.y;
	u.y = v.z*w.x - v.x*w.z;
	u.z = v.x*w.y - v.y*w.x;
	return u;
}

double dot(Vec3 v, Vec3 w)
{
	return v.x*w.x + v.y*w.y + v.z*w.z;
}

double length(Vec3 v)
{
	return sqrt(dot(v, v));
}

Vec3 normalise(Vec3 v)
{
	return v/length(v);
}

/*	Vec4	*/

Vec4 operator+(Vec4 v, Vec4 w)
{
	v.x += w.x;
	v.y += w.y;
	v.z += w.z;
	v.w += w.w;
	return v;
}

Vec4 operator-(Vec4 v, Vec4 w)
{
	v.x -= w.x;
	v.y -= w.y;
	v.z -= w.z;
	v.w -= w.w;
	return v;
}

Vec4 operator-(Vec4 v)
{
	v.x = -v.x;
	v.y = -v.y;
	v.z = -v.z;
	v.w = -v.w;
	return v;
}

Vec4 operator*(Vec4 v, Vec4 w)
{
	v.x *= w.x;
	v.y *= w.y;
	v.z *= w.z;
	v.w *= w.w;
	return v;
}

Vec4 operator*(double d, Vec4 v)
{
	v.x *= d;
	v.y *= d;
	v.z *= d;
	v.w *= d;
	return v;
}

Vec4 operator/(Vec4 v, Vec4 w)
{
	v.x /= w.x;
	v.y /= w.y;
	v.z /= w.z;
	v.w /= w.w;
	return v;
}

Vec4 operator/(Vec4 v, double d)
{
	v.x /= d;
	v.y /= d;
	v.z /= d;
	v.w /= d;
	return v;
}

void operator+=(Vec4& v, Vec4 w)
{
	v = v + w;
}

void operator-=(Vec4& v, Vec4 w)
{
	v = v - w;
}

void operator*=(Vec4& v, Vec4 w)
{
	v = v * w;
}

void operator*=(Vec4& v, double d)
{
	v = d * v;
}

void operator/=(Vec4& v, Vec4 w)
{
	v = v / w;
}

void operator/=(Vec4& v, double d)
{
	v = v / d;
}


Vec4 cross(Vec4 v, Vec4 w)
{
	v.xyz = cross(v.xyz, w.xyz);
	return v;
}

double dot(Vec4 v, Vec4 w)
{
	return dot(v.xyz, w.xyz) + v.w*w.w;
}

double length(Vec4 v)
{
	return sqrt(dot(v, v));
}

Vec4 normalise(Vec4 v)
{
	return v/length(v);
}

/*	Geometry	*/

//Returns true if ray intersects sphere, false otherwise
//t is used to return smallest solution to quadratic equation
bool ray_intersects_sphere(Ray ray, Vec3 sphere_center, double sphere_radius, double* t)
{
	//Sphere: Radius^2 = dot(P-Center, P-Center), where P is point on surface of sphere
	//Ray: R = R_origin + t(R_direction)
	//Derived quadratic equation from two equations above
	
	//a, b, c coefficients in quadratic, solved using quadratic formula
	Vec3 s_center_to_r_origin = ray.origin - sphere_center;
	double a = 1.0;
	double b = 2.0 * dot(ray.direction, s_center_to_r_origin);
	double c = dot(s_center_to_r_origin, s_center_to_r_origin) - sphere_radius*sphere_radius;

	double discriminant = b*b - 4.0 * a * c;

	if(discriminant >= 0.0)
	{//If a solution exists
		double desired_solution = 0.0;
		if(discriminant == 0.0) desired_solution = -b/(2.0*a);
		else
		{//2 solutions exist
			desired_solution = -(b + sqrt(discriminant))/2.0*a;
		}

		*t = desired_solution;
		return desired_solution >= 0.0;
	}
	else return false;
}

//Returns true if ray intersects plane, false otherwise
//t is used to return solution to linear equation
//Plane components consist of: Point on plane p, plane normal n, boundary vectors u and v
//Plane consists of points q which satisfy: 
//	- dot((p - q), n) = 0
//	- 0 <= dot((p - q), u) <= 1
//	- 0 <= dot((p - q), v) <= 1
bool ray_intersects_plane(Ray ray, Vec3 p, Vec3 n, Vec3 u, Vec3 v, double* t)
{
	if(dot(ray.direction, n) != 0.0)
	{
		double desired_solution = dot(p - ray.origin, n)/dot(ray.direction, n);
		*t = desired_solution;

		Vec3 intersection = ray.origin + desired_solution * ray.direction;
		Vec3 u_norm = normalise(u);
		Vec3 v_norm = normalise(v);
		double i_u = dot(intersection - p, u_norm); //How far along u intersection is
		double i_v = dot(intersection - p, v_norm); //How far along v intersection is
		
		bool within_u = i_u >= 0.0 && i_u <= length(u);
		bool within_v = i_v >= 0.0 && i_v <= length(v);
		bool intersection_within_bounds = within_u && within_v;
		return (desired_solution >= 0.0) && intersection_within_bounds;
	}
	else return false;

}
