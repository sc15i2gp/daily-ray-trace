#include "Maths.h"

double uniform_sample()
{
	return (double)(rand() / (RAND_MAX + 1.0));
}

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

double d_max(double a, double b)
{
	return (a < b) ? b : a;
}

double d_min(double a, double b)
{
	return (a < b) ? a : b;
}

double clamp(double d, double low, double high)
{
	if(d < low) return low;
	else if(d > high) return high;
	else return d;
}

double lerp(double x, double x_0, double x_1, double y_0, double y_1)
{
	return y_0 + (x - x_0) * (y_1 - y_0) / (x_1 - x_0);
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

bool operator==(Vec2 v, Vec2 w)
{
	return (v.x == w.x) && (v.y == w.y);
}

bool operator!=(Vec2 v, Vec2 w)
{
	return !(v == w);
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

bool operator==(Vec3 v, Vec3 w)
{
	return (v.x == w.x) && (v.y == w.y) && (v.z == w.z);
}

bool operator!=(Vec3 v, Vec3 w)
{
	return !(v == w);
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

bool operator==(Vec4 v, Vec4 w)
{
	return (v.x == w.x) && (v.y == w.y) && (v.z == w.z) && (v.w == w.w);
}

bool operator!=(Vec4 v, Vec4 w)
{
	return !(v == w);
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

/*	Matrices	*/

Vec3 Mat3x3::row(int i)
{
	Vec3 r = {};
	for(int j = 0; j < 3; ++j) r[j] = this->columns[j][i];
	return r;
}

Vec3 Mat3x3::column(int i)
{
	return (*this)[i];
}

Vec3& Mat3x3::operator[](int i)
{
	return this->columns[i];
}

Mat3x3 operator+(Mat3x3 m, Mat3x3 n)
{
	for(int i = 0; i < 3; ++i) m[i] += n[i];
	return m;
}

Mat3x3 operator-(Mat3x3 m)
{
	for(int i = 0; i < 3; ++i)
	{
		m.columns[i] = -m.columns[i];
	}
	return m;
}

Vec3 operator*(Mat3x3 m, Vec3 v)
{
	Vec3 w = {};
	for(int i = 0; i < 3; ++i) w[i] = dot(m.row(i), v);
	return w;
}

Mat3x3 operator*(Mat3x3 m, Mat3x3 n)
{
	Mat3x3 p = {};
	for(int i = 0; i < 3; ++i)
	{
		for(int j = 0; j < 3; ++j)
		{
			p[i][j] = dot(m.row(i), n.column(j));
		}
	}
	return p;
}

Mat3x3 operator*(double d, Mat3x3 m)
{
	for(int i = 0; i < 3; ++i)
	{
		m[i] *= d;
	}
	return m;
}

Mat3x3 identity3x3()
{
	Mat3x3 i = {};
	for(int j = 0; j < 3; ++j)  i[j][j] = 1.0;
	return i;
}

Mat3x3 find_rotation_between_vectors(Vec3 v, Vec3 w)
{
	Vec3 n = cross(v, w);
	double s = length(n);
	double c = dot(v, w);

	Mat3x3 r;
	if(dot(n, n) == 0)
	{
		if(c > 0) r = identity3x3();
		else r = -identity3x3();
	}
	else
	{
		Mat3x3 m = {};
		m[0][1] = n.z;
		m[0][2] = -n.y;
		m[1][0] = -n.z;
		m[1][2] = n.x;
		m[2][0] = n.y;
		m[2][1] = -n.x;

		r = identity3x3() + m + (1.0/(1.0 + c))*m*m;
	}
	return r;
}

Vec3 reflect_vector(Vec3 v, Vec3 n)
{
	return v - 2.0 * dot(v, n) * n;
}

/*	Geometry	*/

Plane create_plane_from_bounds(Vec3 p, Vec3 u, Vec3 v)
{
	Plane plane = {};
	plane.p = p;
	plane.u = u;
	plane.v = v;
	plane.n = normalise(cross(plane.u, plane.v));

	return plane;
}

//TODO: Make the arguments work ccw instead of cw
Plane create_plane_from_points(Vec3 p, Vec3 u, Vec3 v)
{
	Plane plane = {};
	plane.p = p;
	plane.u = u - p;
	plane.v = v - p;
	plane.n = -normalise(cross(plane.u, plane.v));

	return plane;
}

Vec3 normal(Sphere s, Vec3 p)
{
	return p - s.center;
}

Vec3 normal(Plane p)
{
	return p.n;
}

double area(Sphere s)
{
	return 4.0 * PI * s.radius * s.radius;
}

double area(Plane p)
{
	return length(cross(p.u, p.v));
}

Vec3 uniform_sample_sphere()
{
	double u = 1.0 - 2.0 * uniform_sample();
	double v = uniform_sample();

	double r = sqrt(1.0 - u * u);
	double t = 2.0 * PI * v;

	Vec3 sample = {r * cos(t), r * sin(t), u};
	return sample;
}

Vec3 uniform_sample_sphere_subtended(Sphere s, Vec3 p, double* pdf)
{
	double q = uniform_sample();
	double r = uniform_sample();

	double st_max_sq = s.radius*s.radius / dot(s.center - p, s.center - p);
	double ct_max = sqrt(d_max(0.0, 1.0 - st_max_sq));
	double phi = 2.0 * q * PI;
	double ct = 1.0 - r + r*ct_max;
	double st = sqrt(d_max(0.0, 1.0 - ct * ct));

	double d = length(s.center - p);
	double dc = sqrt(d_max(0.0, s.radius * s.radius - d * d * st * st));
	double ds = d * ct - dc;

	double ca = (d * d + s.radius * s.radius - ds * ds)/(2.0 * d * s.radius);
	double sa = sqrt(d_max(0.0, 1.0 - ca * ca));

	double R = s.radius * sa;

	Vec3 w = {R * cos(phi), R * sin(phi), s.radius * ca};
	Mat3x3 m = find_rotation_between_vectors(Vec3{0.0, 0.0, 1.0}, normalise(s.center - p));
	
	*pdf = 1.0 / (2.0 * PI * (1.0 - ct_max));

	return m * (-w) + s.center;
}

//TODO: Fix this method when I've sorted out handling pdf of 0 in other functions
//	this function needs to do an intersection test to return 0 pdf
//NOTE: Returns pdf value as value over solid angle, not surface area
Vec3 uniform_sample_sphere(Sphere s, Vec3 p, double* pdf)
{
	Vec3 s_p = s.center + s.radius * uniform_sample_sphere();
	Vec3 s_n = normal(s, s_p);
	Vec3 s_p_to_p = p - s_p;
	*pdf = (1.0 / area(s)) * (dot(p - s_p, p - s_p)/abs(dot(normal(s, s_p), p - s_p)));
	return s.center + s.radius * uniform_sample_sphere();
}

Vec3 uniform_sample_hemisphere(Vec3 normal, double* pdf)
{
	for(;;)
	{
		*pdf = 1.0 / (2.0 * PI);
		Vec3 sample = uniform_sample_sphere();
		if(dot(sample, normal) >= 0.0) return sample;
	}
}

Vec3 uniform_sample_plane(Plane p, double* pdf)
{
	double u = uniform_sample();
	double v = uniform_sample();

	*pdf = 1.0 / area(p);
	return p.p + u * p.u + v * p.v;
}

//Returns random point within unit disc with normal (0.0, 0.0, 1.0)
Vec3 uniform_sample_disc()
{
	Vec2 rand_vec = {uniform_sample(), uniform_sample()};
	Vec2 offset = 2.0 * rand_vec - Vec2{1.0, 1.0};
	if(offset.x == 0.0 && offset.y == 0.0) return Vec3{};
	double r = 0.0;
	double t = 0.0;
	if(abs(offset.x) > abs(offset.y))
	{
		r = offset.x;
		t = (PI / 4.0) * (offset.y / offset.x);
	}
	else
	{
		r = offset.y;
		t = (PI / 2.0) - (PI / 4.0) * (offset.x / offset.y);
	}
	return r * Vec3{cos(t), sin(t), 0.0};
}

double cos_weighted_sample_hemisphere_pdf(Vec3 normal, Vec3 v)
{
	return dot(normal, v) / PI;
}

Vec3 cos_weighted_sample_hemisphere(Vec3 normal, double* pdf)
{
	Vec3 p = {};
	for(;;)
	{
		p = uniform_sample_disc();
		if(dot(p, p) < 1.0) break;
	}
	p.z = sqrt(1.0 - dot(p, p));
	//Rotate p by R such that R*(0, 0, 1) = normal
	Mat3x3 r = find_rotation_between_vectors(Vec3{0.0, 0.0, 1.0}, normal);
	Vec3 v = r * p;
	*pdf = cos_weighted_sample_hemisphere_pdf(normal, v);
	return v;
}

//Returns true if ray intersects sphere, false otherwise
//t is used to return smallest solution to quadratic equation
bool ray_intersects_sphere(Ray ray, Sphere s, double* t)
{
	//Sphere: Radius^2 = dot(P-Center, P-Center), where P is point on surface of sphere
	//Ray: R = R_origin + t(R_direction)
	//Derived quadratic equation from two equations above
	
	//a, b, c coefficients in quadratic, solved using quadratic formula
	Vec3 s_center_to_r_origin = ray.origin - s.center;
	double a = 1.0;
	double b = 2.0 * dot(ray.direction, s_center_to_r_origin);
	double c = dot(s_center_to_r_origin, s_center_to_r_origin) - s.radius*s.radius;

	double discriminant = b*b - 4.0 * a * c;

	if(discriminant >= 0.0)
	{//If a solution exists
		double desired_solution = 0.0;
		if(discriminant == 0.0) desired_solution = -b/(2.0*a);
		else
		{//2 solutions exist
			double solution_0 = -(b + sqrt(discriminant))/2.0*a;
			double solution_1 = -(b - sqrt(discriminant))/2.0*a;
			//Choose smallest positive solution
			if(solution_1 < 0.0) desired_solution = solution_0;
			else if(solution_0 < 0.0) desired_solution = solution_1;
			else desired_solution = (solution_0 < solution_1) ? solution_0 : solution_1;
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
bool ray_intersects_plane(Ray ray, Plane p, double* t)
{
	if(dot(ray.direction, p.n) != 0.0)
	{
		double desired_solution = dot(p.p - ray.origin, p.n)/dot(ray.direction, p.n);
		*t = desired_solution;

		Vec3 intersection = ray.origin + desired_solution * ray.direction;
		Vec3 u_norm = normalise(p.u);
		Vec3 v_norm = normalise(p.v);
		double i_u = dot(intersection - p.p, u_norm); //How far along u intersection is
		double i_v = dot(intersection - p.p, v_norm); //How far along v intersection is
		
		bool within_u = i_u >= 0.0 && i_u <= length(p.u);
		bool within_v = i_v >= 0.0 && i_v <= length(p.v);
		bool intersection_within_bounds = within_u && within_v;
		return (desired_solution >= 0.0) && intersection_within_bounds;
	}
	else return false;

}

double point_to_line_distance_sq(Vec3 p, Vec3 l_0, Vec3 l_1)
{
	Vec3 normalised_line = normalise(l_0 - l_1);
	Vec3 q = l_1 + (dot(p - l_1, normalised_line))*normalised_line;
	double distance = dot(p - q, p - q);
	return distance;
}
