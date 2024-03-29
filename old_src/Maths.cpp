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

double point_to_line_distance_sq(Vec3 p, Vec3 l_0, Vec3 l_1)
{
	Vec3 normalised_line = normalise(l_0 - l_1);
	Vec3 q = l_1 + (dot(p - l_1, normalised_line))*normalised_line;
	double distance = dot(p - q, p - q);
	return distance;
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
	plane.u_length = length(u);
	plane.v_length = length(v);
	plane.u = u/plane.u_length;
	plane.v = v/plane.v_length;
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
	plane.u_length = length(plane.u);
	plane.v_length = length(plane.v);
	plane.u /= plane.u_length;
	plane.v /= plane.v_length;
	plane.n = -normalise(cross(plane.u, plane.v));

	return plane;
}

double ray_intersects_sphere(Vec3 o, Vec3 v, Sphere s)
{
	//Sphere: Radius^2 = dot(P-Center, P-Center), where P is point on surface of sphere
	//Ray: R = R_origin + t(R_direction)
	//Derived quadratic equation from two equations above
	
	//a, b, c coefficients in quadratic, solved using quadratic formula
	Vec3 s_center_to_r_origin = o - s.center;
	double a = 1.0;
	double b = 2.0 * dot(v, s_center_to_r_origin);
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

		return desired_solution;
	}
	else return DBL_MAX;

}

double ray_intersects_plane(Vec3 o, Vec3 v, Plane p)
{
	if(dot(v, p.n) != 0.0)
	{
		double desired_solution = dot(p.p - o, p.n)/dot(v, p.n);

		Vec3 intersection = o + desired_solution * v;
		double i_u = dot(intersection - p.p, p.u); //How far along u intersection is
		double i_v = dot(intersection - p.p, p.v); //How far along v intersection is
		
		bool within_u = i_u > -0.0009765625 && i_u <= p.u_length;
		bool within_v = i_v > -0.0009765625 && i_v <= p.v_length;
		bool intersection_within_bounds = within_u && within_v;
		if(intersection_within_bounds) return desired_solution;
	}
	return DBL_MAX;
}

double area(Sphere s)
{
	return 4.0 * PI * s.radius * s.radius;
}

double area(Plane p)
{
	return length(cross(p.u, p.v));
}

// Basic sampling

double uniform_sample()
{
	return (double)(rand() / (RAND_MAX + 1.0));
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

// Geometry sampling

double uniform_sample_sphere_subtended_pdf(Sphere s, Vec3 p)
{
	double st_max_sq = s.radius*s.radius / dot(s.center - p, s.center - p);
	double ct_max = sqrt(d_max(0.0, 1.0 - st_max_sq));
	double pdf = 1.0 / (2.0 * PI * (1.0 - ct_max));
	return pdf;
}

Vec3 uniform_sample_sphere_subtended(Sphere s, Vec3 p)
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

	return m * (-w) + s.center;
}

double uniform_sample_plane_pdf(Plane p)
{
	return 1.0 / area(p);
}

Vec3 uniform_sample_plane(Plane p)
{
	double u = uniform_sample();
	double v = uniform_sample();

	return p.p + u * p.u + v * p.v;
}


//Maths interface 
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

int clamp(int i, int low, int high)
{
	if(i < low) return low;
	else if(i > high) return high;
	else return i;
}

double lerp(double x, double x_0, double x_1, double y_0, double y_1)
{
	return y_0 + (x - x_0) * (y_1 - y_0) / (x_1 - x_0);
}
