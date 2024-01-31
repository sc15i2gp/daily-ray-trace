void print_vector(vec3 v)
{
    printf("%f %f %f", v.x, v.y, v.z);
}

vec3 vec3_sum(vec3 a, vec3 b)
{
    vec3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

vec3 vec3_sub(vec3 a, vec3 b)
{
    vec3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}

f64 vec3_dot(vec3 a, vec3 b)
{
    f64 result = a.x*b.x + a.y*b.y + a.z*b.z;
    return result;
}

vec3 vec3_cross(vec3 a, vec3 b)
{
    vec3 n;
    n.x = a.y*b.z - a.z*b.y;
    n.y = a.z*b.x - a.x*b.z;
    n.z = a.x*b.y - a.y*b.x;
    return n;
}

vec3 vec3_mul_by_f64(vec3 v, f64 f)
{
    vec3 result;
    result.x = f * v.x;
    result.y = f * v.y;
    result.z = f * v.z;
    return result;
}

vec3 vec3_div_by_f64(vec3 v, f64 f)
{
    vec3 result;
    result.x = v.x / f;
    result.y = v.y / f;
    result.z = v.z / f;
    return result;
}

f64 vec3_length(vec3 v)
{
    f64 result = vec3_dot(v, v);
    result = sqrt(result);
    return result;
}

vec3 vec3_normalise(vec3 v)
{
    vec3 result = v;
    result = vec3_div_by_f64(v, vec3_length(v));
    return result;
}

vec3 vec3_reverse(vec3 v)
{
    vec3 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;
    return result;
}

f64 point_to_line_distance(vec3 p, vec3 l_0, vec3 l_1)
{
    vec3 normalised_l = vec3_normalise(vec3_sub(l_0, l_1));
    vec3 q = vec3_sum(l_1, vec3_mul_by_f64(normalised_l, vec3_dot(vec3_sub(p, l_1), normalised_l)));
    vec3 q_to_p = vec3_sub(p, q);
    f64 distance = vec3_dot(q_to_p, q_to_p);
    return distance;
}

//line_o: line origin
//line_d: normalised line direction
//sphere_c: sphere center
//sphere_r: sphere radius
//Returns the smallest positive intersection point from line_o along line_d through the sphere
//Returns NaN if there is no intersection point, or if there are only negative intersection points
f64 line_sphere_intersection(vec3 line_o, vec3 line_d, vec3 sphere_c, f64 sphere_r)
{
    vec3 c_to_o = vec3_sub(line_o, sphere_c);

    f64 a               = 1.0;
    f64 b               = -2.0 * vec3_dot(c_to_o, line_d);
    f64 c               = vec3_dot(c_to_o, c_to_o) - sphere_r*sphere_r;
    f64 discriminant    = b*b - 4.0 * a * c;

    if(discriminant < 0.0) return INFINITY;

    f64 sqrt_discriminant = sqrt(discriminant);

    f64 a_2 = 2.0 * a;

    f64 solution_0 = (b + sqrt_discriminant)/a_2;
    f64 solution_1 = (b - sqrt_discriminant)/a_2;

    if      (solution_0 < 0.0 && solution_1 < 0.0)  return INFINITY;    //No positive solutions
    else if (solution_0 >= 0.0 && solution_1 < 0.0) return solution_0;  //Only solution_0 positive
    else if (solution_1 >= 0.0 && solution_0 < 0.0) return solution_1;  //Only solution_1 positive
    else if (solution_0 <= solution_1)              return solution_0;  //Smallest solution_0
    else                                            return solution_1;  //Smallest solution_1
}

//plane consists of a reference point as origin, a normal and two bounding vectors extending from reference
//line_o: line origin
//line_d: normalised line direction
//plane_p: plane reference point
//plane_n: normalised plane normal
//plane_u: first plane bounds vector
//plane_v: second plane bounds vector
//Returns positive intersection point from line_o along line_d through the plane
//Returns NaN if there is no positive intersection point
f64 line_plane_intersection(vec3 line_o, vec3 line_d, vec3 plane_p, vec3 plane_n, vec3 plane_u, vec3 plane_v)
{
    if(vec3_dot(line_d, plane_n) == 0.0) return INFINITY; //Line parallel to plane
    
    vec3    o_to_p  = vec3_sub(plane_p, line_o);
    f64     l       = vec3_dot(o_to_p, plane_n) / vec3_dot(line_d, plane_n);
    vec3    i       = vec3_sum(line_o, vec3_mul_by_f64(line_d, l));
    vec3    j       = vec3_sub(i, plane_p);

    f64 plane_u_length          = vec3_length(plane_u);
    f64 plane_v_length          = vec3_length(plane_v);
    f64 j_dot_plane_u           = vec3_dot(j, plane_u);
    f64 j_dot_plane_v           = vec3_dot(j, plane_v);

    if( 0.0 <= j_dot_plane_u && j_dot_plane_u <= plane_u_length &&
        0.0 <= j_dot_plane_v && j_dot_plane_v <= plane_v_length)
    {
        return l;
    }

    return INFINITY;
}

//triangle consists of three points in CCW order when viewed down from the normal
//line_o: line origin
//line_d: normalised line direction
f64 line_triangle_intersection(vec3 line_o, vec3 line_d, vec3 triangle_a, vec3 triangle_b, vec3 triangle_c, vec3 triangle_n)
{
    vec3 o_to_a = vec3_sub(triangle_a, line_o);
    f64 l = vec3_dot(o_to_a, triangle_n)/vec3_dot(line_d, triangle_n);
    vec3 i = vec3_sum(line_o, vec3_mul_by_f64(line_d, l));

    f64 a = point_to_line_distance(i, triangle_b, triangle_c)/point_to_line_distance(triangle_a, triangle_b, triangle_c);
    f64 b = point_to_line_distance(i, triangle_c, triangle_a)/point_to_line_distance(triangle_b, triangle_c, triangle_a);
    f64 c = point_to_line_distance(i, triangle_a, triangle_b)/point_to_line_distance(triangle_c, triangle_a, triangle_b);
    f64 d = a + b + c;

    if(d > 1.0) return NAN;

    return l;
}

void create_plane_from_points(vec3 o_point, vec3 u_point, vec3 v_point, vec3 *plane_o, vec3 *plane_u, vec3 *plane_v, vec3 *plane_n)
{
    *plane_o = o_point;
    *plane_u = vec3_sub(u_point, o_point);
    *plane_v = vec3_sub(v_point, o_point);
    *plane_n = vec3_normalise(vec3_cross(*plane_u, *plane_v));
}
