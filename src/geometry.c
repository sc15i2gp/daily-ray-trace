void print_vector(vec3 v)
{
    printf("%f %f %f", v.x, v.y, v.z);
}

u32 vec3_equal(vec3 a, vec3 b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
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

vec3 vec3_reflect(vec3 v, vec3 n)
{
    f64 f = 2.0 * vec3_dot(v, n);
    v = vec3_sub(v, vec3_mul_by_f64(n, f));
    return v;
}

vec3 vec3_transmit(vec3 v, vec3 n, f64 ir, f64 tr)
{
    f64  vn_dot  = vec3_dot(v, n);
    f64  rel_ref = ir/tr;
    
    vec3 m = vec3_mul_by_f64(n, vn_dot);
    v = vec3_sub(m, v);
    
    vec3 perpend = vec3_reverse(vec3_mul_by_f64(v, rel_ref));
    f64  perpend_dot = -sqrt(1.0 - vec3_dot(perpend, perpend));
    vec3 parallel = vec3_mul_by_f64(n, perpend_dot);

    vec3 transmit = vec3_sum(perpend, parallel);
    return transmit;
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

    vec3 u_norm = vec3_normalise(plane_u);
    vec3 v_norm = vec3_normalise(plane_v);
    f64 j_dot_plane_u           = vec3_dot(j, u_norm);
    f64 j_dot_plane_v           = vec3_dot(j, v_norm);

    if( l >= 0.0 &&
        0.0 <= j_dot_plane_u && j_dot_plane_u <= plane_u_length &&
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

vec3 row(mat3x3 m, u32 r)
{
    vec3 v;
    v.x = m.columns[0].xyz[r];
    v.y = m.columns[1].xyz[r];
    v.z = m.columns[2].xyz[r];
    return v;
}

vec3 mat3x3_vec3_mul(mat3x3 m, vec3 v)
{
    vec3 w;
    for(u32 i = 0; i < 3; i += 1)
    {
        w.xyz[i] = vec3_dot(row(m, i), v);
    }
    return w;
}

mat3x3 mat3x3_sum(mat3x3 m, mat3x3 n)
{
    mat3x3 r;
    for(u32 i = 0; i < 3; i += 1)
    {
        r.columns[i] = vec3_sum(m.columns[i], n.columns[i]);
    }
    return r;
}

mat3x3 mat3x3_mul(mat3x3 m, mat3x3 n)
{
    mat3x3 r;
    
    for(u32 i = 0; i < 3; i += 1)
    {
        for(u32 j = 0; j < 3; j += 1)
        {
            r.columns[i].xyz[j] = vec3_dot(row(m, i), n.columns[j]);
        }
    }
    return r;
}

mat3x3 mat3x3_mul_by_f64(mat3x3 m, f64 f)
{
    for(u32 i = 0; i < 3; i += 1)
    {
        m.columns[i] = vec3_mul_by_f64(m.columns[i], f);
    }
    return m;
}

mat3x3 find_rotation_between_vectors(vec3 v, vec3 w)
{
    vec3 n = vec3_cross(v, w);
    f64  s = vec3_length(n);
    f64  c = vec3_dot(v, w);
    
    mat3x3 r = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    if(vec3_dot(n, n) == 0.0 && c <= 0.0)
    {
        r.columns[0].xyz[0] = -1.0;
        r.columns[1].xyz[1] = -1.0;
        r.columns[2].xyz[2] = -1.0;
    }
    else
    {
        mat3x3 m;
        m.columns[0].xyz[0] = 0.0;
        m.columns[0].xyz[1] = n.z;
        m.columns[0].xyz[2] = -n.y;
        m.columns[1].xyz[0] = -n.z;
        m.columns[1].xyz[1] = 0.0;
        m.columns[1].xyz[2] = n.x;
        m.columns[2].xyz[0] = n.y;
        m.columns[2].xyz[1] = -n.x;
        m.columns[2].xyz[2] = 0.0;

        mat3x3 i = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        mat3x3 mm = mat3x3_mul(m, m);
        mm = mat3x3_mul_by_f64(mm, (1.0/(1.0+c)));
        r = mat3x3_sum(mat3x3_sum(i, m), mm);
    }
    return r;
}

mat3x3 rotation_about_axis(vec3 axis, f64 angle_rad)
{
    f64 cos_th = cos(angle_rad);
    f64 sin_th = sin(angle_rad);
    mat3x3 r;
    r.columns[0].x = cos_th+(axis.x*axis.x)*(1-cos_th);
    r.columns[0].y = axis.y*axis.x*(1-cos_th)+axis.z*sin_th;
    r.columns[0].z = axis.z*axis.z*(1-cos_th)-axis.y*sin_th;
    r.columns[1].x = axis.x*axis.y*(1-cos_th)-axis.z*sin_th;
    r.columns[1].y = cos_th+(axis.y*axis.y)*(1-cos_th);
    r.columns[1].z = axis.z*axis.y*(1-cos_th)+axis.x*sin_th;
    r.columns[2].x = axis.x*axis.z*(1-cos_th)+axis.y*sin_th;
    r.columns[2].y = axis.y*axis.z*(1-cos_th)-axis.x*sin_th;
    r.columns[2].z = cos_th+axis.z*axis.z*(1-cos_th);

    return r;
}
