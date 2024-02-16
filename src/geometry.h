typedef struct
{
    union
    {
        struct
        {
            f64 x;
            f64 y;
            f64 z;
        };
        f64 xyz[3];
    };
} vec3;

void print_vector(vec3);
vec3 vec3_sum(vec3, vec3);
vec3 vec3_sub(vec3, vec3);
vec3 vec3_mul_by_f64(vec3, f64);
vec3 vec3_div_by_f64(vec3, f64);
f64  vec3_dot(vec3, vec3);
f64  vec3_length(vec3);
vec3 vec3_normalise(vec3);
vec3 vec3_reverse(vec3);
f64  point_to_line_distance(vec3 p, vec3 l_0, vec3 l_1);
f64  line_sphere_intersection(vec3 line_o, vec3 line_d, vec3 sphere_c, f64 sphere_r);
f64  line_plane_intersection(vec3 line_o, vec3 line_d, vec3 plane_p, vec3 plane_n, vec3 plane_u, vec3 plane_v);
f64  line_triangle_intersection(vec3 line_o, vec3 line_d, vec3 triangle_a, vec3 triangle_b, vec3 triangle_c, vec3 triangle_n);
void create_plane_from_points(vec3 o_point, vec3 u_point, vec3 v_point, vec3 *plane_o, vec3 *plane_u, vec3 *plane_v, vec3 *plane_n);

typedef struct
{
    vec3 columns[3];
} mat3x3;

vec3 row(mat3x3 m, u32 r);
