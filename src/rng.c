//Gens number between 0 and 1
f64 rng()
{
    f64 r = (f64)rand();
    f64 f = r / (f64)RAND_MAX;
    return f;
}

void seed_rng(u32 seed)
{
    srand(seed);
}

vec3 uniform_sample_sphere()
{
    f64  u      = rng();
    f64  v      = rng();
    f64  r      = sqrt(1.0 - u*u);
    f64  t      = 2.0 * PI * v;
    vec3 sample = {r * cos(t), r * sin(t), u};

    return sample;
}

vec3 uniform_sample_disc()
{
    vec3 v = {0.0, 0.0, 0.0};

    f64 r_x = rng();
    f64 r_y = rng();
    f64 o_x = 2.0 * r_x - 1.0;
    f64 o_y = 2.0 * r_y - 1.0;

    if(o_x == 0.0 && o_y == 0.0) return v;

    f64 r, t;
    if(fabs(o_x) > fabs(o_y))
    {
        r = o_x;
        t = (PI/4.0) * (o_y/o_x);
    }
    else
    {
        r = o_y;
        t = (PI/2.0) - (PI/4.0) * (o_x/o_y);
    }
    v.x = r * cos(t);
    v.y = r * sin(t);

    return v;
}

void uniform_sample_hemisphere(vec3 *v, f64 *pdf, vec3 n)
{
    *v = uniform_sample_sphere();
    vec3 i = {0.0, 0.0, 1.0};
    mat3x3 r = find_rotation_between_vectors(i, n);
    *v   = mat3x3_vec3_mul(r, *v);
    *pdf = 1.0 / (2.0 * PI);
}

void cos_weighted_sample_hemisphere(vec3 *v, f64 *pdf, vec3 n)
{
    vec3 p;
    for(;;)
    {
        p = uniform_sample_disc();
        if(vec3_dot(p, p) < 1.0) break;
    }
    p.z = sqrt(1.0 - vec3_dot(p, p));
    vec3 i = {0.0, 0.0, 1.0};
    mat3x3 r = find_rotation_between_vectors(i, n);
    *v   = mat3x3_vec3_mul(r, p);
    *pdf = vec3_dot(n, *v) / PI;
}
