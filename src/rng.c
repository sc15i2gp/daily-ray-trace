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
