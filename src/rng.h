void seed_rng(u32 seed);
f64  rng(); //Gens number between 0 and 1
vec3 uniform_sample_disc();
void uniform_sample_hemisphere(vec3 *dst, f64 *pdf, vec3 n);
