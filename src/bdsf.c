void bp_diffuse_bdsf(spectrum reflectance, scene_point *p, vec3 incoming, vec3 outgoing)
{
    spectral_mul_by_scalar(reflectance, p->material->diffuse_spd, 1.0/PI);
}

void bp_glossy_bdsf(spectrum reflectance, scene_point *p, vec3 incoming, vec3 outgoing)
{
    vec3 bisector         = vec3_normalise(vec3_sum(outgoing, incoming));
    f64  spec_coefficient = pow(f64_max(0.0, vec3_dot(p->normal, bisector)), p->material->shininess);

    spectral_mul_by_scalar(reflectance, p->material->glossy_spd, spec_coefficient);
}

void mirror_bdsf(spectrum reflectance, scene_point *p, vec3 incoming, vec3 outgoing)
{
    outgoing = vec3_reverse(outgoing);
    if(vec3_equal(incoming, vec3_reflect(outgoing, p->normal)))
    {
        copy_spectrum(reflectance, p->material->mirror_spd);
    }
    else
    {
        zero_spectrum(reflectance);
    }
}

void uniform_sample_hemisphere(vec3 *v, f64 *pdf, vec3 n, vec3 w)
{
    *v = uniform_sample_sphere();
    vec3 i = {0.0, 0.0, 1.0};
    mat3x3 r = find_rotation_between_vectors(i, n);
    *v   = mat3x3_vec3_mul(r, *v);
    *pdf = 1.0 / (2.0 * PI);
}

void cos_weighted_sample_hemisphere(vec3 *v, f64 *pdf, vec3 n, vec3 w)
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

void sample_specular_direction(vec3 *v, f64 *pdf, vec3 n, vec3 w)
{
    w = vec3_reverse(w);
    *v = vec3_reflect(w, n);
    *pdf = 1.0;
}
