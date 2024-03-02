//UTIL FUNCTIONS
void fs_dielectric_reflectance(spectrum reflectance, spectrum ir, spectrum tr, f64 inc_cos)
{
    f64 inc_sin_sq = 1.0 - inc_cos * inc_cos;
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        f64 relative_refract = ir.samples[i] / tr.samples[i];
        f64 ts_sin_sq = relative_refract * relative_refract * inc_sin_sq;
        if(ts_sin_sq >= 1.0)
        {
            reflectance.samples[i] = 1.0;
            continue;
        }
        f64 ts_cos = sqrt(1.0 - ts_sin_sq * ts_sin_sq);
        f64 tr_by_on_cos = tr.samples[i] * inc_cos;
        f64 tr_by_ts_cos = tr.samples[i] * ts_cos;
        f64 ir_by_on_cos = ir.samples[i] * inc_cos;
        f64 ir_by_ts_cos = ir.samples[i] * ts_cos;
        f64 parallel_reflectance = (tr_by_on_cos - ir_by_ts_cos)/(tr_by_on_cos + ir_by_ts_cos);
        f64 perpend_reflectance = (ir_by_on_cos - tr_by_ts_cos)/(ir_by_on_cos + tr_by_ts_cos);
        parallel_reflectance *= parallel_reflectance;
        perpend_reflectance  *= perpend_reflectance;
        reflectance.samples[i] = 0.5 * (parallel_reflectance + perpend_reflectance);
    }
}

void fs_dielectric_transmittance(spectrum transmittance, spectrum ir, spectrum tr, f64 inc_cos)
{
    fs_dielectric_reflectance(transmittance, ir, tr, inc_cos);
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        transmittance.samples[i] = 1.0 - transmittance.samples[i];
    }
}

//BDSF FUNCTIONS

void bp_diffuse_bdsf(spectrum reflectance, scene_point *p, vec3 incoming)
{
    spectral_mul_by_scalar(reflectance, p->surface_material->diffuse_spd, 1.0/PI);
    spectral_mul_by_scalar(reflectance, reflectance, fabs(vec3_dot(p->normal, incoming)));
}

void bp_glossy_bdsf(spectrum reflectance, scene_point *p, vec3 incoming)
{
    vec3 outgoing         = p->out;
    vec3 bisector         = vec3_normalise(vec3_sum(outgoing, incoming));
    f64  spec_coefficient = pow(f64_max(0.0, vec3_dot(p->normal, bisector)), p->surface_material->shininess);

    spectral_mul_by_scalar(reflectance, p->surface_material->glossy_spd, spec_coefficient);
    spectral_mul_by_scalar(reflectance, reflectance, fabs(vec3_dot(p->normal, incoming)));
}

void mirror_bdsf(spectrum reflectance, scene_point *p, vec3 incoming)
{
    vec3 outgoing = vec3_reverse(p->out);
    if(vec3_equal(incoming, vec3_reflect(outgoing, p->normal)))
    {
        copy_spectrum(reflectance, p->surface_material->mirror_spd);
    }
    else
    {
        zero_spectrum(reflectance);
    }
}

void fs_conductor_bdsf(spectrum reflectance, scene_point *p, vec3 incoming)
{
    vec3 outgoing = vec3_reverse(p->out);
    if(vec3_equal(incoming, vec3_reflect(outgoing, p->normal)))
    {
        spectrum incident_refract = p->incident_material->refract_spd;
        spectrum transmit_refract = p->transmit_material->refract_spd;
        spectrum transmit_extinct = p->transmit_material->extinct_spd;

        f64 on_dot = p->on_dot;
        f64 in_dot = vec3_dot(p->normal, incoming);
        f64 on_dot_sq = on_dot * on_dot;
        f64 on_sin_sq = 1.0 - on_dot_sq;

        for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
        {
            f64 relative_refract = transmit_refract.samples[i] / incident_refract.samples[i];
            f64 relative_extinct = transmit_extinct.samples[i] / incident_refract.samples[i];
            f64 relative_refract_sq = relative_refract * relative_refract;
            f64 relative_extinct_sq = relative_extinct * relative_extinct;
            f64 r = relative_refract_sq - relative_extinct_sq - on_sin_sq;
            f64 apb_sq = sqrt(r * r + 4.0 * relative_refract_sq * relative_extinct_sq);
            f64 a = sqrt(0.5 * (apb_sq + r));
            f64 s = apb_sq + on_dot_sq;
            f64 t = 2.0 * a * on_dot;
            f64 u = on_dot_sq * apb_sq + on_sin_sq * on_sin_sq;
            f64 v = t * on_sin_sq;
            f64 parallel_reflectance = (s - t) / (s + t);
            f64 perpend_reflectance  = parallel_reflectance * (u - v) / (u + v);

            reflectance.samples[i] = 0.5 * (parallel_reflectance + perpend_reflectance);
        }
    }
}

void fs_dielectric_reflectance_bdsf(spectrum reflectance, scene_point *p, vec3 incoming)
{
    vec3 outgoing = vec3_reverse(p->out);
    if(vec3_equal(incoming, vec3_reflect(outgoing, p->normal)))
    {
        spectrum incident_refract = p->incident_material->refract_spd;
        spectrum transmit_refract = p->transmit_material->refract_spd;

        f64 on_dot = p->on_dot;
        fs_dielectric_reflectance(reflectance, incident_refract, transmit_refract, on_dot);
    }
}

#define TMP_TRANSMIT_WL 50
void fs_dielectric_transmittance_bdsf(spectrum transmittance, scene_point *p, vec3 incoming)
{
    vec3 outgoing = vec3_reverse(p->out);
    spectrum ir_spd = p->incident_material->refract_spd;
    spectrum tr_spd = p->transmit_material->refract_spd;
    f64 ir = ir_spd.samples[TMP_TRANSMIT_WL];
    f64 tr = tr_spd.samples[TMP_TRANSMIT_WL];
    if(vec3_equal(incoming, vec3_transmit(outgoing, p->normal, ir, tr)))
    {
        fs_dielectric_transmittance(transmittance, ir_spd, tr_spd, p->on_dot);
    }
}

//DIRECTION SAMPLING FUNCTIONS

//NOTE: Direction sampling functions should return the reciprocal of their pdf values
void uniform_sample_hemisphere(vec3 *v, f64 *pdf, scene_point *p)
{
    *v = uniform_sample_sphere();
    vec3 i = {0.0, 0.0, 1.0};
    mat3x3 r = find_rotation_between_vectors(i, p->normal);
    *v   = mat3x3_vec3_mul(r, *v);
    *pdf = (2.0 * PI);
}

void cos_weighted_sample_hemisphere(vec3 *v, f64 *pdf, scene_point *p)
{
    vec3 q;
    for(;;)
    {
        q = uniform_sample_disc();
        if(vec3_dot(q, q) < 1.0) break;
    }
    q.z = sqrt(1.0 - vec3_dot(q, q));
    vec3 i = {0.0, 0.0, 1.0};
    mat3x3 r = find_rotation_between_vectors(i, p->normal);
    *v   = mat3x3_vec3_mul(r, q);
    *pdf = PI / vec3_dot(p->normal, *v);
}

void sample_specular_direction(vec3 *v, f64 *pdf, scene_point *p)
{
    vec3 w = vec3_reverse(p->out);
    *v = vec3_reflect(w, p->normal);
    *pdf = 1.0;
}

void sample_transmit_direction(vec3 *v, f64 *pdf, scene_point *p)
{
    spectrum reflectance = alloc_spd();
    spectrum ir_spd = p->incident_material->refract_spd;
    spectrum tr_spd = p->transmit_material->refract_spd;

    f64 ir = ir_spd.samples[TMP_TRANSMIT_WL];
    f64 tr = tr_spd.samples[TMP_TRANSMIT_WL];
    *v   = vec3_transmit(vec3_reverse(p->out), p->normal, ir, tr);
    *pdf = 1.0;

    free_spd(reflectance);
}

void sample_reflect_or_transmit_direction(vec3 *v, f64 *pdf, scene_point *p)
{
    spectrum reflectance = alloc_spd();
    spectrum ir_spd = p->incident_material->refract_spd;
    spectrum tr_spd = p->transmit_material->refract_spd;
    fs_dielectric_reflectance(reflectance, ir_spd, tr_spd, p->on_dot);
    f64 rd = reflectance.samples[TMP_TRANSMIT_WL];
    f64 ir = ir_spd.samples[TMP_TRANSMIT_WL];
    f64 tr = tr_spd.samples[TMP_TRANSMIT_WL];

    f64 f = rng();
    vec3 w = vec3_reverse(p->out);
    if(f < rd)
    {
        *v = vec3_reflect(w, p->normal);
        *pdf = 1.0/rd;
    }
    else
    {
        *v = vec3_transmit(w, p->normal, ir, tr);
        *pdf = 1.0/(1.0 - rd);
    }
    free_spd(reflectance);
}
