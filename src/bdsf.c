//UTIL FUNCTIONS

f64 ggx(vec3 sn, vec3 mn, f64 r)
{
    f64 g;
    f64 sn_mn_dot = vec3_dot(sn, mn);
    f64 r_2 = r * r;
    if(sn_mn_dot <= 0.0)
    {
        g = 0.0;
    }
    else
    {
        f64 sn_mn_dot_2 = sn_mn_dot * sn_mn_dot;
        f64 sn_mn_dot_4 = sn_mn_dot_2 * sn_mn_dot_2;
        f64 sn_mn_tan_sq = (1.0/sn_mn_dot_2) - 1.0;
        g = r_2 / (PI * sn_mn_dot_4 * (r_2 + sn_mn_tan_sq) * (r_2 + sn_mn_tan_sq));
    }
    return g;
}

f64 ggx_att(vec3 v, vec3 sn, vec3 mn, f64 r)
{
    f64 att;
    f64 g = ggx(sn, mn, r);
    f64 v_mn_dot = vec3_dot(v, mn);
    f64 v_sn_dot = vec3_dot(v, sn);
    f64 dot_quot = fabs(v_mn_dot / v_sn_dot);
    f64 r_2 = r * r;

    if(dot_quot <= 0.0)
    {
        att = 0.0;
    }
    else
    {
        f64 vn_tan_sq = (1.0/(v_sn_dot*v_sn_dot)) - 1.0;
        att = 2.0 / (1.0 + sqrt(1.0 + r_2 * vn_tan_sq));
    }

    return g * att;
}

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

void fs_conductor_reflectance(spectrum reflectance, spectrum ir, spectrum tr, spectrum te, f64 inc_cos)
{
    f64 inc_cos_sq = inc_cos * inc_cos;
    f64 inc_sin_sq = 1.0 - inc_cos_sq;

    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        f64 relative_refract = tr.samples[i] / ir.samples[i];
        f64 relative_extinct = te.samples[i] / ir.samples[i];
        f64 relative_refract_sq = relative_refract * relative_refract;
        f64 relative_extinct_sq = relative_extinct * relative_extinct;
        f64 r = relative_refract_sq - relative_extinct_sq - inc_sin_sq;
        f64 apb_sq = sqrt(r * r + 4.0 * relative_refract_sq * relative_extinct_sq);
        f64 a = sqrt(0.5 * (apb_sq + r));
        f64 s = apb_sq + inc_cos_sq;
        f64 t = 2.0 * a * inc_cos;
        f64 u = inc_cos_sq * apb_sq + inc_sin_sq * inc_sin_sq;
        f64 v = t * inc_sin_sq;
        f64 parallel_reflectance = (s - t) / (s + t);
        f64 perpend_reflectance  = parallel_reflectance * (u - v) / (u + v);

        reflectance.samples[i] = 0.5 * (parallel_reflectance + perpend_reflectance);
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
        fs_conductor_reflectance(reflectance, incident_refract, transmit_refract, transmit_extinct, on_dot);
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

void fs_dielectric_transmittance_bdsf(spectrum transmittance, scene_point *p, vec3 incoming)
{
    vec3 outgoing = vec3_reverse(p->out);
    spectrum ir_spd = p->incident_material->refract_spd;
    spectrum tr_spd = p->transmit_material->refract_spd;
    f64 ir = value_at_wl(ir_spd, p->trans_wl);
    f64 tr = value_at_wl(tr_spd, p->trans_wl);
    if(vec3_equal(incoming, vec3_transmit(outgoing, p->normal, ir, tr)))
    {
        fs_dielectric_transmittance(transmittance, ir_spd, tr_spd, p->on_dot);
    }
}

void ct_conductor_bdsf(spectrum reflectance, scene_point *p, vec3 incoming)
{
    vec3 micro_normal = vec3_normalise(vec3_sum(p->out, incoming));
    f64 in_dot = fabs(vec3_dot(p->normal, incoming));
    f64 mn_dot = fabs(vec3_dot(p->normal, micro_normal));

    spectrum incident_refract = p->incident_material->refract_spd;
    spectrum transmit_refract = p->transmit_material->refract_spd;
    spectrum transmit_extinct = p->transmit_material->extinct_spd;
    fs_conductor_reflectance(reflectance, incident_refract, transmit_refract, transmit_extinct, mn_dot);
    f64 reflectance_coefficient = ggx_att(p->out, p->normal, micro_normal, p->surface_material->roughness) * (1.0/(4.0*p->on_dot));
    spectral_mul_by_scalar(reflectance, reflectance, reflectance_coefficient);
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

    f64 ir = value_at_wl(ir_spd, p->trans_wl);
    f64 tr = value_at_wl(tr_spd, p->trans_wl);
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
    f64 rd = value_at_wl(reflectance, p->trans_wl);
    f64 ir = value_at_wl(ir_spd, p->trans_wl);
    f64 tr = value_at_wl(tr_spd, p->trans_wl);

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

void sample_ct_direction(vec3 *v, f64 *pdf, scene_point *p)
{
    do
    {
        f64 f = rng();
        f64 g = rng();

        f64 phi_mn = 2.0 * PI * g;
        f64 tan_mn = (p->surface_material->roughness * sqrt(f)) / sqrt(1.0 - f);
        f64 cos_mn = 1.0 / sqrt(1.0 + tan_mn*tan_mn);
        f64 sin_mn = sqrt(1.0 - cos_mn * cos_mn);

        vec3 i = {0.0, 0.0, 1.0};
        vec3 micro_normal = {sin_mn * cos(phi_mn), sin_mn * sin(phi_mn), cos_mn};
        mat3x3 r = find_rotation_between_vectors(i, p->normal);
        micro_normal = mat3x3_vec3_mul(r, micro_normal);

        f64 sn_mn_dot = vec3_dot(p->normal, micro_normal);
        if(sn_mn_dot < 0.0)
        {
            micro_normal = vec3_reverse(micro_normal);
            sn_mn_dot = -sn_mn_dot;
        }
        f64 o_mn_dot = vec3_dot(p->out, micro_normal);

        vec3 w = vec3_reverse(p->out);
        *v = vec3_reflect(w, micro_normal);
        f64 d = ggx(p->normal, micro_normal, p->surface_material->roughness) * sn_mn_dot;
        *pdf = ((4.0 * o_mn_dot)/d);
    }
    while (vec3_dot(*v, p->normal) < 0.0);
}
