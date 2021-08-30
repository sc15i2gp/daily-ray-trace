
double fresnel_reflectance_dielectric(double incident_refraction_index, double transmit_refraction_index, double incident_cos, double transmit_cos)
{
	double parallel_reflectance = 
		(transmit_refraction_index * incident_cos - incident_refraction_index * transmit_cos) / 
		(transmit_refraction_index * incident_cos + incident_refraction_index * transmit_cos);
	double perpendicular_reflectance = 
		(incident_refraction_index * incident_cos - transmit_refraction_index * transmit_cos) / 
		(incident_refraction_index * incident_cos + transmit_refraction_index * transmit_cos);

	double reflectance = 0.5 * (parallel_reflectance*parallel_reflectance + perpendicular_reflectance*perpendicular_reflectance);
	return reflectance;
}

double fresnel_reflectance_dielectric(double incident_refraction_index, double transmit_refraction_index, double incident_cos)
{
	double relative_refract_index = incident_refraction_index / transmit_refraction_index;
	double incident_sin_sq = d_max(0.0, 1.0 - incident_cos * incident_cos);
	double transmit_sin_sq = relative_refract_index * relative_refract_index * incident_sin_sq;
	
	//Total internal reflection
	if(transmit_sin_sq >= 1.0) 
	{
		return 1.0;
	}

	double transmit_cos = sqrt(d_max(0.0, 1.0 - transmit_sin_sq * transmit_sin_sq));

	return fresnel_reflectance_dielectric(incident_refraction_index, transmit_refraction_index, incident_cos, transmit_cos);
}


void fresnel_reflectance_dielectric(Spectrum& incident_refract_index, Spectrum& transmit_refract_index, double incident_cos, Spectrum& reflectance)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		reflectance.samples[i] = fresnel_reflectance_dielectric(incident_refract_index.samples[i], transmit_refract_index.samples[i], incident_cos);
	}
}

double fresnel_transmittance_dielectric(double incident_refraction_index, double transmit_refraction_index, double incident_cos, double transmit_cos)
{
	double transmittance = 1.0 - fresnel_reflectance_dielectric(incident_refraction_index, transmit_refraction_index, incident_cos, transmit_cos);
	return transmittance;
}

double fresnel_reflectance_conductor(double incident_refract_index, double transmit_refract_index, double transmit_extinct_index, double incident_cos)
{
	double relative_refract_index = transmit_refract_index / incident_refract_index;
	double relative_extinct_index = transmit_extinct_index / incident_refract_index;

	double incident_cos_sq = incident_cos * incident_cos;
	double incident_sin_sq = 1.0 - incident_cos_sq;
	double relative_refract_index_sq = relative_refract_index * relative_refract_index;
	double relative_extinct_index_sq = relative_extinct_index * relative_extinct_index;
	double r = relative_refract_index_sq - relative_extinct_index_sq - incident_sin_sq;
	double a_sq_plus_b_sq = sqrt(r * r + 4.0 * relative_refract_index_sq * relative_extinct_index_sq);
	double a = sqrt(0.5 * (a_sq_plus_b_sq + r));
	double s = a_sq_plus_b_sq + incident_cos_sq;
	double t = 2.0 * a * incident_cos;
	double u = incident_cos_sq * a_sq_plus_b_sq + incident_sin_sq * incident_sin_sq;
	double v = t * incident_sin_sq;
	double parallel_reflectance = (s - t) / (s + t);
	double perpendicular_reflectance = parallel_reflectance * (u - v) / (u + v);

	return 0.5 * (parallel_reflectance + perpendicular_reflectance);
}

void fresnel_reflectance_conductor(Spectrum& incident_refract_index, Spectrum& transmit_refract_index, Spectrum& transmit_extinct_index, double incident_cos, Spectrum& reflectance)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		reflectance.samples[i] = fresnel_reflectance_conductor(incident_refract_index.samples[i], transmit_refract_index.samples[i], transmit_extinct_index.samples[i], incident_cos);
	}
}

void fresnel_transmittance_dielectric(Spectrum& incident_refract_index, Spectrum& transmit_refract_index, double incident_cos, double transmit_cos, Spectrum& transmittance)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		transmittance.samples[i] = fresnel_transmittance_dielectric(incident_refract_index.samples[i], transmit_refract_index.samples[i], incident_cos, transmit_cos);
	}
}

//NOTE: The index of refraction chosen is the spectrum's value at 630 nm
//	Most materials' given refraction indices are the single values at 633 nm
//	630 is the closest value which can be used to this which is likely to work for any interval length
//	This is temporary and a better solution will be implemented as necessary.
#define TRANSMISSION_WAVELENGTH 630

Vec3 transmit_vector(double incident_refract_index, double transmit_refract_index, Vec3 normal, Vec3 outgoing)
{
	//Handles entering and leaving medium
	if(dot(outgoing, normal) < 0.0) normal = -normal;

	//Parallel and perpendicular to p.normal
	Vec3 transmit_perpendicular_component = (incident_refract_index/transmit_refract_index) * (dot(outgoing, normal)*normal -outgoing);
	Vec3 transmit_parallel_component = -sqrt(1.0 - dot(transmit_perpendicular_component, transmit_perpendicular_component))*normal;
	
	Vec3 transmit_direction = transmit_perpendicular_component + transmit_parallel_component;
	return transmit_direction;
}

Vec3 transmit_vector(Spectrum& incident_refract_index_spd, Spectrum& transmit_refract_index_spd, Vec3 normal, Vec3 outgoing)
{

	//Handles entering and leaving medium
	if(dot(outgoing, normal) < 0.0) normal = -normal;

	double incident_refract_index = spd_value_at_wavelength(incident_refract_index_spd, TRANSMISSION_WAVELENGTH);
	double transmit_refract_index = spd_value_at_wavelength(transmit_refract_index_spd, TRANSMISSION_WAVELENGTH);
	return transmit_vector(incident_refract_index, transmit_refract_index, normal, outgoing);
}

//NOTE: Isotropic roughness
double microfacet_distribution(Surface_Point& p, Vec3 v)
{
	double cos_th_sq = dot(p.normal, v) * dot(p.normal, v);
	double sin_th_sq = 1.0 - cos_th_sq;
	double sin_th = sqrt(sin_th_sq);
	double tan_th_sq = sin_th_sq / cos_th_sq;
	double cos_ph = (sin_th == 0.0) ? 1.0 : clamp(v.x / sin_th, -1.0, 1.0);
	double sin_ph = (sin_th == 0.0) ? 0.0 : clamp(v.y / sin_th, -1.0, 1.0);
	double cos_ph_sq = cos_ph * cos_ph;
	double sin_ph_sq = sin_ph * sin_ph;

	double r = *TEXTURE_SAMPLE(double, p.surface_material->roughness_texture, p.texture_coordinates);
	double r_sq = r * r;
	double d = exp(-tan_th_sq / r_sq) / (PI * r_sq * cos_th_sq * cos_th_sq);

	return d;
}

double lambda(Surface_Point& p, Vec3 v)
{
	double cos_th = dot(p.normal, v);
	double sin_th = sqrt(1.0 - cos_th * cos_th);
	double abs_tan_th = abs(sin_th/cos_th);
	double cos_ph = (sin_th == 0.0) ? 1.0 : clamp(v.x / sin_th, -1.0, 1.0);
	double sin_ph = (sin_th == 0.0) ? 0.0 : clamp(v.y / sin_th, -1.0, 1.0);
	double cos_ph_sq = cos_ph * cos_ph;
	double sin_ph_sq = sin_ph * sin_ph;

	double r = *TEXTURE_SAMPLE(double, p.surface_material->roughness_texture, p.texture_coordinates);
	double r_sq = r * r; 
	double alpha = sqrt(cos_ph_sq * r_sq + sin_ph_sq * r_sq); //alpha = sqrt(2) * sigma
	double a = 1.0 / (alpha * abs_tan_th);
	
	if(a >= 1.6) return 1.0;

	double d = (1.0 - 1.259*a + 0.396*a*a) / (3.535*a + 2.181*a*a);
	return d;
}

double geometric_attenuation(Surface_Point& p, Vec3 incoming, Vec3 outgoing)
{
	return 1.0 / (1.0 + lambda(p, outgoing) + lambda(p, incoming));
}

double geometric_attenuation(Surface_Point& p, Vec3 v)
{
	return 1.0 / (1.0 + lambda(p, v));
}

//Cite: Microfacet models for refraction through rough surfaces
//TODO: Compare performance of working out trig values from identities vs inverse functions
//	e.g. sin(acos(n_h_dot) = sqrt(1 - n_h_dot*n_h_dot) 
double beckmann_distribution(Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double n_h_dot = dot(surface_normal, microfacet_normal);
	if(n_h_dot <= 0.0) return 0.0;
	else
	{
		double a_sq = lobe_width * lobe_width;
		double cos_h_sq = n_h_dot * n_h_dot;
		double cos_h_quart = cos_h_sq * cos_h_sq;
		double tan_h_sq = (1.0/cos_h_sq) - 1.0;
		double d = 1.0/(PI * a_sq * cos_h_quart);
		double e = exp(-tan_h_sq/a_sq);

		return d * e;
		/*
		double nh_dot_sq = n_h_dot * n_h_dot;
		double e_quot = (nh_dot_sq - 1.0) / (lobe_width*lobe_width*nh_dot_sq);
		double denominator = PI * lobe_width * lobe_width * nh_dot_sq * nh_dot_sq;
		double d = exp(e_quot) / denominator;
		return d;
		*/
		
	}

}

double ggx_distribution(Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double nh_dot = dot(surface_normal, microfacet_normal);
	if(nh_dot <= 0.0) return 0.0;
	else
	{
		double a_sq = lobe_width * lobe_width;
		double cos_h_sq = nh_dot * nh_dot;
		double cos_h_quart = cos_h_sq * cos_h_sq;
		double tan_h_sq = (1.0/cos_h_sq) - 1.0;

		double d = a_sq / (PI * cos_h_quart * (a_sq + tan_h_sq) * (a_sq + tan_h_sq));
		return d;
	}
}

double ggx_g_1(Vec3 v, Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double vh_dot = dot(v, microfacet_normal);
	double vn_dot = dot(v, surface_normal);
	double vnh_dot_quot = abs(vh_dot / vn_dot);

	if(vnh_dot_quot <= 0.0) return 0.0;
	else
	{
		double a_sq = lobe_width * lobe_width;
		double tan_vn_sq = (1.0/(vn_dot*vn_dot)) - 1.0;
		double d = 2.0 / (1.0 + sqrt(1.0 + a_sq*tan_vn_sq));
		return d;
	}
}

//Cite: Microfacet models for refraction through rough surfaces
//TODO: Make these maths variable names more consistent/better thought out
double g_1(Vec3 v, Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double v_h_dot = dot(v, microfacet_normal);
	double v_n_dot = dot(v, surface_normal);
	double v_h_n_dot_quot = abs(v_h_dot / v_n_dot);
	double tan_v_n = sqrt((1.0/(v_n_dot*v_n_dot)) - 1.0);
	double a = 1.0/(lobe_width * tan_v_n);
	
	if(v_h_n_dot_quot <= 0.0) return 0.0;
	else if(a < 1.6)
	{
		double d = v_h_n_dot_quot * (3.535*a + 2.181*a*a)/(1.0 + 2.276*a + 2.577*a*a);
		return d;
	}
	else return 1.0;
}

double geometric_attenuation(Vec3 incoming, Vec3 outgoing, Vec3 surface_normal, Vec3 microfacet_normal, double lobe_width)
{
	double g = ggx_g_1(outgoing, surface_normal, microfacet_normal, lobe_width) * ggx_g_1(incoming, surface_normal, microfacet_normal, lobe_width);
	return g;
}


/*	BDSFs	*/


void diffuse_phong_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& reflectance, Spectrum&)
{
	Spectrum* diffuse_spd = TEXTURE_SAMPLE(Spectrum, p.surface_material->diffuse_spd_texture, p.texture_coordinates);
	spectral_multiply(*diffuse_spd, 1.0/PI, reflectance);
}

void glossy_phong_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& reflectance, Spectrum&)
{
	Vec3 bisector = (incoming + outgoing)/2.0;
	double* shininess = TEXTURE_SAMPLE(double, p.surface_material->shininess_texture, p.texture_coordinates);
	double specular_coefficient = pow(d_max(0.0, dot(p.normal, bisector)), *shininess);
	Spectrum* glossy_spd = TEXTURE_SAMPLE(Spectrum, p.surface_material->glossy_spd_texture, p.texture_coordinates);
	spectral_multiply(*glossy_spd, specular_coefficient, reflectance);
}

void perfect_specular_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& reflectance, Spectrum&)
{
	if(incoming == reflect_vector(-outgoing, p.normal))
	{
		set_spectrum_to_value(reflectance, 1.0);
	}
	else set_spectrum_to_value(reflectance, 0.0);
}

//TODO: Make it so specular bdsfs don't have to check if incoming is specular reflection vector
void fresnel_specular_reflection_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& reflectance, Spectrum&)
{
	if(incoming == reflect_vector(-outgoing, p.normal))
	{
		double incident_cos = abs(dot(outgoing, p.normal));
		double reflectance_cos = abs(dot(incoming, p.normal));
		Spectrum& incident_refract_index = *TEXTURE_SAMPLE(Spectrum, p.incident_material->refract_index_texture, p.texture_coordinates);
		Spectrum& transmit_refract_index = *TEXTURE_SAMPLE(Spectrum, p.transmit_material->refract_index_texture, p.texture_coordinates);

		if(p.surface_material->type == MAT_TYPE_CONDUCTOR) 
		{
			Spectrum& transmit_extinct_index = *TEXTURE_SAMPLE(Spectrum, p.transmit_material->extinct_index_texture, p.texture_coordinates);
			fresnel_reflectance_conductor(incident_refract_index, transmit_refract_index, transmit_extinct_index, incident_cos, reflectance);
			spectral_multiply(reflectance, 1.0 / reflectance_cos, reflectance);
		}
		else if(p.surface_material->type == MAT_TYPE_DIELECTRIC)
		{
			fresnel_reflectance_dielectric(incident_refract_index, transmit_refract_index, incident_cos, reflectance);
			spectral_multiply(reflectance, 1.0 / reflectance_cos, reflectance);
		}
	}
}


void fresnel_transmission_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& transmittance, Spectrum&)
{
	Spectrum& incident_refract_index = *TEXTURE_SAMPLE(Spectrum, p.incident_material->refract_index_texture, p.texture_coordinates);
	Spectrum& transmit_refract_index = *TEXTURE_SAMPLE(Spectrum, p.transmit_material->refract_index_texture, p.texture_coordinates);
	if(incoming == transmit_vector(incident_refract_index, transmit_refract_index, p.normal, outgoing))
	{
		double incident_cos = abs(dot(outgoing, p.normal));
		double transmit_cos = abs(dot(incoming, p.normal));
		fresnel_transmittance_dielectric(incident_refract_index, transmit_refract_index, incident_cos, transmit_cos, transmittance);
		spectral_multiply(transmittance, 1.0 / transmit_cos, transmittance);
	}
}

void fresnel_reflection_transmission_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& f, Spectrum& temp_result)
{
	fresnel_specular_reflection_bdsf(p, incoming, outgoing, f, f);
	fresnel_transmission_bdsf(p, incoming, outgoing, temp_result, temp_result);
	spectral_sum(f, temp_result, f);
}

//MATERIAL BDSFs
void plastic_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& reflectance, Spectrum& temp_result)
{
	if(dot(p.normal, incoming) > 0.0)
	{
		diffuse_phong_bdsf(p, incoming, outgoing, reflectance, temp_result);
		glossy_phong_bdsf(p, incoming, outgoing, temp_result, temp_result);
		spectral_sum(reflectance, temp_result, reflectance);
	}
}

void mirror_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& reflectance, Spectrum&)
{
	if(dot(p.normal, incoming) > 0.0)
	{
		perfect_specular_bdsf(p, incoming, outgoing, reflectance, reflectance);
	}
}

//Cite: Microfacet models for refraction through rough surfaces
void cook_torrance_reflectance_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& fr, Spectrum&)
{
	Vec3 microfacet_normal = normalise(outgoing + incoming);
	Vec3 surface_normal = p.normal;
	double cos_sn_out = abs(dot(outgoing, surface_normal));
	double cos_sn_in = abs(dot(incoming, surface_normal));

	double cos_mn_out = abs(dot(outgoing, microfacet_normal));
	double cos_mn_in = abs(dot(incoming, microfacet_normal));

	Spectrum& incident_refract_index = *TEXTURE_SAMPLE(Spectrum, p.incident_material->refract_index_texture, p.texture_coordinates);
	Spectrum& transmit_refract_index = *TEXTURE_SAMPLE(Spectrum, p.transmit_material->refract_index_texture, p.texture_coordinates);
	Spectrum& transmit_extinct_index = *TEXTURE_SAMPLE(Spectrum, p.transmit_material->extinct_index_texture, p.texture_coordinates);
	if(p.surface_material->type == MAT_TYPE_CONDUCTOR)
	{
		fresnel_reflectance_conductor(incident_refract_index, transmit_refract_index, transmit_extinct_index, cos_mn_out, fr);
		spectral_multiply(fr, 1.0 / cos_mn_in, fr);
	}
	else if(p.surface_material->type == MAT_TYPE_DIELECTRIC)
	{
		fresnel_reflectance_dielectric(incident_refract_index, transmit_refract_index, cos_mn_out, fr);
		spectral_multiply(fr, 1.0 / cos_mn_in, fr);
	}

	double roughness = *TEXTURE_SAMPLE(double, p.surface_material->roughness_texture, p.texture_coordinates);
	double reflectance_coefficient = ggx_distribution(surface_normal, microfacet_normal, roughness) * geometric_attenuation(incoming, outgoing, surface_normal, microfacet_normal, roughness);
	spectral_multiply(fr, reflectance_coefficient, fr);
	reflectance_coefficient = 1.0/(4.0 * cos_sn_out * cos_sn_in);
	spectral_multiply(fr, reflectance_coefficient, fr);
}

void torrance_sparrow_bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& fr, Spectrum&)
{
	double cos_th_out = dot(outgoing, p.normal);
	double cos_th_in = dot(incoming, p.normal);

	double incident_cos = abs(cos_th_out);
	double reflectance_cos = abs(cos_th_in);

	Spectrum& incident_refract_index = *TEXTURE_SAMPLE(Spectrum, p.incident_material->refract_index_texture, p.texture_coordinates);
	Spectrum& transmit_refract_index = *TEXTURE_SAMPLE(Spectrum, p.transmit_material->refract_index_texture, p.texture_coordinates);
	Spectrum& transmit_extinct_index = *TEXTURE_SAMPLE(Spectrum, p.transmit_material->extinct_index_texture, p.texture_coordinates);
	if(p.surface_material->type == MAT_TYPE_CONDUCTOR)
	{
		fresnel_reflectance_conductor(incident_refract_index, transmit_refract_index, transmit_extinct_index, incident_cos, fr);
		spectral_multiply(fr, 1.0 / reflectance_cos, fr);
	}
	else if(p.surface_material->type == MAT_TYPE_DIELECTRIC)
	{
		fresnel_reflectance_dielectric(incident_refract_index, transmit_refract_index, incident_cos, fr);
		spectral_multiply(fr, 1.0 / reflectance_cos, fr);
	}

	Vec3 halfway = (outgoing + incoming)/2.0;
	double reflectance_coefficient = microfacet_distribution(p, halfway) * geometric_attenuation(p, incoming, outgoing);
	spectral_multiply(fr, reflectance_coefficient, fr);
	reflectance_coefficient = 1.0/(4.0 * cos_th_in, cos_th_out);
	spectral_multiply(fr, reflectance_coefficient, fr);
}

//GENERAL BDSF METHOD
void bdsf(Surface_Point& p, Vec3 incoming, Vec3 outgoing, Spectrum& reflectance)
{
	TIMED_FUNCTION;
	set_spectrum_to_zero(reflectance);
	Spectrum spare_spd;
	Spectrum bdsf_result;
	BDSF* material_bdsfs = p.surface_material->bdsfs;
	for(int i = 0; i < p.surface_material->number_of_bdsfs; ++i)
	{
		set_spectrum_to_zero(bdsf_result);
		set_spectrum_to_zero(spare_spd);
		material_bdsfs[i].bdsf(p, incoming, outgoing, bdsf_result, spare_spd);
		spectral_sum(reflectance, bdsf_result, reflectance);
	}
}


/*	Direction sampling	*/


double diffuse_pdf(Surface_Point& p, Vec3 outgoing, Vec3 sampled)
{
	return cos_weighted_sample_hemisphere_pdf(p.normal, sampled);
}

Vec3 sample_diffuse_direction(Surface_Point& p, Vec3 outgoing, double* pdf_value)
{
	return cos_weighted_sample_hemisphere(p.normal, pdf_value);
}

//NOTE: Although these two sample functions are the same now, the glossy one may be changed soon
Vec3 sample_glossy_direction(Surface_Point& p, Vec3 outgoing, double* pdf_value)
{
	/*
	*pdf_value = 1.0;
	Vec3 dir = reflect_vector(-outgoing, p.normal);
	return dir;
	*/
	return sample_diffuse_direction(p, outgoing, pdf_value);
}

Vec3 sample_cook_torrance_reflection_direction(Surface_Point& p, Vec3 outgoing, double* pdf_value)
{
	double f = uniform_sample();
	double g = uniform_sample();

	double a = *TEXTURE_SAMPLE(double, p.surface_material->roughness_texture, p.texture_coordinates);
	double a_sq = a * a;
#if 0
	double l = log(1.0 - f);
	if(1.0 - f <= 0.0) l = 0.0;

	double tan_th_mn_sq = -a_sq * l;
#else
	double tan_th_mn = (a * sqrt(f)) / sqrt(1.0 - f);
	double tan_th_mn_sq = tan_th_mn * tan_th_mn;
#endif
	double cos_th_mn = 1.0 / sqrt(1.0 + tan_th_mn_sq);
	double sin_th_mn = sqrt(1.0 - cos_th_mn * cos_th_mn);
	double phi_mn = 2.0 * PI * g;

	Vec3 microfacet_normal = {sin_th_mn * cos(phi_mn), sin_th_mn * sin(phi_mn), cos_th_mn};
	Mat3x3 r = find_rotation_between_vectors(Vec3{0.0, 0.0, 1.0}, p.normal);
	microfacet_normal = r * microfacet_normal;
	double mn_sn_dot = dot(microfacet_normal, p.normal);
	if(mn_sn_dot < 0.0) 
	{
		microfacet_normal = -microfacet_normal;
		mn_sn_dot = -mn_sn_dot;
	}
	
	Vec3 reflection_direction = reflect_vector(-outgoing, microfacet_normal);
	//*pdf_value = dot(outgoing, microfacet_normal) * geometric_attenuation(reflection_direction, outgoing, p.normal, microfacet_normal, p.material.roughness) / (dot(outgoing, p.normal) * dot(microfacet_normal, p.normal));
	double d = ggx_distribution(p.normal, microfacet_normal, a) * mn_sn_dot;

	*pdf_value = d * (1.0/(4.0 * dot(outgoing, microfacet_normal)));

	return reflection_direction;
}

Vec3 sample_torrance_sparrow_direction(Surface_Point& p, Vec3 outgoing, double* pdf_value)
{
	double f = uniform_sample();
	double g = uniform_sample();

	double l = log(1.0 - f);
	if(1.0 - f <= 0.0) l = 0.0;
	
	double r = *TEXTURE_SAMPLE(double, p.surface_material->roughness_texture, p.texture_coordinates);
	double tan_th_sq = - r * r * l;
	double phi = 2.0 * PI * g;

	double cos_th = 1.0 / sqrt(1.0 + tan_th_sq);
	double sin_th = sqrt(1.0 - cos_th * cos_th);
	
	Vec3 halfway = {sin_th * cos(phi), sin_th * sin(phi), cos_th};
	if(dot(halfway, p.normal) < 0.0) halfway = -halfway;
	*pdf_value = microfacet_distribution(p, halfway) * geometric_attenuation(p, outgoing) * abs(dot(outgoing, halfway)) / abs(dot(outgoing, p.normal));

	return reflect_vector(-outgoing, halfway);
}

Vec3 sample_specular_direction(Surface_Point& p, Vec3 outgoing, double* pdf_value)
{
	*pdf_value = 1.0;
	return reflect_vector(-outgoing, p.normal);
}

Vec3 sample_specular_transmission_direction(Surface_Point& p, Vec3 outgoing, double* pdf_value)
{
	*pdf_value = 1.0;
	Spectrum& incident_refract_index = *TEXTURE_SAMPLE(Spectrum, p.incident_material->refract_index_texture, p.texture_coordinates);
	Spectrum& transmit_refract_index = *TEXTURE_SAMPLE(Spectrum, p.transmit_material->refract_index_texture, p.texture_coordinates);
	return transmit_vector(incident_refract_index, transmit_refract_index, p.normal, outgoing);
}


//TODO: Total internal reflection
//TODO: Fix this function
//	- fresnel_reflectance should return 1.0 with total internal reflection
Vec3 sample_specular_reflection_or_transmission_direction(Surface_Point& p, Vec3 outgoing, double* pdf_value)
{
	Spectrum& incident_refract_index = *TEXTURE_SAMPLE(Spectrum, p.incident_material->refract_index_texture, p.texture_coordinates);
	Spectrum& transmit_refract_index = *TEXTURE_SAMPLE(Spectrum, p.transmit_material->refract_index_texture, p.texture_coordinates);
	double incident_refract_index_d = spd_value_at_wavelength(incident_refract_index, TRANSMISSION_WAVELENGTH);
	double transmit_refract_index_d = spd_value_at_wavelength(transmit_refract_index, TRANSMISSION_WAVELENGTH);
	double incident_cos = abs(dot(outgoing, p.normal));

	double reflectance = fresnel_reflectance_dielectric(incident_refract_index_d, transmit_refract_index_d, incident_cos);

	double f = uniform_sample();

	Vec3 sampled_direction = {};
	if(f < reflectance)
	{
		sampled_direction = reflect_vector(-outgoing, p.normal);
		*pdf_value = reflectance;
	}
	else
	{
		sampled_direction = transmit_vector(incident_refract_index_d, transmit_refract_index_d, p.normal, outgoing);
		*pdf_value = 1.0 - reflectance;
	}
	
	return sampled_direction;
}


/*	Create materials	*/


Material create_plastic(Spectrum diffuse_spd, Spectrum glossy_spd, double shininess)
{
	Material plastic = {};
	//Reflectances
	plastic.type = MAT_TYPE_DIELECTRIC;
	plastic.diffuse_spd_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(plastic.diffuse_spd_texture, diffuse_spd);
	plastic.glossy_spd_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(plastic.glossy_spd_texture, glossy_spd);
	plastic.shininess_texture = TEXTURE_CREATE(double, 1, 1);
	TEXTURE_CLEAR(plastic.shininess_texture, shininess);
	
	//BDSFs
	plastic.number_of_bdsfs = 2;
	BDSF diffuse_reflection = {};
	diffuse_reflection.type = BDSF_TYPE_DIFFUSE;
	diffuse_reflection.pdf = diffuse_pdf;
	diffuse_reflection.bdsf = diffuse_phong_bdsf;
	diffuse_reflection.sample_direction = sample_diffuse_direction;

	BDSF glossy_reflection = {};
	glossy_reflection.type = BDSF_TYPE_DIFFUSE;
	glossy_reflection.pdf = diffuse_pdf;
	glossy_reflection.bdsf = glossy_phong_bdsf;
	glossy_reflection.sample_direction = sample_glossy_direction;
	plastic.bdsfs[0] = diffuse_reflection;
	plastic.bdsfs[1] = glossy_reflection;

	return plastic;
}

Material create_textured_plastic(Texture diffuse_spd_texture, Spectrum glossy_spd, double shininess)
{
	Material plastic = {};
	plastic.type = MAT_TYPE_DIELECTRIC;
	plastic.diffuse_spd_texture = diffuse_spd_texture;
	plastic.glossy_spd_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(plastic.glossy_spd_texture, glossy_spd);
	plastic.shininess_texture = TEXTURE_CREATE(double, 1, 1);
	TEXTURE_CLEAR(plastic.shininess_texture, shininess);

	plastic.number_of_bdsfs = 2;
	BDSF diffuse_reflection = {};
	diffuse_reflection.type = BDSF_TYPE_DIFFUSE;
	diffuse_reflection.pdf = diffuse_pdf;
	diffuse_reflection.bdsf = diffuse_phong_bdsf;
	diffuse_reflection.sample_direction = sample_diffuse_direction;

	BDSF glossy_reflection = {};
	glossy_reflection.type = BDSF_TYPE_DIFFUSE;
	glossy_reflection.pdf = diffuse_pdf;
	glossy_reflection.bdsf = glossy_phong_bdsf;
	glossy_reflection.sample_direction = sample_glossy_direction;
	plastic.bdsfs[0] = diffuse_reflection;
	plastic.bdsfs[1] = glossy_reflection;

	return plastic;
}

Material create_mirror()
{
	Material mirror = {};
	mirror.type = MAT_TYPE_CONDUCTOR;
	mirror.number_of_bdsfs = 1;
	mirror.bdsfs[0].type = BDSF_TYPE_SPECULAR;
	mirror.bdsfs[0].bdsf = perfect_specular_bdsf;
	mirror.bdsfs[0].sample_direction = sample_specular_direction;

	return mirror;
}

Material create_conductor(Spectrum refract_index, Spectrum extinct_index, double roughness)
{
	Material conductor = {};
	conductor.type = MAT_TYPE_CONDUCTOR;
	conductor.number_of_bdsfs = 1;

	if(roughness >= 0.001)
	{
		conductor.bdsfs[0].type = BDSF_TYPE_SPECULAR;
		conductor.bdsfs[0].bdsf = cook_torrance_reflectance_bdsf;
		conductor.bdsfs[0].sample_direction = sample_cook_torrance_reflection_direction;
	}
	else
	{
		conductor.bdsfs[0].type = BDSF_TYPE_SPECULAR;
		conductor.bdsfs[0].bdsf = fresnel_specular_reflection_bdsf;
		conductor.bdsfs[0].sample_direction = sample_specular_direction;
	}

	conductor.refract_index_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(conductor.refract_index_texture, refract_index);
	conductor.extinct_index_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(conductor.extinct_index_texture, extinct_index);

	conductor.roughness_texture = TEXTURE_CREATE(double, 1, 1);
	TEXTURE_CLEAR(conductor.roughness_texture, roughness);

	return conductor;
}

Material create_dielectric(Spectrum refract_index) 
{
	Material dielectric = {};
	dielectric.type = MAT_TYPE_DIELECTRIC;
	dielectric.number_of_bdsfs = 1;
	dielectric.bdsfs[0].type = BDSF_TYPE_SPECULAR;
	dielectric.bdsfs[0].bdsf = fresnel_reflection_transmission_bdsf;
	dielectric.bdsfs[0].sample_direction = sample_specular_reflection_or_transmission_direction;

	dielectric.bdsfs[1].type = BDSF_TYPE_DIFFUSE;
	dielectric.bdsfs[1].pdf = diffuse_pdf;
	dielectric.bdsfs[1].bdsf = glossy_phong_bdsf;
	dielectric.bdsfs[1].sample_direction = sample_glossy_direction;

	dielectric.glossy_spd_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	Spectrum glossy_spd = generate_constant_spd(1.0);
	TEXTURE_CLEAR(dielectric.glossy_spd_texture, glossy_spd);
	dielectric.shininess_texture = TEXTURE_CREATE(double, 1, 1);
	double shininess = 16.0;
	TEXTURE_CLEAR(dielectric.shininess_texture, shininess);

	dielectric.refract_index_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(dielectric.refract_index_texture, refract_index);

	return dielectric;
}
