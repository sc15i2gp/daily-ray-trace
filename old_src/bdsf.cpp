//reflectance functions
void mirror_reflectance(Surface_Point* p, Vec3 in_direction, Vec3 out_direction, Spectrum* reflectance)
{
	if(out_direction == reflect_vector(-in_direction, p->normal)) set_spectrum_to_value(reflectance, 1.0);
	else set_spectrum_to_value(reflectance, 0.0);
}

void diffuse_phong_reflectance(Surface_Point* p, Vec3 in_direction, Vec3 out_direction, Spectrum* reflectance)
{
	spectral_multiply(p->diffuse_spd, 1.0/PI, reflectance);
}

void glossy_phong_reflectance(Surface_Point* p, Vec3 in_direction, Vec3 out_direction, Spectrum* reflectance)
{
	Vec3 bisector = (in_direction + out_direction)/2.0;
	double shininess = p->shininess;
	double specular_coefficient = pow(d_max(0.0, dot(p->normal, bisector)), p->shininess);
	spectral_multiply(p->glossy_spd, specular_coefficient, reflectance);
}

//NEW bdsf pdf values
double const_1_pdf(Surface_Point* p, Vec3 in_direction, Vec3 out_direction)
{
	return 1.0;
}

double cos_weighted_hemisphere_pdf(Surface_Point* p, Vec3 in_direction, Vec3 out_direction)
{
	return dot(p->normal, out_direction) / PI;
}

// bdsf direction sampling

Vec3 sample_camera_lens_direction(Surface_Point* p, Vec3 in_direction)
{
	Vec3 plane_of_focus_point = p->position + p->camera->focal_depth * normalise(p->camera->pinhole_position - p->position);
	Mat3x3 r = find_rotation_between_vectors(p->camera->film_normal, Vec3{0.0, 0.0, 1.0});
	Vec3 sampled_lens_point = p->camera->pinhole_position + r * ((p->camera->focal_length / p->camera->aperture_radius) * uniform_sample_disc());
	return normalise(p->camera->pinhole_position - p->position);
}

Vec3 sample_camera_pinhole_direction(Surface_Point* p, Vec3 in_direction)
{
	return normalise(p->camera->pinhole_position - p->position);
}

Vec3 sample_diffuse_direction(Surface_Point* p, Vec3 in_direction)
{
	Vec3 v = {};
	for(;;)
	{
		v = uniform_sample_disc();
		if(dot(v, v) < 1.0) break;
	}
	v.z = sqrt(1.0 - dot(v, v));
	//Rotate v by R such that R*(0, 0, 1) = normal
	Mat3x3 r = find_rotation_between_vectors(Vec3{0.0, 0.0, 1.0}, p->normal);
	return r * v;
}

Vec3 sample_glossy_direction(Surface_Point* p, Vec3 in_direction)
{
	return sample_diffuse_direction(p, in_direction);
}

Vec3 sample_specular_reflection_direction(Surface_Point* p, Vec3 in_direction)
{
	return reflect_vector(-in_direction, p->normal);
}
