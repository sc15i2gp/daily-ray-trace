
/*
bool ray_intersects_model(Ray ray, Model& m, double* t, int* ret_intersecting_triangle_index, Vec3& ret_vec)
{
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;

	double min_distance_to_triangle = DBL_MAX;
	int intersecting_triangle_index = -1;
	for(int i = 0; i < m.number_of_vertices; i += 3)
	{
		Vec3 triangle_normal = (m.vertices[i].normal + m.vertices[i+1].normal + m.vertices[i+2].normal)/3.0;
		if(dot(ray.direction, triangle_normal) != 0.0)
		{
			double distance_to_triangle = dot(m.vertices[i].position - ray.origin, triangle_normal)/dot(ray.direction, triangle_normal);
			if(distance_to_triangle < min_distance_to_triangle && distance_to_triangle > 0.0)
			{
				Vec3 p = ray.origin + distance_to_triangle * ray.direction;

				double temp_a = sqrt(point_to_line_distance_sq(p, m.vertices[i+1].position, m.vertices[i+2].position)/point_to_line_distance_sq(m.vertices[i].position, m.vertices[i+1].position,m.vertices[i+2].position));
				double temp_b = sqrt(point_to_line_distance_sq(p, m.vertices[i+2].position, m.vertices[i].position)/point_to_line_distance_sq(m.vertices[i+1].position, m.vertices[i+2].position, m.vertices[i].position));
				double temp_c = sqrt(point_to_line_distance_sq(p, m.vertices[i].position, m.vertices[i+1].position)/point_to_line_distance_sq(m.vertices[i+2].position, m.vertices[i].position, m.vertices[i+1].position));

				double d = temp_a + temp_b + temp_c;
				if(d <= 1.00001)
				{
					intersecting_triangle_index = i;
					min_distance_to_triangle = distance_to_triangle;
					a = temp_a;
					b = temp_b;
					c = temp_c;
				}
			}
		}
	}
	ret_vec = {a, b, c};

	*t = min_distance_to_triangle;
	*ret_intersecting_triangle_index = intersecting_triangle_index;
	return intersecting_triangle_index >= 0;
}

Model create_model()
{
	Model model = {};
	model.number_of_vertices = 12;
	model.vertices = (Model_Vertex*)alloc(model.number_of_vertices * sizeof(Model_Vertex));

	Vec3 position_0 = {1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0)};
	Vec3 position_1 = {-1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0)};
	Vec3 position_2 = {0.0, 2.0/sqrt(3.0), -1.0/sqrt(6.0)};
	Vec3 position_3 = {0.0, 0.0, 3.0/sqrt(6.0)};

	Vec3 normal_0 = normalise(cross(position_2 - position_0, position_1 - position_0));
	Vec3 normal_1 = normalise(cross(position_3 - position_1, position_2 - position_1));
	Vec3 normal_2 = normalise(cross(position_2 - position_3, position_0 - position_3));
	Vec3 normal_3 = normalise(cross(position_1 - position_0, position_3 - position_0));

	int v = 0;
	//Back triangle
	model.vertices[v++] = {position_0, normal_0, {0.5f, 1.0f}};
	model.vertices[v++] = {position_1, normal_0, {0.0f, 0.0f}};
	model.vertices[v++] = {position_2, normal_0, {1.0f, 0.0f}};

	//Left triangle
	model.vertices[v++] = {position_1, normal_1, {0.5f, 1.0f}};
	model.vertices[v++] = {position_3, normal_1, {0.0f, 0.0f}};
	model.vertices[v++] = {position_2, normal_1, {1.0f, 0.0f}};

	//Right triangle
	model.vertices[v++] = {position_3, normal_2, {0.5f, 1.0f}};
	model.vertices[v++] = {position_0, normal_2, {0.0f, 0.0f}};
	model.vertices[v++] = {position_2, normal_2, {1.0f, 0.0f}};

	//Bottom triangle
	model.vertices[v++] = {position_0, normal_3, {0.5f, 1.0f}};
	model.vertices[v++] = {position_3, normal_3, {0.0f, 0.0f}};
	model.vertices[v++] = {position_1, normal_3, {1.0f, 0.0f}};


	return model;

}

Texture create_default_spd_texture()
{
	Texture texture = TEXTURE_CREATE(Spectrum, 32, 32);
	Spectrum purple = RGB64_to_spectrum(RGB64{0.875, 0.0, 0.996});
	Spectrum black = {};
	
	for(int y = 0; y < texture.height; ++y)
	{
		for(int x = 0; x < texture.width; ++x)
		{
			Spectrum* pixel = (Spectrum*)(get_pixel(texture, x, y));
			*pixel = ((x/8) % 2 == (y/8) % 2) ? purple : black;
		}
	}

	return texture;
}

void load_scene(Scene* scene)
{
	//Spectrum light_spd = generate_black_body_spd(4000.0);
	//normalise(light_spd);
	Material mirror = create_mirror();
	Spectrum gold_refract_index = load_spd("spectra/au_spec_n.csv");
	Spectrum gold_extinct_index = load_spd("spectra/au_spec_k.csv");
	Material gold = create_conductor(gold_refract_index, gold_extinct_index, 0.344);
	Spectrum glass_refract_index = load_spd("spectra/glass.csv");
	Material glass = create_dielectric(glass_refract_index);

}
*/

void NEW_RGB64_to_spectrum(RGB64 c, Spectrum* dst, NEW_Spectrum_Buffer* spectrum_buffer)
{
	NEW_RGB64_to_spectrum(c, dst, &spectrum_buffer->spectra[0], 
			&spectrum_buffer->spectra[1], &spectrum_buffer->spectra[2], &spectrum_buffer->spectra[3], 
			&spectrum_buffer->spectra[4], &spectrum_buffer->spectra[5], &spectrum_buffer->spectra[6]);
}

void NEW_load_scene(NEW_Scene* scene, NEW_Spectrum_Buffer* spectrum_buffer)
{
	NEW_load_spd("spectra/white_rgb_to_spd.csv", &spectrum_buffer->spectra[0]);
	NEW_load_spd("spectra/red_rgb_to_spd.csv", &spectrum_buffer->spectra[1]);
	NEW_load_spd("spectra/blue_rgb_to_spd.csv", &spectrum_buffer->spectra[2]);
	NEW_load_spd("spectra/green_rgb_to_spd.csv", &spectrum_buffer->spectra[3]);
	NEW_load_spd("spectra/cyan_rgb_to_spd.csv", &spectrum_buffer->spectra[4]);
	NEW_load_spd("spectra/magenta_rgb_to_spd.csv", &spectrum_buffer->spectra[5]);
	NEW_load_spd("spectra/yellow_rgb_to_spd.csv", &spectrum_buffer->spectra[6]);

	//Geometries
	double h = 3.0;
	scene->planes[0] = create_plane_from_points(Vec3{-h, h, -h}, Vec3{h, h, -h}, Vec3{-h, -h, -h}); //back wall
	scene->planes[1] = create_plane_from_points(Vec3{-h, h, h}, Vec3{-h, h, -h}, Vec3{-h, -h, h}); //left wall
	scene->planes[2] = create_plane_from_points(Vec3{h, h, -h}, Vec3{h, h, h}, Vec3{h, -h, -h}); //right wall
	scene->planes[3] = create_plane_from_points(Vec3{-h, -h, -h}, Vec3{h, -h, -h}, Vec3{-h, -h, h}); //floor
	scene->planes[4] = create_plane_from_points(Vec3{-h, h, h}, Vec3{h, h, h}, Vec3{-h, h, -h});  //ceiling
	scene->planes[5] = create_plane_from_points(Vec3{-0.5, 2.9, 0.5}, Vec3{0.5, 2.9, 0.5}, Vec3{-0.5, 2.9, -0.5}); //Light

	scene->number_of_geometries = 6;
	scene->geometries[0].type = GEO_TYPE_PLANE;
	scene->geometries[0].plane = &scene->planes[0];
	scene->geometries[1].type = GEO_TYPE_PLANE;
	scene->geometries[1].plane = &scene->planes[1];
	scene->geometries[2].type = GEO_TYPE_PLANE;
	scene->geometries[2].plane = &scene->planes[2];
	scene->geometries[3].type = GEO_TYPE_PLANE;
	scene->geometries[3].plane = &scene->planes[3];
	scene->geometries[4].type = GEO_TYPE_PLANE;
	scene->geometries[4].plane = &scene->planes[4];
	scene->geometries[5].type = GEO_TYPE_PLANE;
	scene->geometries[5].plane = &scene->planes[5];

	//Materials
	scene->number_of_materials = 5;
	//back wall (blue)
	scene->materials[0] = {};
	scene->materials[0].shininess = 32.0;
	NEW_RGB64_to_spectrum(RGB64{0.2, 0.2, 0.8}, &scene->materials[0].diffuse_spd, spectrum_buffer);
	NEW_RGB64_to_spectrum(RGB64{0.8, 0.8, 0.9}, &scene->materials[0].glossy_spd, spectrum_buffer);
	scene->materials[0].number_of_bdsfs = 2;
	scene->materials[0].bdsfs[0].reflectance = NEW_diffuse_phong_reflectance;
	scene->materials[0].bdsfs[0].sample_direction = NEW_sample_diffuse_direction;
	scene->materials[0].bdsfs[0].pdf = NEW_cos_weighted_hemisphere_pdf;
	scene->materials[0].bdsfs[1].reflectance = NEW_glossy_phong_reflectance;
	scene->materials[0].bdsfs[1].sample_direction = NEW_sample_glossy_direction;
	scene->materials[0].bdsfs[1].pdf = NEW_cos_weighted_hemisphere_pdf;

	//left wall (green)
	scene->materials[1] = {};
	scene->materials[1].shininess = 32.0;
	NEW_RGB64_to_spectrum(RGB64{0.1, 0.35, 0.1}, &scene->materials[1].diffuse_spd, spectrum_buffer);
	NEW_RGB64_to_spectrum(RGB64{0.45, 0.55, 0.45}, &scene->materials[1].glossy_spd, spectrum_buffer);
	scene->materials[1].number_of_bdsfs = 2;
	scene->materials[1].bdsfs[0].reflectance = NEW_diffuse_phong_reflectance;
	scene->materials[1].bdsfs[0].sample_direction = NEW_sample_diffuse_direction;
	scene->materials[1].bdsfs[0].pdf = NEW_cos_weighted_hemisphere_pdf;
	scene->materials[1].bdsfs[1].reflectance = NEW_glossy_phong_reflectance;
	scene->materials[1].bdsfs[1].sample_direction = NEW_sample_glossy_direction;
	scene->materials[1].bdsfs[1].pdf = NEW_cos_weighted_hemisphere_pdf;
	
	//right wall (red)
	scene->materials[2] = {};
	scene->materials[2].shininess = 32.0;
	NEW_RGB64_to_spectrum(RGB64{0.5, 0.0, 0.0}, &scene->materials[2].diffuse_spd, spectrum_buffer);
	NEW_RGB64_to_spectrum(RGB64{0.7, 0.6, 0.6}, &scene->materials[2].glossy_spd, spectrum_buffer);
	scene->materials[2].number_of_bdsfs = 2;
	scene->materials[2].bdsfs[0].reflectance = NEW_diffuse_phong_reflectance;
	scene->materials[2].bdsfs[0].sample_direction = NEW_sample_diffuse_direction;
	scene->materials[2].bdsfs[0].pdf = NEW_cos_weighted_hemisphere_pdf;
	scene->materials[2].bdsfs[1].reflectance = NEW_glossy_phong_reflectance;
	scene->materials[2].bdsfs[1].sample_direction = NEW_sample_glossy_direction;
	scene->materials[2].bdsfs[1].pdf = NEW_cos_weighted_hemisphere_pdf;

	//floor + ceiling (white)
	scene->materials[3] = {};
	scene->materials[3].shininess = 32.0;
	NEW_RGB64_to_spectrum(RGB64{0.55, 0.55, 0.55}, &scene->materials[3].diffuse_spd, spectrum_buffer);
	NEW_RGB64_to_spectrum(RGB64{0.7, 0.7, 0.7}, &scene->materials[3].glossy_spd, spectrum_buffer);
	scene->materials[3].number_of_bdsfs = 2;
	scene->materials[3].bdsfs[0].reflectance = NEW_diffuse_phong_reflectance;
	scene->materials[3].bdsfs[0].sample_direction = NEW_sample_diffuse_direction;
	scene->materials[3].bdsfs[0].pdf = NEW_cos_weighted_hemisphere_pdf;
	scene->materials[3].bdsfs[1].reflectance = NEW_glossy_phong_reflectance;
	scene->materials[3].bdsfs[1].sample_direction = NEW_sample_glossy_direction;
	scene->materials[3].bdsfs[1].pdf = NEW_cos_weighted_hemisphere_pdf;

	//light
	scene->materials[4] = {};
	scene->materials[4].is_emissive = true;
	NEW_set_spectrum_to_value(&scene->materials[4].emission_spd, 1.0);

	//Set scene objects and light sources
	scene->number_of_objects = 6;

	//back wall
	scene->object_geometries[0] = &scene->geometries[0];
	scene->object_materials[0] = &scene->materials[0];

	//left wall
	scene->object_geometries[1] = &scene->geometries[1];
	scene->object_materials[1] = &scene->materials[1];
	
	//right wall
	scene->object_geometries[2] = &scene->geometries[2];
	scene->object_materials[2] = &scene->materials[2];
	
	//floor and ceiling
	scene->object_geometries[3] = &scene->geometries[3];
	scene->object_geometries[4] = &scene->geometries[4];
	scene->object_materials[3] = &scene->materials[3];
	scene->object_materials[4] = &scene->materials[3];

	//Light
	scene->object_geometries[5] = &scene->geometries[5];
	scene->object_materials[5] = &scene->materials[4];
	scene->number_of_light_sources = 1;
	scene->light_source_geometries[0] = &scene->geometries[5];
	scene->light_source_materials[0] = &scene->materials[4];

	//Load camera data
	scene->camera_position = {0.0, 0.0, 8.0};
	scene->camera_fov = 90.0;
	scene->camera_up = {0.0, 1.0, 0.0};
	scene->camera_right = {1.0, 0.0, 0.0};
	scene->camera.film_normal = {0.0, 0.0, -1.0};	
	scene->camera.focal_depth = 8.0;
	scene->camera.focal_length = 0.5;
	scene->camera.aperture_radius = 0.0; //NOTE: change for lens

	//film material
	scene->film_material = {};
	scene->film_material.number_of_bdsfs = 1;
	scene->film_material.bdsfs[0].sample_direction = (scene->camera.aperture_radius > 0.0) ? NEW_sample_camera_lens_direction : NEW_sample_camera_pinhole_direction;
	scene->film_material.bdsfs[0].reflectance = NEW_const_1_reflectance;
	scene->film_material.bdsfs[0].pdf = NEW_const_1_pdf; //TODO: Change this for lens (probably shouldn't be 1)
}


void NEW_sample_geometry(NEW_Scene_Geometry* geometry, Vec3 origin, NEW_Geometry_Sample* geometry_sample)
{
	switch(geometry->type)
	{
		case GEO_TYPE_SPHERE: geometry_sample->position = uniform_sample_sphere_subtended(*geometry->sphere, origin, &geometry_sample->pdf); break;
		case GEO_TYPE_PLANE: geometry_sample->position = uniform_sample_plane(*geometry->plane, &geometry_sample->pdf); break;
		//TODO: case GEO_TYPE_MODEL
	}
}

Vec3 NEW_find_ray_scene_intersection_point(NEW_Scene* scene, Ray ray, int* ret_intersecting_object = NULL)
{
	double min_length_along_ray = DBL_MAX;
	double length_along_ray = DBL_MAX;
	Vec3 min_intersection_along_ray = {DBL_MAX, DBL_MAX, DBL_MAX};
	int intersecting_object = -1;
	
	for(int i = 0; i < scene->number_of_geometries; ++i)
	{
		NEW_Scene_Geometry* geometry = scene->object_geometries[i];
		switch(geometry->type)
		{
			case GEO_TYPE_SPHERE: length_along_ray = NEW_ray_intersects_sphere(ray, *geometry->sphere); break;
			case GEO_TYPE_PLANE: length_along_ray = NEW_ray_intersects_plane(ray, *geometry->plane); break;
			//TODO: case GEO_TYPE_MODEL
		}
		if(length_along_ray > -0.0009765625 && length_along_ray < min_length_along_ray)
		{
			min_length_along_ray = length_along_ray;
			min_intersection_along_ray = ray.origin + min_length_along_ray * ray.direction;
			intersecting_object = i;
		}
	}
	if(ret_intersecting_object) *ret_intersecting_object = intersecting_object;

	return min_intersection_along_ray;
}

void NEW_find_ray_scene_intersection(NEW_Scene* scene, Ray ray, NEW_Scene_Point* intersection_point)
{
	int intersection_object;
	intersection_point->position = NEW_find_ray_scene_intersection_point(scene, ray, &intersection_object);
	intersection_point->object = intersection_object;
	if(intersection_object >= 0)
	{
		intersection_point->in_ray = ray;
		intersection_point->in_ray.direction = -intersection_point->in_ray.direction;

		switch(scene->object_geometries[intersection_object]->type)
		{
			case GEO_TYPE_SPHERE:
			{
				intersection_point->normal = normalise(intersection_point->position - scene->object_geometries[intersection_object]->sphere->center);
				break;
			}
			case GEO_TYPE_PLANE:
			{
				intersection_point->normal = scene->object_geometries[intersection_object]->plane->n;
				break;
			}
		}
		intersection_point->roughness = scene->object_materials[intersection_object]->roughness;
		//Set incident and transmission materials
		//For now, if entering an object air is the incident material, otherwise the object's
		double in_dot = dot(ray.direction, intersection_point->normal);
		intersection_point->incident_material = (in_dot < 0.0) ? &scene->air_material : scene->object_materials[intersection_object];
		intersection_point->transmit_material = (in_dot < 0.0) ? scene->object_materials[intersection_object] : &scene->air_material;
		intersection_point->material = scene->object_materials[intersection_object];
	}
	else
	{
		intersection_point->material = &scene->null_material;
	}
}

int NEW_sample_scene_path(NEW_Scene* scene, Vec3 film_point, NEW_Scene_Point* scene_path, int max_depth)
{
	int depth = 0;
	bool is_there_next_point = true;
	NEW_Scene_Point* current_point = &scene_path[0];
	NEW_Scene_Point* next_point = &scene_path[1];
	current_point->position = film_point;
	current_point->normal = scene->camera.film_normal;
	current_point->in_ray.direction = {1.0, 1.0, 0.0};
	current_point->in_ray.origin = {};
	current_point->material = &scene->film_material;
	current_point->incident_material = &scene->film_material;
	current_point->transmit_material = &scene->film_material;
	current_point->pdf = 1.0;
	current_point->object = -1;

	NEW_Surface_Point surface_point = {};

	while(is_there_next_point && depth < max_depth)
	{
		//Sample incoming light rays
		for(int i = 0; depth != 0 && i < scene->number_of_light_sources; ++i)
		{
			NEW_Geometry_Sample* light_sample = &current_point->light_contributions[current_point->number_of_light_contributions];
			NEW_sample_geometry(scene->light_source_geometries[i], current_point->position, light_sample);
			light_sample->object = i;
			
			Ray shadow_ray = {};
			shadow_ray.direction = normalise(light_sample->position - current_point->position);
			shadow_ray.origin = current_point->position + 0.03125*shadow_ray.direction;
			
			//Test if an object in the scene obscures the shadow ray
			int shadow_test_object;
			Vec3 shadow_test_point = NEW_find_ray_scene_intersection_point(scene, shadow_ray, &shadow_test_object);
			bool light_point_visible = true;
			if(shadow_test_object >= 0)
			{
				Vec3 shadow_test_ray = shadow_test_point - current_point->position;
				double shadow_test_ray_length_sq = dot(shadow_test_ray, shadow_test_ray);
				Vec3 light_test_ray = light_sample->position - current_point->position;
				double light_test_ray_length_sq = dot(light_test_ray, light_test_ray);
				light_point_visible = !(shadow_test_ray_length_sq - light_test_ray_length_sq < -0.0009765625);
			}

			if(light_point_visible) ++current_point->number_of_light_contributions;
		}

		//Sample next direction in the scene path
		double r = uniform_sample();
		double n = (double)current_point->material->number_of_bdsfs;
		int sampled_bdsf = (int)floor(r * n);

		surface_point.position = current_point->position;
		surface_point.normal = current_point->normal;
		surface_point.roughness = current_point->roughness;
		surface_point.incident_refract_index_spd = &current_point->incident_material->refract_index_spd;
		surface_point.transmit_refract_index_spd = &current_point->transmit_material->refract_index_spd;
		surface_point.incident_extinct_index_spd = &current_point->incident_material->extinct_index_spd;
		surface_point.transmit_extinct_index_spd = &current_point->transmit_material->extinct_index_spd;
		surface_point.camera = &scene->camera;

		current_point->out_ray.direction = current_point->material->bdsfs[sampled_bdsf].sample_direction(&surface_point, current_point->in_ray.direction);
		for(int i = 0; i < current_point->material->number_of_bdsfs; ++i)
		{
			next_point->pdf += current_point->material->bdsfs[i].pdf(&surface_point, current_point->in_ray.direction, current_point->out_ray.direction);
		}
		current_point->out_ray.origin = current_point->position + 0.03125*current_point->out_ray.direction;
		NEW_find_ray_scene_intersection(scene, current_point->out_ray, next_point);
		++depth;

		//Store and iterate points
		current_point = next_point;
		is_there_next_point = !next_point->material->is_emissive && next_point->object != -1;
		++next_point;
	}
	if(!is_there_next_point && depth < max_depth) 
	{
		//TODO: Fix the need for this due to light_sample object ints referring to light sources
		if(current_point->object != -1) current_point->object = 0;
		scene_path[depth-1].light_contributions[scene_path[depth-1].number_of_light_contributions++] = *((NEW_Geometry_Sample*)(&current_point->position));
	}
	
	return depth;
}

void NEW_load_surface_and_emission_materials(NEW_Scene* scene, NEW_Scene_Point* p, NEW_Surface_Point* surface_point, NEW_Spectrum_Buffer* spectrum_buffer)
{
	NEW_copy_spectrum(&p->incident_material->refract_index_spd, &spectrum_buffer->spectra[4]);
	NEW_copy_spectrum(&p->incident_material->extinct_index_spd, &spectrum_buffer->spectra[5]);
	NEW_copy_spectrum(&p->transmit_material->refract_index_spd, &spectrum_buffer->spectra[6]);
	NEW_copy_spectrum(&p->transmit_material->extinct_index_spd, &spectrum_buffer->spectra[7]);
	NEW_copy_spectrum(&p->material->glossy_spd, &spectrum_buffer->spectra[8]);
	NEW_copy_spectrum(&p->material->diffuse_spd, &spectrum_buffer->spectra[9]);

	for(int i = 0; i < p->number_of_light_contributions; ++i)
	{
		if(p->light_contributions[i].object >= 0)
		NEW_copy_spectrum(&scene->light_source_materials[p->light_contributions[i].object]->emission_spd, &spectrum_buffer->spectra[10+i]);
		else NEW_copy_spectrum(&scene->null_material.emission_spd, &spectrum_buffer->spectra[10+i]);
	}
	
	
	//Set surface point
	surface_point->normal = p->normal;
	surface_point->incident_refract_index_spd = &spectrum_buffer->spectra[4];
	surface_point->incident_extinct_index_spd = &spectrum_buffer->spectra[5];
	surface_point->transmit_refract_index_spd = &spectrum_buffer->spectra[6];
	surface_point->transmit_extinct_index_spd = &spectrum_buffer->spectra[7];
	surface_point->glossy_spd = &spectrum_buffer->spectra[8];
	surface_point->diffuse_spd = &spectrum_buffer->spectra[9];
	surface_point->roughness = p->material->roughness;
	surface_point->shininess = p->material->shininess;
}

void NEW_bdsf(NEW_Surface_Point* p, NEW_BDSF* bdsfs, int number_of_bdsfs, Vec3 incoming, Vec3 outgoing, Spectrum* reflectance, Spectrum* bdsf_result)
{
	NEW_set_spectrum_to_value(reflectance, 0.0);
	NEW_set_spectrum_to_value(bdsf_result, 0.0);
	for(int i = 0; i < number_of_bdsfs; ++i)
	{
		bdsfs[i].reflectance(p, incoming, outgoing, reflectance);
		NEW_spectral_sum(reflectance, bdsf_result, bdsf_result);
	}
}

void NEW_compute_path_radiance(NEW_Scene* scene, NEW_Spectrum_Buffer* spectrum_buffer, NEW_Scene_Point* scene_path, int depth)
{
	Spectrum* eye_ray_radiance = spectrum_buffer->spectra;
	Spectrum* f = &spectrum_buffer->spectra[1];
	Spectrum* bdsf_result = &spectrum_buffer->spectra[2];
	Spectrum* reflectance = &spectrum_buffer->spectra[3];
	Spectrum* initial_emission_spd = &spectrum_buffer->spectra[4 + 6];
	NEW_set_spectrum_to_value(eye_ray_radiance, 0.0);
	NEW_set_spectrum_to_value(f, 1.0);

	NEW_Surface_Point prev_surface_point = {};
	NEW_Surface_Point current_surface_point = {};

	//Set previous scene point data
	NEW_Scene_Material prev_scene_point_material_buffer = {};
	NEW_Scene_Point prev_scene_point_buffer = {};
	NEW_Scene_Point* prev_scene_point = &prev_scene_point_buffer;
	prev_scene_point->material = &prev_scene_point_material_buffer;
	prev_scene_point->material->number_of_bdsfs = 1;
	prev_scene_point->material->bdsfs[0].reflectance = NEW_const_1_reflectance;
	prev_scene_point->pdf = 1.0;
	prev_scene_point->out_ray.direction = {1.0, 0.0};
	prev_scene_point->normal = prev_scene_point->out_ray.direction;

	for(int i = 0; i < depth; ++i)
	{
		//Load surface point and emission spectra of point's contributing light sources
		NEW_Scene_Point* current_scene_point = &scene_path[i];
		NEW_load_surface_and_emission_materials(scene, current_scene_point, &current_surface_point, spectrum_buffer);
		double pdf = prev_scene_point->pdf/(double)(prev_scene_point->material->number_of_bdsfs);
		double theta = abs(dot(prev_scene_point->normal, prev_scene_point->out_ray.direction));
		NEW_bdsf(&prev_surface_point, prev_scene_point->material->bdsfs, prev_scene_point->material->number_of_bdsfs, prev_scene_point->in_ray.direction, prev_scene_point->out_ray.direction, reflectance, bdsf_result);
		NEW_spectral_multiply(f, bdsf_result, theta/pdf, f);

		for(int j = 0; j < current_scene_point->number_of_light_contributions; ++j)
		{
			Spectrum* emission_spd = initial_emission_spd + j;
			NEW_Geometry_Sample* light = &current_scene_point->light_contributions[j];
			Vec3 light_direction = normalise(light->position - current_scene_point->position);
			NEW_bdsf(&current_surface_point, current_scene_point->material->bdsfs, current_scene_point->material->number_of_bdsfs, current_scene_point->in_ray.direction, light_direction, reflectance, bdsf_result);
			double light_theta = abs(dot(light_direction, current_scene_point->normal));
			NEW_spectral_multiply(emission_spd, bdsf_result, light_theta/light->pdf, emission_spd);
			NEW_spectral_sum_and_multiply(eye_ray_radiance, emission_spd, f, eye_ray_radiance);
		}

		prev_scene_point = current_scene_point;
		prev_surface_point = current_surface_point;
	}
}


void NEW_render_image(RGB64* render_target, int render_target_width, int render_target_height, int number_of_render_samples, bool invert_final_image)
{
	int number_of_pixels = render_target_width * render_target_height;
	int max_path_depth = 4;

	//Load scene data
	srand(100000);
	NEW_Spectrum_Buffer* spectrum_buffer = (NEW_Spectrum_Buffer*)alloc(sizeof(NEW_Spectrum_Buffer));
	Spectrum* spectral_render_target = (Spectrum*)alloc(number_of_pixels * sizeof(Spectrum));
	NEW_Scene* scene = (NEW_Scene*)alloc(sizeof(NEW_Scene));
	int scene_path_size = (max_path_depth + 1)*sizeof(NEW_Scene_Point);
	NEW_Scene_Point* scene_path = (NEW_Scene_Point*)alloc(scene_path_size);
	NEW_load_scene(scene, spectrum_buffer);

	//Compute camera data
	double aspect_ratio = (double)(render_target_width) / (double)(render_target_height);
	double image_plane_aperture_distance = (scene->camera.focal_length * scene->camera.focal_depth) / (scene->camera.focal_length + scene->camera.focal_depth);
	double image_plane_width = 2.0 * image_plane_aperture_distance * tan_deg(scene->camera_fov/2.0);
	double image_plane_height = image_plane_width / aspect_ratio;
	Vec3 film_top_left = scene->camera_position - 0.5 * image_plane_width * scene->camera_right + 0.5 * image_plane_height * scene->camera_up;
	double pixel_width = image_plane_width/(double)(render_target_width);
	double pixel_height = image_plane_height/(double)(render_target_height);

	scene->camera.pinhole_position = scene->camera_position + image_plane_aperture_distance * scene->camera.film_normal;

	printf("Starting render...\n");
	//Render image
	for(int pass = 0; pass < number_of_render_samples; ++pass)
	{
		printf("Pass %d/%d\n", pass+1, number_of_render_samples);
		for(int pixel = 0; pixel < number_of_pixels; ++pixel)
		{
			int x = pixel % render_target_width;
			int y = pixel / render_target_width;

			//Sample (center) point on the pixel
			Vec3 pixel_top_left = film_top_left + ((double)x)*pixel_width*scene->camera_right - ((double)y)*pixel_height*scene->camera_up;
			Vec3 pixel_sample = Vec3{0.5*pixel_width, 0.5*pixel_height, 0.0};
			Vec3 sampled_pixel_point = pixel_top_left + pixel_sample; 

			//Cast ray to get scene path
			zero_mem(scene_path, scene_path_size);
			int path_depth = NEW_sample_scene_path(scene, sampled_pixel_point, scene_path, max_path_depth);

			//Compute path radiance
			zero_mem(spectrum_buffer, sizeof(NEW_Spectrum_Buffer));
			NEW_compute_path_radiance(scene, spectrum_buffer, scene_path, path_depth);
			/*
			double depth = (double)path_depth/4.0;
			NEW_set_spectrum_to_value(&spectral_render_target[pixel], depth);
			*/
			double pass_d = (double)(pass);
			double pass_inc_d = pass_d + 1.0;
			NEW_spectral_sum_and_multiply(&spectrum_buffer->spectra[0], &spectral_render_target[pixel], pass_d, &spectrum_buffer->spectra[0]);
			NEW_spectral_multiply(&spectrum_buffer->spectra[0], 1.0/pass_inc_d, &spectral_render_target[pixel]);
		}
	}

	//Convert spectral image to RGB64
	NEW_set_spectrum_to_value(&spectrum_buffer->spectra[0], 1.0);
	NEW_load_spd("spectra/cmf_x.csv", &spectrum_buffer->spectra[1]);
	NEW_load_spd("spectra/cmf_y.csv", &spectrum_buffer->spectra[2]);
	NEW_load_spd("spectra/cmf_z.csv", &spectrum_buffer->spectra[3]);
	for(int i = 0; i < number_of_pixels; ++i)
	{
		render_target[(invert_final_image) ? (number_of_pixels - 1 - i) : i] = NEW_spectrum_to_RGB64(&spectral_render_target[i], &spectrum_buffer->spectra[0], &spectrum_buffer->spectra[1], &spectrum_buffer->spectra[4]);
	}
}

/*
void NEW_bdsf(NEW_Surface_Point* p, NEW_BDSF* bdsfs, int number_of_bdsfs, Vec3 incoming, Vec3 outgoing, Spectrum* reflectance, Spectrum* bdsf_result)
{
	//assume reflectance is 0
	
	for(int i = 0; i < number_of_bdsfs; ++i)
	{
		//so long as bdsf functions do not have any initial assumptions about contents of the result spectrum
		//	then bdsf result spectrum does not need to be set to zero each iteration
		bdsfs[i].reflectance(p, incoming, outgoing, bdsf_result);
		spectral_sum(*reflectance, bdsf_result, *reflectance);
	}
}

void NEW_load_surface_and_emission_materials(NEW_Scene* scene, NEW_Scene_Point* p, NEW_Surface_Point* surface_point, NEW_Spectrum_Buffer* spectrum_buffer)
{
	//NOTE: When texture support is readded, this will be one of the places where it will need to happen
	//Load into spectrum buffer (could involve texture sampling):
	//	- Scene point's surface material spds
	//		- Incident + transmit extinct + refract indices
	//		- Glossy and diffuse spds
	//	- For each light contribution
	//		- Light point's emission spd
	
	copy_spectrum(&p->incident_material->refract_index_spd, &spectrum_buffer->spectra[4]);
	copy_spectrum(&p->incident_material->extinct_index_spd, &spectrum_buffer->spectra[5]);
	copy_spectrum(&p->transmit_material->refract_index_spd, &spectrum_buffer->spectra[6]);
	copy_spectrum(&p->transmit_material->extinct_index_spd, &spectrum_buffer->spectra[7]);
	copy_spectrum(&p->material->glossy_spd, &spectrum_buffer->spectra[8]);
	copy_spectrum(&p->material->diffuse_spd, &spectrum_buffer->spectra[9]);

	for(int i = 0; i < p->number_of_light_contributions; ++i)
	{
		copy_spectrum(&scene->light_source_materials[p->light_contributions[i].object]->emission_spd, spectrum_buffer->spectra[10+i]);
	}
	
	
	//Set surface point
	surface_point->normal = scene_point->normal;
	surface_point->incident_refract_index_spd = &spectrum_buffer->spectra[4];
	surface_point->incident_extinct_index_spd = &spectrum_buffer->spectra[5];
	surface_point->transmit_refract_index_spd = &spectrum_buffer->spectra[6];
	surface_point->transmit_extinct_index_spd = &spectrum_buffer->spectra[7];
	surface_point->glossy_spd = &spectrum_buffer->spectra[8];
	surface_point->diffuse_spd = &spectrum_buffer->spectra[9];
	surface_point->roughness = p->material->roughness;
	surface_point->shininess = p->material->shininess;
}

void NEW_compute_path_radiance(NEW_Scene* scene, NEW_Spectrum_Buffer* spectrum_buffer, NEW_Scene_Point* scene_path, int depth)
{
	Spectrum* eye_ray_radiance = spectrum_buffer->spectra;
	Spectrum* f = &spectrum_buffer->spectra[1];
	Spectrum* bdsf_result = &spectrum_buffer->spectra[2];
	Spectrum* reflectance = &spectrum_buffer->spectra[3];
	Spectrum* initial_emission_spd = &spectrum_buffer[4 + 6]; //+6 for surface point's spectra
	set_spectrum_to_value(*f, 1.0);

	NEW_Surface_Point prev_surface_point = {};
	NEW_Scene_Point prev_scene_point_buffer = {};
	NEW_Surface_Point surface_point = {};

	New_Scene_Point* prev_scene_point = &surface_point;
	prev_surface_point->normal = {1.0, 1.0};
	prev_scene_point->number_of_bdsfs = 1;
	prev_scene_point->bdsfs[0].reflectance = &NEW_const_1_bdsf;
	prev_scene_point->number_of_bdsfs = 1;
	prev_scene_point->pdf = 1.0;

	for(int i = 0; i < depth; ++i)
	{//For each point in the scene path
		//Load surface point and emission spectra of point's contributing light sources
		NEW_Scene_Point* scene_point = &scene_path[i];
		NEW_load_surface_and_emission_materials(scene, scene_point, surface_point, spectrum_buffer);

		//Find coefficient for lambert's law and pdf
		double pdf = prev_scene_point->pdf/(double)(prev_scene_point->number_of_bdsfs);
		double theta = abs(dot(prev_scene_point->normal, prev_scene_point->out_ray.direction));
		//Compute previous point's bdsf and compute new f value
		NEW_bdsf(prev_surface_point, prev_scene_point->bdsfs, prev_scene_point->number_of_bdsfs, prev_scene_point->in_direction, prev_scene_point->out_direction, reflectance, bdsf_result);
		spectral_multiply(f, bdsf_result, theta/pdf, f);

		for(int j = 0; j < scene_point.number_of_light_contributions; ++j)
		{
			Spectrum* emission_spd = initial_emission_spd + j;
			NEW_Geometry_Sample* light = &scene_point.light_contributions[j];
			Vec3 light_direction = normalise(light->position - scene_point->position);
			//do bdsf of current point, in direction and light direction
			NEW_bdsf(surface_point, scene_point->bdsfs, scene_point->number_of_bdsfs, scene_point->in_ray.direction, light_direction, reflectance, bdsf_result);
			//Find lambert cosine of light direction to scene point surface
			double light_theta = abs(dot(light_direction, scene_point->normal));
			//Multiply light point's emission spd with bdsf result and cosine
			spectral_multiply(emission_spd, bdsf_result, light_theta/light->pdf, emission_spd);
			//Eye ray radiance = eye_ray_radiance + light_theta * previous spd result
			spectral_sum_and_multiply(eye_ray_radiance, light_theta, emission_spd, eye_ray_radiance);
		}
	}
}
*/
