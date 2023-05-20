
void RGB64_to_spectrum(RGB64 c, Spectrum* dst, Spectrum* spectra)
{
	RGB64_to_spectrum(c, dst, &spectra[0], 
			&spectra[1], &spectra[2], &spectra[3], 
			&spectra[4], &spectra[5], &spectra[6]);
}

void load_scene(Scene* scene, Spectrum* spectra)
{
	load_spd("spectra/white_rgb_to_spd.csv", &spectra[0]);
	load_spd("spectra/red_rgb_to_spd.csv", &spectra[1]);
	load_spd("spectra/blue_rgb_to_spd.csv", &spectra[2]);
	load_spd("spectra/green_rgb_to_spd.csv", &spectra[3]);
	load_spd("spectra/cyan_rgb_to_spd.csv", &spectra[4]);
	load_spd("spectra/magenta_rgb_to_spd.csv", &spectra[5]);
	load_spd("spectra/yellow_rgb_to_spd.csv", &spectra[6]);

	//Geometries
	double h = 3.0;
	scene->planes[0] = create_plane_from_points(Vec3{-h, h, -h}, Vec3{h, h, -h}, Vec3{-h, -h, -h}); //back wall
	scene->planes[1] = create_plane_from_points(Vec3{-h, h, h}, Vec3{-h, h, -h}, Vec3{-h, -h, h}); //left wall
	scene->planes[2] = create_plane_from_points(Vec3{h, h, -h}, Vec3{h, h, h}, Vec3{h, -h, -h}); //right wall
	scene->planes[3] = create_plane_from_points(Vec3{-h, -h, -h}, Vec3{h, -h, -h}, Vec3{-h, -h, h}); //floor
	scene->planes[4] = create_plane_from_points(Vec3{-h, h, h}, Vec3{h, h, h}, Vec3{-h, h, -h});  //ceiling
	scene->planes[5] = create_plane_from_points(Vec3{-0.5, 2.9, 0.5}, Vec3{0.5, 2.9, 0.5}, Vec3{-0.5, 2.9, -0.5}); //Light
	scene->planes[6] = create_plane_from_points(Vec3{-1.0, 1.0, -2.4}, Vec3{1.0, 1.0, -2.4}, Vec3{-1.0, -1.0, -2.9});//Mirror
	scene->spheres[0] = {{-2.0, -2.2, -0.5}, 0.75};

	scene->number_of_geometries = 7;
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
	scene->geometries[6].type = GEO_TYPE_PLANE;
	scene->geometries[6].plane = &scene->planes[6];
	scene->geometries[7].type = GEO_TYPE_SPHERE;
	scene->geometries[7].sphere = &scene->spheres[0];

	//Materials
	load_spd("spectra/air_spec.csv", &scene->air_material.refract_index_spd);
	scene->number_of_materials = 7;
	//back wall (blue)
	scene->materials[0] = {};
	scene->materials[0].shininess = 32.0;
	RGB64_to_spectrum(RGB64{0.2, 0.2, 0.8}, &scene->materials[0].diffuse_spd, spectra);
	RGB64_to_spectrum(RGB64{0.8, 0.8, 0.9}, &scene->materials[0].glossy_spd, spectra);
	scene->materials[0].number_of_bdsfs = 2;
	scene->materials[0].bdsfs[0].reflectance = diffuse_phong_reflectance;
	scene->materials[0].bdsfs[0].sample_direction = sample_diffuse_direction;
	scene->materials[0].bdsfs[0].pdf = cos_weighted_hemisphere_pdf;
	scene->materials[0].bdsfs[1].reflectance = glossy_phong_reflectance;
	scene->materials[0].bdsfs[1].sample_direction = sample_glossy_direction;
	scene->materials[0].bdsfs[1].pdf = cos_weighted_hemisphere_pdf;

	//left wall (green)
	scene->materials[1] = {};
	scene->materials[1].shininess = 32.0;
	RGB64_to_spectrum(RGB64{0.1, 0.35, 0.1}, &scene->materials[1].diffuse_spd, spectra);
	RGB64_to_spectrum(RGB64{0.45, 0.55, 0.45}, &scene->materials[1].glossy_spd, spectra);
	scene->materials[1].number_of_bdsfs = 2;
	scene->materials[1].bdsfs[0].reflectance = diffuse_phong_reflectance;
	scene->materials[1].bdsfs[0].sample_direction = sample_diffuse_direction;
	scene->materials[1].bdsfs[0].pdf = cos_weighted_hemisphere_pdf;
	scene->materials[1].bdsfs[1].reflectance = glossy_phong_reflectance;
	scene->materials[1].bdsfs[1].sample_direction = sample_glossy_direction;
	scene->materials[1].bdsfs[1].pdf = cos_weighted_hemisphere_pdf;
	
	//right wall (red)
	scene->materials[2] = {};
	scene->materials[2].shininess = 32.0;
	RGB64_to_spectrum(RGB64{0.5, 0.0, 0.0}, &scene->materials[2].diffuse_spd, spectra);
	RGB64_to_spectrum(RGB64{0.7, 0.6, 0.6}, &scene->materials[2].glossy_spd, spectra);
	scene->materials[2].number_of_bdsfs = 2;
	scene->materials[2].bdsfs[0].reflectance = diffuse_phong_reflectance;
	scene->materials[2].bdsfs[0].sample_direction = sample_diffuse_direction;
	scene->materials[2].bdsfs[0].pdf = cos_weighted_hemisphere_pdf;
	scene->materials[2].bdsfs[1].reflectance = glossy_phong_reflectance;
	scene->materials[2].bdsfs[1].sample_direction = sample_glossy_direction;
	scene->materials[2].bdsfs[1].pdf = cos_weighted_hemisphere_pdf;

	//floor + ceiling (white)
	scene->materials[3] = {};
	scene->materials[3].shininess = 32.0;
	RGB64_to_spectrum(RGB64{0.55, 0.55, 0.55}, &scene->materials[3].diffuse_spd, spectra);
	RGB64_to_spectrum(RGB64{0.7, 0.7, 0.7}, &scene->materials[3].glossy_spd, spectra);
	scene->materials[3].number_of_bdsfs = 2;
	scene->materials[3].bdsfs[0].reflectance = diffuse_phong_reflectance;
	scene->materials[3].bdsfs[0].sample_direction = sample_diffuse_direction;
	scene->materials[3].bdsfs[0].pdf = cos_weighted_hemisphere_pdf;
	scene->materials[3].bdsfs[1].reflectance = glossy_phong_reflectance;
	scene->materials[3].bdsfs[1].sample_direction = sample_glossy_direction;
	scene->materials[3].bdsfs[1].pdf = cos_weighted_hemisphere_pdf;

	//light
	scene->materials[4] = {};
	scene->materials[4].is_emissive = true;
	set_spectrum_to_value(&scene->materials[4].emission_spd, 1.0);

	//Mirror
	scene->materials[5] = {};
	scene->materials[5].number_of_bdsfs = 1;
	scene->materials[5].bdsfs[0].reflectance = mirror_reflectance;
	scene->materials[5].bdsfs[0].sample_direction = sample_specular_reflection_direction;
	scene->materials[5].bdsfs[0].pdf = const_1_pdf;

	//Gold
	/*
	double roughness = 0.344;
	scene->materials[6] = {};
	scene->materials[6].number_of_bdsfs = 1;
	scene->materials[6].bdsfs[0].reflectance = (roughness >= 0.001) ? cook_torrance_conductor_reflectance : fresnel_conductor_specular_reflectance;
	scene->materials[6].bdsfs[0].sample_direction = (roughness >= 0.001) ? sample_cook_torrance_reflection_direction : sample_specular_reflection_direction;
	scene->materials[6].bdsfs[0].pdf = (roughness >= 0.001) ? cook_torrance_reflection_pdf : const_1_pdf;
	scene->materials[6].roughness = roughness;
	load_spd("spectra/au_spec_n.csv", &scene->materials[6].refract_index_spd);
	load_spd("spectra/au_spec_k.csv", &scene->materials[6].extinct_index_spd);
	*/

	//Set scene objects and light sources
	scene->number_of_objects = 7;

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

	//Mirror
	scene->object_geometries[6] = &scene->geometries[6];
	scene->object_materials[6] = &scene->materials[5];

	//Gold
	scene->object_geometries[7] = &scene->geometries[7];
	scene->object_materials[7] = &scene->materials[6];

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
	scene->film_material.bdsfs[0].sample_direction = (scene->camera.aperture_radius > 0.0) ? sample_camera_lens_direction : sample_camera_pinhole_direction;
	scene->film_material.bdsfs[0].reflectance = const_1_reflectance;
	scene->film_material.bdsfs[0].pdf = const_1_pdf; //TODO: Change this for lens (probably shouldn't be 1)
}


Vec3 find_ray_scene_intersection_point(Scene* scene, Vec3 ray_origin, Vec3 ray_direction, int* ret_intersecting_object = NULL)
{
	double min_length_along_ray = DBL_MAX;
	double length_along_ray = DBL_MAX;
	Vec3 min_intersection_along_ray = {DBL_MAX, DBL_MAX, DBL_MAX};
	int intersecting_object = -1;
	
	for(int i = 0; i < scene->number_of_geometries; ++i)
	{
		Scene_Geometry* geometry = scene->object_geometries[i];
		switch(geometry->type)
		{
			case GEO_TYPE_SPHERE: length_along_ray = ray_intersects_sphere(ray_origin, ray_direction, *geometry->sphere); break;
			case GEO_TYPE_PLANE: length_along_ray = ray_intersects_plane(ray_origin, ray_direction, *geometry->plane); break;
			//TODO: case GEO_TYPE_MODEL
		}
		if(length_along_ray > -0.0009765625 && length_along_ray < min_length_along_ray)
		{
			min_length_along_ray = length_along_ray;
			min_intersection_along_ray = ray_origin + min_length_along_ray * ray_direction;
			intersecting_object = i;
		}
	}
	if(ret_intersecting_object) *ret_intersecting_object = intersecting_object;

	return min_intersection_along_ray;
}

int sample_scene_path(Scene* scene, Vec3 film_point, Scene_Path_Point* scene_path, int max_depth)
{
	bool is_there_next_point = true;
	Scene_Path_Point* current_point = &scene_path[0];
	Scene_Path_Point* next_point = &scene_path[1];
	current_point->properties.position = film_point;
	current_point->properties.normal = scene->camera.film_normal;
	current_point->properties.camera = &scene->camera;
	current_point->in_direction = {1.0, 1.0, 0.0};
	current_point->number_of_bdsfs = scene->film_material.number_of_bdsfs;
	current_point->bdsfs = scene->film_material.bdsfs;

	Scene_Material* current_point_material = &scene->film_material;
	Scene_Material* next_point_material = NULL;

	int depth = 0;
	for(; depth < max_depth && is_there_next_point; ++depth)
	{
		//Sample direct lighting contributions for the current point
		for(int i = 0; depth != 0 & i < scene->number_of_light_sources; ++i)
		{ // If the current point isn't on the camera's film surface and for each light source in the scene
			Light_Contribution* potential_light_contribution = &current_point->light_contributions[current_point->number_of_light_contributions];
			Scene_Geometry* light_geometry = scene->light_source_geometries[i];
			//Sample the current light's geometry for a point on its surface
			switch(geometry->type)
			{
				case GEO_TYPE_SPHERE: light_sample->position = uniform_sample_sphere_subtended(*geometry->sphere, current_point->position, &light_sample->pdf); break;
				case GEO_TYPE_PLANE: light_sample->position = uniform_sample_plane(*geometry->plane, &light_sample->pdf); break;
				//TODO: case GEO_TYPE_MODEL
			}
			
			Ray shadow_ray = {};
			shadow_ray.direction = normalise(light_sample->position - current_point->position);
			shadow_ray.origin = current_point->position + 0.03125*shadow_ray.direction;
			
			//Test if an object in the scene obscures the light position from the current point's position
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

		//Sample next point in scene path
		double r = uniform_sample();
		double n = (double)current_path_point->number_of_bdsfs;
		int sampled_bdsf = (int)floor(r * n);
		current_path_point->out_direction = current_path_point->bdsfs[sampled_bdsf].sample_direction(&current_path_point->properties, current_path_point->in_direction);
		Vec3 out_ray_origin = current_path_point->properties.position + 0.03125*current_path_point->out_direction;
		
		//Compute next point in path
		int next_object = -1;
		next_path_point->properties.position = find_ray_scene_intersection_point(scene, out_ray_origin, current_path_point->out_direction, &next_object);
		bool is_next_point_emissive = false;
		//Set next point properties
		if(next_object >= 0)
		{
			next_path_point->in_direction = -current_path_point->out_direction;
			switch(scene->object_geometries[next_object]->type)
			{
				case GEO_TYPE_SPHERE:
				{
					next_path_point->properties.normal = normalise(next_path_point->properties.position - scene->object_geometries[next_object]->sphere->center);
					break;
				}
				case GEO_TYPE_PLANE:
				{
					next_path_point->properties.normal = scene->object_geometries[next_object]->plane->n;
					break;
				}
			}
			next_path_point->properties.roughness = scene->object_materials[next_object]->roughness;
			next_path_point->properties.shininess = scene->object_materials[next_object]->shininess;
			next_path_point->properties.glossy_spd = &scene->object_materials[next_object]->glossy_spd;
			next_path_point->properties.diffuse_spd = &scene->object_materials[next_object]->diffuse_spd;
			next_path_point->number_of_bdsfs = scene->object_materials[next_object]->number_of_bdsfs;
			next_path_point->bdsfs = scene->object_materials[next_object]->bdsfs;
			next_point_is_emissive = scene->object_materials[next_object]->is_emissive;
			//Set incident and transmission materials
			//For now, if entering an object air is the incident material, otherwise the object's
			double in_dot = dot(current_path_point->out_direction, next_path_point->properties.normal);
			if(in_dot < 0.0)
			{
				next_path_point->properties.incident_refract_index_spd = &scene->air_material.refract_index_spd;
				next_path_point->properties.incident_extinct_index_spd = &scene->air_material.extinct_index_spd;
				next_path_point->properties.transmit_refract_index_spd = &scene->object_materials[next_object]->refract_index_spd;
				next_path_point->properties.transmit_extinct_index_spd = &scene->object_materials[next_object]->extinct_index_spd;
			}
			else
			{
				next_path_point->properties.incident_refract_index_spd = &scene->object_materials[next_object].refract_index_spd;
				next_path_point->properties.incident_extinct_index_spd = &scene->object_materials[next_object].extinct_index_spd;
				next_path_point->properties.transmit_refract_index_spd = &scene->air_material.refract_index_spd;
				next_path_point->properties.transmit_extinct_index_spd = &scene->air_material.extinct_index_spd;
			}
		}
		else
		{
			next_point_is_emissive = false;
		}

		//Iterate points
		current_point = next_point;
		is_there_next_point = !next_point_is_emissive && next_object != -1;
		++next_point;
	}
	if(!is_there_next_point && depth <= max_depth) 
	{
		scene_path[depth-1].light_contributions[scene_path[depth-1].number_of_light_contributions++] = *((NEW_Geometry_Sample*)(&current_point->position));
	}
	
	return depth;
}

void bdsf(Scene_Path_Point* path_point, Vec3 out_direction, Spectrum* reflectance, Spectrum* bdsf_result)
{
	set_spectrum_to_value(reflectance, 0.0);
	set_spectrum_to_value(bdsf_result, 0.0);
	for(int i = 0; i < path_point->number_of_bdsfs; ++i)
	{
		path_point->bdsfs[i].reflectance(&path_point->properties,  path_point->in_direction, out_direction, reflectance);
		spectral_sum(reflectance, bdsf_result, bdsf_result);
	}
}

void compute_path_radiance(Scene_Path_Point* scene_path, int depth, Spectrum* eye_ray_radiance, Spectrum* f, Spectrum* bdsf_result, Spectrum* reflectance)
{
	set_spectrum_to_value(eye_ray_radiance, 0.0);
	set_spectrum_to_value(f, 1.0);

	//Set previous scene point data
	Scene_Path_Point prev_point_buffer = {};
	prev_point_buffer->number_of_bdsfs = 1;
	prev_point_buffer->bdsfs[0].reflectance = const_1_reflectance;
	prev_point_buffer->bdsfs[0].pdf = const_1_pdf;
	prev_point_buffer->out_direction = {1.0, 0.0, 0.0};
	prev_point_buffer->normal = prev_point_buffer->out_direction;

	Scene_Path_Point* prev_path_point = &prev_point_buffer;

	for(int i = 0; i < depth; ++i)
	{
		//Load surface point and emission spectra of point's contributing light sources
		Scene_Path_Point* current_path_point = &scene_path[i];

		bdsf(prev_path_point, prev_path_point->out_direction, reflectance, bdsf_result);
		double pdf = 0.0;
		for(int j = 0; j < prev_path_point->number_of_bdsfs; ++j)
		{
			pdf += current_point->material->bdsfs[j].pdf(&prev_path_point.properties, prev_path_point->in_direction, current_point->out_direction);
		}
		pdf /= (double)(prev_path_point->number_of_bdsfs);
		double f_coefficient = abs(dot(prev_path_point->properties.normal, prev_path_point->out_direction)) / pdf;
		spectral_multiply(f, bdsf_result, f_coefficient, f);

		for(int j = 0; j < current_path_point->number_of_light_contributions; ++j)
		{
			Light_Contribution* light = &current_scene_point->light_contributions[j];
			bdsf(current_path_point, light->direction, reflectance, bdsf_result);
			double light_pdf = 0.0;
			switch(light->geometry.type)
			{
				case GEO_TYPE_SPHERE:
				{
					light_pdf = uniform_sample_sphere_subtended_pdf(light->geometry->sphere, current_path_point->properties.position);
					break;
				}
				case GEO_TYPE_PLANE:
				{
					light_pdf = uniform_sample_plane_pdf(light->geometry->plane);
					break;
				}
			}
			double light_coefficient = abs(dot(light->direction, current_path_point->properties.normal)) / light_pdf;
			spectral_multiply(emission_spd, bdsf_result, light_coefficient, emission_spd);
			spectral_sum_and_multiply(eye_ray_radiance, light->emission_spd, f, eye_ray_radiance);
		}

		prev_path_point = current_path_point;
	}
}


void render_image(RGB64* render_target, int render_target_width, int render_target_height, int number_of_render_samples, bool invert_final_image)
{
	int number_of_pixels = render_target_width * render_target_height;
	int max_path_depth = 4;

	//Load scene data
	srand(100000);
	int scene_path_size = (max_path_depth + 1)*sizeof(Scene_Point);
	Spectrum* spectral_render_target = (Spectrum*)alloc(number_of_pixels * sizeof(Spectrum));
	Scene* scene = (Scene*)alloc(sizeof(Scene));
	Spectrum* spectrum_buffer = (Spectrum_Buffer*)alloc(16*sizeof(Spectrum));
	Scene_Path_Point* scene_path = (Scene_Path_Point*)alloc(scene_path_size);
	load_scene(scene, spectrum_buffer);

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
	double max_render_time = 0.0;
	double min_render_time = DBL_MAX;
	double max_path_time = 0.0;
	double min_path_time = DBL_MAX;
	double max_radiance_time = 0.0;
	double min_radiance_time = DBL_MAX;
	double total_render_time = 0.0;
	double total_path_time = 0.0;
	double total_radiance_time = 0.0;
	Timer path_t = {};
	Timer radiance_t = {};
	Timer t = {};

	DEBUG
	(
		debug_info.scene = scene;
		debug_info.scene_path = scene_path;
		debug_info.spectrum_buffer = spectrum_buffer;
		debug_info.render_target = spectral_render_target;
	)

	//Render image
	for(int pass = 0; pass < number_of_render_samples; ++pass)
	{
		printf("Pass %d/%d\n", pass+1, number_of_render_samples);
		start_timer(&t);
		for(int pixel = 0; pixel < number_of_pixels; ++pixel)
		{
			int x = pixel % render_target_width;
			int y = pixel / render_target_width;

			DEBUG
			(
				debug_info.pixel_x = x;
				debug_info.pixel_y = y;
				debug_info.pass = pass;
			);
			zero_mem(scene_path, scene_path_size);


			//Sample (center) point on the pixel
			Vec3 pixel_top_left = film_top_left + ((double)x)*pixel_width*scene->camera_right - ((double)y)*pixel_height*scene->camera_up;
			Vec3 pixel_sample = Vec3{0.5*pixel_width, 0.5*pixel_height, 0.0};
			Vec3 sampled_pixel_point = pixel_top_left + pixel_sample; 

			//Cast ray to get scene path
			start_timer(&path_t);
			int path_depth = sample_scene_path(scene, sampled_pixel_point, scene_path, max_path_depth);
			stop_timer(&path_t);

			//Compute path radiance
			start_timer(&radiance_t);
			compute_path_radiance(scene_path, path_depth, &spectrum_buffer[0], &spectrum_buffer[1], &spectrum_buffer[2], &spectrum_buffer[3]);
			stop_timer(&radiance_t);
			/*
			double depth = (double)path_depth/4.0;
			NEW_set_spectrum_to_value(&spectral_render_target[pixel], depth);
			*/
			double pass_d = (double)(pass);
			double pass_inc_d = pass_d + 1.0;
			spectral_sum_and_multiply(&spectrum_buffer[0], &spectral_render_target[pixel], pass_d, &spectrum_buffer[0]);
			spectral_multiply(&spectrum_buffer[0], 1.0/pass_inc_d, &spectral_render_target[pixel]);
			double path_time = elapsed_time_in_ms(&path_t);
			double radiance_time = elapsed_time_in_ms(&radiance_t);
			if(path_time > max_path_time) max_path_time = path_time;
			if(path_time < min_path_time) min_path_time = path_time;
			if(radiance_time > max_radiance_time) max_radiance_time = radiance_time;
			if(radiance_time < min_radiance_time) min_radiance_time = radiance_time;
			total_path_time += path_time;
			total_radiance_time += radiance_time;
		}
		stop_timer(&t);
		double render_time = elapsed_time_in_ms(&t);
		if(render_time > max_render_time) max_render_time = render_time;
		if(render_time < min_render_time) min_render_time = render_time;
		total_render_time += render_time;
	}
	printf("Max render time = %fms\n", max_render_time);
	printf("Min render time = %fms\n", min_render_time);
	printf("Total render time = %fms\n", total_render_time);
	printf("Max path time = %fms\n", max_path_time);
	printf("Min path time = %fms\n", min_path_time);
	printf("Total path time = %fms\n", total_path_time);
	printf("Max radiance time = %fms\n", max_radiance_time);
	printf("Min radiance time = %fms\n", min_radiance_time);
	printf("Total radiance time = %fms\n", total_radiance_time);

	//Convert spectral image to RGB64
	set_spectrum_to_value(&spectrum_buffer[0], 1.0);
	load_spd("spectra/cmf_x.csv", &spectrum_buffer[1]);
	load_spd("spectra/cmf_y.csv", &spectrum_buffer[2]);
	load_spd("spectra/cmf_z.csv", &spectrum_buffer[3]);
	for(int i = 0; i < number_of_pixels; ++i)
	{
		render_target[(invert_final_image) ? (number_of_pixels - 1 - i) : i] = spectrum_to_RGB64(&spectral_render_target[i], &spectrum_buffer[0], &spectrum_buffer[1], &spectrum_buffer[4]);
	}
}

