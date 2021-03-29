
void add_object_to_scene(Scene* scene, Scene_Object object)
{
	scene->objects[scene->number_of_objects++] = object;
}

void add_point_light_to_scene(Scene* scene, Point p, Spectrum emission_spd)
{
	Scene_Object point_light = {};
	Scene_Geometry point_geometry = {};
	point_geometry.type = GEO_TYPE_POINT;
	point_geometry.point = p;
	
	point_light.geometry = point_geometry;
	point_light.material.emission_spd = emission_spd;
	point_light.light_type = LIGHT_TYPE_POINT;
	point_light.is_emissive = true;

	add_object_to_scene(scene, point_light);
}

void add_sphere_light_to_scene(Scene* scene, Sphere s, Spectrum emission_spd)
{
	Scene_Object sphere_light = {};
	Scene_Geometry sphere_geometry = {};
	sphere_geometry.type = GEO_TYPE_SPHERE;
	sphere_geometry.sphere = s;

	sphere_light.geometry = sphere_geometry;
	sphere_light.material.emission_spd = emission_spd;
	sphere_light.light_type = LIGHT_TYPE_AREA;
	sphere_light.is_emissive = true;

	add_object_to_scene(scene, sphere_light);
}

void add_plane_light_to_scene(Scene* scene, Plane p, Spectrum emission_spd, char* name)
{
	Scene_Object plane_light = {};
	Scene_Geometry plane_geometry = {};
	plane_geometry.type = GEO_TYPE_PLANE;
	plane_geometry.plane = p;

	plane_light.name = name;
	plane_light.geometry = plane_geometry;
	plane_light.material.emission_spd = emission_spd;
	plane_light.light_type = LIGHT_TYPE_AREA;
	plane_light.is_emissive = true;

	add_object_to_scene(scene, plane_light);
}

void add_sphere_to_scene(Scene* scene, Sphere s, Material material, char* name)
{
	Scene_Object sphere = {};
	Scene_Geometry sphere_geometry = {};
	sphere_geometry.type = GEO_TYPE_SPHERE;
	sphere_geometry.sphere = s;

	sphere.name = name;
	sphere.geometry = sphere_geometry;
	sphere.material = material;

	add_object_to_scene(scene, sphere);
}

void add_plane_to_scene(Scene* scene, Plane p, Material material, char* name)
{
	Scene_Object plane = {};
	Scene_Geometry plane_geometry = {};
	plane_geometry.type = GEO_TYPE_PLANE;
	plane_geometry.plane = p;

	plane.name = name;
	plane.geometry = plane_geometry;
	plane.material = material;

	add_object_to_scene(scene, plane);
}

void add_model_to_scene(Scene* scene, Model m, Material material, char* name)
{
	Scene_Object model = {};
	Scene_Geometry model_geometry = {};
	model_geometry.type = GEO_TYPE_MODEL;
	model_geometry.model = m;

	model.name = name;
	model.geometry = model_geometry;
	model.material = material;

	add_object_to_scene(scene, model);
}

bool ray_intersects_model(Ray ray, Model m, double* t, int* ret_intersecting_triangle_index, double* ret_a, double* ret_b, double* ret_c)
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
			Vec3 p = ray.origin + distance_to_triangle * ray.direction;

			double temp_a = sqrt(point_to_line_distance_sq(p, m.vertices[i+1].position, m.vertices[i+2].position)/point_to_line_distance_sq(m.vertices[i].position, m.vertices[i+1].position,m.vertices[i+2].position));
			double temp_b = sqrt(point_to_line_distance_sq(p, m.vertices[i+2].position, m.vertices[i].position)/point_to_line_distance_sq(m.vertices[i+1].position, m.vertices[i+2].position, m.vertices[i].position));
			double temp_c = sqrt(point_to_line_distance_sq(p, m.vertices[i].position, m.vertices[i+1].position)/point_to_line_distance_sq(m.vertices[i+2].position, m.vertices[i].position, m.vertices[i+1].position));

			double d = temp_a + temp_b + temp_c;
			if(d <= 1.00001 && distance_to_triangle < min_distance_to_triangle && distance_to_triangle > 0.0)
			{
				intersecting_triangle_index = i;
				min_distance_to_triangle = distance_to_triangle;
				a = temp_a;
				b = temp_b;
				c = temp_c;
			}
		}
	}
	
	*ret_a = a;
	*ret_b = b;
	*ret_c = c;

	*t = min_distance_to_triangle;
	*ret_intersecting_triangle_index = intersecting_triangle_index;
	return intersecting_triangle_index >= 0;
}

Geometry_Intersection_Point find_ray_scene_intersection(Scene* scene, Ray ray)
{
	Geometry_Intersection_Point intersection_point;
	intersection_point.scene_object = -1;

	double length_along_ray = DBL_MAX;
	double ith_length_along_ray = DBL_MAX;
	for(int i = 0; i < scene->number_of_objects; ++i)
	{
		Scene_Geometry* ith_object_geometry = &(scene->objects[i].geometry);
		if(ith_object_geometry->type)
		{
			bool ray_intersects_ith_object = false;
			double ith_length_along_ray = DBL_MAX;
			switch(ith_object_geometry->type)
			{
				case GEO_TYPE_SPHERE:
				{
					ray_intersects_ith_object = ray_intersects_sphere(ray, ith_object_geometry->sphere, &ith_length_along_ray);
					break;
				}
				case GEO_TYPE_PLANE:
				{
					ray_intersects_ith_object = ray_intersects_plane(ray, ith_object_geometry->plane, &ith_length_along_ray);
					break;
				}
				case GEO_TYPE_MODEL:
				{
					ray_intersects_ith_object = ray_intersects_model(ray, ith_object_geometry->model, &ith_length_along_ray, &intersection_point.model_triangle_index, &intersection_point.bc_a, &intersection_point.bc_b, &intersection_point.bc_c);
					break;
				}
			}

			if(ray_intersects_ith_object && ith_length_along_ray < length_along_ray)
			{
				intersection_point.scene_object = i;
				length_along_ray = ith_length_along_ray;
			}
		}

	}
	intersection_point.position = ray.origin + length_along_ray * ray.direction;

	return intersection_point;
}

Radiance direct_light_contribution(Scene* scene, Surface_Point p, Ray outgoing)
{
	Radiance contribution = {};
	if(!p.is_emissive)
	{
		for(int i = 0; i < scene->number_of_objects; ++i)
		{
			if(scene->objects[i].is_emissive)
			{
				double light_pdf = 1.0;
				Scene_Object light = scene->objects[i];
				Vec3 light_point = {};
				switch(light.geometry.type)
				{
					case GEO_TYPE_POINT: 
					{
						light_pdf = 1.0;
						light_point = light.geometry.point.position;
					}
					case GEO_TYPE_SPHERE: 
					{
						light_point = uniform_sample_sphere_subtended(light.geometry.sphere, p.position, &light_pdf);
					}
					case GEO_TYPE_PLANE: 
					{
						light_point = uniform_sample_plane(light.geometry.plane, &light_pdf);
					}
				}
				
				Ray shadow_ray = {};
				shadow_ray.direction = normalise(light_point - p.position);
				shadow_ray.origin = p.position + 0.001*shadow_ray.direction;

				Geometry_Intersection_Point shadow_test_point = find_ray_scene_intersection(scene, shadow_ray);

				double shadow_ray_length = length(shadow_test_point.position - p.position);
				double ray_length = length(light_point - p.position);
				double difference_between_ray_lengths = shadow_ray_length - ray_length;
				bool light_point_visible = !(difference_between_ray_lengths < -5e-12);

				if(light_point_visible && light_pdf > 0.0)
				{
					Vec3 incoming = light_point - p.position;
					double dist = length(incoming);
					incoming = normalise(incoming);

					double attenuation_factor = (light.light_type == LIGHT_TYPE_POINT) ? 1.0/(4.0 * PI * dist * dist) : 1.0;
					double d = abs(dot(incoming, p.normal));
					double f = attenuation_factor * d / light_pdf;

					contribution += f * bsdf(p, incoming, outgoing.direction) * light.material.emission_spd;
				}
			}
		}
	}

	return contribution;
}

//MIRROR CHANGE: Probability of sampling specular direction depends on material
Vec3 choose_incoming_direction(Surface_Point p, Vec3 outgoing, double* pdf_value, bool* consider_emissive)
{
	double r = uniform_sample();
	double n = (double)p.material.number_of_bsdfs;
	int chosen_bsdf = (int)floor(r * n);
	BSDF_TYPE chosen_bsdf_type = p.material.bsdfs[chosen_bsdf].type;
	*consider_emissive = chosen_bsdf_type & BSDF_TYPE_SPECULAR;
	Vec3 direction = p.material.bsdfs[chosen_bsdf].sample_direction(p, outgoing, pdf_value);

	//Sum probabilities for BSDFs with same distributions (that aren't specular)
	for(int i = 0; i < p.material.number_of_bsdfs; ++i)
	{
		if(i != chosen_bsdf && chosen_bsdf_type != BSDF_TYPE_SPECULAR && p.material.bsdfs[i].type == chosen_bsdf_type)
		{
			*pdf_value += p.material.bsdfs[i].pdf(p, outgoing, direction);
		}
	}
	*pdf_value /= n;
	return direction;
}

int max_depth = 8; //NOTE: Arbitrarily chosen

Radiance cast_ray(Scene* scene, Ray eye_ray)
{
	Radiance eye_ray_radiance = {};
	Radiance direct_contribution = {};
	Ray outgoing = eye_ray;
	Ray incoming = {};
	Spectrum f = generate_constant_spd(1.0);
	double dir_pdf = 0.0;
	bool consider_emissive = true;
	for(int depth = 0; depth < max_depth; ++depth)
	{
		Geometry_Intersection_Point gp = find_ray_scene_intersection(scene, outgoing);
		Surface_Point p = {};
		if(gp.scene_object >= 0)
		{
			Vec3 surface_normal = {};
			Vec2 texture_coordinates = {};
			Scene_Object* object = scene->objects + gp.scene_object;
			Scene_Geometry* geometry = &(object->geometry);
			switch(geometry->type)
			{
				case GEO_TYPE_SPHERE: 
				{
					surface_normal = normalise(gp.position - geometry->sphere.center);
					//TODO: texture_coords = ?
					break;
				}
				case GEO_TYPE_PLANE: 
				{
					surface_normal = geometry->plane.n;
					//TODO: texture_coords = ?
					break;
				}
				case GEO_TYPE_MODEL:
				{
					Model m = geometry->model;
					int i = gp.model_triangle_index;
					surface_normal = gp.bc_a * m.vertices[i].normal + gp.bc_b * m.vertices[i+1].normal + gp.bc_c * m.vertices[i+2].normal;
					texture_coordinates = gp.bc_a * m.vertices[i].texture_coords + gp.bc_b * m.vertices[i+1].texture_coords + gp.bc_c * m.vertices[i+2].texture_coords;
					break;
				}
			}

			p.name = object->name;
			p.material = object->material;
			if(object->material.emission_spd_texture.in_use)
			{
				p.material.emission_spd = *(Spectrum*)(sample_texture(object->material.emission_spd_texture, texture_coordinates));
			}
			if(object->material.diffuse_spd_texture.in_use)
			{
				p.material.diffuse_spd = *TEXTURE_SAMPLE(Spectrum, object->material.diffuse_spd_texture, texture_coordinates);
			}
			if(object->material.glossy_spd_texture.in_use)
			{
				p.material.glossy_spd = *(Spectrum*)(sample_texture(object->material.glossy_spd_texture, texture_coordinates));
			}
			if(object->material.shininess_texture.in_use)
			{
				p.material.shininess = *(double*)(sample_texture(object->material.shininess_texture, texture_coordinates));
			}
			if(object->material.refract_index_texture.in_use)
			{
				p.material.refract_index = *(Spectrum*)(sample_texture(object->material.refract_index_texture, texture_coordinates));
			}
			if(object->material.extinct_index_texture.in_use)
			{
				p.material.extinct_index = *(Spectrum*)(sample_texture(object->material.extinct_index_texture, texture_coordinates));
			}
			if(object->material.roughness_texture.in_use)
			{
				p.material.roughness = *(double*)(sample_texture(object->material.roughness_texture, texture_coordinates));
			}

			//If object transmits and the light is incident to the object
			if(dot(outgoing.direction, surface_normal) < 0.0)
			{
				p.incident_refract_index = generate_constant_spd(1.0);
				p.transmit_refract_index = p.material.refract_index;
			}
			//Else if object transmits and the light is transmitting through the object
			else
			{
				p.incident_refract_index = p.material.refract_index;
				p.transmit_refract_index = generate_constant_spd(1.0);
			}

			p.normal = surface_normal;
			p.position = gp.position;
			p.exists = true;
			p.is_emissive = object->is_emissive;
		}

		if(p.exists && p.is_emissive && consider_emissive)
		{
			eye_ray_radiance += f * p.material.emission_spd;
			break;
		}
		else if(p.exists && !p.is_emissive)
		{
			outgoing.direction = -outgoing.direction; //Reverse for bsdf computation, needs to start other way round for intersection test
			eye_ray_radiance += f * direct_light_contribution(scene, p, outgoing);
			
			//Choose new incoming direction
			incoming.direction = choose_incoming_direction(p, outgoing.direction, &dir_pdf, &consider_emissive);
			incoming.origin = p.position + 0.001*incoming.direction;
			
			if(dir_pdf <= 0.0) break;

			//Compute new direction pdf value
			f *= (abs(dot(p.normal, incoming.direction)/dir_pdf)) * bsdf(p, incoming.direction, outgoing.direction);
			
			outgoing = incoming;
			outgoing.origin += 0.001*outgoing.direction;
		}
		else break;
	}

	if(eye_ray_radiance.samples[0] < 0.0) 
	{
		printf("RADIANCE IS NEGATIVE\n");
	}
	else if(isnan(eye_ray_radiance.samples[0])) 
	{
		printf("RADIANCE IS NAN\n");
	}
	return eye_ray_radiance;
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
	Spectrum light_spd = generate_constant_spd(1.0);
	Spectrum white_diffuse_spd = RGB64_to_spectrum(RGB64{0.55, 0.55, 0.55});
	Spectrum white_glossy_spd = RGB64_to_spectrum(RGB64{0.7, 0.7, 0.7});
	Spectrum red_diffuse_spd = RGB64_to_spectrum(RGB64{0.5, 0.0, 0.0});
	Spectrum red_glossy_spd = RGB64_to_spectrum(RGB64{0.7, 0.6, 0.6});
	Spectrum green_diffuse_spd = RGB64_to_spectrum(RGB64{0.1, 0.35, 0.1});
	Spectrum green_glossy_spd = RGB64_to_spectrum(RGB64{0.45, 0.55, 0.45});
	Spectrum blue_diffuse_spd = RGB64_to_spectrum(RGB64{0.2, 0.2, 0.8});
	Spectrum blue_glossy_spd = RGB64_to_spectrum(RGB64{0.8, 0.8, 0.9});
	Spectrum orange_diffuse_spd = RGB64_to_spectrum(RGB64{1.0, 0.64, 0.0});
	Spectrum orange_glossy_spd = RGB64_to_spectrum(RGB64{1.0, 0.8, 0.0});

	Spectrum sphere_diffuse_spd = RGB64_to_spectrum(RGB64{0.2, 0.8, 0.8});
	Spectrum sphere_glossy_spd = RGB64_to_spectrum(RGB64{0.8, 0.9, 0.9});

	Material white_material = create_plastic(white_diffuse_spd, white_glossy_spd, 32.0);
	Material red_material = create_plastic(red_diffuse_spd, red_glossy_spd, 32.0);
	Material green_material = create_plastic(green_diffuse_spd, green_glossy_spd, 32.0);
	Material blue_material = create_plastic(blue_diffuse_spd, blue_glossy_spd, 32.0);
	Material sphere_material = create_plastic(sphere_diffuse_spd, sphere_glossy_spd, 50.0);
	Material mirror = create_mirror();
	Spectrum gold_refract_index = load_spd("au_spec_n.csv");
	Spectrum gold_extinct_index = load_spd("au_spec_k.csv");
	Material gold = create_conductor(gold_refract_index, gold_extinct_index, 0.344);
	Spectrum glass_refract_index = load_spd("glass.csv");
	Material glass = create_dielectric(glass_refract_index);


	Sphere glass_sphere =
	{
		{1.5, -1.8, 2.0}, 1.0
	};

	Sphere gold_sphere = 
	{
		{-2.0, -2.2, -0.5}, 0.75
	};

	Sphere mirror_sphere =
	{
		{0.0, 1.0, -1.0}, 1.0
	};

	Sphere plastic_sphere =
	{
		{2.0, -1.0, -2.0}, 1.0
	};

	Plane mirror_plane = create_plane_from_points(Vec3{-1.0, 1.0, -2.4}, Vec3{1.0, 1.0, -2.4}, Vec3{-1.0, -1.0, -2.9});

	Point light_p = {{0.0, 2.0, 2.0}};
	Sphere light_s = 
	{
		{0.0, 2.5, 2.0}, 0.5
	};

	double h = 3.0;
	Plane back_wall = create_plane_from_points(Vec3{-h, h, -h}, Vec3{h, h, -h}, Vec3{-h, -h, -h});
	Plane left_wall = create_plane_from_points(Vec3{-h, h, h}, Vec3{-h, h, -h}, Vec3{-h, -h, h});
	Plane right_wall = create_plane_from_points(Vec3{h, h, -h}, Vec3{h, h, h}, Vec3{h, -h, -h});
	Plane front_wall = create_plane_from_points(Vec3{h, h, h}, Vec3{-h, h, h}, Vec3{h, -h, h});
	Plane floor = create_plane_from_points(Vec3{-h, -h, -h}, Vec3{h, -h, -h}, Vec3{-h, -h, h});
	Plane ceiling = create_plane_from_points(Vec3{-h, h, h}, Vec3{h, h, h}, Vec3{-h, h, -h});

	Plane light_plane = create_plane_from_points(Vec3{-0.5, 2.9, 0.5}, Vec3{0.5, 2.9, 0.5}, Vec3{-0.5, 2.9, -0.5});

	Model triangle = create_model();
	Material triangle_material = create_textured_plastic(create_default_spd_texture(), white_glossy_spd, 32.0);
	//add_sphere_light_to_scene(scene, light_s, 10.0 * light_spd);
	add_plane_light_to_scene(scene, light_plane, light_spd, "Scene light");
	//add_point_light_to_scene(scene, light_p, 64.0 * light_spd);
	//add_sphere_to_scene(scene, mirror_sphere, mirror);
	//add_sphere_to_scene(scene, glass_sphere, glass);
	//add_sphere_to_scene(scene, gold_sphere, gold);
	//add_sphere_to_scene(scene, plastic_sphere, sphere_material);
	//add_plane_to_scene(scene, mirror_plane, mirror);

	add_plane_to_scene(scene, back_wall, blue_material, "Back wall");
	add_plane_to_scene(scene, left_wall, red_material, "Left wall");
	add_plane_to_scene(scene, right_wall, green_material, "Right wall");
	add_plane_to_scene(scene, floor, white_material, "Floor");
	add_plane_to_scene(scene, ceiling, white_material, "Ceiling");

	add_model_to_scene(scene, triangle, triangle_material, "Model");
}

int number_of_render_samples = 256;
double total_render_time = 0.0;
double average_sample_render_time = 0.0;
double max_sample_render_time = 0.0;
double min_sample_render_time = DBL_MAX;

void print_render_profile()
{
	printf("Time to render image: %fms\n", total_render_time);
	printf("Average sample render time: %fms\n", average_sample_render_time);
	printf("Max sample render time: %fms\n", max_sample_render_time);
	printf("Min sample render time: %fms\n", min_sample_render_time);
}

void render_image(Texture* render_target)
{
	printf("Starting render...\n");
	Timer timer = {};

	Spectrum clear_spectrum = {};
	Texture spectrum_buffer = TEXTURE_CREATE(Spectrum, render_target->width, render_target->height);
	
	TEXTURE_CLEAR(spectrum_buffer, clear_spectrum);

	Scene scene = {};
	srand(100000);
	load_colour_data();
	load_scene(&scene);

	double fov = 90.0;
	double focal_length = 0.5;
	double focal_depth = 8.0;
	double aperture_radius = 0.0;

	Vec3 image_plane_position = {0.0, 0.0, 8.0};
	Vec3 forward = {0.0, 0.0, -1.0};
	Vec3 right = {1.0, 0.0, 0.0};
	Vec3 up = {0.0, 1.0, 0.0};

	int image_plane_width_px = spectrum_buffer.width;
	int image_plane_height_px = spectrum_buffer.height;

	double aspect_ratio = (double)(image_plane_width_px)/(double)(image_plane_height_px);
	
	double image_plane_aperture_distance = (focal_length * focal_depth) / (focal_length + focal_depth);
	double image_plane_width = 2.0 * image_plane_aperture_distance * tan_deg(fov/2.0);
	double image_plane_height = image_plane_width / aspect_ratio;
	
	double pixel_width = image_plane_width/(double)(image_plane_width_px);
	double pixel_height = image_plane_height/(double)(image_plane_height_px);

	Vec3 image_plane_top_left = image_plane_position - 0.5 * image_plane_width * right + 0.5 * image_plane_height * up;

	Vec3 pinhole_position = image_plane_position + image_plane_aperture_distance * forward;

	for(int pass = 0; pass < number_of_render_samples; ++pass)
	{
		int i = pass % 4;

		printf("Pass %d/%d\r", pass+1, number_of_render_samples);
		start_timer(&timer);

		Spectrum pixel_spectrum = {};
		Vec3 pixel_top_left = {};
		Vec3 sampled_pixel_point = {};
		Ray eye_ray = {};

		for(int y = 0; y < image_plane_height_px; ++y)
		{
			for(int x = 0; x < image_plane_width_px; ++x)
			{
				DEBUG(debug_set_pixel(x, y));
				pixel_top_left = image_plane_top_left + ((double)x)*pixel_width*right - ((double)y)*pixel_height*up;
				Vec3 x_pixel_sample = (0.5 * pixel_width + (double)(i%2) * pixel_width) * Vec3{1.0, 0.0, 0.0};
				Vec3 y_pixel_sample = (0.5 * pixel_height + (double)(i / 2) * pixel_height) * Vec3{0.0, 1.0, 0.0};
				sampled_pixel_point = pixel_top_left + x_pixel_sample + y_pixel_sample;
				if(aperture_radius > 0)
				{
					Vec3 plane_of_focus_point = sampled_pixel_point + focal_depth * normalise(pinhole_position - sampled_pixel_point);
					Mat3x3 r = find_rotation_between_vectors(forward, Vec3{0.0, 0.0, 1.0});
					Vec3 sampled_lens_point = pinhole_position + r * ((focal_length / aperture_radius) * uniform_sample_disc());
					eye_ray.origin = sampled_lens_point;
					eye_ray.direction = normalise(plane_of_focus_point - eye_ray.origin);
				}
				else
				{
					eye_ray.origin = sampled_pixel_point;
					eye_ray.direction = normalise(pinhole_position - eye_ray.origin);
				}
				pixel_spectrum = cast_ray(&scene, eye_ray);			
				Spectrum* target_pixel = TEXTURE_READ(Spectrum, spectrum_buffer, x, y);
				*target_pixel = ((double)pass * *target_pixel + pixel_spectrum)/(double)(pass+1);
				RGB8 pixel_colour = rgb64_to_rgb8(spectrum_to_RGB64(*target_pixel));
				TEXTURE_WRITE(*render_target, x, y, pixel_colour);
			}
		}

		stop_timer(&timer);

		double elapsed = elapsed_time_in_ms(&timer);
		total_render_time += elapsed;
		average_sample_render_time = ((double)pass * average_sample_render_time + elapsed)/(double)(pass + 1);
		max_sample_render_time = d_max(max_sample_render_time, elapsed);
		min_sample_render_time = d_min(min_sample_render_time, elapsed);
	}
	printf("Render completed\n");
}

