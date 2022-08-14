
void add_object_to_scene(Scene* scene, Scene_Object object)
{
	scene->objects[scene->number_of_objects++] = object;
}

void add_sphere_light_to_scene(Scene* scene, Sphere s, Spectrum emission_spd)
{
	Scene_Object sphere_light = {};
	Scene_Geometry sphere_geometry = {};
	sphere_geometry.type = GEO_TYPE_SPHERE;
	sphere_geometry.sphere = s;

	sphere_light.geometry = sphere_geometry;
	sphere_light.material.emission_spd_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(sphere_light.material.emission_spd_texture, emission_spd);
	sphere_light.material.is_emissive = true;

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
	plane_light.material.emission_spd_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(plane_light.material.emission_spd_texture, emission_spd);
	plane_light.material.is_emissive = true;

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

Geometry_Intersection_Point find_ray_scene_intersection(Scene* scene, Ray ray)
{
	TIMED_FUNCTION;
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
					ray_intersects_ith_object = ray_intersects_model(ray, ith_object_geometry->model, &ith_length_along_ray, &intersection_point.model_triangle_index, intersection_point.barycentric_coordinates);
					break;
				}
			}
			
			DEBUG(if(ray_intersects_ith_object)
			{
				intersection_point.position = ray.origin + length_along_ray * ray.direction;
				debug_set_last_intersection_computed(intersection_point);
			})

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

void direct_light_contribution(Scene* scene, Surface_Point& p, Ray outgoing, Radiance& contribution)
{
	TIMED_FUNCTION;
	Spectrum bdsf_result;
	set_spectrum_to_zero(bdsf_result);
	set_spectrum_to_zero(contribution);
	if(!p.surface_material->is_emissive)
	{
		for(int i = 0; i < scene->number_of_objects; ++i)
		{
			if(scene->objects[i].material.is_emissive)
			{
				double light_pdf = 1.0;
				Scene_Object& light = scene->objects[i];
				Vec3 light_point = {};
				switch(light.geometry.type)
				{
					case GEO_TYPE_SPHERE: 
					{
						light_point = uniform_sample_sphere_subtended(light.geometry.sphere, p.position, &light_pdf);
						break;
					}
					case GEO_TYPE_PLANE: 
					{
						light_point = uniform_sample_plane(light.geometry.plane, &light_pdf);
						break;
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

				if(light_point_visible)
				{
					Vec3 incoming = light_point - p.position;
					incoming = normalise(incoming);

					double d = abs(dot(incoming, p.normal));
					double f = d / light_pdf;

					bdsf(p, incoming, outgoing.direction, bdsf_result);
					Vec2 t = {0.5, 0.5};
					Spectrum& emission_spd = *TEXTURE_SAMPLE(Spectrum, light.material.emission_spd_texture, t);
					spectral_multiply(bdsf_result, emission_spd, bdsf_result);
					spectral_sum_and_multiply(contribution, bdsf_result, f, contribution);
				}
			}
		}
	}
}

//MIRROR CHANGE: Probability of sampling specular direction depends on material
Vec3 choose_incoming_direction(Surface_Point& p, Vec3 outgoing, double* pdf_value)
{
	TIMED_FUNCTION;
	double r = uniform_sample();
	double n = (double)p.surface_material->number_of_bdsfs;
	int chosen_bdsf = (int)floor(r * n);
	BDSF_TYPE chosen_bdsf_type = p.surface_material->bdsfs[chosen_bdsf].type;
	Vec3 direction = p.surface_material->bdsfs[chosen_bdsf].sample_direction(p, outgoing, pdf_value);

	//Sum probabilities for BDSFs with same distributions (that aren't specular)
	for(int i = 0; i < p.surface_material->number_of_bdsfs; ++i)
	{
		if(i != chosen_bdsf && chosen_bdsf_type != BDSF_TYPE_SPECULAR && p.surface_material->bdsfs[i].type == chosen_bdsf_type)
		{
			*pdf_value += p.surface_material->bdsfs[i].pdf(p, outgoing, direction);
		}
	}
	*pdf_value /= n;
	return direction;
}



uint8_t* sample_texture_or_default(uint8_t* default_value, Texture texture, Vec2 texture_coords)
{
	if(texture.in_use) return sample_texture(texture, texture_coords);
	else return default_value;
}

#define TEXTURE_SAMPLE_OR_DEFAULT(type, default_value, texture, texture_coords) *(type*)sample_texture_or_default((uint8_t*)&default_value, texture, texture_coords)

void find_intersection_surface_point(Scene* scene, Ray outgoing, Surface_Point& p)
{
	TIMED_FUNCTION;
	Geometry_Intersection_Point gp = find_ray_scene_intersection(scene, outgoing);
	p.exists = false;
	p.name = NULL;
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
				texture_coordinates = {0.5, 0.5};
				break;
			}
			case GEO_TYPE_PLANE: 
			{
				surface_normal = geometry->plane.n;
				//TODO: texture_coords = ?
				texture_coordinates = {0.5, 0.5};
				break;
			}
			case GEO_TYPE_MODEL:
			{
				Model m = geometry->model;
				int i = gp.model_triangle_index;
				surface_normal = gp.barycentric_coordinates.x * m.vertices[i].normal + gp.barycentric_coordinates.y * m.vertices[i+1].normal + gp.barycentric_coordinates.z * m.vertices[i+2].normal;
				texture_coordinates = gp.barycentric_coordinates.x * m.vertices[i].texture_coords + gp.barycentric_coordinates.y * m.vertices[i+1].texture_coords + gp.barycentric_coordinates.z * m.vertices[i+2].texture_coords;
				break;
			}
		}

		p.name = object->name;
		p.surface_material = &(object->material);
		if(dot(outgoing.direction, surface_normal) < 0.0)
		{
			p.incident_material = &(scene->air_material);
			p.transmit_material = &(object->material);
		}
		else
		{
			p.incident_material = &(object->material);
			p.transmit_material = &(scene->air_material);
		}

		p.normal = surface_normal;
		p.position = gp.position;
		p.texture_coordinates = texture_coordinates;
		p.exists = true;
	}
}

void cast_ray(Scene* scene, Ray eye_ray, Radiance& eye_ray_radiance)
{
int max_depth = 4; //NOTE: Arbitrarily chosen
	TIMED_FUNCTION;
	Radiance direct_contribution;
	Radiance bdsf_result;
	set_spectrum_to_zero(direct_contribution);
	set_spectrum_to_zero(bdsf_result);
	Ray outgoing = eye_ray;
	Ray incoming = {};
	Spectrum f;
	set_spectrum_to_value(f, 1.0);

	Surface_Point p = {};
	for(int depth = 0; depth < max_depth; ++depth)
	{
		double dir_pdf = 0.0;
		DEBUG(debug_set_current_eye_radiance(eye_ray_radiance);)

		find_intersection_surface_point(scene, outgoing, p);

		if(p.exists && p.surface_material->is_emissive)
		{
			Spectrum& emission_spd = *TEXTURE_SAMPLE(Spectrum, p.surface_material->emission_spd_texture, p.texture_coordinates);
			spectral_sum_and_multiply(eye_ray_radiance, f, emission_spd, eye_ray_radiance);
			break;
		}
		else if(p.exists && !p.surface_material->is_emissive)
		{
			outgoing.direction = -outgoing.direction; //Reverse for bdsf computation, needs to start other way round for intersection test
			direct_light_contribution(scene, p, outgoing, direct_contribution);
			spectral_sum_and_multiply(eye_ray_radiance, f, direct_contribution, eye_ray_radiance);
			
			//Choose new incoming direction
			incoming.direction = choose_incoming_direction(p, outgoing.direction, &dir_pdf);
			incoming.origin = p.position + 0.001*incoming.direction;

			//Compute new direction pdf value
			double f_coefficient = abs(dot(p.normal, incoming.direction)) * (1.0/dir_pdf);
			bdsf(p, incoming.direction, outgoing.direction, bdsf_result);
			spectral_multiply(f, bdsf_result, f_coefficient, f);
			
			outgoing = incoming;
		}
		else break;
	}
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

	scene->air_material.refract_index_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(scene->air_material.refract_index_texture, light_spd);
	scene->air_material.extinct_index_texture = TEXTURE_CREATE(Spectrum, 1, 1);
	TEXTURE_CLEAR(scene->air_material.extinct_index_texture, light_spd);

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
	Spectrum gold_refract_index = load_spd("spectra/au_spec_n.csv");
	Spectrum gold_extinct_index = load_spd("spectra/au_spec_k.csv");
	Material gold = create_conductor(gold_refract_index, gold_extinct_index, 0.344);
	Spectrum glass_refract_index = load_spd("spectra/glass.csv");
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
	add_sphere_to_scene(scene, glass_sphere, glass, "Glass sphere");
	add_sphere_to_scene(scene, gold_sphere, gold, "Gold sphere");
	//add_sphere_to_scene(scene, plastic_sphere, sphere_material);
	add_plane_to_scene(scene, mirror_plane, mirror, "Mirror");

	add_plane_to_scene(scene, back_wall, blue_material, "Back wall");
	add_plane_to_scene(scene, left_wall, red_material, "Left wall");
	add_plane_to_scene(scene, right_wall, green_material, "Right wall");
	add_plane_to_scene(scene, floor, white_material, "Floor");
	add_plane_to_scene(scene, ceiling, white_material, "Ceiling");

	add_model_to_scene(scene, triangle, triangle_material, "Model");

	bool not_sorted = true;
	while(not_sorted)
	{
		not_sorted = false;
		for(int i = 0; i < scene->number_of_objects-1; ++i)
		{
			if(scene->objects[i].geometry.type < scene->objects[i+1].geometry.type)
			{
				Scene_Object temp = scene->objects[i];
				scene->objects[i] = scene->objects[i+1];
				scene->objects[i+1] = temp;
				not_sorted = true;
			}
		}
	}
}

void render_image(Texture* render_target, int number_of_render_samples)
{
	TIMED_FUNCTION;
	double total_render_time = 0.0;
	double average_sample_render_time = 0.0;
	double max_sample_render_time = 0.0;
	double min_sample_render_time = DBL_MAX;
	printf("Starting render...\n");
	Timer timer = {};

	Spectrum clear_spectrum = {};
	Texture spectrum_buffer = TEXTURE_CREATE(Spectrum, render_target->width, render_target->height);
	
	TEXTURE_CLEAR(spectrum_buffer, clear_spectrum);

	Scene* scene = (Scene*)alloc(sizeof(Scene));
	srand(100000);
	load_colour_data();
	load_scene(scene);

	double fov = 90.0;

	Vec3 image_plane_position = {0.0, 0.0, 8.0};
	Camera camera;
	camera.focal_length = 0.5;
	camera.focal_depth = 8.0;
	camera.aperture_radius = 0.0;
	camera.forward = {0.0, 0.0, -1.0};
	camera.right = {1.0, 0.0, 0.0};
	camera.up = {0.0, 1.0, 0.0};

	camera.film_width_px = spectrum_buffer.width;
	camera.film_height_px = spectrum_buffer.height;

	double aspect_ratio = (double)(camera.film_width_px)/(double)(camera.film_height_px);
	
	double image_plane_aperture_distance = (camera.focal_length * camera.focal_depth) / (camera.focal_length + camera.focal_depth);
	double image_plane_width = 2.0 * image_plane_aperture_distance * tan_deg(fov/2.0);
	double image_plane_height = image_plane_width / aspect_ratio;
	
	camera.film_pixel_width = image_plane_width/(double)(camera.film_width_px);
	camera.film_pixel_height = image_plane_height/(double)(camera.film_height_px);

	camera.film_top_left = image_plane_position - 0.5 * image_plane_width * camera.right + 0.5 * image_plane_height * camera.up;

	camera.pinhole_position = image_plane_position + image_plane_aperture_distance * camera.forward;

	Spectrum pixel_spectrum;
	Vec3 pixel_top_left = {};
	Vec3 sampled_pixel_point = {};
	Ray eye_ray = {};

	for(int pass = 0; pass < number_of_render_samples; ++pass)
	{
		printf("Pass %d/%d\n", pass+1, number_of_render_samples);
		start_timer(&timer);


		for(int pixel = 0; pixel < spectrum_buffer.width * spectrum_buffer.height; ++pixel)
		{
			int y = pixel / spectrum_buffer.width;
			int x = pixel % spectrum_buffer.width;
			DEBUG(debug_set_pixel(x, y));
			set_spectrum_to_zero(pixel_spectrum);

			//Sample image plane
			pixel_top_left = camera.film_top_left + ((double)x)*camera.film_pixel_width*camera.right - ((double)y)*camera.film_pixel_height*camera.up;
			Vec3 x_pixel_sample = (0.5 * camera.film_pixel_width) * Vec3{1.0, 0.0, 0.0};
			Vec3 y_pixel_sample = (0.5 * camera.film_pixel_height) * Vec3{0.0, 1.0, 0.0};
			sampled_pixel_point = pixel_top_left + x_pixel_sample + y_pixel_sample;

			//Sample camera lens
			if(camera.aperture_radius > 0)
			{
				Vec3 plane_of_focus_point = sampled_pixel_point + camera.focal_depth * normalise(camera.pinhole_position - sampled_pixel_point);
				Mat3x3 r = find_rotation_between_vectors(camera.forward, Vec3{0.0, 0.0, 1.0});
				Vec3 sampled_lens_point = camera.pinhole_position + r * ((camera.focal_length / camera.aperture_radius) * uniform_sample_disc());
				eye_ray.origin = sampled_lens_point;
				eye_ray.direction = normalise(plane_of_focus_point - eye_ray.origin);
			}
			else
			{
				eye_ray.origin = sampled_pixel_point;
				eye_ray.direction = normalise(camera.pinhole_position - eye_ray.origin);
			}

			//Sample scene
			cast_ray(scene, eye_ray, pixel_spectrum);

			//Filter scene sample with pixel buffer contents
			Spectrum* target_pixel = TEXTURE_READ(Spectrum, spectrum_buffer, x, y);
			double pass_d = (double)pass;
			double pass_inc_d = (double)(pass+1);
			spectral_sum_and_multiply(pixel_spectrum, *target_pixel, pass_d, *target_pixel);
			spectral_multiply(*target_pixel, 1.0/pass_inc_d, *target_pixel);
			//*target_pixel = ((double)pass * *target_pixel + pixel_spectrum)/(double)(pass+1);
		}

		stop_timer(&timer);

		double elapsed = elapsed_time_in_ms(&timer);
		total_render_time += elapsed;
		average_sample_render_time = ((double)pass * average_sample_render_time + elapsed)/(double)(pass + 1);
		max_sample_render_time = d_max(max_sample_render_time, elapsed);
		min_sample_render_time = d_min(min_sample_render_time, elapsed);
	}

	for(int y = 0; y < camera.film_height_px; ++y)
	{
		for(int x = 0; x < camera.film_width_px; ++x)
		{
			Spectrum* spectrum_pixel = TEXTURE_READ(Spectrum, spectrum_buffer, x, y);
			RGB8 pixel_colour = rgb64_to_rgb8(spectrum_to_RGB64(*spectrum_pixel));
			TEXTURE_WRITE(*render_target, x, y, pixel_colour);
		}
	}

	printf("Render completed\n");
	printf("Time to render image: %fms\n", total_render_time);
	printf("Average sample render time: %fms\n", average_sample_render_time);
	printf("Max sample render time: %fms\n", max_sample_render_time);
	printf("Min sample render time: %fms\n", min_sample_render_time);
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

	NEW_Surface_Point surface_point = {};

	while(is_there_next_point && depth <= max_depth)
	{
		//Sample incoming light rays
		for(int i = 0; depth != 0 && i < scene->number_of_light_sources; ++i)
		{
			NEW_Geometry_Sample* light_sample = &current_point->light_contributions[current_point->number_of_light_contributions];
			NEW_sample_geometry(scene->light_source_geometries[i], current_point->position, light_sample);
			light_sample->object = i;
			
			Ray shadow_ray = {};
			shadow_ray.direction = normalise(light_sample->position - current_point->position);
			shadow_ray.origin = current_point->position + 0.0009765625*shadow_ray.direction;
			
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
				light_point_visible = !(shadow_test_ray_length_sq - light_test_ray_length_sq < 0.0);
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

		current_point->out_ray.direction = current_point->material->bdsfs[sampled_bdsf].sample_direction(&surface_point, current_point->in_ray.direction, &current_point->pdf);
		current_point->out_ray.origin = current_point->position + 0.0009765625*current_point->out_ray.direction;
		NEW_find_ray_scene_intersection(scene, current_point->out_ray, next_point);
		++depth;

		//Store and iterate points
		current_point = next_point;
		is_there_next_point = !next_point->material->is_emissive && next_point->object != -1;
		++next_point;
	}
	--depth;
	if(!is_there_next_point && depth <= max_depth) scene_path[depth].light_contributions[scene_path[depth].number_of_light_contributions++] = *((NEW_Geometry_Sample*)(&current_point->position));
	
	return depth;
}

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
	scene->materials[0].bdsfs[1].reflectance = NEW_glossy_phong_reflectance;
	scene->materials[0].bdsfs[1].sample_direction = NEW_sample_glossy_direction;

	//left wall (green)
	scene->materials[1] = {};
	scene->materials[1].shininess = 32.0;
	NEW_RGB64_to_spectrum(RGB64{0.1, 0.35, 0.1}, &scene->materials[1].diffuse_spd, spectrum_buffer);
	NEW_RGB64_to_spectrum(RGB64{0.45, 0.55, 0.45}, &scene->materials[1].glossy_spd, spectrum_buffer);
	scene->materials[1].number_of_bdsfs = 2;
	scene->materials[1].bdsfs[0].reflectance = NEW_diffuse_phong_reflectance;
	scene->materials[1].bdsfs[0].sample_direction = NEW_sample_diffuse_direction;
	scene->materials[1].bdsfs[1].reflectance = NEW_glossy_phong_reflectance;
	scene->materials[1].bdsfs[1].sample_direction = NEW_sample_glossy_direction;
	
	//right wall (red)
	scene->materials[2] = {};
	scene->materials[2].shininess = 32.0;
	NEW_RGB64_to_spectrum(RGB64{0.5, 0.0, 0.0}, &scene->materials[2].diffuse_spd, spectrum_buffer);
	NEW_RGB64_to_spectrum(RGB64{0.7, 0.6, 0.6}, &scene->materials[2].glossy_spd, spectrum_buffer);
	scene->materials[2].number_of_bdsfs = 2;
	scene->materials[2].bdsfs[0].reflectance = NEW_diffuse_phong_reflectance;
	scene->materials[2].bdsfs[0].sample_direction = NEW_sample_diffuse_direction;
	scene->materials[2].bdsfs[1].reflectance = NEW_glossy_phong_reflectance;
	scene->materials[2].bdsfs[1].sample_direction = NEW_sample_glossy_direction;

	//floor + ceiling (white)
	scene->materials[3] = {};
	scene->materials[3].shininess = 32.0;
	NEW_RGB64_to_spectrum(RGB64{0.55, 0.55, 0.55}, &scene->materials[3].diffuse_spd, spectrum_buffer);
	NEW_RGB64_to_spectrum(RGB64{0.7, 0.7, 0.7}, &scene->materials[3].glossy_spd, spectrum_buffer);
	scene->materials[3].number_of_bdsfs = 2;
	scene->materials[3].bdsfs[0].reflectance = NEW_diffuse_phong_reflectance;
	scene->materials[3].bdsfs[0].sample_direction = NEW_sample_diffuse_direction;
	scene->materials[3].bdsfs[1].reflectance = NEW_glossy_phong_reflectance;
	scene->materials[3].bdsfs[1].sample_direction = NEW_sample_glossy_direction;

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

	//Render image
	for(int pass = 0; pass < number_of_render_samples; ++pass)
	{
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
			double first_point_depth = length(scene_path[1].position - scene->camera.pinhole_position)/12.0;	

			NEW_set_spectrum_to_value(&spectral_render_target[pixel], first_point_depth);
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
