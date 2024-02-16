void load_scene(const char *scene_path, camera_data *camera, scene_data *scene, u32 width_px, u32 height_px)
{
    HANDLE scene_file_handle = CreateFile(scene_path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);

    DWORD bytes_read;
    DWORD scene_file_size    = GetFileSize(scene_file_handle, NULL);
    char  *scene_file_buffer = (char*)VirtualAlloc(0, scene_file_size, MEM_COMMIT, PAGE_READWRITE);
    ReadFile(scene_file_handle, scene_file_buffer, scene_file_size, &bytes_read, NULL);
    CloseHandle(scene_file_handle);

    camera_input_data camera_input;
    scene_input_data  scene_input;
    memset(&scene_input, 0, sizeof(scene_input));
    memset(&camera_input, 0, sizeof(camera_input));
    camera_input.width_px  = width_px;
    camera_input.height_px = height_px;
    parse_scene(scene_file_buffer, scene_file_size, &camera_input, &scene_input);
    init_camera(camera, &camera_input);
    init_scene(scene, &scene_input);
}

void init_camera(camera_data *camera, camera_input_data *input)
{
    camera->forward = input->forward;
    camera->right   = input->right;
    camera->up      = input->up;

    f64 aperture_distance     = (input->flength * input->fdepth) / (input->flength + input->fdepth);
    camera->aperture_position = vec3_sum(input->position, vec3_mul_by_f64(input->forward, aperture_distance));
    camera->aperture_radius   = input->aperture;

    f64  fov_rad             = input->fov * (PI/180.0);
    f64  aspect_ratio        = (f64)input->width_px / (f64)input->height_px;
    f64  film_width          = 2.0 * aperture_distance * tan(fov_rad/2.0);
    f64  film_height         = film_width / aspect_ratio;
    vec3 film_right          = vec3_mul_by_f64(input->right, 0.5 * film_width);
    vec3 film_top            = vec3_mul_by_f64(input->up,    0.5 * film_height);
    camera->film_bottom_left = vec3_sub(vec3_sub(input->position, film_right), film_top);

    camera->pixel_width  = film_width  / (f64)input->width_px;
    camera->pixel_height = film_height / (f64)input->height_px;
}

void init_spd(spectrum *dst, spd_input_data *input, spectrum white, spectrum rgb_red, spectrum rgb_green, spectrum rgb_blue, spectrum rgb_cyan, spectrum rgb_magenta, spectrum rgb_yellow)
{
    *dst = alloc_spd();
    switch(input->method)
    {
        case SPD_METHOD_RGB:
        {
            rgb_f64_to_spectrum(input->rgb, *dst, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
            break;
        }
        case SPD_METHOD_CSV:
        {
            load_csv_file_to_spectrum(*dst, input->csv);
            break;
        }
        case SPD_METHOD_BLACKBODY:
        {
            generate_blackbody_spectrum(*dst, input->blackbody_temp);
            spectrum_normalise(*dst);
            break;
        }
        default:
        {
            free_spd(*dst);
            dst->samples = NULL;
            return;
        }
    }
    if(input->has_scale_factor)
    {
        spectral_mul_by_scalar(*dst, *dst, input->scale_factor);
    }
}

void init_scene(scene_data* scene, scene_input_data *scene_input)
{
    scene->num_surfaces        = scene_input->num_surfaces;
    scene->num_scene_materials = scene_input->num_scene_materials;

    u32 surfaces_size         = scene->num_surfaces * sizeof(object_geometry);
    u32 material_indices_size = scene->num_surfaces * sizeof(u32);
    u32 surfaces_buffer_size  = surfaces_size + material_indices_size;
    u32 materials_buffer_size = scene->num_scene_materials * sizeof(object_material);
    char *surfaces_buffer  = VirtualAlloc(NULL, surfaces_buffer_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    char *materials_buffer = VirtualAlloc(NULL, materials_buffer_size, MEM_COMMIT  | MEM_RESERVE, PAGE_READWRITE);

    scene->surfaces = (object_geometry*)surfaces_buffer;
    scene->surface_material_indices = (u32*)(surfaces_buffer + surfaces_size);
    scene->scene_materials = (object_material*)materials_buffer;

    spectrum white       = alloc_spd();
    spectrum rgb_red     = alloc_spd();
    spectrum rgb_green   = alloc_spd();
    spectrum rgb_blue    = alloc_spd();
    spectrum rgb_cyan    = alloc_spd();
    spectrum rgb_magenta = alloc_spd();
    spectrum rgb_yellow  = alloc_spd();
    const_spectrum(white, 1.0); //TODO: Parameterise
    load_csv_file_to_spectrum(rgb_red, "spectra\\red_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_green, "spectra\\green_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_blue, "spectra\\blue_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_cyan, "spectra\\cyan_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_magenta, "spectra\\magenta_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_yellow, "spectra\\yellow_rgb_to_spd.csv");
    for(u32 i = 0; i < scene->num_scene_materials; i += 1)
    {
        material_input_data *input = &scene_input->scene_materials[i];
        object_material     *dst   = &scene->scene_materials[i];

        u32 name_len = strlen(input->name);
        memcpy(dst->name, input->name, name_len);
        dst->is_black_body = input->is_black_body;
        dst->is_emissive   = input->is_emissive;
        dst->shininess     = input->shininess;

        init_spd(&dst->emission_spd, &input->emission_input, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
        init_spd(&dst->diffuse_spd, &input->diffuse_input, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
        init_spd(&dst->glossy_spd, &input->glossy_input, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
    }
    free_spd(white);
    free_spd(rgb_red);
    free_spd(rgb_green);
    free_spd(rgb_blue);
    free_spd(rgb_cyan);
    free_spd(rgb_magenta);
    free_spd(rgb_yellow);

    for(u32 i = 0; i < scene->num_surfaces; i += 1)
    {
        surface_input_data *input = &scene_input->surfaces[i];
        object_geometry    *dst   = &scene->surfaces[i];

        u32 name_len = strlen(input->name);
        memcpy(dst->name, input->name, name_len);
        dst->type     = input->type;
        dst->position = input->position;
        
        switch(dst->type)
        {
            case GEO_TYPE_SPHERE:
            {
                dst->radius = input->radius;
                break;
            }
            case GEO_TYPE_PLANE:
            {
                create_plane_from_points(input->position, input->u, input->v, &dst->position, &dst->u, &dst->v, &dst->normal);
                break;
            }
        }

        for(u32 j = 0; j < scene->num_scene_materials; j += 1)
        {
            if(strcmp(input->material_name, scene->scene_materials[j].name) == 0)
            {
                scene->surface_material_indices[i] = j;
                break;
            }
        }
    }
}

f64 f64_max(f64 f0, f64 f1)
{
    return (f0 > f1) ? f0 : f1;
}

void blinn_phong_diffuse_bdsf(spectrum reflectance, scene_point *p, vec3 incoming, vec3 outgoing)
{
    spectral_mul_by_scalar(reflectance, p->material->diffuse_spd, 1.0/PI);
}

void blinn_phong_glossy_bdsf(spectrum reflectance, scene_point *p, vec3 incoming, vec3 outgoing)
{
    vec3 bisector         = vec3_normalise(vec3_sum(outgoing, incoming));
    f64  spec_coefficient = pow(f64_max(0.0, vec3_dot(p->normal, bisector)), p->material->shininess);

    spectral_mul_by_scalar(reflectance, p->material->glossy_spd, spec_coefficient);
}

//Out = back towards camera
//In  = towards light source/next intersection
void bdsf(spectrum reflectance, scene_point *p, vec3 in, vec3 out)
{
    spectrum bdsf_result = alloc_spd();
    zero_spectrum(reflectance);
    blinn_phong_diffuse_bdsf(bdsf_result, p, in, out);
    spectral_sum(reflectance, bdsf_result, reflectance);
    blinn_phong_glossy_bdsf(bdsf_result, p, in, out);
    spectral_sum(reflectance, bdsf_result, reflectance);
    free_spd(bdsf_result);
}

//Used in cases where one point is on a light source (sampled)
//Instead of moving non-light point outside of object:
//  - Test object point's normal to see if it points away from the light
//  - Do the normal test for everything else
//Maybe even don't call this function if the light source is on the other side of the object
//Assumes p0 is from scene_point
u32 points_mutually_visible(vec3 p0, vec3 p1, scene_data *scene)
{
    u32 visible = 1;
    vec3 ray_dir_full  = vec3_sub(p1, p0);
    vec3 ray_direction = vec3_normalise(ray_dir_full);
    vec3 ray_origin    = vec3_sum(p0, vec3_mul_by_f64(ray_direction, 0.001));
    f64 vis_dist = vec3_length(vec3_sub(p1, ray_origin)) - 0.001;
    f64 dist     = INFINITY;
    for(u32 i = 0; i < scene->num_surfaces; i += 1)
    {
        object_geometry *surface = &scene->surfaces[i];
        switch(surface->type)
        {
            case GEO_TYPE_POINT: continue;
            case GEO_TYPE_SPHERE:
            {
                dist = line_sphere_intersection(ray_origin, ray_direction, surface->position, surface->radius);
                break;
            };
            case GEO_TYPE_PLANE:
            {
                dist = line_plane_intersection(ray_origin, ray_direction, surface->position, surface->normal, surface->u, surface->v);
                break;
            };
        }
        if(dist < vis_dist)
        {
            visible = 0;
            break;
        }
    }
    return visible;
}

void direct_light_contribution(spectrum contribution, scene_point *intersection, scene_data *scene, vec3 ray_direction)
{
    spectrum reflectance = alloc_spd();
    zero_spectrum(reflectance);
    zero_spectrum(contribution);
    for(u32 i = 0; i < scene->num_surfaces; i += 1)
    {
        object_material *light_material = &scene->scene_materials[scene->surface_material_indices[i]];
        if(light_material->is_emissive)
        {
            f64 light_pdf;
            f64 attenuation_factor = 1.0;
            object_geometry *light_surface = &scene->surfaces[i];
            vec3 light_position; 
            switch(light_surface->type)
            {
                case GEO_TYPE_POINT:
                {
                    light_position = light_surface->position;
                    f64 dist = vec3_length(vec3_sub(light_position, intersection->position));
                    light_pdf = 1.0;
                    attenuation_factor = 1.0 / (4.0 * PI * dist * dist);
                    break;
                }
                case GEO_TYPE_SPHERE:
                {
                    f64 u = rng();
                    f64 v = rng();
                    f64 r = sqrt(1.0 - u * u);
                    f64 t = 2.0 * PI * v;
                    vec3 sphere_point = {r * cos(t), r * sin(t), u};
                    light_position = vec3_sum(light_surface->position, vec3_mul_by_f64(sphere_point, light_surface->radius));
                    light_pdf = 1.0 / (4.0 * PI * light_surface->radius * light_surface->radius);
                    break;
                }
                case GEO_TYPE_PLANE:
                {
                    f64 u = rng();
                    f64 v = rng();
                    vec3 u_pos = vec3_mul_by_f64(light_surface->u, u);
                    vec3 v_pos = vec3_mul_by_f64(light_surface->v, v);
                    light_position = vec3_sum(vec3_sum(light_surface->position, u_pos), v_pos);
                    light_pdf = 1.0/vec3_length(vec3_cross(light_surface->u, light_surface->v));
                    break;
                }
            }
            if(points_mutually_visible(intersection->position, light_position, scene))
            {
                vec3 outgoing = vec3_reverse(ray_direction);
                vec3 incoming = vec3_normalise(vec3_sub(light_position, intersection->position));

                bdsf(reflectance, intersection, incoming, outgoing);
                spectral_sum(contribution, contribution, reflectance);
                spectral_mul_by_spectrum(contribution, contribution, light_material->emission_spd);

                f64 d = vec3_dot(incoming, intersection->normal);
                f64 c = fabs(d) * attenuation_factor * (1.0/light_pdf);
                spectral_mul_by_scalar(contribution, contribution, c);
            }
        }
    }
    free_spd(reflectance);
}

void find_ray_intersection(scene_point *intersection, scene_data *scene, vec3 ray_origin, vec3 ray_direction)
{
    f64 min_dist = INFINITY;
    object_geometry *intersection_surface = NULL;
    u32 intersection_index = -1;
    ray_origin = vec3_sum(ray_origin, vec3_mul_by_f64(ray_direction, 0.001));
    for(u32 i = 0; i < scene->num_surfaces; i += 1)
    {
        f64 dist = INFINITY;
        object_geometry *surface = &scene->surfaces[i];
        switch(surface->type)
        {
            case GEO_TYPE_POINT: continue;
            case GEO_TYPE_SPHERE:
            {
                dist = line_sphere_intersection(ray_origin, ray_direction, surface->position, surface->radius);
                break;
            }
            case GEO_TYPE_PLANE:
            {
                dist = line_plane_intersection(ray_origin, ray_direction, surface->position, surface->normal, surface->u, surface->v);
                break;
            }
        }
        if(dist < min_dist)
        {
            min_dist = dist;
            intersection_surface = surface;
            intersection_index = i;
        }
    }
    if(intersection_surface)
    {
        intersection->position = vec3_sum(ray_origin, vec3_mul_by_f64(ray_direction, min_dist));
        switch(intersection_surface->type)
        {
            case GEO_TYPE_SPHERE:
            {
                intersection->normal = vec3_normalise(vec3_sub(intersection->position, intersection_surface->position));
                break;
            }
            case GEO_TYPE_PLANE:
            {
                intersection->normal = intersection_surface->normal;
                break;
            }
        }
        intersection->material = &scene->scene_materials[scene->surface_material_indices[intersection_index]];
        intersection->is_black_body = intersection->material->is_black_body;
        intersection->is_emissive   = intersection->material->is_emissive;
        intersection->surface = &scene->surfaces[intersection_index];
    }
    else
    {
        intersection->is_black_body = 1;
        intersection->is_emissive   = 0;
    }
}

u32 estimate_indirect_contribution(spectrum contribution, scene_point *intersection, scene_data *scene, vec3 ray_origin, vec3 ray_direction)
{
    u32 ret;
    find_ray_intersection(intersection, scene, ray_origin, ray_direction);
    if(intersection->is_black_body)
    {
        if(intersection->is_emissive)
        {
            copy_spectrum(contribution, intersection->material->emission_spd);
        }
        else
        {
            const_spectrum(contribution, 0.0);
        }
        ret = 0;
    }
    else
    {
        direct_light_contribution(contribution, intersection, scene, ray_direction);
        ret = 1;
    }
    return ret;
}

vec3 cos_weighted_sample_hemisphere(vec3 n)
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
    vec3 v = mat3x3_vec3_mul(r, p);
    return v;
}

void cast_ray(spectrum dst, scene_data *scene, vec3 ray_origin, vec3 ray_direction, u32 max_depth)
{
    spectrum contribution = alloc_spd();
    spectrum throughput   = alloc_spd();
    spectrum reflectance  = alloc_spd();
    spectrum tmp_spectrum = alloc_spd();
    zero_spectrum(contribution);
    zero_spectrum(reflectance);
    const_spectrum(throughput, 1.0);

    vec3 out = ray_direction;
    vec3 in;

    scene_point intersection;
    memset(&intersection, 0, sizeof(scene_point));
    for(u32 depth = 0; depth < max_depth; depth += 1)
    {
        find_ray_intersection(&intersection, scene, ray_origin, out);

        object_material *mat = intersection.material;
        if(intersection.is_black_body && !intersection.is_emissive) break;
        else if(mat->is_black_body && mat->is_emissive)
        {
            spectral_mul_by_spectrum(tmp_spectrum, throughput, mat->emission_spd);
            spectral_sum(dst, dst, tmp_spectrum);
            break;
        }
        else
        {
            out = vec3_reverse(out);
            direct_light_contribution(contribution, &intersection, scene, out);
            spectral_mul_by_spectrum(tmp_spectrum, throughput, contribution);
            spectral_sum(dst, dst, tmp_spectrum);

#if 0
            do
            {
                in = uniform_sample_sphere();
            }
            while(vec3_dot(in, intersection.normal) <= 0.0);
            f64 dir_pdf = 1.0/(2.0*PI);
#else
            in = cos_weighted_sample_hemisphere(intersection.normal);
            f64 dir_pdf = vec3_dot(intersection.normal, in) / PI;
#endif
            f64 throughput_coefficient = fabs(vec3_dot(intersection.normal, in)) * (1.0/dir_pdf);
            bdsf(reflectance, &intersection, in, out);
            spectral_mul_by_scalar(reflectance, reflectance, throughput_coefficient);
            spectral_mul_by_spectrum(throughput, throughput, reflectance);

            out = in;
            ray_origin = intersection.position;
        }
    }
    free_spd(tmp_spectrum);
    free_spd(reflectance);
    free_spd(throughput);
    free_spd(contribution);
}

void print_camera(camera_data *camera)
{
    printf("CAMERA:\n");
    printf("Forward: "); print_vector(camera->forward);  printf("\n");
    printf("Right:   "); print_vector(camera->right);    printf("\n");
    printf("Up:      "); print_vector(camera->up);       printf("\n");
    printf("Aperture position: "); print_vector(camera->aperture_position); printf("\n");
    printf("Aperture radius: %f\n", camera->aperture_radius);
    printf("Film bottom left: "); print_vector(camera->film_bottom_left); printf("\n");
    printf("Pixel width: %f Pixel height: %f\n", camera->pixel_width, camera->pixel_height);
}

void print_material(object_material *mat)
{
    printf("Material: %s\n", mat->name);
    printf("Emission: %p\n", mat->emission_spd.samples);
    printf("Diffuse:  %p\n", mat->diffuse_spd.samples);
    printf("Glossy:   %p\n", mat->glossy_spd.samples);
    printf("Shininess:%f\n", mat->shininess);
    printf("Is black body: %u\n", mat->is_black_body);
}

void print_surface(object_geometry *surface, object_material *mat)
{
    const char *surface_type_name;
    switch(surface->type)
    {
        case GEO_TYPE_POINT:  surface_type_name = "point";  break;
        case GEO_TYPE_SPHERE: surface_type_name = "sphere"; break;
        case GEO_TYPE_PLANE:  surface_type_name = "plane";  break;
        default:              surface_type_name = "???";    break;
    }
    printf("Surface: %s\n", surface->name);
    printf("Type:    %s\n", surface_type_name);
    printf("Material:%s\n", mat->name);
    printf("Position:%f %f %f\n", surface->position.x, surface->position.y, surface->position.z);
    switch(surface->type)
    {
        case GEO_TYPE_SPHERE: printf("Radius:  %f\n", surface->radius); break;
        case GEO_TYPE_PLANE:
        {
            printf("U:       %f %f %f\n", surface->u.x, surface->u.y, surface->u.z);
            printf("V:       %f %f %f\n", surface->v.x, surface->v.y, surface->v.z);
            printf("N:       %f %f %f\n", surface->normal.x, surface->normal.y, surface->normal.z);
            break;
        }
    }
}

void print_scene(scene_data *scene)
{
    printf("SCENE:\n");
    printf("Num materials: %u Num surfaces: %u\n", scene->num_scene_materials, scene->num_surfaces);
    for(u32 i = 0; i < scene->num_scene_materials; i += 1)
    {
        print_material(&scene->scene_materials[i]); printf("\n");
    }
    for(u32 i = 0; i < scene->num_surfaces; i += 1)
    {
        print_surface(&scene->surfaces[i], &scene->scene_materials[scene->surface_material_indices[i]]);
        printf("\n");
    }
}

//TODO: Make work with non-pinhole camera
vec3 sample_camera(camera_data *camera)
{
    vec3 p = camera->aperture_position;
    return p;
}

void render_image(f64 *dst_pixels, u32 dst_width, u32 dst_height, scene_data *scene, camera_data *camera, u32 samples_per_pixel)
{
    const u32 max_cast_depth = 4;
    u32 num_pixels = dst_width * dst_height;

    u32 pixel_filter_sums_size = num_pixels * sizeof(f64);
    f64 *pixel_filter_sums     = (f64*)VirtualAlloc(NULL, pixel_filter_sums_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    spectrum contribution  = alloc_spd();

    for(u32 sample = 0; sample < samples_per_pixel; sample += 1)
    {
        printf("Sample %u / %u\n", sample+1, samples_per_pixel);
        //Filter final contribution and write to dst
        //pixel value = sum(filter * weight * radiance)/sum(filter)
        for(u32 y = 0; y < dst_height; y += 1)
        {
            for(u32 x = 0; x < dst_width; x += 1)
            {
                //Sample point on film plane
                //vec3 sampled_pixel_point = camera->film_bottom_left + film_x * camera->right + film_y * camera->up;
                f64 p_sample_x = 0.5; //TODO: Replace this with some scheme (e.g. completely random, stratified etc.)
                f64 p_sample_y = 0.5;
                f64 film_x     = ((f64)x + p_sample_x) * camera->pixel_width;
                f64 film_y     = ((f64)y + p_sample_y) * camera->pixel_height;

                vec3 sampled_pixel_point_bottom = vec3_mul_by_f64(camera->up,    film_y);
                vec3 sampled_pixel_point_left   = vec3_mul_by_f64(camera->right, film_x);
                vec3 sampled_pixel_point_bl     = vec3_sum(sampled_pixel_point_left, sampled_pixel_point_bottom);
                vec3 sampled_pixel_point        = vec3_sum(sampled_pixel_point_bl,   camera->film_bottom_left);

                //Sample camera lens
                vec3 sampled_camera_point = sample_camera(camera);
                vec3 ray_origin           = sampled_pixel_point;
                vec3 ray_direction        = vec3_normalise(vec3_sub(sampled_camera_point, ray_origin));
                
                //Only one dst pixel for now (EDIT: What does this comment mean?)
                zero_spectrum(contribution);
                cast_ray(contribution, scene, ray_origin, ray_direction, max_cast_depth);

                f64 pixel_filter_value = 1.0;
                f64 pixel_weight_value = 1.0;

                spectrum dst_pixel;
                dst_pixel.samples = dst_pixels + number_of_spectrum_samples * ((y*dst_width) + x);
                spectral_mul_by_scalar(contribution, contribution, pixel_filter_value);
                spectral_mul_by_scalar(contribution, contribution, pixel_weight_value);
                spectral_sum(dst_pixel, dst_pixel, contribution);

                f64 *dst_filter = pixel_filter_sums + (y*dst_width) + x;
                *dst_filter += pixel_filter_value;
            }
        }
    }
    for(u32 pixel = 0; pixel < num_pixels; pixel += 1)
    {
        f64 pixel_filter = 1.0 / pixel_filter_sums[pixel];
        spectrum dst_pixel;
        dst_pixel.samples = dst_pixels + pixel * number_of_spectrum_samples;
        spectral_mul_by_scalar(dst_pixel, dst_pixel, pixel_filter);
    }
}
