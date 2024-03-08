void load_scene(const char *scene_path, camera_data *camera, scene_data *scene, u32 width_px, u32 height_px)
{
    file_handle scene_file_handle = open_file(scene_path, ACCESS_READ, FILE_EXISTS);
    u32 scene_file_size = get_file_size(scene_file_handle);
    char * scene_file_buffer = alloc(scene_file_size);
    read_file(scene_file_handle, scene_file_size, scene_file_buffer);
    close_file(scene_file_handle);

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
    vec3 ref_forward = {0.0, 0.0, -1.0};
    vec3 ref_up      = {0.0, 1.0, 0.0};
    camera->forward = vec3_normalise(vec3_sub(input->target, input->position));
    mat3x3 or = find_rotation_between_vectors(ref_forward, camera->forward);
    f64 roll_rad = input->roll * (PI/180.0);
    mat3x3 rr = rotation_about_axis(camera->forward, roll_rad);
    camera->up = mat3x3_vec3_mul(rr, mat3x3_vec3_mul(or, ref_up));
    camera->right = vec3_normalise(vec3_cross(camera->forward, camera->up));

    f64 aperture_distance     = (input->flength * input->fdepth) / (input->flength + input->fdepth);
    camera->aperture_position = vec3_sum(input->position, vec3_mul_by_f64(camera->forward, aperture_distance));
    camera->aperture_radius   = input->aperture;

    f64  fov_rad             = input->fov * (PI/180.0);
    f64  aspect_ratio        = (f64)input->width_px / (f64)input->height_px;
    f64  film_width          = 2.0 * aperture_distance * tan(fov_rad/2.0);
    f64  film_height         = film_width / aspect_ratio;
    vec3 film_right          = vec3_mul_by_f64(camera->right, 0.5 * film_width);
    vec3 film_top            = vec3_mul_by_f64(camera->up,    0.5 * film_height);
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
            char csv_name[256];
            memset(csv_name, 0, 256);
            const char *prefix = "spectra\\";
            u32 prefix_size = strlen(prefix);
            memcpy(csv_name, prefix, prefix_size);
            memcpy(csv_name + prefix_size, input->csv, strlen(input->csv));
            load_csv_file_to_spectrum(*dst, csv_name);
            //spectrum_normalise(*dst); //TODO: Allow this normalise as an option
            break;
        }
        case SPD_METHOD_BLACKBODY:
        {
            generate_blackbody_spectrum(*dst, input->blackbody_temp);
            spectrum_normalise(*dst);
            break;
        }
        case SPD_METHOD_CONST:
        {
            const_spectrum(*dst, input->constant);
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
    scene->num_scene_materials = scene_input->num_scene_materials + 1;

    u32 surfaces_size         = scene->num_surfaces * sizeof(object_geometry);
    u32 material_indices_size = scene->num_surfaces * sizeof(u32);
    u32 surfaces_buffer_size  = surfaces_size + material_indices_size;
    u32 materials_buffer_size = scene->num_scene_materials * sizeof(object_material);
    char *surfaces_buffer = alloc(surfaces_buffer_size);
    char *materials_buffer = alloc(materials_buffer_size);

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
        dst->is_black_body    = (input->is_escape_material) ? 1 : input->is_black_body;
        dst->is_emissive      = input->is_emissive;
        dst->shininess        = input->shininess;
        dst->roughness        = input->roughness;

        dst->sample_direction = input->sample_direction_function;
        dst->num_bdsfs        = input->num_bdsfs;
        for(u32 j = 0; j < dst->num_bdsfs; j += 1)
        {
            dst->bdsfs[j] = input->bdsfs[j];
        }

        init_spd(&dst->emission_spd, &input->emission_input, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
        init_spd(&dst->diffuse_spd, &input->diffuse_input, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
        init_spd(&dst->glossy_spd, &input->glossy_input, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
        init_spd(&dst->mirror_spd, &input->mirror_input, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
        init_spd(&dst->refract_spd, &input->refract_input, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
        init_spd(&dst->extinct_spd, &input->extinct_input, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
        if(input->is_escape_material) scene->escape_material = dst;
        if(input->is_base_material)   scene->base_material   = dst;
    }
    if(!scene->escape_material)
    {
        //TODO
    }
    if(!scene->base_material)
    {
        //TODO
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

//Out = back towards camera
//In  = towards light source/next intersection
void bdsf(spectrum reflectance, scene_point *p, vec3 in)
{
    spectrum bdsf_result = alloc_spd();
    zero_spectrum(bdsf_result);
    zero_spectrum(reflectance);
    
    object_material *mat = p->surface_material;
    for(u32 i = 0; i < mat->num_bdsfs; i += 1)
    {
        bdsf_func f = mat->bdsfs[i];
        f(bdsf_result, p, in);
        spectral_sum(reflectance, bdsf_result, reflectance);
    }
    free_spd(bdsf_result);
}

//Used in cases where one point is on a light source (sampled)
//Instead of moving non-light point outside of object:
//  - Test object point's normal to see if it points away from the light
//  - Do the normal test for everything else
//Maybe even don't call this function if the light source is on the other side of the object
//Assumes p0 is from scene_point
#define vis_fudge 0.0001
u32 points_mutually_visible(vec3 p0, vec3 p1, scene_data *scene)
{
    u32 visible = 1;
    vec3 ray_dir_full  = vec3_sub(p1, p0);
    vec3 ray_direction = vec3_normalise(ray_dir_full);
    vec3 ray_origin    = vec3_sum(p0, vec3_mul_by_f64(ray_direction, vis_fudge));
    f64 vis_dist = vec3_length(vec3_sub(p1, ray_origin)) - vis_fudge;
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

void direct_light_contribution(spectrum contribution, scene_point *intersection, scene_data *scene)
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
                    attenuation_factor = (4.0 * PI * dist * dist);
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
                    light_pdf = (4.0 * PI * light_surface->radius * light_surface->radius);
                    break;
                }
                case GEO_TYPE_PLANE:
                {
                    f64 u = rng();
                    f64 v = rng();
                    vec3 u_pos = vec3_mul_by_f64(light_surface->u, u);
                    vec3 v_pos = vec3_mul_by_f64(light_surface->v, v);
                    light_position = vec3_sum(vec3_sum(light_surface->position, u_pos), v_pos);
                    light_pdf = vec3_length(vec3_cross(light_surface->u, light_surface->v));
                    break;
                }
            }
            if(points_mutually_visible(intersection->position, light_position, scene))
            {
                vec3 incoming = vec3_normalise(vec3_sub(light_position, intersection->position));

                bdsf(reflectance, intersection, incoming);
                spectral_sum(contribution, contribution, reflectance);
                spectral_mul_by_spectrum(contribution, contribution, light_material->emission_spd);

                f64 c = attenuation_factor * (light_pdf);
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
    ray_origin = vec3_sum(ray_origin, vec3_mul_by_f64(ray_direction, vis_fudge));
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
        intersection->out = vec3_reverse(ray_direction);
        intersection->on_dot = vec3_dot(intersection->normal, intersection->out);
        intersection->transmit_material = &scene->scene_materials[scene->surface_material_indices[intersection_index]];
        intersection->incident_material = scene->base_material;
        if(intersection->on_dot < 0.0) 
        {
            if(intersection_surface->type != GEO_TYPE_PLANE)
            {
                intersection->transmit_material = scene->base_material;
                intersection->incident_material = &scene->scene_materials[scene->surface_material_indices[intersection_index]];
            }
            intersection->normal = vec3_reverse(intersection->normal);
            intersection->on_dot = vec3_dot(intersection->normal, intersection->out);
        }
        intersection->surface_material = &scene->scene_materials[scene->surface_material_indices[intersection_index]];
        intersection->surface = &scene->surfaces[intersection_index];
    }
    else
    {
        intersection->surface_material = scene->escape_material;
    }
}

/*
u32 estimate_indirect_contribution(spectrum contribution, scene_point *intersection, scene_data *scene, vec3 ray_origin, vec3 ray_direction)
{
    u32 ret;
    find_ray_intersection(intersection, scene, ray_origin, ray_direction);
    object_material *mat = intersection->material;
    if(mat->is_black_body)
    {
        if(mat->is_emissive)
        {
            copy_spectrum(contribution, mat->emission_spd);
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
*/

void cast_ray(spectrum dst, scene_data *scene, vec3 ray_origin, vec3 ray_direction, u32 max_depth)
{
    spectrum contribution = alloc_spd();
    spectrum throughput   = alloc_spd();
    spectrum reflectance  = alloc_spd();
    spectrum tmp_spectrum = alloc_spd();
    zero_spectrum(contribution);
    zero_spectrum(reflectance);
    const_spectrum(throughput, 1.0);

    vec3 in;

    scene_point intersection;
    memset(&intersection, 0, sizeof(scene_point));
    for(u32 depth = 0; depth < max_depth; depth += 1)
    {
        find_ray_intersection(&intersection, scene, ray_origin, ray_direction);

        object_material *mat = intersection.surface_material;
        if(mat->is_black_body && !mat->is_emissive) break;
        else if(mat->is_black_body && mat->is_emissive)
        {
            spectral_mul_by_spectrum(tmp_spectrum, throughput, mat->emission_spd);
            spectral_sum(dst, dst, tmp_spectrum);
            break;
        }
        else
        {
            direct_light_contribution(contribution, &intersection, scene);
            spectral_mul_by_spectrum(tmp_spectrum, throughput, contribution);
            spectral_sum(dst, dst, tmp_spectrum);

            f64 dir_pdf;
            mat->sample_direction(&in, &dir_pdf, &intersection);

            bdsf(reflectance, &intersection, in);
            spectral_mul_by_scalar(reflectance, reflectance, dir_pdf);
            spectral_mul_by_spectrum(throughput, throughput, reflectance);

            ray_direction = in;
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

//Data to output to file(s):
//  - Pixel spectra: sum(filter * weight * incident radiance)
//  - Variance spectra
//  - Filter sums

void sample_scene(spectrum contribution, f64 *filter, u32 x, u32 y, scene_data *scene, camera_data *camera, u32 max_cast_depth)
{
    zero_spectrum(contribution);

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
    cast_ray(contribution, scene, ray_origin, ray_direction, max_cast_depth);

    f64 pixel_filter_value = 1.0;
    f64 pixel_weight_value = 1.0;
    f64 vignette_factor    = vec3_dot(ray_direction, camera->forward);
    spectral_mul_by_scalar(contribution, contribution, vignette_factor * pixel_filter_value);

    *filter = pixel_filter_value;
}

//Render image TODO:
//  - Open spd output file
//  - Open variance output file
//  - Allocate spd output buffer
//  - Allocate variance output buffer
//  - Assume output pixels size is a multiple of output buffer size for simplicity
//  - For each sample to take:
//      - For output pixels size / output buffer size
//          - Read dst pixels into output buffer
//          - Start pixel = ...
//          - For output buffer size
//              - Current pixel = ...
//              - (Sample scene)
//          - write file

void render_image(const char *output_spd_path, const char *output_avg_path, const char *output_var_path, u32 dst_width, u32 dst_height, scene_data *scene, camera_data *camera, u32 samples_per_pixel)
{
    file_handle output_spd_file = open_file(output_spd_path, ACCESS_READWRITE, FILE_NEW);
    file_handle output_var_file = open_file(output_var_path, ACCESS_READWRITE, FILE_NEW);
    file_handle output_avg_file = open_file(output_avg_path, ACCESS_READWRITE, FILE_NEW);

    spd_file_header header;
    memset(&header, 0, sizeof(header));
    header.id = 0xedfeefbe;
    header.has_filter_values = 1;
    header.width_in_pixels = dst_width;
    header.height_in_pixels = dst_height;
    header.number_of_wavelengths = number_of_spectrum_samples;
    header.min_wavelength = smallest_wavelength;
    header.wavelength_interval = sample_interval;
    write_file(output_spd_file, sizeof(header), &header);
    
    header.has_filter_values = 0;
    write_file(output_var_file, sizeof(header), &header);
    write_file(output_avg_file, sizeof(header), &header);

    u32 num_pixels = dst_width * dst_height;

    u32 spd_pixel_num_f64 = number_of_spectrum_samples + 1;
    u32 spd_pixel_size = spd_pixel_num_f64 * sizeof(f64);
    u32 spd_pixel_data_size = num_pixels * spd_pixel_size;
    u32 spd_avg_data_size   = num_pixels * spectrum_size;
    u32 spd_var_data_size   = num_pixels * spectrum_size;
    f64 *dst_pixels = alloc(spd_pixel_data_size);
    f64 *dst_avgs   = alloc(spd_avg_data_size);
    f64 *dst_vars   = alloc(spd_var_data_size);
    printf("Size = %u\n", spd_pixel_data_size);

    f64 *contribution_buffer = alloc(spd_pixel_size);
    spectrum contribution;
    contribution.samples = contribution_buffer;
    f64 *filter = &contribution_buffer[number_of_spectrum_samples];
    spectrum tmp_0_spd = alloc_spd();
    spectrum tmp_1_spd = alloc_spd();
    const u32 max_cast_depth = 4;
    for(u32 sample = 0; sample < samples_per_pixel; sample += 1)
    {
        printf("Sample %u / %u\n", sample+1, samples_per_pixel);
        //Filter final contribution and write to dst
        //pixel value = sum(filter * weight * radiance)/sum(filter)
        for(u32 y = 0; y < dst_height; y += 1)
        {
            for(u32 x = 0; x < dst_width; x += 1)
            {
                u32 pixel_offset = (y * dst_width) + x;
                f64 *dst_pixel = dst_pixels + spd_pixel_num_f64 * pixel_offset;
                f64 *dst_avg   = dst_avgs + number_of_spectrum_samples * pixel_offset;
                f64 *dst_var   = dst_vars + number_of_spectrum_samples * pixel_offset;
                spectrum dst_pixel_spd, dst_pixel_avg, dst_pixel_var;
                dst_pixel_spd.samples = dst_pixel;
                dst_pixel_avg.samples = dst_avg;
                dst_pixel_var.samples = dst_var;

                sample_scene(contribution, filter, x, y, scene, camera, max_cast_depth);

                //Accumulate contribution and filter values
                spectral_sum(dst_pixel_spd, dst_pixel_spd, contribution);
                dst_pixel[number_of_spectrum_samples] += *filter;

                //Compute new mean and variance
                copy_spectrum(tmp_0_spd, dst_pixel_avg);
                spectral_sub(tmp_0_spd, contribution, tmp_0_spd);
                copy_spectrum(tmp_1_spd, tmp_0_spd);
                spectral_div_by_scalar(tmp_0_spd, tmp_0_spd, (f64)(sample+1));
                spectral_sum(dst_pixel_avg, dst_pixel_avg, tmp_0_spd);
                spectral_sub(tmp_0_spd, contribution, dst_pixel_avg);
                spectral_mul_by_spectrum(tmp_0_spd, tmp_1_spd, tmp_0_spd);
                spectral_sum(dst_pixel_var, dst_pixel_var, tmp_0_spd);
            }
        }
    }

    u32 dst_pixel_size = spectrum_size + sizeof(f64);
    for(u32 pixel = 0; pixel < num_pixels; pixel += 1)
    {
        f64 *dst_pixel = dst_pixels + spd_pixel_num_f64 * pixel;
        f64 *dst_avg   = dst_avgs + number_of_spectrum_samples * pixel;
        f64 *dst_var   = dst_vars + number_of_spectrum_samples * pixel;
        write_file(output_spd_file, dst_pixel_size, dst_pixel);
        write_file(output_avg_file, spectrum_size, dst_avg);
        spectrum var_spd;
        var_spd.samples = dst_var;
        spectrum_normalise(var_spd);
        write_file(output_var_file, spectrum_size, dst_var);
    }
    close_file(output_avg_file);
    close_file(output_var_file);
    close_file(output_spd_file);
    unalloc(dst_vars, 0);
    unalloc(dst_avgs, 0);
    unalloc(dst_pixels, 0);
}
