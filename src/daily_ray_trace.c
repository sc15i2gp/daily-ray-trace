void init_camera(camera_data *camera, u32 width_px, u32 height_px, vec3 position, vec3 up, vec3 right, vec3 forward, f64 fov, f64 focal_depth, f64 focal_length, f64 aperture_radius)
{
    camera->forward = forward;
    camera->right   = right;
    camera->up      = up;

    f64 aperture_distance     = (focal_length * focal_depth) / (focal_length + focal_depth);
    camera->aperture_position = vec3_sum(position, vec3_mul_by_f64(forward, aperture_distance));
    camera->aperture_radius   = aperture_radius;

    f64  fov_rad             = fov * (PI/180.0);
    f64  aspect_ratio        = (f64)width_px / (f64)height_px;
    f64  film_width          = 2.0 * aperture_distance * tan(fov_rad/2.0);
    f64  film_height         = film_width / aspect_ratio;
    vec3 film_right          = vec3_mul_by_f64(right, 0.5 * film_width);
    vec3 film_top            = vec3_mul_by_f64(up,    0.5 * film_height);
    camera->film_bottom_left = vec3_sub(vec3_sub(position, film_right), film_top);

    camera->pixel_width  = film_width  / (f64)width_px;
    camera->pixel_height = film_height / (f64)height_px;
}

void init_scene(scene_data *scene)
{
    scene->num_surfaces = 2;
    scene->surfaces = VirtualAlloc(NULL, sizeof(object_geometry), MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);

    scene->surfaces[0].type   = GEO_TYPE_SPHERE;
    scene->surfaces[0].center.x = 0.0;
    scene->surfaces[0].center.y = 0.0;
    scene->surfaces[0].center.z = 0.0;
    scene->surfaces[0].radius = 0.3;
    scene->surfaces[1].type = GEO_TYPE_POINT;
    scene->surfaces[1].position.x = 0.0; 
    scene->surfaces[1].position.y = 1.0; 
    scene->surfaces[1].position.z = 1.0; 
    /*
    vec3 pp = {-0.5, -0.5, -0.5};
    vec3 pv = {-0.5, 0.5, -0.5};
    vec3 pu = {0.5, -0.5, -0.5};
    scene->surfaces[0].type = GEO_TYPE_PLANE;
    scene->surfaces[0].origin = pp;
    scene->surfaces[0].u = vec3_sub(pu, pp);
    scene->surfaces[0].v = vec3_sub(pv, pp);
    scene->surfaces[0].normal = vec3_normalise(vec3_cross(scene->surfaces[0].u, scene->surfaces[0].v));
    */

    scene->num_surface_materials = 2;
    scene->surface_materials = VirtualAlloc(NULL, scene->num_surface_materials * sizeof(object_material), MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    scene->surface_materials[0].diffuse_spd = alloc_spd();
    scene->surface_materials[0].glossy_spd  = alloc_spd();
    scene->surface_materials[0].shininess   = 1.0;
    scene->surface_materials[1].emission_spd = alloc_spd();
    //const_spectrum(scene->surface_materials[1].emission_spd, 1.0);
    generate_blackbody_spectrum(scene->surface_materials[1].emission_spd, 4000.0L);
    spectrum_normalise(scene->surface_materials[1].emission_spd);
    spectral_mul_by_scalar(scene->surface_materials[1].emission_spd, scene->surface_materials[1].emission_spd, 16);
    scene->surface_materials[1].is_emissive = 1;

    spectrum white       = alloc_spd();
    spectrum rgb_red     = alloc_spd();
    spectrum rgb_green   = alloc_spd();
    spectrum rgb_blue    = alloc_spd();
    spectrum rgb_cyan    = alloc_spd();
    spectrum rgb_magenta = alloc_spd();
    spectrum rgb_yellow  = alloc_spd();
    const_spectrum(white, 1.0);
    load_csv_file_to_spectrum(rgb_red, "spectra\\red_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_green, "spectra\\green_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_blue, "spectra\\blue_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_cyan, "spectra\\cyan_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_magenta, "spectra\\magenta_rgb_to_spd.csv");
    load_csv_file_to_spectrum(rgb_yellow, "spectra\\yellow_rgb_to_spd.csv");

    rgb_f64 diffuse_rgb = {0.25, 1.0, 0.25};
    rgb_f64_to_spectrum(diffuse_rgb, scene->surface_materials[0].diffuse_spd, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);
    rgb_f64 glossy_rgb  = {0.5, 1.0, 0.5};
    rgb_f64_to_spectrum(glossy_rgb, scene->surface_materials[0].glossy_spd, white, rgb_red, rgb_green, rgb_blue, rgb_cyan, rgb_magenta, rgb_yellow);

    free_spd(white);
    free_spd(rgb_red);
    free_spd(rgb_green);
    free_spd(rgb_blue);
    free_spd(rgb_cyan);
    free_spd(rgb_magenta);
    free_spd(rgb_yellow);
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
    f64 vis_dist = vec3_length(ray_dir_full);
    f64 dist     = INFINITY;
    for(u32 i = 0; i < scene->num_surfaces; ++i)
    {
        object_geometry *surface = &scene->surfaces[i];
        switch(surface->type)
        {
            case GEO_TYPE_POINT: continue;
            case GEO_TYPE_SPHERE:
            {
                dist = line_sphere_intersection(ray_origin, ray_direction, surface->center, surface->radius);
                break;
            };
            case GEO_TYPE_PLANE:
            {
                dist = line_plane_intersection(ray_origin, ray_direction, surface->origin, surface->normal, surface->u, surface->v);
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

u32 estimate_indirect_contribution(spectrum contribution, scene_point *intersection, scene_data *scene, vec3 ray_origin, vec3 ray_direction)
{
    u32 ret;
    spectrum reflectance = alloc_spd();
    intersection->is_black_body = 0;

    //Find scene intersection
    f64 min_dist = INFINITY;
    object_geometry *intersection_surface = NULL;
    u32 intersection_index = -1;
    for(u32 i = 0; i < scene->num_surfaces; ++i)
    {
        f64 dist = INFINITY;
        object_geometry *surface = &scene->surfaces[i];
        switch(surface->type)
        {
            case GEO_TYPE_POINT: continue;
            case GEO_TYPE_SPHERE:
            {
                dist = line_sphere_intersection(ray_origin, ray_direction, surface->center, surface->radius);
                break;
            }
            case GEO_TYPE_PLANE:
            {
                dist = line_plane_intersection(ray_origin, ray_direction, surface->origin, surface->normal, surface->u, surface->v);
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
                intersection->normal = vec3_normalise(vec3_sub(intersection->position, intersection_surface->center));
                break;
            }
            case GEO_TYPE_PLANE:
            {
                intersection->normal = intersection_surface->normal;
                break;
            }
        }
        intersection->material = &scene->surface_materials[intersection_index];
    }
    else
    {
        intersection->is_black_body = 1;
    }

    //Compute contribution
    if(intersection->is_black_body)
    {
        const_spectrum(contribution, 0.0);
        ret = 0;
    }
    else
    {
        for(u32 i = 0; i < scene->num_surfaces; ++i)
        {
            object_material *light_material = &scene->surface_materials[i];
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
                }
                if(points_mutually_visible(intersection->position, light_position, scene))
                {
                    vec3 outgoing = vec3_reverse(ray_direction);
                    vec3 incoming = vec3_normalise(vec3_sub(light_position, intersection->position));

                    bdsf(reflectance, intersection, incoming, outgoing);
                    spectral_sum(contribution, contribution, reflectance);
                    spectral_mul_by_spectrum(contribution, contribution, light_material->emission_spd);

                    f64  c = vec3_dot(incoming, intersection->normal) * attenuation_factor * (1.0/light_pdf);
                    spectral_mul_by_scalar(contribution, contribution, c);
                }
            }
        }
        ret = 1;
    }
    free_spd(reflectance);
    return ret;
}

//In  direction = towards light/away from camera
//Out direction = away from light/towards camera
void cast_ray(spectrum dst, scene_data* scene, vec3 ray_origin, vec3 ray_direction, u32 max_depth)
{
    f64  throughput_coefficient = 1.0;
    vec3 out_direction;
    vec3 in_direction = ray_direction;
    scene_point intersection;
    spectrum contribution = alloc_spd();
    spectrum throughput   = alloc_spd();
    spectrum reflectance  = alloc_spd();
    const_spectrum(throughput, 1.0);
    u32 ray_continues = estimate_indirect_contribution(contribution, &intersection, scene, ray_origin, in_direction);
    spectral_sum(dst, contribution, dst);
    for(u32 depth = 0; ray_continues && depth < max_depth; ++depth)
    {
        //Sample next direction
        out_direction = vec3_reverse(in_direction);
        do
        {
            in_direction = uniform_sample_sphere();
        }
        while(vec3_dot(in_direction, intersection.normal) <= 0.0);
        f64 dir_pdf = 1.0/(2.0*PI);
        
        //Update throughput
        throughput_coefficient = abs(vec3_dot(in_direction, intersection.normal)) * 1.0/dir_pdf; 
        bdsf(reflectance, &intersection, in_direction, out_direction);
        spectral_mul_by_spectrum(throughput, reflectance, throughput);
        spectral_mul_by_scalar(throughput, throughput, throughput_coefficient);

        //Estimate contribution
        ray_continues = estimate_indirect_contribution(contribution, &intersection, scene, intersection.position, in_direction);

        //Multiply contribution by throughput and acc
        spectral_mul_by_spectrum(contribution, contribution, throughput);
        spectral_sum(dst, contribution, dst);
    }
    free_spd(reflectance);
    free_spd(throughput);
    free_spd(contribution);
}

void print_camera(camera_data *camera)
{
    printf("CAMERA:\n");
    printf("Forward: "); print_vector(camera->forward); printf("\n");
    printf("Right:   "); print_vector(camera->right);   printf("\n");
    printf("Up:      "); print_vector(camera->up);      printf("\n");
    printf("Aperture position: "); print_vector(camera->aperture_position); printf("\n");
    printf("Aperture radius: %f\n", camera->aperture_radius);
    printf("Film bottom left: "); print_vector(camera->film_bottom_left); printf("\n");
    printf("Pixel width: %f Pixel height: %f\n", camera->pixel_width, camera->pixel_height);
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

    for(u32 sample = 0; sample < samples_per_pixel; ++sample)
    {
        //Filter final contribution and write to dst
        //pixel value = sum(filter * weight * radiance)/sum(filter)
        for(u32 y = 0; y < dst_height; ++y)
        {
            for(u32 x = 0; x < dst_width; ++x)
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
    for(u32 pixel = 0; pixel < num_pixels; ++pixel)
    {
        f64 pixel_filter = 1.0 / pixel_filter_sums[pixel];
        spectrum dst_pixel;
        dst_pixel.samples = dst_pixels + pixel * number_of_spectrum_samples;
        spectral_mul_by_scalar(dst_pixel, dst_pixel, pixel_filter);
    }
}
