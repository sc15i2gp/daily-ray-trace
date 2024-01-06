/*
vec3 cos_weighted_sample_hemisphere()
{
}

void blinn_phong_bsdf(spectrum *reflectance, scene_point *p, vec3 in, vec3 out)
{
}

void bsdf(spectrum *reflectance, scene_point *p, vec3 in, vec3 out)
{
    blinn_phong_bsdf(reflectance, p, in, out);
}

vec3 sample_point_on_surface(object_geometry *geometry, f64 *pdf)
{

}

u32 points_mutually_visible(vec3 p0, vec3 p1, object_geometry *surfaces)
{

}

void find_scene_intersection(scene_point *intersection, scene_data *scene, vec3 ray_origin, vec3 ray_direction)
{
    //Go through scene objects to see if ray intersects them
    //Get the nearest intersection point's local material data
}

u32 estimate_indirect_contribution(spectrum *contribution, scene_point *intersection, scene_data *scene, vec3 ray_origin, vec3 ray_direction)
{
    spectrum *reflectance;
    zero_spectrum(contribution);

    //Find scene intersection
    find_scene_intersection(intersection, scene, ray_origin, ray_direction);

    //If scene intersection is a black body
    if(intersection->is_black_body)
    {
        //If the black body is emissive
        if(intersection->is_emissive)
        {
            copy_spectrum(contribution, intersection->emission_spd);
        }
        return 0; //Escaped or absorbed
    }

    //For each emissive surface in the scene
    for(u32 i = 0; i < scene->num_emissive_surfaces; ++i)
    {
        f64 light_pdf;
        object_geometry *light_geometry = &scene->emissive_surfaces[i];
        object_material *light_material = &scene->emissive_materials[i];

        //Sample a point on its surface
        vec3 light_pos = sample_point_on_surface(light_geometry, &light_pdf);
        vec3 light_dir = vec3_normalise(vec3_sub(light_pos, intersection->position));

        //If intersection->position and light_pos are not mutually visible, continue
        if(!(points_mutually_visible(intersection->position, light_pos, scene->surfaces))) continue;

        //contribution += emission_spd * 1/light_pdf * dost(ray_direction, light_direction) * bsdf(reflectance, intersection, ray_direction, light_direction)
        bsdf(reflectance, intersection, ray_direction, light_dir);
        spectral_mul_by_scalar(reflectance, reflectance, 1.0/light_pdf);
        spectral_mul_by_spectrum(reflectance, reflectance, &light_material->emission_spd);
    }

    //return 1 for bounce
    return 1;
}
*/

f64 f64_max(f64 f0, f64 f1)
{
    return (f0 > f1) ? f0 : f1;
}

void cast_ray(spectrum *dst, scene_data* scene, vec3 ray_origin, vec3 ray_direction, u32 max_depth)
{
    //Render blinn-phong sphere with point light
    spectrum tmp_spd;
    zero_spectrum(&tmp_spd);
    f64 dist = line_sphere_intersection(ray_origin, ray_direction, scene->sphere_center, scene->sphere_radius);
    if(!isnan(dist))
    {
        vec3 outgoing      = vec3_reverse(ray_direction);
        vec3 intersection  = vec3_sum(ray_origin, vec3_mul_by_f64(ray_direction, dist));
        vec3 normal        = vec3_normalise(vec3_sub(intersection, scene->sphere_center));
        vec3 incoming      = vec3_normalise(vec3_sub(scene->light_position, intersection));
        f64  c = f64_max(0.0, vec3_dot(incoming, normal));
        //Diffuse
        spectral_mul_by_scalar(&tmp_spd, &scene->diffuse_spd, 1.0/PI);
        //Glossy
        vec3 bisector         = vec3_normalise(vec3_sum(outgoing, incoming));
        f64  spec_coefficient = pow(f64_max(0.0, vec3_dot(normal, bisector)), scene->shininess);
        spectral_mul_by_scalar(dst, &scene->glossy_spd, spec_coefficient);
        spectral_sum(&tmp_spd, &tmp_spd, dst);

        spectral_mul_by_spectrum(&tmp_spd, &tmp_spd, &scene->light_spd);
        spectral_mul_by_scalar(dst, &tmp_spd, c);
        /*
        vec3 light_dir = vec3_normalise(vec3_sub(scene->light_position, sphere_point));
        vec3 bisector = vec3_div_by_f64(vec3_sum(ray_direction, light_dir), 2.0);
        f64 spec_coefficient = pow(f64_max(0.0, vec3_dot(sphere_normal, bisector)), scene->shininess);
        spectral_mul_by_scalar(dst, &scene->glossy_spd, spec_coefficient);
        spectral_sum(dst, dst, &tmp_spd);
        */
    }
    else
    {
        const_spectrum(dst, 0.0);
    }
    /*
    spectrum *contribution;
    spectrum *throughput;
    spectrum *reflectance;
    scene_point *intersection;
    vec3 ray_outgoing;
    f64 dir_pdf;
    f64 throughput_coefficient = 1.0;
    //Estimate contribution
    u32 ray_continues = estimate_indirect_contribution(contribution, intersection, scene, ray_origin, ray_direction);
    spectral_sum(dst, contribution, dst);
    //for each point until max depth
    for(u32 depth = 0; depth < max_depth && ray_continues; ++depth)
    {
        //sample next direction
        ray_outgoing = cos_weighted_sample_hemisphere();
        throughput_coefficient = abs(vec3_dot(intersection->normal, ray_outgoing)) * (1.0/dir_pdf);

        //update throughput
        //f = f * dot(p normal, next_out) * (1/pdf) * bdsf(p, current_in, next_out)
        bsdf(reflectance, intersection, ray_direction, ray_outgoing);
        spectral_mul_by_spectrum(throughput, reflectance, throughput);
        spectral_mul_by_scalar(throughput, throughput, throughput_coefficient);

        //Estimate contribution
        ray_continues = estimate_indirect_contribution(contribution, intersection, scene, intersection->position, ray_outgoing);

        //Multiply contribution by throughput
        spectral_mul_by_spectrum(contribution, contribution, throughput);
        spectral_sum(dst, contribution, dst);
    }
    */
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

//TODO: Make work with non-pinhole camera
vec3 sample_camera(camera_data *camera)
{
    vec3 p = camera->aperture_position;
    return p;
}

void render_image(spectrum *dst_pixels, u32 dst_width, u32 dst_height, scene_data *scene, camera_data *camera, u32 samples_per_pixel)
{
    const u32 max_cast_depth = 4;
    u32 num_pixels = dst_width * dst_height;

    u32 pixel_filter_sums_size = num_pixels * sizeof(f64);
    u32 contribution_size = sizeof(spectrum);

    u8 *mem = (u8*)VirtualAlloc(NULL, pixel_filter_sums_size + contribution_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    f64 *pixel_filter_sums = (f64*)(mem);
    spectrum *contribution = (spectrum*)(mem + pixel_filter_sums_size);

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
                cast_ray(contribution, scene, ray_origin, ray_direction, max_cast_depth);

                f64 pixel_filter_value = 1.0;
                f64 pixel_weight_value = 1.0;

                spectrum *dst_pixel = dst_pixels + (y*dst_width) + x;
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
        spectrum *dst_pixel = &dst_pixels[pixel];
        spectral_mul_by_scalar(dst_pixel, dst_pixel, pixel_filter);
    }
}
