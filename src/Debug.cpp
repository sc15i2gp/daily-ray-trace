
Debug_Info debug_info;

void debug_set_pixel(int x, int y)
{
	debug_info.pixel_x = x;
	debug_info.pixel_y = y;
}

void debug_set_sample_number(int n)
{
	debug_info.sample_number = n;
}

void debug_set_current_eye_radiance(Spectrum radiance)
{
	debug_info.current_eye_ray_radiance = radiance;
}

void debug_set_last_intersection_computed(Geometry_Intersection_Point p)
{
	debug_info.last_intersection_computed = p;
}

void debug_set_last_sampled_texture(Texture t)
{
	debug_info.last_sampled_texture = t;
}
