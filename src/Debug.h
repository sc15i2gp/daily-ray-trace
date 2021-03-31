struct Debug_Info
{
	char* name = "Debug";
	int pixel_x;
	int pixel_y;
	int sample_number;
	Spectrum current_eye_ray_radiance;
	Geometry_Intersection_Point last_intersection_computed;
	Texture last_sampled_texture;
};

#ifdef DEBUG_BUILD
#define DEBUG(statements) statements
#else
#define DEBUG(statements)
#endif

void debug_set_pixel(int x, int y);
void debug_set_sample_number(int n);
void debug_set_current_eye_radiance(Spectrum radiance);
void debug_set_last_intersection_computed(Geometry_Intersection_Point);
void debug_set_last_sampled_texture(Texture);
