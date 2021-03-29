struct Debug_Info
{
	char* name = "Debug";
	int pixel_x;
	int pixel_y;
	Spectrum current_eye_ray_radiance;
	Geometry_Intersection_Point last_intersection_computed;
};

#ifdef DEBUG_BUILD
#define DEBUG(statements) statements
#else
#define DEBUG(statements)
#endif

void debug_set_pixel(int x, int y);
