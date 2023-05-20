
struct Debug_Info
{
	Scene* scene;
	Scene_Point* scene_path;
	Spectrum* render_target;
	int scene_path_depth;
	int pixel_x;
	int pixel_y;
	int pass;
	char* scene_func;
};

Debug_Info debug_info;

#ifdef DEBUG_BUILD
#define DEBUG(statement) statement
#else
#define DEBUG(statement)
#endif
