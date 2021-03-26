#pragma once

struct Debug_Info
{
	char* name = "Debug";
	int pixel_x;
	int pixel_y;
};

#ifdef DEBUG_BUILD
#define DEBUG(statements) statements
#else
#define DEBUG(statements)
#endif
extern Debug_Info debug_info;

void debug_set_pixel(int x, int y);
