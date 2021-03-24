#pragma once
#include "Maths.h"

struct Texture
{
	bool in_use;
	int width;
	int height;
	int pixel_size;
	uint8_t* pixels;
};

uint8_t* get_pixel(Texture, int x, int y);uint8_t* sample_texture(Texture, Vec2 coordinates);
