#pragma once
#include <string.h>
#include "Maths.h"
#include "Platform.h"

struct Texture
{
	bool in_use;
	int width;
	int height;
	int pixel_size;
	uint8_t* pixels;
};

uint8_t* get_pixel(Texture, int x, int y);
void set_pixel(Texture, int x, int y, uint8_t* value);
void clear_texture(Texture, uint8_t* clear_value);
uint8_t* sample_texture(Texture, Vec2 coordinates);
Texture create_texture(int width, int height, int pixel_size);

#define TEXTURE_SAMPLE(type, texture, coordinates) (type*)sample_texture(texture, coordinates)
#define TEXTURE_READ(type, texture, x, y) (type*)get_pixel(texture, x, y)
#define TEXTURE_WRITE(texture, x, y, value) set_pixel(texture, x, y, (uint8_t*)(&value))
#define TEXTURE_CLEAR(texture, value) clear_texture(texture, (uint8_t*)(&value))
#define TEXTURE_CREATE(type, width, height) create_texture(width, height, sizeof(type))
