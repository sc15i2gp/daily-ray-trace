#include "Texture.h"

uint8_t* get_pixel(Texture texture, int x, int y)
{
	//Clamp for now
	x = clamp(x, 0, texture.width-1);
	y = clamp(y, 0, texture.height -1);
	uint8_t* pixel = texture.pixels + ((y * texture.width) + x)*texture.pixel_size;
	return pixel;
}

void set_pixel(Texture texture, int x, int y, uint8_t* value)
{
	uint8_t* pixel = get_pixel(texture, x, y);
	memcpy(pixel, value, texture.pixel_size);
}

void clear_texture(Texture texture, uint8_t* clear_value)
{
	for(int y = 0; y < texture.height; ++y)
	{
		for(int x = 0; x < texture.width; ++x)
		{
			set_pixel(texture, x, y, clear_value);
		}
	}
}

uint8_t* sample_texture(Texture texture, Vec2 texture_coordinates)
{
	Vec2 texture_dimensions = {(double)texture.width, (double)texture.height};
	Vec2 actual_texture_coordinates = texture_dimensions * texture_coordinates;
	
	//Nearest neighbour for now
	int pixel_x = round(actual_texture_coordinates.x); 
	int pixel_y = round(actual_texture_coordinates.y);

	uint8_t* pixel = get_pixel(texture, pixel_x, pixel_y);

	return pixel;
}

Texture create_texture(int width, int height, int pixel_size)
{
	Texture t = {};
	t.in_use = true;
	t.width = width;
	t.height = height;
	t.pixel_size = pixel_size;
	t.pixels = (uint8_t*)alloc(width * height * pixel_size);
	return t;
}
