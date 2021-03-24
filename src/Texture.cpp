#include "Texture.h"

uint8_t* get_pixel(Texture texture, int x, int y)
{
	//Clamp for now
	x = clamp(x, 0, texture.width-1);
	y = clamp(y, 0, texture.height -1);
	uint8_t* pixel = texture.pixels + ((y * texture.width) + x)*texture.pixel_size;
	return pixel;
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
