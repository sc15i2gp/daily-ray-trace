#pragma once
#include "Maths.h"
#include "Platform.h"

struct RGB8
{
	uint8_t B;
	uint8_t G;
	uint8_t R;
};
typedef Vec3 RGB64;

RGB8 rgb64_to_rgb8(RGB64);

#define SPECTRUM_SAMPLE_MAX 512

struct Spectrum
{
	double samples[SPECTRUM_SAMPLE_MAX];
};
Spectrum operator+(Spectrum, Spectrum);
Spectrum operator*(Spectrum, Spectrum);
Spectrum operator*(double, Spectrum);
Spectrum operator/(Spectrum, double);
void operator+=(Spectrum&, Spectrum);
void operator*=(Spectrum&, Spectrum);
void operator/=(Spectrum&, double);
void normalise(Spectrum& spd);

Spectrum generate_black_body_spd(double temperature);
Spectrum generate_constant_spd(double constant);

RGB64 spectrum_to_RGB64(Spectrum);
Spectrum RGB64_to_spectrum(RGB64);

Spectrum load_spd(const char* spd_path);
void load_colour_data();

void set_reference_white(Spectrum);
