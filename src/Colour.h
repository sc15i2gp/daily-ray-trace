#include "Maths.h"

struct RGB8
{
	uint8_t B;
	uint8_t G;
	uint8_t R;
};
typedef Vec3 RGB64;

RGB8 rgb64_to_rgb8(RGB64);

#define SPECTRUM_SAMPLE_MAX 1024

struct Spectrum
{
	double samples[SPECTRUM_SAMPLE_MAX];
	double start_wavelength; //In nm
	double end_wavelength; //In nm
	int number_of_samples;
};
Spectrum operator+(Spectrum, Spectrum);
Spectrum operator*(Spectrum, Spectrum);
Spectrum operator*(double, Spectrum);
Spectrum operator/(Spectrum, double);
void operator+=(Spectrum&, Spectrum);
void normalise(Spectrum& spd);
Spectrum generate_black_body_spd(double temperature, double start_wavelength = 300.0, double end_wavelength = 780.0, double sample_interval = 5.0);
