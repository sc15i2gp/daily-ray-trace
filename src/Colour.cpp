#include "Colour.h"

RGB8 rgb64_to_rgb8(RGB64 rgb64)
{
	RGB8 rgb8 = {};

	if(rgb64.R < 0.0) rgb64.R = 0.0;
	if(rgb64.R > 1.0) rgb64.R = 1.0;
	if(rgb64.G < 0.0) rgb64.G = 0.0;
	if(rgb64.G > 1.0) rgb64.G = 1.0;
	if(rgb64.B < 0.0) rgb64.B = 0.0;
	if(rgb64.B > 1.0) rgb64.B = 1.0;

	rgb8.R = (uint8_t)(rgb64.R * 255.0);
	rgb8.G = (uint8_t)(rgb64.G * 255.0);
	rgb8.B = (uint8_t)(rgb64.B * 255.0);

	return rgb8;
}

Spectrum operator+(Spectrum spd_0, Spectrum spd_1)
{
	Spectrum spd = {};
	spd.start_wavelength = spd_0.start_wavelength;
	spd.end_wavelength = spd_0.end_wavelength;
	spd.number_of_samples = spd_0.number_of_samples;

	for(int i = 0; i < spd.number_of_samples; ++i)
	{
		spd.samples[i] = spd_0.samples[i] + spd_1.samples[i];
	}
	return spd;
}

//NOTE: Currently assumes both spectra have same wavelength ranges and intervals
Spectrum operator*(Spectrum spd_0, Spectrum spd_1)
{
	Spectrum spd = {};
	spd.start_wavelength = spd_0.start_wavelength;
	spd.end_wavelength = spd_0.end_wavelength;
	spd.number_of_samples = spd_0.number_of_samples;

	for(int i = 0; i < spd.number_of_samples; ++i)
	{
		spd.samples[i] = spd_0.samples[i] * spd_1.samples[i];
	}

	return spd;
}

Spectrum operator*(double d, Spectrum spd_0)
{
	Spectrum spd = {};
	spd.start_wavelength = spd_0.start_wavelength;
	spd.end_wavelength = spd_0.end_wavelength;
	spd.number_of_samples = spd_0.number_of_samples;

	for(int i = 0; i < spd.number_of_samples; ++i)
	{
		spd.samples[i] = d * spd_0.samples[i];
	}

	return spd;
}

Spectrum operator/(Spectrum spd_0, double d)
{
	Spectrum spd = {};
	spd.start_wavelength = spd_0.start_wavelength;
	spd.end_wavelength = spd_0.end_wavelength;
	spd.number_of_samples = spd_0.number_of_samples;

	for(int i = 0; i < spd.number_of_samples; ++i)
	{
		spd.samples[i] = spd_0.samples[i] / d;
	}

	return spd;
}

void operator+=(Spectrum& spd_0, Spectrum spd_1)
{
	spd_0 = spd_0 + spd_1;
}

void normalise(Spectrum& spd)
{
	double highest_value = 0.0;
	for(int i = 0; i < spd.number_of_samples; ++i) if(spd.samples[i] > highest_value) highest_value = spd.samples[i];

	for(int i = 0; i < spd.number_of_samples; ++i) spd.samples[i] /= highest_value;
}

long double c = 2.99792458e8L; //Speed of light
long double h = 6.626176e-34L; //Planck constant
long double k = 1.380662e-23L; //Boltzmann constant
//Uses Planck's formula to compute the power of a black body radiator's emission at a given temperature and wavelength
//Temperature in kelvin and wavelength in meters
long double compute_black_body_power(long double temperature, long double wavelength)
{
	long double numerator = 2.0L * PI * h * c * c;

	long double lambda_5 = powl(wavelength, 5.0L);
	long double e_power_numerator = (h * c) / k;
	long double e_power_denominator = temperature * wavelength;
	long double e_power = e_power_numerator / e_power_denominator;
	long double e_term = expl(e_power);

	long double denominator = lambda_5 * (e_term - 1.0L);

	return numerator/denominator;
}

inline
long double m_to_nm(long double m)
{
	long double nm = m * 1e9;
	return nm;
}

inline
long double nm_to_m(long double nm)
{
	long double m = nm * 1e-9;
	return m;
}

long double sample_interval(Spectrum spd)
{
	return (spd.end_wavelength - spd.start_wavelength)/(long double)(spd.number_of_samples - 1);
}

Spectrum generate_black_body_spd(double temperature, double start_wavelength, double end_wavelength, double sample_interval)
{
	Spectrum spd = {};
	spd.number_of_samples = (int)((end_wavelength - start_wavelength)/sample_interval) + 1;
	for(int i = 0; i < spd.number_of_samples; ++i)
	{
		double wavelength = nm_to_m(start_wavelength + i * sample_interval);
		spd.samples[i] = (double)nm_to_m(compute_black_body_power(temperature, wavelength));
	}
	spd.start_wavelength = start_wavelength;
	spd.end_wavelength = end_wavelength;
	return spd;
}
