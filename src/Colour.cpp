#include "Colour.h"

double start_wavelength = 380.0; //In nm
double end_wavelength = 720.0; //In nm
double sample_interval = 5.0; //In nm
int number_of_samples = 69;

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
	for(int i = 0; i < number_of_samples; ++i)
	{
		spd.samples[i] = spd_0.samples[i] + spd_1.samples[i];
	}
	return spd;
}

//NOTE: Currently assumes both spectra have same wavelength ranges and intervals
Spectrum operator*(Spectrum spd_0, Spectrum spd_1)
{
	Spectrum spd = {};
	for(int i = 0; i < number_of_samples; ++i)
	{
		spd.samples[i] = spd_0.samples[i] * spd_1.samples[i];
	}
	return spd;
}

Spectrum operator*(double d, Spectrum spd_0)
{
	Spectrum spd = {};
	for(int i = 0; i < number_of_samples; ++i)
	{
		spd.samples[i] = d * spd_0.samples[i];
	}
	return spd;
}

Spectrum operator/(Spectrum spd_0, double d)
{
	Spectrum spd = {};
	for(int i = 0; i < number_of_samples; ++i)
	{
		spd.samples[i] = spd_0.samples[i] / d;
	}
	return spd;
}

void operator+=(Spectrum& spd_0, Spectrum spd_1)
{
	spd_0 = spd_0 + spd_1;
}

void operator*=(Spectrum& spd_0, Spectrum spd_1)
{
	spd_0 = spd_0 * spd_1;
}

void operator/=(Spectrum& spd, double d)
{
	spd = spd/d;
}

void normalise(Spectrum& spd)
{
	double highest_value = 0.0;
	for(int i = 0; i < number_of_samples; ++i) if(spd.samples[i] > highest_value) highest_value = spd.samples[i];

	for(int i = 0; i < number_of_samples; ++i) spd.samples[i] /= highest_value;
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

//WAVELENGTH DEPENDENT
Spectrum generate_black_body_spd(double temperature)
{
	Spectrum spd = {};
	for(int i = 0; i < number_of_samples; ++i)
	{
		double wavelength = nm_to_m(start_wavelength + i * sample_interval);
		spd.samples[i] = (double)nm_to_m(compute_black_body_power(temperature, wavelength));
	}
	return spd;
}

Spectrum generate_constant_spd(double constant)
{
	Spectrum spd = {};
	for(int i = 0; i < number_of_samples; ++i) spd.samples[i] = 1.0L;
	return spd;
}

char* find_next_line(char* c)
{
	for(; *c != '\n' && *c != 0; ++c);
	++c;
	return c;
}

char* find_next_character(char* c, char character)
{
	for(; *c != character && *c != 0; ++c);
	return c;
}

bool is_number_character(char c)
{
	return (c >= '0' && c <= '9') || c == '-';
}

char* find_next_number(char* c)
{
	for(; !is_number_character(*c) && *c != 0; ++c);
	return c;
}

Spectrum reference_white = {};
Spectrum colour_matching_functions[3] = {};

Spectrum white_rgb_to_spd = {};
Spectrum cyan_rgb_to_spd = {};
Spectrum magenta_rgb_to_spd = {};
Spectrum yellow_rgb_to_spd = {};
Spectrum red_rgb_to_spd = {};
Spectrum green_rgb_to_spd = {};
Spectrum blue_rgb_to_spd = {};

Spectrum load_spd(const char* spd_path)
{
	Spectrum spd = {};
	char* spd_contents = read_file_contents(spd_path);
	
	char* c = spd_contents;
	for(int i = 0; i < number_of_samples; ++i)
	{
		c = find_next_number(find_next_character(c, ','));
		spd.samples[i] = (long double) atof(c);

		c = find_next_line(c);
	}

	dealloc(spd_contents);

	return spd;
}

void load_colour_matching_functions()
{
	colour_matching_functions[0] = load_spd("cmf_x.csv");
	colour_matching_functions[1] = load_spd("cmf_y.csv");
	colour_matching_functions[2] = load_spd("cmf_z.csv");
}

void load_d65_illuminant()
{
	reference_white = load_spd("d65.csv");
	//normalise(reference_white);
	for(int i = 0; i < number_of_samples; ++i) reference_white.samples[i] /= 100.0L;
}

void load_e_illuminant()
{
	reference_white = generate_constant_spd(1.0);
	for(int i = 0; i < number_of_samples; ++i) reference_white.samples[i] = 1.0L;
}

void load_rgb_to_spd_functions()
{
	white_rgb_to_spd = load_spd("white_rgb_to_spd.csv");
	cyan_rgb_to_spd = load_spd("cyan_rgb_to_spd.csv");
	magenta_rgb_to_spd = load_spd("magenta_rgb_to_spd.csv");
	yellow_rgb_to_spd = load_spd("yellow_rgb_to_spd.csv");
	red_rgb_to_spd = load_spd("red_rgb_to_spd.csv");
	green_rgb_to_spd = load_spd("green_rgb_to_spd.csv");
	blue_rgb_to_spd = load_spd("blue_rgb_to_spd.csv");
}

void load_colour_data()
{
	load_rgb_to_spd_functions();
	load_e_illuminant();
	//load_d65_illuminant();
	load_colour_matching_functions();
}

//WAVELENGTH DEPENDENT
Vec3 spectrum_to_xyz(Spectrum spd)
{
	//SPD -> XYZ
	Spectrum spd_xyz[3] = {};
	spd_xyz[0] = spd * colour_matching_functions[0];
	spd_xyz[1] = spd * colour_matching_functions[1];
	spd_xyz[2] = spd * colour_matching_functions[2];

	Vec3 xyz = {};
	double l = (end_wavelength - start_wavelength)/(double)(number_of_samples);
	for(int i = 0; i < 3; ++i)
	{
		for(int j = 0; j < number_of_samples; ++j)
		{
			xyz[i] += spd_xyz[i].samples[j];
		}
	}
	xyz *= l;
	
	return xyz;
}

RGB64 gamma_linear_rgb(RGB64 linear_rgb)
{
	RGB64 rgb = {};
	for(int i = 0; i < 3; ++i)
	{
		rgb[i] = (linear_rgb[i] <= 0.0031308) ? 12.92 * linear_rgb[i] : (1.055 * pow(linear_rgb[i], 1.0/2.4)) - 0.055;
	}
	return rgb;
}

RGB64 inverse_gamma_rgb(RGB64 rgb)
{
	RGB64 linear_rgb = {};
	for(int i = 0; i < 3; ++i)
	{
		linear_rgb[i] = (rgb[i] <= 0.04045) ? rgb[i]/12.92 : pow(((rgb[i] + 0.055)/1.055), 2.4);
	}
	return linear_rgb;
}

RGB64 spectrum_to_RGB64(Spectrum spd)
{
	Vec3 xyz = spectrum_to_xyz(spd);
	Vec3 white_xyz = spectrum_to_xyz(reference_white);

	xyz /= white_xyz.y;

	//XYZ -> RGB
	RGB64 linear_rgb = {};
	linear_rgb.R = 3.24097 * xyz.x - 1.53738 * xyz.y - 0.49831 * xyz.z;
	linear_rgb.G = -0.96934 * xyz.x + 1.87596 * xyz.y + 0.04156 * xyz.z;
	linear_rgb.B = 0.05563 * xyz.x - 0.20398 * xyz.y + 1.05697 * xyz.z;

	RGB64 rgb = gamma_linear_rgb(linear_rgb);
	return rgb;
}

Spectrum RGB64_to_spectrum(RGB64 rgb)
{
	Spectrum spd = {};
	if(rgb.R <= rgb.G && rgb.R <= rgb.B)
	{
		spd += rgb.R * white_rgb_to_spd;
		if(rgb.G <= rgb.B)
		{
			spd += (rgb.G - rgb.R)*cyan_rgb_to_spd;
			spd += (rgb.B - rgb.G)*blue_rgb_to_spd;
		}
		else
		{
			spd += (rgb.B - rgb.R)*cyan_rgb_to_spd;
			spd += (rgb.G - rgb.B)*green_rgb_to_spd;
		}
	}
	else if(rgb.G <= rgb.R && rgb.G <= rgb.B)
	{
		spd += rgb.G * white_rgb_to_spd;
		if(rgb.R <= rgb.B)
		{
			spd += (rgb.R - rgb.G)*magenta_rgb_to_spd;
			spd += (rgb.B - rgb.R)*blue_rgb_to_spd;
		}
		else
		{
			spd += (rgb.B - rgb.G)*magenta_rgb_to_spd;
			spd += (rgb.R - rgb.B)*red_rgb_to_spd;
		}
	}
	else
	{
		spd += rgb.B * white_rgb_to_spd;
		if(rgb.R <= rgb.G)
		{
			spd += (rgb.R - rgb.B)*yellow_rgb_to_spd;
			spd += (rgb.G - rgb.R)*green_rgb_to_spd;
		}
		else
		{
			spd += (rgb.G - rgb.B)*yellow_rgb_to_spd;
			spd += (rgb.R - rgb.G)*red_rgb_to_spd;
		}
	}

	return spd;
}
