
double start_wavelength = 380.0; //In nm
double end_wavelength = 720.0; //In nm
double sample_interval = 5.0; //In nm
int number_of_samples = 69;

//Spectral operations and example locations:
//	- spectrum = 0
//	- spectrum = double | cast_ray in Scene.cpp
//	- spectrum_0 = spectrum_1

//	- spectrum_0 = (double_0 * spectrum_0 + spectrum_1) | render_image in Scene.cpp
//	- spectrum_0 = spectrum_0 + spectrum_1 * spectrum_2 | cast_ray in Scene.cpp
//	- spectrum_0 = spectrum_0 + double * spectrum_1 * spectrum_2 | direct_light_contribution in Scene.cpp

//	- spectrum_0 = spectrum_0 * double * spectrum_1
//	- spectrum_0 = (double * double) * spectrum | torrance_sparrow_bsdf in bsdf.cpp
//	- spectrum_0 = double * spectrum | glossy_phong_bsdf in bsdf.cpp

//	- spectrum_0 = spectrum_1 / double | diffuse_phong_bsdf in bsdf.cpp
//	- spectrum_0 = spectrum_1 / double | cook_torrance_reflectance_bsdf in bsdf.cpp

//	- spectrum_0 = spectrum_1 + spectrum_2 | plastic_bsdf in bsdf.cpp
//	- spectrum_0 = spectrum_0 + spectrum_1 | bsdf in bsdf.cpp

//	- spectrum_0 = spectrum_0 + double * spectrum_1

void NEW_normalise(Spectrum* spd)
{
	double highest_value = 0.0;
	for(int i = 0; i < number_of_samples; ++i) if(spd->samples[i] > highest_value) highest_value = spd->samples[i];

	for(int i = 0; i < number_of_samples; ++i) spd->samples[i] /= highest_value;
}

double NEW_spectrum_value_at_wavelength(Spectrum* spd, double wavelength)
{
	int index = (int)((wavelength - start_wavelength)/sample_interval);
	return spd->samples[index];
}

bool NEW_zero_spectrum(Spectrum* spd)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		if(spd->samples[i] > 0.0) return false;
	}
	return true;
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
void NEW_generate_black_body_spd(Spectrum* spd, double temperature)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		double wavelength = nm_to_m(start_wavelength + i * sample_interval);
		spd->samples[i] = (double)nm_to_m(compute_black_body_power(temperature, wavelength));
	}
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

void NEW_write_spd(const char* spd_path, Spectrum* spd)
{
	char* spd_contents = (char*)alloc(MEGABYTES(4));

	int spd_size = sprintf(spd_contents, "Wavelength,Power\n");

	for(int i = 0; i < number_of_samples; ++i)
	{
		double ith_wavelength = start_wavelength + (double)i * sample_interval;
		spd_size += sprintf(spd_contents+spd_size, "%f,%f\n", ith_wavelength, spd->samples[i]);
	}
	write_file_contents(spd_path, spd_contents, spd_size);
	dealloc(spd_contents);
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

void NEW_set_spectrum_to_value(Spectrum* spd, double d)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		spd->samples[i] = d;
	}
}

void NEW_copy_spectrum(Spectrum* src, Spectrum* dst)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		dst->samples[i] = src->samples[i];
	}
}

void NEW_spectral_sum(Spectrum* spd_0, Spectrum* spd_1, Spectrum* dst)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		dst->samples[i] = spd_0->samples[i] + spd_1->samples[i];
	}
}

void NEW_spectral_multiply(Spectrum* spd_0, Spectrum* spd_1, Spectrum* dst)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		dst->samples[i] = spd_0->samples[i] * spd_1->samples[i];
	}
}

void NEW_spectral_multiply(Spectrum* spd, double d, Spectrum* dst)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		dst->samples[i] = d * spd->samples[i];
	}
}

void NEW_spectral_multiply(Spectrum* spd_0, Spectrum* spd_1, double d, Spectrum* dst)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		dst->samples[i] = spd_0->samples[i] * spd_1->samples[i] * d;
	}
}

void NEW_spectral_sum_and_multiply(Spectrum* spd_0, Spectrum* spd_1, double d, Spectrum* dst)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		dst->samples[i] = spd_0->samples[i] + (d * spd_1->samples[i]);
	}
}

void NEW_spectral_sum_and_multiply(Spectrum* spd_0, Spectrum* spd_1, Spectrum* spd_2, Spectrum* dst)
{
	for(int i = 0; i < number_of_samples; ++i)
	{
		dst->samples[i] = spd_0->samples[i] + (spd_1->samples[i] * spd_2->samples[i]);
	}
}

void NEW_load_spd(const char* spd_path, Spectrum* spd)
{
	//Store the spd file's contents in these 2 arrays
	double file_wavelengths[SPECTRUM_FILE_SAMPLE_MAX] = {};
	double file_samples[SPECTRUM_FILE_SAMPLE_MAX] = {};
	int file_number_of_samples = 0;
	
	char* spd_contents = read_file_contents(spd_path);
	
	//Find the number of samples given in the file
	char* c = find_next_line(spd_contents);
	for(char* d = c; d != NULL && *d != 0; d = find_next_line(d)) if(find_next_number(d)) ++file_number_of_samples;
	//Read the file's samples and wavelengths
	for(int i = 0; i < file_number_of_samples; ++i)
	{
		c = find_next_number(c);
		file_wavelengths[i] = (double) atof(c);
		c = find_next_number(find_next_character(c, ','));
		file_samples[i] = (double) atof(c);
		c = find_next_line(c);
	}
	
	//This is so that an spd file can be given in nm or um
	//If the first wavelength given in the file is less than an arbitrary number (10 in this case), then it is assumed that file is given in um
	//Linearly interpolate the desired wavelength values from the file's contents
	double um_threshold = 10.0;
	if(file_wavelengths[0] < um_threshold)
	{
		for(int i = 0; i < file_number_of_samples; ++i) file_wavelengths[i] *= 1000.0;
	}
	int t = 0;
	for(int i = 0; i < number_of_samples; ++i)
	{
		double ith_wavelength = start_wavelength + (double)i * sample_interval;
		for(; file_wavelengths[t+1] < ith_wavelength; ++t);
		spd->samples[i] = lerp(ith_wavelength, file_wavelengths[t], file_wavelengths[t+1], file_samples[t], file_samples[t+1]);
	}
}

void NEW_RGB64_to_spectrum(RGB64 rgb, Spectrum* spd, Spectrum* white_rgb_to_spd, Spectrum* red_rgb_to_spd, Spectrum* blue_rgb_to_spd, Spectrum* green_rgb_to_spd, Spectrum* cyan_rgb_to_spd, Spectrum* magenta_rgb_to_spd, Spectrum* yellow_rgb_to_spd)
{
	if(rgb.R <= rgb.G && rgb.R <= rgb.B)
	{
		NEW_spectral_sum_and_multiply(spd, white_rgb_to_spd, rgb.R, spd);
		if(rgb.G <= rgb.B)
		{
			NEW_spectral_sum_and_multiply(spd, cyan_rgb_to_spd, rgb.G - rgb.R, spd);
			NEW_spectral_sum_and_multiply(spd, blue_rgb_to_spd, rgb.B - rgb.G, spd);
		}
		else
		{
			NEW_spectral_sum_and_multiply(spd, cyan_rgb_to_spd, rgb.B - rgb.R, spd);
			NEW_spectral_sum_and_multiply(spd, green_rgb_to_spd, rgb.G - rgb.B, spd);
		}
	}
	else if(rgb.G <= rgb.R && rgb.G <= rgb.B)
	{
		NEW_spectral_sum_and_multiply(spd, white_rgb_to_spd, rgb.G, spd);
		if(rgb.R <= rgb.B)
		{
			NEW_spectral_sum_and_multiply(spd, magenta_rgb_to_spd, rgb.R - rgb.G, spd);
			NEW_spectral_sum_and_multiply(spd, blue_rgb_to_spd, rgb.B - rgb.R, spd);
		}
		else
		{
			NEW_spectral_sum_and_multiply(spd, magenta_rgb_to_spd, rgb.B - rgb.G, spd);
			NEW_spectral_sum_and_multiply(spd, red_rgb_to_spd, rgb.R - rgb.B, spd);
		}
	}
	else
	{
		NEW_spectral_sum_and_multiply(spd, white_rgb_to_spd, rgb.B, spd);
		if(rgb.R <= rgb.G)
		{
			NEW_spectral_sum_and_multiply(spd, yellow_rgb_to_spd, rgb.R - rgb.B, spd);
			NEW_spectral_sum_and_multiply(spd,green_rgb_to_spd, rgb.G - rgb.R, spd);
		}
		else
		{
			NEW_spectral_sum_and_multiply(spd, yellow_rgb_to_spd, rgb.G - rgb.B, spd); 
			NEW_spectral_sum_and_multiply(spd, red_rgb_to_spd, rgb.R - rgb.G, spd); 
		}
	}
}

Vec3 NEW_spectrum_to_xyz(Spectrum* spectrum, Spectrum* colour_matching_spectra, Spectrum* xyz_spectra)
{
	for(int i = 0; i < 3; ++i)
	{
		NEW_spectral_multiply(spectrum, &colour_matching_spectra[i], &xyz_spectra[i]);
	}

	Vec3 xyz = {};
	double l = (end_wavelength - start_wavelength)/(double)(number_of_samples);
	for(int i = 0; i < 3; ++i)
	{
		for(int j = 0; j < number_of_samples; ++j) xyz[i] += xyz_spectra[i].samples[j];
	}
	xyz *= l;
	return xyz;
}

RGB64 NEW_spectrum_to_RGB64(Spectrum* spectrum, Spectrum* reference_white, Spectrum* colour_matching_spectra, Spectrum* xyz_spectra)
{
	Vec3 xyz = NEW_spectrum_to_xyz(spectrum, colour_matching_spectra, xyz_spectra);
	Vec3 white_xyz = NEW_spectrum_to_xyz(reference_white, colour_matching_spectra, xyz_spectra);

	xyz /= white_xyz.y;

	//XYZ -> RGB
	RGB64 linear_rgb = {};
	linear_rgb.R = 3.24097 * xyz.x - 1.53738 * xyz.y - 0.49831 * xyz.z;
	linear_rgb.G = -0.96934 * xyz.x + 1.87596 * xyz.y + 0.04156 * xyz.z;
	linear_rgb.B = 0.05563 * xyz.x - 0.20398 * xyz.y + 1.05697 * xyz.z;

	RGB64 rgb = gamma_linear_rgb(linear_rgb);
	return rgb;
}

