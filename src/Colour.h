
extern int number_of_samples;
typedef Vec3 RGB64;

#define SPECTRUM_SAMPLE_MAX 128
#define SPECTRUM_FILE_SAMPLE_MAX 2048

struct Spectrum
{
	double samples[SPECTRUM_SAMPLE_MAX];
};

void NEW_normalise_spectrum(Spectrum*);
double NEW_spectrum_value_at_wavelength(Spectrum*, double);
void NEW_generate_black_body_spd(double temperature);
bool NEW_zero_spectrum(Spectrum*);

void NEW_load_spd(const char* spd_path, Spectrum* dst);
void NEW_RGB64_to_spectrum(RGB64, Spectrum* result, Spectrum* white_rgb_to_spd, Spectrum* red_rgb_to_spd, Spectrum* blue_rgb_to_spd, Spectrum* green_rgb_to_spd, Spectrum* cyan_rgb_to_spd, Spectrum* magenta_rgb_to_spd, Spectrum* yellow_rgb_to_spd);
RGB64 NEW_spectrum_to_RGB64(Spectrum* spectrum, Spectrum* reference_white, Spectrum* colour_matching_spectra, Spectrum* xyz_spectra);

void NEW_set_spectrum_to_value(Spectrum*, double d);
void NEW_copy_spectrum(Spectrum* src, Spectrum* dst);
void NEW_spectral_sum(Spectrum* spd_0, Spectrum* spd_1, Spectrum* dst);
void NEW_spectral_sum_and_multiply(Spectrum* spd_0, Spectrum* spd_1, double d, Spectrum* dst); //dst = spd0 + d*spd1
void NEW_spectral_sum_and_multiply(Spectrum* spd_0, Spectrum* spd_1, Spectrum* spd_2, Spectrum* dst); //dst = spd_0 + spd_1 * spd_2
void NEW_spectral_multiply(Spectrum* spd_0, Spectrum* spd_1, double d, Spectrum* dst);
void NEW_spectral_multiply(Spectrum* spd_0, Spectrum* spd_1, Spectrum* dst); //dst = spd_0 * spd_1
void NEW_spectral_multiply(Spectrum* spd_0, double d, Spectrum* dst); //dst = d * spd_0
