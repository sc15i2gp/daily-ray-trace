
extern const int number_of_samples;
typedef Vec3 RGB64;

#define SPECTRUM_SAMPLE_MAX 128
#define SPECTRUM_FILE_SAMPLE_MAX 2048

struct Spectrum
{
	double samples[SPECTRUM_SAMPLE_MAX];
};

void normalise_spectrum(Spectrum*);
double spectrum_value_at_wavelength(Spectrum*, double);
void generate_black_body_spd(double temperature);
bool zero_spectrum(Spectrum*);

void load_spd(const char* spd_path, Spectrum* dst);
void RGB64_to_spectrum(RGB64, Spectrum* result, Spectrum* white_rgb_to_spd, Spectrum* red_rgb_to_spd, Spectrum* blue_rgb_to_spd, Spectrum* green_rgb_to_spd, Spectrum* cyan_rgb_to_spd, Spectrum* magenta_rgb_to_spd, Spectrum* yellow_rgb_to_spd);
RGB64 spectrum_to_RGB64(Spectrum* spectrum, Spectrum* reference_white, Spectrum* colour_matching_spectra, Spectrum* xyz_spectra);

void set_spectrum_to_value(Spectrum*, double d);
void copy_spectrum(Spectrum* src, Spectrum* dst);
void spectral_sum(Spectrum* spd_0, Spectrum* spd_1, Spectrum* dst);
void spectral_sum_and_multiply(Spectrum* spd_0, Spectrum* spd_1, double d, Spectrum* dst); //dst = spd0 + d*spd1
void spectral_sum_and_multiply(Spectrum* spd_0, Spectrum* spd_1, Spectrum* spd_2, Spectrum* dst); //dst = spd_0 + spd_1 * spd_2
void spectral_multiply(Spectrum* spd_0, Spectrum* spd_1, double d, Spectrum* dst);
void spectral_multiply(Spectrum* spd_0, Spectrum* spd_1, Spectrum* dst); //dst = spd_0 * spd_1
void spectral_multiply(Spectrum* spd_0, double d, Spectrum* dst); //dst = d * spd_0
