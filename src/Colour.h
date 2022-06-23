
extern int number_of_samples;
struct RGB8
{
	uint8_t B;
	uint8_t G;
	uint8_t R;
	uint8_t A;
};
typedef Vec3 RGB64;

RGB8 rgb64_to_rgb8(RGB64);

#define SPECTRUM_SAMPLE_MAX 512
#define SPECTRUM_FILE_SAMPLE_MAX 2048

struct Spectrum
{
	double samples[SPECTRUM_SAMPLE_MAX];
};


void set_spectrum_to_value(Spectrum& spd, double d);
void set_spectrum_to_zero(Spectrum& spd);
void copy_spectrum(Spectrum& src_spd, Spectrum& dst_spd);
void spectral_sum(Spectrum& spd_0, Spectrum& spd_1, Spectrum& dst_spd);
void spectral_multiply(Spectrum& spd_0, Spectrum& spd_1, Spectrum& dst_spd);
void spectral_multiply(Spectrum& spd_0, double d, Spectrum& dst_spd);
void spectral_multiply(Spectrum& spd_0, Spectrum& spd_1, double d, Spectrum& dst_spd);
void spectral_sum_and_multiply(Spectrum& spd_0, Spectrum& spd_1, Spectrum& spd_2, Spectrum& dst_spd);
void spectral_sum_and_multiply(Spectrum& spd_0, Spectrum& spd_1, double d, Spectrum& dst_spd);
void normalise(Spectrum& spd);

double spd_value_at_wavelength(Spectrum&, double wavelength);

bool zero_spectrum(Spectrum&);

Spectrum generate_black_body_spd(double temperature);
Spectrum generate_constant_spd(double constant);

RGB64 spectrum_to_RGB64(Spectrum);
Spectrum RGB64_to_spectrum(RGB64);

Spectrum load_spd(const char* spd_path);
void write_spd(const char* spd_path, Spectrum);
void load_colour_data();

void set_reference_white(Spectrum);

void NEW_load_spd(const char* spd_path, Spectrum* dst);
void NEW_RGB64_to_spectrum(RGB64, Spectrum* result, Spectrum* white_rgb_to_spd, Spectrum* red_rgb_to_spd, Spectrum* blue_rgb_to_spd, Spectrum* green_rgb_to_spd, Spectrum* cyan_rgb_to_spd, Spectrum* magenta_rgb_to_spd, Spectrum* yellow_rgb_to_spd);
RGB64 NEW_spectrum_to_RGB64(Spectrum* spectrum, Spectrum* reference_white, Spectrum* colour_matching_spectra, Spectrum* xyz_spectra);

void NEW_set_spectrum_to_value(Spectrum*, double d);
void NEW_copy_spectrum(Spectrum* src, Spectrum* dst);
void NEW_spectral_sum_and_multiply(Spectrum* spd_0, Spectrum* spd_1, double d, Spectrum* dst); //dst = spd0 + d*spd1
void NEW_spectral_multiply(Spectrum* spd_0, Spectrum* spd_1, Spectrum* dst); //dst = spd_0 * spd_1
