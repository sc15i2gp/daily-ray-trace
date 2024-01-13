#define MAX_NUM_SPECTRUM_VALUES 128
//TODO: Is this info nececssary outside of i/o?
u32 number_of_spectrum_samples;
f64 smallest_wavelength; //In nm
f64 largest_wavelength; //In nm
f64 sample_interval;   //In nm
u32 spectrum_size;

//Spectrum as point samples over a range of wavelengths
//All spectra have the same number of samples
//All samples "sample[i]" are at the same wavelength (smallest_wavelength + i*sample_interval)
typedef struct
{
    f64 samples[MAX_NUM_SPECTRUM_VALUES];
}   spectrum;

typedef struct
{
    f64 *samples;
}   new_spectrum;

typedef struct
{
    u32      capacity;
    u32      allocated;
    u32      *is_allocated;
    f64      *spectrum_buffer;
} spectrum_allocation_table;

spectrum_allocation_table spd_table;

new_spectrum alloc_spd();
void     free_spd(new_spectrum);

void zero_spectrum(spectrum *dst);
void const_spectrum(spectrum *dst, f64 f);
void copy_spectrum(spectrum *dst, spectrum *src);
void spectrum_normalise(spectrum *dst);
void spectral_mul_by_spectrum(spectrum *dst, spectrum *src0, spectrum *src1); //dst = src0 * src1
void spectral_mul_by_scalar(spectrum *dst, spectrum *src0, f64 src1); //dst = src0 * src1
void spectral_acc(spectrum *dst, spectrum *src); //dst += src
void spectral_acc_mul_by_scalar(spectrum *dst, spectrum *src0, f64 src1); //dst += src0 * src1
void spectral_sum(spectrum *dst, spectrum *src0, spectrum *src1); //dst = src0 + src1

rgb_f64 spectrum_to_xyz(spectrum *dst, spectrum *cmf_x, spectrum *cmf_y, spectrum *cmf_z, spectrum *ref_white);
rgb_f64 spectrum_to_rgb_f64(spectrum *spd, spectrum *cmf_x, spectrum *cmf_y, spectrum *cmf_z, spectrum *ref_white);
void rgb_f64_to_spectrum(rgb_f64 rgb, spectrum *dst, spectrum *white, spectrum *red, spectrum *green, spectrum *blue, spectrum *cyan, spectrum *magenta, spectrum *yellow);

void generate_blackbody_spectrum(spectrum *dst, f64 temperature);
u32 load_csv_file_to_spectrum(spectrum *dst, const char *path);

void new_zero_spectrum(new_spectrum dst);
void new_const_spectrum(new_spectrum dst, f64 f);
void new_copy_spectrum(new_spectrum dst, new_spectrum src);
void new_spectrum_normalise(new_spectrum dst);
void new_spectral_mul_by_spectrum(new_spectrum dst, new_spectrum src0, new_spectrum src1); //dst = src0 * src1
void new_spectral_mul_by_scalar(new_spectrum dst, new_spectrum src0, f64 src1); //dst = src0 * src1
void new_spectral_acc(new_spectrum dst, new_spectrum src); //dst += src
void new_spectral_acc_mul_by_scalar(new_spectrum dst, new_spectrum src0, f64 src1); //dst += src0 * src1
void new_spectral_sum(new_spectrum dst, new_spectrum src0, new_spectrum src1); //dst = src0 + src1

rgb_f64 new_spectrum_to_xyz(new_spectrum dst, new_spectrum cmf_x, new_spectrum cmf_y, new_spectrum cmf_z, new_spectrum ref_white);
rgb_f64 new_spectrum_to_rgb_f64(new_spectrum spd, new_spectrum cmf_x, new_spectrum cmf_y, new_spectrum cmf_z, new_spectrum ref_white);
void new_rgb_f64_to_spectrum(rgb_f64 rgb, new_spectrum dst, new_spectrum white, new_spectrum red, new_spectrum green, new_spectrum blue, new_spectrum cyan, new_spectrum magenta, new_spectrum yellow);

void new_generate_blackbody_spectrum(new_spectrum dst, f64 temperature);
u32 new_load_csv_file_to_spectrum(new_spectrum dst, const char *path);
