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
    f64 *samples;
}   spectrum;

struct
{
    spectrum rw;
    spectrum x;
    spectrum y;
    spectrum z;
}   cmfs;

struct
{
    spectrum white;
    spectrum red;
    spectrum green;
    spectrum blue;
    spectrum cyan;
    spectrum magenta;
    spectrum yellow;
}   rgb_spds;

struct
{
    u32 capacity;
    u32 allocated;
    f64 *base;
    f64 *next;
}   spd_stack;

struct
{
    u32 size;
    f64 *base;
    f64 *cmf_base;
    f64 *rgb_base;
    f64 *stack_base;
}   spd_alloc_table;

typedef struct
{
    const char *white;
    const char *cmf_x;
    const char *cmf_y;
    const char *cmf_z;
    const char *rgb_red;
    const char *rgb_green;
    const char *rgb_blue;
    const char *rgb_cyan;
    const char *rgb_magenta;
    const char *rgb_yellow;
} spd_tables_csvs;

spectrum alloc_spd();
void free_spd(spectrum);
void init_spd_tables(spd_tables_csvs csvs, u32 stack_capacity, f64 low_wl, f64 hi_wl, f64 wl_interval);
rgb_f64 spectrum_to_xyz(spectrum spd);
rgb_f64 spectrum_to_rgb_f64(spectrum spd);
void rgb_f64_to_spectrum(rgb_f64 rgb, spectrum dst);

f64  value_at_wl(spectrum, f64);
void zero_spectrum(spectrum dst);
void const_spectrum(spectrum dst, f64 f);
void copy_spectrum(spectrum dst, spectrum src);
void spectrum_normalise(spectrum dst);
void spectral_mul_by_spectrum(spectrum dst, spectrum src0, spectrum src1); //dst = src0 * src1
void spectral_mul_by_scalar(spectrum dst, spectrum src0, f64 src1); //dst = src0 * src1
void spectral_acc(spectrum dst, spectrum src); //dst += src
void spectral_acc_mul_by_scalar(spectrum dst, spectrum src0, f64 src1); //dst += src0 * src1
void spectral_sum(spectrum dst, spectrum src0, spectrum src1); //dst = src0 + src1
void spectral_div_by_scalar(spectrum dst, spectrum src0, f64 src1); //dst = src0/src1
void spectral_sub(spectrum dst, spectrum src0, spectrum src1); //dst = src0 - src1

void generate_blackbody_spectrum(spectrum dst, f64 temperature);
