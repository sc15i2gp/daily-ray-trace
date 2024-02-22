void init_spd_table(u32 capacity, f64 low_wl, f64 hi_wl, f64 wl_interval)
{
    u32 num_samples = (u32)(((hi_wl - low_wl)/wl_interval) + 1.0);
    number_of_spectrum_samples = num_samples;
    smallest_wavelength        = low_wl;
    largest_wavelength         = hi_wl;
    sample_interval            = wl_interval;
    spectrum_size = number_of_spectrum_samples * sizeof(f64);

    spd_table.capacity        = capacity;
    spd_table.allocated       = 0;
    spd_table.is_allocated    = alloc(spd_table.capacity * sizeof(u32));
    spd_table.spectrum_buffer = alloc(spectrum_size * spd_table.capacity);
}

spectrum alloc_spd()
{
    spectrum spd = {NULL};
    if(spd_table.allocated == spd_table.capacity)
    {
        printf("WARN: RUNNING OUT OF SPECTRA %u\n", spd_table.allocated);
    }
    for(u32 i = 0; i < spd_table.capacity; i += 1)
    {
        if(!spd_table.is_allocated[i])
        {
            spd_table.is_allocated[i] = 1;
            spd_table.allocated += 1;
            spd.samples = spd_table.spectrum_buffer + i * number_of_spectrum_samples; 
            break;
        }
    }
    return spd;
}

void free_spd(spectrum spd)
{
    u32 spd_offset = (u32)(spd.samples - spd_table.spectrum_buffer);
    u32 spd_index  = spd_offset / number_of_spectrum_samples;
    spd_table.is_allocated[spd_index] = 0;
    spd_table.allocated -= 1;
}

void zero_spectrum(spectrum dst)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        dst.samples[i] = 0.0;
    }
}

void const_spectrum(spectrum dst, f64 value)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1) dst.samples[i] = value;
}

void copy_spectrum(spectrum dst, spectrum src)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1) dst.samples[i] = src.samples[i];
}

void spectrum_normalise(spectrum dst)
{
    f64 highest_value = 0.0;
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1) if(dst.samples[i] > highest_value) highest_value = dst.samples[i];
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1) dst.samples[i] /= highest_value;
}

void spectral_mul_by_spectrum(spectrum dst, spectrum src0, spectrum src1)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        dst.samples[i] = src0.samples[i] * src1.samples[i];
    }
}

void spectral_mul_by_scalar(spectrum dst, spectrum src, f64 d)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        dst.samples[i] = src.samples[i] * d;
    }
}

void spectral_acc(spectrum dst, spectrum src)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        dst.samples[i] += src.samples[i];
    }
}

void spectral_acc_mul_by_scalar(spectrum dst, spectrum src, f64 d)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        dst.samples[i] += src.samples[i] * d;
    }
}

void spectral_sum(spectrum dst, spectrum src0, spectrum src1)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        dst.samples[i] = src0.samples[i] + src1.samples[i];
    }
}

void spectral_div_by_scalar(spectrum dst, spectrum src0, f64 d)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        dst.samples[i] = src0.samples[i] / d;
    }
}

void spectral_sub(spectrum dst, spectrum src0, spectrum src1)
{
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        dst.samples[i] = src0.samples[i] - src1.samples[i];
    }
}

rgb_f64 spectrum_to_xyz(spectrum spd, spectrum cmf_x, spectrum cmf_y, spectrum cmf_z, spectrum ref_white)
{
    rgb_f64 xyz = {0.0, 0.0, 0.0};
    f64 n = 0.0;
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        n += (cmf_y.samples[i] * ref_white.samples[i]);
    }
    n *= sample_interval;

    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        xyz.x += (cmf_x.samples[i] * spd.samples[i] * ref_white.samples[i]);
        xyz.y += (cmf_y.samples[i] * spd.samples[i] * ref_white.samples[i]);
        xyz.z += (cmf_z.samples[i] * spd.samples[i] * ref_white.samples[i]);
    }
    xyz.x *= (sample_interval/n);
    xyz.y *= (sample_interval/n);
    xyz.z *= (sample_interval/n);

    return xyz;
}

rgb_f64 spectrum_to_rgb_f64(spectrum spd, spectrum cmf_x, spectrum cmf_y, spectrum cmf_z, spectrum ref_white)
{
    rgb_f64 spd_xyz = spectrum_to_xyz(spd, cmf_x, cmf_y, cmf_z, ref_white);
    
    rgb_f64 linear_rgb;
    linear_rgb.r = (2.3706743 * spd_xyz.x) - (0.9000405 * spd_xyz.y) - (0.4706338 * spd_xyz.z);
    linear_rgb.g = (-0.5138850 * spd_xyz.x) + (1.4253036 * spd_xyz.y) + (0.0885814 * spd_xyz.z);
    linear_rgb.b = (0.0052982 * spd_xyz.x) - (0.0146949 * spd_xyz.y) + (1.0093968 * spd_xyz.z);

    return linear_rgb;
}

//RGB to spectrum
//  - smallest_rgb
//  - mid_rgb
//  - largest_rgb
//  - rgb_to_spd is spectrum of largest_rgb to spd
//  - cym_to_spd is:
//      - If R is smallest_rgb then cyan
//      - If G is smallest_rgb then magenta
//      - If B is smallest_rgb then yellow
//  - spectrum = smallest_rgb * white_rgb_spd + cym_spd * (mid_rgb - smallest_rgb) + rgb_spd * (largest_rgb - mid_rgb)

void rgb_f64_to_spectrum(   rgb_f64 rgb,         spectrum dst, 
                            spectrum white_spd, spectrum red_spd, 
                            spectrum green_spd, spectrum blue_spd,
                            spectrum cyan_spd,  spectrum magenta_spd,
                            spectrum yellow_spd)
{
    spectrum rgb_spectra[] = {red_spd,  green_spd,   blue_spd};
    spectrum cmy_spectra[] = {cyan_spd, magenta_spd, yellow_spd};
    f64 rgb_values[] =        {rgb.r,    rgb.g,       rgb.b};
    u32 indices[] = {0, 1, 2};
    u32 tmp_val;

    if(rgb_values[indices[0]] > rgb_values[indices[1]])
    {
        tmp_val = indices[1];
        indices[1] = indices[0];
        indices[0] = tmp_val;
    }
    if(rgb_values[indices[1]] > rgb_values[indices[2]])
    {
        tmp_val = indices[2];
        indices[2] = indices[1];
        indices[1] = tmp_val;
    }
    if(rgb_values[indices[0]] > rgb_values[indices[1]])
    {
        tmp_val = indices[1];
        indices[1] = indices[0];
        indices[0] = tmp_val;
    }

    u32 small_index = indices[0];
    u32 mid_index = indices[1];
    u32 large_index = indices[2];
    f64 rgb_diff_mid_small = rgb_values[mid_index] - rgb_values[small_index];
    f64 rgb_diff_large_mid = rgb_values[large_index] - rgb_values[mid_index];
    spectral_mul_by_scalar(dst, white_spd, rgb_values[small_index]); //dst = white_rgb_spd * smallest_rgb
    spectral_acc_mul_by_scalar(dst, cmy_spectra[small_index], rgb_diff_mid_small); //dst += cym_spd * (mid - small)
    spectral_acc_mul_by_scalar(dst, rgb_spectra[large_index], rgb_diff_large_mid); //dst += rgb_spd * (large - mid)
}

f128 compute_blackbody_power(f128 temperature, f128 wavelength)
{
    const f128 c = 2.99792458e8L; //Speed of light
    const f128 h = 6.626176e-34L; //Planck's constant
    const f128 k = 1.380662e-23L; //Boltzmann constant
    f128 numerator = 2.0L * PI * h * c * c;

    f128 lambda_5 = powl(wavelength, 5.0L);
    f128 e_power_numerator = (h * c) / k;
    f128 e_power_denominator = temperature * wavelength;
    f128 e_power = e_power_numerator / e_power_denominator;
    f128 e_term = expl(e_power);
    f128 denominator = lambda_5 * (e_term - 1.0L);
    f128 power = numerator / denominator;

    return power;
}

void generate_blackbody_spectrum(spectrum spd, f64 temperature)
{
    f128 temp_long = (f128)temperature;
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        f128 wavelength_in_nm = (f128)(smallest_wavelength + (i * sample_interval));
        f128 wavelength_in_m = wavelength_in_nm * 1e-9L;
        f128 blackbody_power = compute_blackbody_power(temp_long, wavelength_in_m);
        spd.samples[i] = (f64)(blackbody_power * 1e9L);
    }
}
