void init_spd_tables(spd_tables_csvs csvs, u32 stack_capacity, f64 low_wl, f64 hi_wl, f64 wl_interval)
{
    u32 num_samples = (u32)(((hi_wl - low_wl)/wl_interval) + 1.0);
    number_of_spectrum_samples = num_samples;
    smallest_wavelength        = low_wl;
    largest_wavelength         = hi_wl;
    sample_interval            = wl_interval;
    spectrum_size = number_of_spectrum_samples * sizeof(f64);

    u32 num_cmfs = 4;
    u32 num_rgbs = 7;
    spd_alloc_table.size = (num_cmfs + num_rgbs + stack_capacity) * spectrum_size;
    spd_alloc_table.base = alloc(spd_alloc_table.size);
    spd_alloc_table.cmf_base = spd_alloc_table.base;
    spd_alloc_table.rgb_base = spd_alloc_table.cmf_base + num_cmfs * number_of_spectrum_samples;
    spd_alloc_table.stack_base = spd_alloc_table.rgb_base + num_rgbs * number_of_spectrum_samples;

    cmfs.rw.samples = spd_alloc_table.cmf_base + 0 * number_of_spectrum_samples;
    cmfs.x.samples  = spd_alloc_table.cmf_base + 1 * number_of_spectrum_samples;
    cmfs.y.samples  = spd_alloc_table.cmf_base + 2 * number_of_spectrum_samples;
    cmfs.z.samples  = spd_alloc_table.cmf_base + 3 * number_of_spectrum_samples;

    rgb_spds.white.samples   = spd_alloc_table.rgb_base + 0 * number_of_spectrum_samples;
    rgb_spds.red.samples     = spd_alloc_table.rgb_base + 1 * number_of_spectrum_samples;
    rgb_spds.green.samples   = spd_alloc_table.rgb_base + 2 * number_of_spectrum_samples;
    rgb_spds.blue.samples    = spd_alloc_table.rgb_base + 3 * number_of_spectrum_samples;
    rgb_spds.cyan.samples    = spd_alloc_table.rgb_base + 4 * number_of_spectrum_samples;
    rgb_spds.magenta.samples = spd_alloc_table.rgb_base + 5 * number_of_spectrum_samples;
    rgb_spds.yellow.samples  = spd_alloc_table.rgb_base + 6 * number_of_spectrum_samples;

    spd_stack.capacity  = stack_capacity;
    spd_stack.allocated = 0;
    spd_stack.base      = spd_alloc_table.stack_base;
    spd_stack.next      = spd_stack.base;

    load_csv_file_to_spectrum(cmfs.rw, csvs.white);
    load_csv_file_to_spectrum(cmfs.x, csvs.cmf_x);
    load_csv_file_to_spectrum(cmfs.y, csvs.cmf_y);
    load_csv_file_to_spectrum(cmfs.z, csvs.cmf_z);
    load_csv_file_to_spectrum(rgb_spds.white, csvs.white);
    load_csv_file_to_spectrum(rgb_spds.red, csvs.rgb_red);
    load_csv_file_to_spectrum(rgb_spds.blue, csvs.rgb_blue);
    load_csv_file_to_spectrum(rgb_spds.green, csvs.rgb_green);
    load_csv_file_to_spectrum(rgb_spds.cyan, csvs.rgb_cyan);
    load_csv_file_to_spectrum(rgb_spds.magenta, csvs.rgb_magenta);
    load_csv_file_to_spectrum(rgb_spds.yellow, csvs.rgb_yellow);
}

rgb_f64 spectrum_to_xyz(spectrum spd)
{
    rgb_f64 xyz = {0.0, 0.0, 0.0};
    f64 n = 0.0;
    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        n += (cmfs.y.samples[i] * cmfs.rw.samples[i]);
    }
    n *= sample_interval;

    for(u32 i = 0; i < number_of_spectrum_samples; i += 1)
    {
        xyz.x += (cmfs.x.samples[i] * spd.samples[i] * cmfs.rw.samples[i]);
        xyz.y += (cmfs.y.samples[i] * spd.samples[i] * cmfs.rw.samples[i]);
        xyz.z += (cmfs.z.samples[i] * spd.samples[i] * cmfs.rw.samples[i]);
    }
    xyz.x *= (sample_interval/n);
    xyz.y *= (sample_interval/n);
    xyz.z *= (sample_interval/n);

    return xyz;
}

rgb_f64 spectrum_to_rgb_f64(spectrum spd)
{
    rgb_f64 spd_xyz = spectrum_to_xyz(spd);
    
    rgb_f64 linear_rgb;
    linear_rgb.r = (2.3706743 * spd_xyz.x) - (0.9000405 * spd_xyz.y) - (0.4706338 * spd_xyz.z);
    linear_rgb.g = (-0.5138850 * spd_xyz.x) + (1.4253036 * spd_xyz.y) + (0.0885814 * spd_xyz.z);
    linear_rgb.b = (0.0052982 * spd_xyz.x) - (0.0146949 * spd_xyz.y) + (1.0093968 * spd_xyz.z);

    return linear_rgb;
}

void rgb_f64_to_spectrum(rgb_f64 rgb, spectrum dst)
{
    spectrum rgb_spectra[] = {rgb_spds.red,  rgb_spds.green,   rgb_spds.blue};
    spectrum cmy_spectra[] = {rgb_spds.cyan, rgb_spds.magenta, rgb_spds.yellow};
    f64 rgb_values[] =       {rgb.r, rgb.g, rgb.b};
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
    spectral_mul_by_scalar(dst, rgb_spds.white, rgb_values[small_index]); //dst = white_rgb_spd * smallest_rgb
    spectral_acc_mul_by_scalar(dst, cmy_spectra[small_index], rgb_diff_mid_small); //dst += cym_spd * (mid - small)
    spectral_acc_mul_by_scalar(dst, rgb_spectra[large_index], rgb_diff_large_mid); //dst += rgb_spd * (large - mid)
}

spectrum alloc_spd()
{
    if(spd_stack.allocated == spd_stack.capacity)
    {
        printf("OUT OF STACK SPECTRA\n");
        exit(-1);
    }

    spectrum dst;
    dst.samples = spd_stack.next;

    spd_stack.allocated += 1;
    spd_stack.next += number_of_spectrum_samples;

    return dst;
}

void free_spd(spectrum spd)
{
    if(spd_stack.allocated == 0)
    {
        printf("FREEING TOO MANY SPECTRA\n");
        exit(-1);
    }

    spd_stack.allocated -= 1;
    spd_stack.next -= number_of_spectrum_samples;
}

f64 value_at_wl(spectrum spd, f64 wl)
{
    u32 i_0 = (u32)((wl - smallest_wavelength) / sample_interval);
    u32 i_1 = i_0 + 1;

    f64 w_0 = smallest_wavelength + i_0*sample_interval;
    f64 w_1 = smallest_wavelength + i_1*sample_interval;
    f64 s_0 = spd.samples[i_0];
    f64 s_1 = spd.samples[i_1];
    
    f64 v = lerp(wl, w_0, w_1, s_0, s_1);
    return v;
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
