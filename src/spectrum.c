char *find_next_newline(char *c)
{
    while(1)
    {
        if(*c == '\n') return c;
        if(*c == EOF || *c == 0) return NULL;
        c += 1;
    }
}

char *find_next_number(char *c)
{
    while(1)
    {
        if(*c >= '0' && *c <= '9') return c;
        if(*c == EOF || *c == 0) return NULL;
        c += 1;
    }
}

char *find_next_char(char *c, char desired)
{
    while(1)
    {
        if(*c == desired) return c;
        if(*c == EOF || *c == 0) return NULL;
        c += 1;
    }
}

f64 lerp(f64 x, f64 x0, f64 x1, f64 y0, f64 y1)
{
    return y0 + ((x - x0) * ((y1 - y0) / (x1 - x0)));
}

void init_spd_table(u32 capacity, u32 num_samples, f64 low_wl, f64 hi_wl, f64 wl_interval)
{
    number_of_spectrum_samples = 69;
    smallest_wavelength        = 380.0;
    largest_wavelength         = 720.0;
    sample_interval            = 5.0;
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

//CSV file structure:
//  - First line ignored (for headings)
//  - Each line is a single sample
//  - A sample consists of two csvs: "wavelength, value"
u32 load_csv_file_to_spectrum(spectrum dst, const char *csv_path)
{
    printf("Opening csv file %s...\n", csv_path);
    file_handle csv_file_handle = open_file(csv_path, ACCESS_READ, FILE_EXISTS);
    u32 csv_file_size = get_file_size(csv_file_handle) + 1;
    char *csv_file_buffer = alloc(csv_file_size);
    printf("Reading csv file %s...\n", csv_path);
    read_file(csv_file_handle, csv_file_size, csv_file_buffer);
    csv_file_buffer[csv_file_size] = 0;
    close_file(csv_file_handle);
    printf("Closed csv file %s.\n", csv_path);

    f64 file_wavelengths[MAX_NUM_SPECTRUM_VALUES];
    f64 file_wavelength_values[MAX_NUM_SPECTRUM_VALUES];
    memset(file_wavelengths, 0.0, sizeof(file_wavelengths));
    memset(file_wavelength_values, 0.0, sizeof(file_wavelength_values));
    u32 file_number_of_samples = 0;

    //TODO: Add support for files with wavelengths in um
    //NOTE: Currently only supporting wavelengths in nm
    
    //Count number of samples in csv file
    printf("Counting number of samples in char buffer...\n");
    for(char *c = find_next_newline(csv_file_buffer); c != NULL; c = find_next_newline(c)) 
    {
        c += 1;
        if(find_next_number(c)) file_number_of_samples += 1;
    }
    //Read wavelengths and sample values
    printf("Reading sample data in char buffer...\n");
    char *c = find_next_newline(csv_file_buffer) + 1;
    for(u32 i = 0; i < file_number_of_samples; i += 1)
    {
        c = find_next_number(c);
        file_wavelengths[i] = atof(c);
        c = find_next_char(c, ',');
        c = find_next_number(c);
        file_wavelength_values[i] = atof(c);
        c = find_next_newline(c);
    }
    printf("Freeing char buffer...\n");
    unalloc(csv_file_buffer, csv_file_size);

    printf("Number of samples read = %u\n", file_number_of_samples);

    //Lerp values and write them to spectrum
    u32 file_wl_index = 0;
    for(u32 sample = 0; sample < number_of_spectrum_samples; sample += 1)
    {
        f64 sample_wavelength = smallest_wavelength + ((f64)sample) * sample_interval;
        for(; file_wavelengths[file_wl_index+1] < sample_wavelength; file_wl_index += 1);
        dst.samples[sample] = lerp(sample_wavelength, file_wavelengths[file_wl_index], file_wavelengths[file_wl_index + 1], file_wavelength_values[file_wl_index], file_wavelength_values[file_wl_index+1]);
    }

    /*
    printf("Wavelength | Value\n");
    for(u32 i = 0; i < file_number_of_samples; ++i)
    {
        f64 sample_wl = smallest_wavelength + ((f64)i) * sample_interval;
        printf("%f | %f\n", sample_wl, dst[i]);
    }
*/
    //return success;
    return 1;
}
