#define KEYWORD(word) #word
const char *keyword_table[] = 
{
    #include "Keywords.h"
};
#undef KEYWORD
u32 num_keywords = sizeof(keyword_table)/sizeof(const char*);

#define KEYWORD(word) TOKEN_##word
typedef enum
{
    TOKEN_NONE,
    TOKEN_NEWLINE,
    TOKEN_SPACE,
    TOKEN_FLOAT,
    TOKEN_WORD,
    TOKEN_END,
    TOKEN_KEYWORD_START,
    #include "keywords.h"
    TOKEN_COUNT
} scene_token_type;
#undef KEYWORD

typedef struct
{
    scene_token_type type;
    char             *loc;
    u32              length;
    f64              value;
    u32              uint;
} scene_token;

void print_token(scene_token *t)
{
    printf("TYPE:%u,LOC:%p,LEN:%u,", t->type, t->loc, t->length);
    switch(t->type)
    {
        case TOKEN_NEWLINE:
        {
            printf("<NEW>");
            break;
        }
        case TOKEN_SPACE:
        {
            printf("<SPACE>");
            break;
        }
        case TOKEN_FLOAT:
        {
            printf("%f", t->value);
            break;
        }
        case TOKEN_END:
        {
            printf("<END>");
            break;
        }
        default:
        {
            printf("%.*s", t->length, t->loc);
            break;
        }
    }
}

u32 is_number(char c)
{
    return (c >= '0' && c <= '9');
}

u32 is_letter(char c)
{
    return ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'));
}

u32 is_word_char(char c)
{
    return is_letter(c) || is_number(c) || c == '_' || c == '.' || c == '\\';
}

u32 is_number_char(char c)
{
    return is_number(c) || c == '.' || c == '-';
}

u32 is_word_keyword(const char *word, const char *kw_candidate)
{
    u32 word_length = 0;
    for(const char *c = word; is_word_char(*c); c += 1, word_length += 1);
    u32 kw_candidate_length = 0;
    for(const char *c = kw_candidate; *c != 0; c += 1, kw_candidate_length += 1);

    if(word_length == kw_candidate_length)
    {
        u32 is_kw = 1;
        for(u32 i = 0; i < word_length; i += 1)
        {
            if(word[i] != kw_candidate[i])
            {
                is_kw = 0;
                break;
            }
        }
        return is_kw;
    }
    else
    {
        return 0;
    }
}

void read_token(scene_token *dst, char *src, char *src_end)
{
    dst->length = 0;
    dst->loc    = src;
    if(src >= src_end)
    {
        dst->type   = TOKEN_END;
        dst->loc    = src_end;
    }
    else if(src[0] == '\n' || src[0] == '\r')
    {
        dst->type   = TOKEN_NEWLINE;
        for(char *c = src; c < src_end && (*c == '\n' || *c == '\r'); c += 1, dst->length += 1);
    }
    else if(is_number_char(src[0]))
    {
        dst->type  = TOKEN_FLOAT;
        dst->value = atof(src);
        for(dst->length = 0; src < src_end && is_number_char(*src); src += 1, dst->length += 1);
    }
    else if(is_letter(src[0]))
    {
        dst->type   = TOKEN_NONE;
        for(char *c = src; c < src_end && is_word_char(*c); c += 1, dst->length += 1);
        for(u32 i = 0; i < num_keywords; i += 1)
        {
            if(is_word_keyword(src, keyword_table[i]))
            {
                dst->type = TOKEN_KEYWORD_START + 1 + i;
                break;
            }
        }
        if(dst->type == TOKEN_NONE)
        {
            dst->type = TOKEN_WORD;
        }
    }
    else
    {
        for(src; src < src_end && (!(is_letter(*src) || is_number_char(*src))); src += 1);
        dst->type   = TOKEN_SPACE;
        dst->length = src - dst->loc;
    }
}

typedef struct
{
    u32  scene_contents_length;
    char *scene_contents;
    char *scene_contents_end;
    char *scene_loc;

    scene_token current_token;
    scene_token lookahead_token;
} scene_tokeniser;

scene_tokeniser tokeniser;

scene_token *next_token()
{
    char *loc = tokeniser.scene_loc;
    do
    {
        read_token(&tokeniser.current_token, loc, tokeniser.scene_contents_end);
        loc += tokeniser.current_token.length;
    }   while(tokeniser.current_token.type == TOKEN_SPACE || tokeniser.current_token.type == TOKEN_NEWLINE);
    tokeniser.scene_loc = loc;

    return &tokeniser.current_token;
}

scene_token *lookahead_token()
{
    char *loc = tokeniser.scene_loc;
    do
    {
        read_token(&tokeniser.lookahead_token, loc, tokeniser.scene_contents_end);
        loc += tokeniser.lookahead_token.length;
    }   while(tokeniser.lookahead_token.type == TOKEN_SPACE || tokeniser.lookahead_token.type == TOKEN_NEWLINE);

    return &tokeniser.lookahead_token;
}

//TODO: Give more error info
void parse_error()
{
    printf("ERROR: Parse error at: %p out of %p\n", tokeniser.scene_loc, tokeniser.scene_loc, tokeniser.scene_contents_end);
    printf("AT TOKEN: ");
    print_token(&tokeniser.current_token);
    printf("\n");
    exit(-1);
}

void parse_float(f64 *dst)
{
    scene_token *t = next_token();
    if(t->type != TOKEN_FLOAT) parse_error();
    *dst = t->value;
}

//TODO: Get rid of this shit
void parse_uint(u32 *dst)
{
    f64 f;
    parse_float(&f);
    *dst = (u32)(f);
}

void parse_bool(u32 *dst)
{
    scene_token *t = next_token();
    switch(t->type)
    {
        case TOKEN_true:
        {
            *dst = 1;
            break;
        }
        case TOKEN_false:
        {
            *dst = 0;
            break;
        }
        default:
        {
            parse_error();
            break;
        }
    }
}

void parse_vec3(vec3 *dst)
{
    parse_float(&dst->x);
    parse_float(&dst->y);
    parse_float(&dst->z);
}

void parse_rgb(rgb_f64 *dst)
{
    parse_float(&dst->r);
    parse_float(&dst->g);
    parse_float(&dst->b);
}

void parse_word(char *dst)
{
    scene_token *t = next_token();
    memcpy(dst, t->loc, t->length);
}

void parse_spd_method(spd_input_data *dst)
{
    scene_token *t = next_token();
    switch(t->type)
    {
        case TOKEN_rgb:
        {
            dst->method = SPD_METHOD_RGB;
            parse_rgb(&dst->rgb);
            break;
        }
        case TOKEN_csv:
        {
            dst->method = SPD_METHOD_CSV;
            parse_word(dst->csv);
            break;
        }
        case TOKEN_blackbody:
        {
            dst->method = SPD_METHOD_BLACKBODY;
            parse_float(&dst->blackbody_temp);
            break;
        }
        case TOKEN_constant:
        {
            dst->method = SPD_METHOD_CONST;
            parse_float(&dst->constant);
            break;
        }
        default:
        {
            parse_error();
            break;
        }
    }
    scene_token *l = lookahead_token();
    if(l->type == TOKEN_scale)
    {
        dst->has_scale_factor = 1;
        t = next_token();
        parse_float(&dst->scale_factor);
    }
}

//TODO: Change this garbage
void parse_bdsfs(u32 *num_bdsfs, bdsf_func *bdsfs)
{
    char bdsf_name[64];
    *num_bdsfs = 0;
    scene_token *l = lookahead_token();
    while(l->type != TOKEN_dir_func)
    {
        memset(bdsf_name, 0, 64);
        parse_word(bdsf_name);
        for(u32 i = 0; i < num_bdsfs_defined; i += 1)
        {
            if(strcmp(bdsf_name, bdsf_name_list[i]) == 0)
            {
                bdsfs[*num_bdsfs] = bdsf_list[i];
                break;
            }
        }
        *num_bdsfs += 1;
        l = lookahead_token();
    }
}

void parse_dir_func(dir_func *sample_direction_function)
{
    char dir_func_name[64];
    memset(dir_func_name, 0, 64);
    parse_word(dir_func_name);
    for(u32 i = 0; i < num_dir_funcs_defined; i += 1)
    {
        if(strcmp(dir_func_name, dir_func_name_list[i]) == 0)
        {
            *sample_direction_function = dir_func_list[i];
            break;
        }
    }
}

void parse_camera(camera_input_data *camera)
{
    scene_token *l = lookahead_token();
    while(l->type != TOKEN_Camera && l->type != TOKEN_Material && l->type != TOKEN_Surface && l->type != TOKEN_END)
    {
        scene_token *t = next_token();
        switch(t->type)
        {
            case TOKEN_position:
            {
                parse_vec3(&camera->position);
                break;
            }
            case TOKEN_up:
            {
                parse_vec3(&camera->up);
                break;
            }
            case TOKEN_right:
            {
                parse_vec3(&camera->right);
                break;
            }
            case TOKEN_forward:
            {
                parse_vec3(&camera->forward);
                break;
            }
            case TOKEN_fov:
            {
                parse_float(&camera->fov);
                break;
            }
            case TOKEN_fdepth:
            {
                parse_float(&camera->fdepth);
                break;
            }
            case TOKEN_flength:
            {
                parse_float(&camera->flength);
                break;
            }
            case TOKEN_aperture:
            {
                parse_float(&camera->aperture);
                break;
            }
            default:
            {
                parse_error();
                break;
            }
        }
        l = lookahead_token();
    }
}

void parse_material(scene_input_data *scene)
{
    material_input_data *dst_material = &scene->scene_materials[scene->num_scene_materials];
    scene->num_scene_materials += 1;

    scene_token *l = lookahead_token();
    while(l->type != TOKEN_Camera && l->type != TOKEN_Material && l->type != TOKEN_Surface && l->type != TOKEN_END)
    {
        scene_token *t = next_token();
        switch(t->type)
        {
            case TOKEN_name:
            {
                parse_word(dst_material->name);
                break;
            }
            case TOKEN_diffuse:
            {
                parse_spd_method(&dst_material->diffuse_input);
                break;
            }
            case TOKEN_glossy:
            {
                parse_spd_method(&dst_material->glossy_input);
                break;
            }
            case TOKEN_emission:
            {
                parse_spd_method(&dst_material->emission_input);
                dst_material->is_emissive = 1;
                break;
            }
            case TOKEN_mirror:
            {
                parse_spd_method(&dst_material->mirror_input);
                break;
            }
            case TOKEN_is_black_body:
            {
                parse_bool(&dst_material->is_black_body);
                break;
            }
            case TOKEN_shininess:
            {
                parse_float(&dst_material->shininess);
                break;
            }
            case TOKEN_bdsfs:
            {
                parse_bdsfs(&dst_material->num_bdsfs, dst_material->bdsfs);
                break;
            }
            case TOKEN_dir_func:
            {
                parse_dir_func(&dst_material->sample_direction_function);
                break;
            }
            case TOKEN_refract:
            {
                parse_spd_method(&dst_material->refract_input);
                break;
            }
            case TOKEN_extinct:
            {
                parse_spd_method(&dst_material->extinct_input);
                break;
            }
            case TOKEN_base_material:
            {
                dst_material->is_base_material = 1;
                break;
            }
            case TOKEN_escape_material:
            {
                dst_material->is_escape_material = 1;
                break;
            }
            default:
            {
                parse_error();
                break;
            }
        }
        l = lookahead_token();
    }
}

void parse_surface(scene_input_data *scene)
{
    surface_input_data *dst_surface = &scene->surfaces[scene->num_surfaces];
    scene->num_surfaces += 1;

    scene_token *l = lookahead_token();
    while(l->type != TOKEN_Camera && l->type != TOKEN_Material && l->type != TOKEN_Surface && l->type != TOKEN_END)
    {
        scene_token *t = next_token();
        switch(t->type)
        {
            case TOKEN_name:
            {
                parse_word(dst_surface->name);
                break;
            }
            case TOKEN_type:
            {
                t = next_token();
                switch(t->type)
                {
                    case TOKEN_point:
                    {
                        dst_surface->type = GEO_TYPE_POINT;
                        break;
                    }
                    case TOKEN_sphere:
                    {
                        dst_surface->type = GEO_TYPE_SPHERE;
                        break;
                    }
                    case TOKEN_plane:
                    {
                        dst_surface->type = GEO_TYPE_PLANE;
                        break;
                    }
                    default:
                    {
                        parse_error();
                        break;
                    }
                }
                break;
            }
            case TOKEN_position:
            {
                parse_vec3(&dst_surface->position);
                break;
            }
            case TOKEN_radius:
            {
                parse_float(&dst_surface->radius);
                break;
            }
            case TOKEN_pointu:
            {
                parse_vec3(&dst_surface->u);
                break;
            }
            case TOKEN_pointv:
            {
                parse_vec3(&dst_surface->v);
                break;
            }
            case TOKEN_material:
            {
                parse_word(dst_surface->material_name);
                break;
            }
            default:
            {
                parse_error();
                break;
            }
        }
        l = lookahead_token();
    }
}

void parse_scene(char *scene_file_contents, u32 scene_file_size, camera_input_data *camera, scene_input_data *scene)
{
    tokeniser.scene_contents        = scene_file_contents;
    tokeniser.scene_loc             = scene_file_contents;
    tokeniser.scene_contents_length = scene_file_size;
    tokeniser.scene_contents_end    = scene_file_contents + scene_file_size;

    for(scene_token *t = next_token(); t->type != TOKEN_END; t = next_token())
    {
        switch(t->type)
        {
            case TOKEN_Camera:
            {
                parse_camera(camera);
                break;
            }
            case TOKEN_Material:
            {
                parse_material(scene);
                break;
            }
            case TOKEN_Surface:
            {
                parse_surface(scene);
                break;
            }
            default:
            {
                parse_error();
                break;
            }
        }
    }
}

void parse_config(char *config_contents, u32 config_contents_size, config_arguments *config)
{
    tokeniser.scene_contents        = config_contents;
    tokeniser.scene_loc             = config_contents;
    tokeniser.scene_contents_length = config_contents_size;
    tokeniser.scene_contents_end    = config_contents + config_contents_size;

    for(scene_token *t = next_token(); t->type != TOKEN_END; t = next_token())
    {
        switch(t->type)
        {
            case TOKEN_num_pixel_samples:
            {
                parse_uint(&config->num_pixel_samples);
                break;
            }
            case TOKEN_output_width:
            {
                parse_uint(&config->output_width);
                break;
            }
            case TOKEN_output_height:
            {
                parse_uint(&config->output_height);
                break;
            }
            case TOKEN_min_wl:
            {
                parse_float(&config->min_wl);
                break;
            }
            case TOKEN_max_wl:
            {
                parse_float(&config->max_wl);
                break;
            }
            case TOKEN_wl_interval:
            {
                parse_float(&config->wl_interval);
                break;
            }
            case TOKEN_input_scene:
            {
                parse_word(config->input_scene);
                break;
            }
            case TOKEN_output_spd:
            {
                parse_word(config->output_spd);
                break;
            }
            case TOKEN_average_spd:
            {
                parse_word(config->average_spd);
                break;
            }
            case TOKEN_variance_spd:
            {
                parse_word(config->variance_spd);
                break;
            }
            case TOKEN_output_bmp:
            {
                parse_word(config->output_bmp);
                break;
            }
            case TOKEN_average_bmp:
            {
                parse_word(config->average_bmp);
                break;
            }
            case TOKEN_variance_bmp:
            {
                parse_word(config->variance_bmp);
                break;
            }
            default:
            {
                parse_error();
                break;
            }
        }
    }
}

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
