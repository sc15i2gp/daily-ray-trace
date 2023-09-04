#define PI 3.1415926535897932385L

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float       f32;
typedef double      f64;
typedef long double f128;
//64 bit rgb values 0.0-1.0
typedef struct
{
    union
    {
        struct
        {
            f64 r;
            f64 g;
            f64 b;
        };
        struct
        {
            f64 x;
            f64 y;
            f64 z;
        };
        f64 rgb[3];
        f64 xyz[3];
    };
}   rgb_f64;

//8 bit rgb values 0-255
typedef struct
{
    union
    {
        u32 value;
        struct
        {
            u8 b;
            u8 g;
            u8 r;
            u8 a;
        };
        u8  bgra[4];
    };
}   rgb_u8;
