u32 urng()
{
    u32 u = (u32)rand();
    return u;
}

u32 range_urng(u32 low, u32 high)
{
    u32 u = urng();
    u32 r = u % (high + 1 - low) + low;
    return r;
}

f64 range_frng(f64 low, f64 high)
{
    u32 u = urng();
    f64 fu = (f64)u;
    f64 fr = low + fu/(((f64)RAND_MAX)/high);
    return fr;
}
