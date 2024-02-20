f64 lerp(f64 x, f64 x0, f64 x1, f64 y0, f64 y1)
{
    return y0 + ((x - x0) * ((y1 - y0) / (x1 - x0)));
}

f64 clamp(f64 f, f64 min, f64 max)
{
    f = (f < min) ? min : f;
    f = (f > max) ? max : f;
    return f;
}
