//Gens number between 0 and 1
f64 rng()
{
    f64 r = (f64)rand();
    f64 f = r / (f64)RAND_MAX;
    return f;
}

void seed_rng(u32 seed)
{
    srand(seed);
}
