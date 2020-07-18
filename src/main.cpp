#include <stdio.h>
#include <immintrin.h>
#include "Platform.h"
//NOTE: Alignment is not preserved when calling a function (nor when declaring on stack w/ alignment > 16)
//NOTE: Possible way of speeding up spectra (after verifying spectra are where most of the program is spent):
//	- Preallocate array of spectra to copy into
//	- Functions which call sum, mult etc, receive a pointer to the dst in this array

#define SPECTRUM_SAMPLE_MAX 512

int number_of_samples = 512;
struct alignas(32) Spectrum
{
	double samples[SPECTRUM_SAMPLE_MAX];
};

Spectrum* __aligned_spectra__ = NULL;

Spectrum gen_spectrum(double d)
{
	Spectrum s = {};
	for(int i = 0; i < number_of_samples; ++i)
	{
		s.samples[i] = d;
	}
	return s;
}

#if 0
Spectrum operator+(Spectrum s, Spectrum t)
{
	Spectrum u = {};
	for(int i = 0; i < number_of_samples; ++i)
	{
		u.samples[i] = s.samples[i] + t.samples[i];
	}
	return u;
}

void sum(Spectrum* s, Spectrum* t, Spectrum* u)
{
	for(int i = 0; i < number_of_samples; ++i) s->samples[i] = t->samples[i] + u->samples[i];
}
#else
Spectrum operator+(Spectrum s, Spectrum  t)
{
	Spectrum* _s = __aligned_spectra__;
	Spectrum* _t = __aligned_spectra__ + 1;
	Spectrum* u = __aligned_spectra__ + 2;
	/*
	memcpy(_s, &s, sizeof(Spectrum));
	memcpy(_t, &t, sizeof(Spectrum));
	*/
	*_s = s;
	*_t = t;
	__m256d m_s;
	__m256d m_t;
	__m256d m_u;
	unsigned int m = (number_of_samples + 3) & ~3;
	for(int i = 0; i < m; i += 4)
	{
		m_s = _mm256_load_pd(_s->samples + i);
		m_t = _mm256_load_pd(_t->samples + i);
		m_u = _mm256_add_pd(m_s, m_t);
		_mm256_store_pd(u->samples + i, m_u);
	}
	return *u;
}

//FASTEST (SIGNIFICANTLY)
void sum(Spectrum* s, Spectrum* t, Spectrum* u)
{
	__m256d m_s;
	__m256d m_t;
	__m256d m_u;
	unsigned int m = (number_of_samples + 3) & ~3;
	for(int i = 0; i < m; i += 4)
	{
		m_s = _mm256_load_pd(s->samples + i);
		m_t = _mm256_load_pd(t->samples + i);
		m_u = _mm256_add_pd(m_s, m_t);
		_mm256_store_pd(u->samples + i, m_u);
	}
}
#endif

void mult(Spectrum* s, double d)
{
	for(int i = 0; i < number_of_samples; ++i) s->samples[i] *= d;
}

int main()
{
	__aligned_spectra__ = (Spectrum*)alloc(3 * sizeof(Spectrum));

	Spectrum a = gen_spectrum(1.0);
	Spectrum b = gen_spectrum(2.0);
	query_pc_frequency();
	Timer t = {};
	start_timer(&t);
	Spectrum c = {}; 
	int m = 150000;
	for(int i = 0; i < m; ++i) 
	{
		c = a + b;
		//sum(&c, &a, &b);
	}
	stop_timer(&t);
	double elapsed = elapsed_time_in_ms(&t);
	printf("Elapsed = %f\n", elapsed);
	return 0;
}
