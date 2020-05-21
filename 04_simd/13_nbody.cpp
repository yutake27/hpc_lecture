#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>


float sum_vec(__m256 vec, int N) {
    __m256 sumvec = _mm256_permute2f128_ps(vec, vec, 1);
    sumvec = _mm256_add_ps(sumvec, vec);
    sumvec = _mm256_hadd_ps(sumvec, sumvec);
    sumvec = _mm256_hadd_ps(sumvec, sumvec);
    float sum[N];
    _mm256_store_ps(sum, sumvec);    
    return sum[0];
}


int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  for(int i=0; i<N; i++) {
    /*for(int j=0; j<N; j++) {
      if(i != j) {
        float rx = x[i] - x[j];
        float ry = y[i] - y[j];
        float r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
      }
    }*/
    // rx = x[i]-x[j]
    __m256 xvec = _mm256_load_ps(x);
    __m256 xivec = _mm256_set1_ps(x[i]);
    __m256 rxvec = _mm256_sub_ps(xivec, xvec);
    // ry = y[i]-y[j]
    __m256 yvec = _mm256_load_ps(y);
    __m256 yivec = _mm256_set1_ps(y[i]);
    __m256 ryvec = _mm256_sub_ps(yivec, yvec);
    // r = rsqrt(rx*rx+ry*ry)
    __m256 rvec = _mm256_rsqrt_ps(_mm256_add_ps(_mm256_mul_ps(rxvec, rxvec), _mm256_mul_ps(ryvec, ryvec)));
    //  r[i] = 0 if i==j 
    __m256 mask = _mm256_cmp_ps(rvec, _mm256_set1_ps(INFINITY), _CMP_NEQ_OQ);
    rvec = _mm256_blendv_ps(_mm256_setzero_ps(), rvec, mask);
    // r3 = r**3
    __m256 r3vec = _mm256_mul_ps(rvec, _mm256_mul_ps(rvec, rvec));
    // make mask
    __m256 mvec = _mm256_load_ps(m);
    mask = _mm256_set1_ps(1);
    mask[i] = 0;
    //  m[i] = 0 if i==j
    mvec = _mm256_mul_ps(mvec, mask);
    // rx*m[j]*r**3
    __m256 fxivec = _mm256_mul_ps(_mm256_mul_ps(rxvec, mvec), r3vec);
    // ry*m[j]*r**3
    __m256 fyivec = _mm256_mul_ps(_mm256_mul_ps(ryvec, mvec), r3vec);
    // caluculate sum of vec
    fx[i] -= sum_vec(fxivec, N);
    fy[i] -= sum_vec(fyivec, N);

    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
