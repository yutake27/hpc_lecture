#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void init_bucket(int *bucket) {
    bucket[blockIdx.x*blockDim.x+threadIdx.x] = 0;
}

__global__ void add_bucket(int *bucket, int *key) {
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    atomicAdd(&bucket[key[i]], 1);
}

__global__ void sort_key(int *bucket, int *key) {
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = 0;
    __syncthreads();
    for (int k=0;k<i;k++) {
        j+=bucket[k];
    }
    __syncthreads();
    for (; bucket[i]>0; bucket[i]--) {
        key[j++] = i;
    }
}

int main() {
  int n = 50;
  int range = 5;
  //std::vector<int> key(n);
  int *key;
  cudaMallocManaged(&key, n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  int *bucket;
  cudaMallocManaged(&bucket, range*sizeof(int));
  init_bucket<<<1, range>>>(bucket);
  add_bucket<<<1, n>>>(bucket, key);
  sort_key<<<1, range>>>(bucket, key);
  cudaDeviceSynchronize();

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  cudaFree(bucket);
  cudaFree(key);
}
