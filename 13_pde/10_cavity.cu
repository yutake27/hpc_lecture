#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;
const int nx = 41;
const int ny = 41;
const int size = nx * ny * sizeof(float);
const int N = 1024;
const dim3 grid((nx+N-1)/N, ny);


void save_matrix(float *u, float *v, float *p) {
    ofstream u_file("u.txt");
    ofstream v_file("v.txt");
    ofstream p_file("p.txt");
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            u_file << u[ny*i+j] << ' ';
            v_file << v[ny*i+j] << ' ';
            p_file << p[ny*i+j] << ' ';
        }
        u_file << endl;
        v_file << endl;
        p_file << endl;
    }
    u_file.close();
    v_file.close();
    p_file.close();
}


__global__ void build_up_b(float *b, float rho, float dt, float *u, float *v, float dx, float dy) {
    int i = blockIdx.y;
    int j = threadIdx.x + blockDim.x * blockIdx.x;
    if (i!=0 && j!=0 && i<nx-1 && j<ny-1) {
        b[ny*i+j] = (rho * (1 / dt *
                      ((u[ny*i+j+1] - u[ny*i+j-1]) /
                       (2 * dx) + (v[ny*(i+1)+j] - v[ny*(i-1)+j]) / (2 * dy)) -
                      ((u[ny*i+j+1] - u[ny*i+j-1]) / (2 * dx))*((u[ny*i+j+1] - u[ny*i+j-1]) / (2 * dx)) -
                        2 * ((u[ny*(i+1)+j] - u[ny*(i-1)+j]) / (2 * dy) *
                             (v[ny*i+j+1] - v[ny*i+j-1]) / (2 * dx)) -
                            ((v[ny*(i+1)+j] - v[ny*(i-1)+j]) / (2 * dy))*((v[ny*(i+1)+j] - v[ny*(i-1)+j]) / (2 * dy))));

    }
}

__global__ void pressure_poisson_main(float *p, float *pn, float *b, float dx, float dy) {
    int i = blockIdx.y;
    int j = threadIdx.x + blockDim.x * blockIdx.x;
    if (i!=0 && j!=0 && i<nx-1 && j<ny-1) {
        p[ny*i+j] = (((pn[ny*i+j+1] + pn[ny*i+j-1]) * dy * dy +
                      (pn[ny*(i+1)+j] + pn[ny*(i-1)+j]) * dx * dx) /
                      (2 * (dx * dx + dy * dy)) -
                      (dx * dx * dy * dy) / (2 * (dx * dx + dy * dy)) *
                       b[ny*i+j]);
    }
}


__global__ void pressure_poisson_border_x(float *p) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i<nx) {
        p[ny*i+ny-1] = p[ny*i+ny-2];
        p[ny*i] = p[ny*i+1];
    }
}


__global__ void pressure_poisson_border_y(float *p) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j<ny) {
        p[ny*0+j] = p[ny*1+j];
        p[ny*(nx-1)+j] = 0;
    }
}


void pressure_poisson(float *p, float dx, float dy, float *b, int nit) {
    float *pn;
    cudaMallocManaged(&pn, size);

    for (int n=0; n<nit; n++) {
        cudaMemcpy(pn, p, size, cudaMemcpyDeviceToDevice);
        pressure_poisson_main<<<grid, N>>>(p, pn, b, dx, dy);
        cudaDeviceSynchronize();
    
        pressure_poisson_border_x<<<(nx+N-1)/N, N>>>(p);
        pressure_poisson_border_y<<<(ny+N-1)/N, N>>>(p);
        cudaDeviceSynchronize();
       }
    cudaFree(pn);
}


__global__ void cavity_flow_main(float *u, float *v, float *un, float *vn, float *p, float dt, float dx, float dy, float rho, float nu) {
    int i = blockIdx.y;
    int j = threadIdx.x + blockDim.x * blockIdx.x;
    if (i!=0 && j!=0 && i<nx-1 && j<ny-1) {
        u[ny*i+j] =(un[ny*i+j]-
                    un[ny*i+j] * dt / dx *
                   (un[ny*i+j] - un[ny*i+j-1]) -
                    vn[ny*i+j] * dt / dy *
                   (un[ny*i+j] - un[ny*(i-1)+j]) -
                    dt / (2 * rho * dx) * (p[ny*i+j+1] - p[ny*i+j-1]) +
                    nu * (dt / (dx * dx) *
                   (un[ny*i+j+1] - 2 * un[ny*i+j] + un[ny*i+j-1]) +
                    dt / (dy * dy) *
                   (un[ny*(i+1)+j] - 2 * un[ny*i+j] + un[ny*(i-1)+j])));
                
        v[ny*i+j] = (vn[ny*i+j] -
                 un[ny*i+j] * dt / dx *
                (vn[ny*i+j] - vn[ny*i+j-1]) -
                 vn[ny*i+j] * dt / dy *
                (vn[ny*i+j] - vn[ny*(i-1)+j]) -
                 dt / (2 * rho * dy) * (p[ny*(i+1)+j] - p[ny*(i-1)+j]) +
                 nu * (dt / (dx * dx) *
                (vn[ny*i+j+1] - 2 * vn[ny*i+j] + vn[ny*i+j-1]) +
                 dt / (dy * dy) *
                (vn[ny*(i+1)+j] - 2 * vn[ny*i+j] + vn[ny*(i-1)+j])));
    }
}


__global__ void cavity_flow_border_x(float *u, float *v) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i<nx) {
        u[ny*i+0] = 0;
        u[ny*i+ny-1] = 0;
        v[ny*i+0] = 0;
        v[ny*i+ny-1] = 0;
    }
}


__global__ void cavity_flow_border_y(float *u, float *v) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j<ny) {
        u[ny*0+j] = 0;
        u[ny*(nx-1)+j] = 1;
        v[ny*0+j] = 0;
        v[ny*(nx-1)+j] = 0;
    }
}


void cavity_flow(int nt, float *u, float *v, float dt, float dx, float dy, float *p, float rho, float nu, int nit) {
    float *un;
    float *vn;
    float *b;
    cudaMallocManaged(&un, size);
    cudaMallocManaged(&vn, size);
    cudaMallocManaged(&b, size);

    for (int n=0; n<nt; n++) {
        cudaMemcpy(un, u, size, cudaMemcpyDeviceToDevice);
        cudaMemcpy(vn, v, size, cudaMemcpyDeviceToDevice);

        build_up_b<<<grid, N>>>(b, rho, dt, u, v, dx, dy);
        cudaDeviceSynchronize();
        pressure_poisson(p, dx, dy, b, nit);

        cavity_flow_main<<<grid, N>>>(u, v, un, vn, p, dt, dx, dy, rho, nu);
        cudaDeviceSynchronize();
        cavity_flow_border_x<<<(nx+N-1)/N, N>>>(u, v);
        cavity_flow_border_y<<<(ny+N-1)/N, N>>>(u, v);
        cudaDeviceSynchronize();
    }
    cudaFree(un);
    cudaFree(vn);
    cudaFree(b);
}


int main() {
    const int nit = 50;
    const float dx = 2.0 / (nx - 1);
    const float dy = 2.0 / (ny - 1);

    float rho = 1.0;
    float nu = 0.1;
    float dt = 0.001;

    float *u;
    float *v;
    float *p;
    cudaMallocManaged(&u, size);
    cudaMallocManaged(&v, size);
    cudaMallocManaged(&p, size);

    auto tic = chrono::steady_clock::now();
    cavity_flow(100, u, v, dt, dx, dy, p, rho, nu, nit);
    auto toc = chrono::steady_clock::now();
    double time = chrono::duration<double>(toc-tic).count();
    cout << time << endl;
    save_matrix(u, v, p);

    cudaFree(u);
    cudaFree(v);
    cudaFree(p);
    return 0;
}
