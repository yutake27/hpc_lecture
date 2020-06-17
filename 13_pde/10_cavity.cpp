#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
typedef vector<vector<float> > matrix;


void copy_matrix(matrix &from, matrix &to, int nx, int ny) {
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            to[i][j] = from[i][j];
        }
    }
}

void print_matrix(matrix &a, int nx, int ny) {
    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++) {
            cout << a[i][j] <<' ';
        }
        cout << endl;
    }
}

void save_matrix(matrix &u, matrix &v, matrix &p, int nx, int ny) {
    ofstream u_file("u.txt");
    ofstream v_file("v.txt");
    ofstream p_file("p.txt");
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            u_file << u[i][j] << ' ';
            v_file << v[i][j] << ' ';
            p_file << p[i][j] << ' ';
        }
        u_file << endl;
        v_file << endl;
        p_file << endl;
    }
    u_file.close();
    v_file.close();
    p_file.close();
}

void build_up_b(matrix &b, float rho, float dt, matrix &u, matrix &v, float dx, float dy, int nx, int ny) {
    matrix bn(nx, vector<float>(ny));
    for (int i=1; i<nx-1; i++){
        for (int j=1; j<ny-1; j++){
            b[i][j] = (rho * (1 / dt *
                      ((u[i][j+1] - u[i][j-1]) /
                       (2 * dx) + (v[i+1][j] - v[i-1][j]) / (2 * dy)) -
                      ((u[i][j+1] - u[i][j-1]) / (2 * dx))*((u[i][j+1] - u[i][j-1]) / (2 * dx)) -
                        2 * ((u[i+1][j] - u[i-1][j]) / (2 * dy) *
                             (v[i][j+1] - v[i][j-1]) / (2 * dx)) -
                            ((v[i+1][j] - v[i-1][j]) / (2 * dy))*((v[i+1][j] - v[i-1][j]) / (2 * dy))));
        }
    }
}


void pressure_poisson(matrix &p, float dx, float dy, matrix &b, int nx, int ny, int nit) {
    matrix pn(nx, vector<float>(ny));

    for (int n=0; n<nit; n++) {
        copy_matrix(p, pn, nx, ny);
        for (int i=1; i<nx-1; i++) {
            for (int j=1; j<ny-1; j++) {
                p[i][j] = (((pn[i][j+1] + pn[i][j-1]) * dy * dy +
                            (pn[i+1][j] + pn[i-1][j]) * dx * dx) /
                            (2 * (dx * dx + dy * dy)) -
                            (dx * dx * dy * dy) / (2 * (dx * dx + dy * dy)) *
                            b[i][j]);
            }
        }

        for (int i=0; i<nx; i++) {
            p[i][ny-1] = p[i][ny-2];
            p[i][0] = p[i][1];
        }

        for (int j=0; j< ny; j++) { 
            p[0][j] = p[1][j];
            p[nx-1][j] = 0;
        }
    }
}


void cavity_flow(int nt, matrix &u, matrix &v, float dt, float dx, float dy, matrix &p, float rho, float nu, int nx, int ny, int nit) {
    matrix un(nx, vector<float>(ny));
    matrix vn(nx, vector<float>(ny));
    matrix b(nx, vector<float>(ny));

    for (int n=0; n<nt; n++) {
        copy_matrix(u, un, nx, ny);
        copy_matrix(v, vn, nx, ny);

        build_up_b(b, rho, dt, u, v, dx, dy, nx, ny);
        pressure_poisson(p, dx, dy, b, nx, ny, nit);

        for (int i=1; i<nx-1; i++) {
            for (int j=1; j<ny-1; j++) {
                u[i][j] =(un[i][j]-
                          un[i][j] * dt / dx *
                         (un[i][j] - un[i][j-1]) -
                          vn[i][j] * dt / dy *
                         (un[i][j] - un[i-1][j]) -
                          dt / (2 * rho * dx) * (p[i][j+1] - p[i][j-1]) +
                          nu * (dt / (dx * dx) *
                         (un[i][j+1] - 2 * un[i][j] + un[i][j-1]) +
                          dt / (dy * dy) *
                         (un[i+1][j] - 2 * un[i][j] + un[i-1][j])));
                
                v[i][j] = (vn[i][j] -
                           un[i][j] * dt / dx *
                          (vn[i][j] - vn[i][j-1]) -
                           vn[i][j] * dt / dy *
                          (vn[i][j] - vn[i-1][j]) -
                           dt / (2 * rho * dy) * (p[i+1][j] - p[i-1][j]) +
                           nu * (dt / (dx * dx) *
                          (vn[i][j+1] - 2 * vn[i][j] + vn[i][j-1]) +
                           dt / (dy * dy) *
                          (vn[i+1][j] - 2 * vn[i][j] + vn[i-1][j])));
            }
        }

        for (int i=0; i<nx; i++) {
            u[i][0] = 0;
            u[i][ny-1] = 0;
            v[i][0] = 0;
            v[i][ny-1] = 0;
        }

        for (int j=0; j<ny; j++) {
            u[0][j] = 0;
            u[nx-1][j] = 1;
            v[0][j] = 0;
            v[nx-1][j] = 0;
        }
    }
}


int main() {
    const int nx = 41;
    const int ny = 41;
    const int nt = 500;
    const int nit = 50;
    const int c = 1;
    const float dx = 2.0 / (nx - 1);
    const float dy = 2.0 / (ny - 1);

    float rho = 1.0;
    float nu = 0.1;
    float dt = 0.001;

    matrix u(nx, vector<float>(ny));
    matrix v(nx, vector<float>(ny));
    matrix p(nx, vector<float>(ny));
    matrix b(nx, vector<float>(ny));
    cavity_flow(100, u, v, dt, dx, dy, p, rho, nu, nx, ny, nit);
    save_matrix(u, v, p, nx, ny);
    return 0;
}
