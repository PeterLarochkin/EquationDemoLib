#include <stdio.h>          /* printf */
#include <string.h>         /* strcspn */
#include <stdlib.h>         /* strtol */
#include <math.h>
#include <time.h>

double q(double x, double y) {
    return x + y;
}

double u(double x, double y) {
    return sqrt(4.0+x*y);
}
double F(double x, double y) {
    double u = sqrt(4+x*y);
    return 1/(4*u*u*u)*(x*x + y*y) + (x + y)*u;
}

double psi(double x, double y, double A1, double A2, double B1, double B2, double h1, double h2) {
    
    if (x == A2 && B1<y && y<B2) {
        // 1
        double u_ = sqrt(4.0+x*y);
        return y/(2*u_)+u_;
    } else if (x == A1 && B1<y && y<B2) {
        // 2
        return -y/4+2;
    } else if (y==B2 && A1<x && x<A2) {
        // 3
        double u_ = sqrt(4.0+x*y);
        return x/(2*u_)+u_;
    } else if (y==B1 && A1<x && x<A2) {
        // 4
        return -x/4+2;

    } else if (x==A1 && y==B1 ) {
        
        return (h1*(-x/4+2) + h2*(-y/4+2)) / (h1 + h2);
    } else if (x==A1 && y==B2 ) {
        
        double u = sqrt(4.0+x*y);
        return (h1*(x/(2*u)+u) + h2*(-y/4+2)) / (h1 + h2);
    } else if (x==A2 && y==B1 ) {
        double u = sqrt(4.0+x*y);
        return (h1*(-x/4+2) + h2*(y/(2*u)+u)) / (h1 + h2);
    } else if (x==A2 && y==B2 ) {
        double u = sqrt(4.0+x*y);
        return (h1*(x/(2*u)+u) + h2*(y/(2*u)+u)) / (h1 + h2);
    } else {
        printf("ERROR:(%.10f, %.10f)", x, y);
        return 0.0;
    } 
}

void applyA(double** whatApplyTo, double** whatWriteTo, int M, int N, double h1, double h2, double A1, double B1) {

    //with padding, inside "picture"
    for (int i = 1; i < M; ++i) {
        for (int j = 1; j < N; ++j) {
            // here is (7) equation works
            whatWriteTo[i][j] = whatApplyTo[i][j] * (2/(h1*h1) + 2/(h2*h2) + q(A1+ i*h1, B1+ j*h2)) + 
                                whatApplyTo[i-1][j] * (-1/(h1*h1)) +
                                whatApplyTo[i+1][j] * (-1/(h1*h1)) +
                                whatApplyTo[i][j-1] * (-1/(h2*h2)) +
                                whatApplyTo[i][j+1] * (-1/(h2*h2));
        }
    }
    for (int i=1; i<M; ++i) {
        // it's (10) equations
        // i=1,M-1
        // top applying
        whatWriteTo[i][N] = 2/(h2*h2) * (whatApplyTo[i][N] - whatApplyTo[i][N-1]) +
                            ( q(A1+ i*h1,B1+ N*h2) + 2/h2 ) * whatApplyTo[i][N] -
                            1/(h1*h1)*(whatApplyTo[i+1][N] - whatApplyTo[i][N] - whatApplyTo[i][N]+ whatApplyTo[i-1][N]); 
        // bottom applying
        whatWriteTo[i][0] = -2/(h2*h2) * (whatApplyTo[i][1] - whatApplyTo[i][0]) +
                            ( q(A1+ i*h1,B1+ 0*h2) + 2/h2 ) * whatApplyTo[i][0] -
                            1/(h1*h1)*(whatApplyTo[i+1][0] - whatApplyTo[i][0] - whatApplyTo[i][0]+ whatApplyTo[i-1][0]);
    }

    for (int j=1; j<N; ++j) {
        // it's (9) equations
        // j=1,N-1
        // right applying
        whatWriteTo[M][j] = 2/(h1*h1) * (whatApplyTo[M][j] - whatApplyTo[M-1][j]) + 
                            (q(A1+ M*h1,B1+ j*h2) + 2/h1) * whatApplyTo[M][j] - 
                            1/(h2*h2)*(whatApplyTo[M][j+1] - whatApplyTo[M][j] - whatApplyTo[M][j]+ whatApplyTo[M][j-1]);
        // left applying
        whatWriteTo[0][j] = -2/(h1*h1) * (whatApplyTo[1][j] - whatApplyTo[0][j]) + 
                            (q(A1+ 0*h1,B1+ j*h2) + 2/h1) * whatApplyTo[0][j] - 
                            1/(h2*h2)*(whatApplyTo[0][j+1] - whatApplyTo[0][j] - whatApplyTo[0][j]+ whatApplyTo[0][j-1]);
    }

    // remaining corner points
    // bottom left
    // it's (11) equation
    whatWriteTo[0][0] = -2/(h1*h1)*(whatApplyTo[1][0] - whatApplyTo[0][0]) - 
                         2/(h2*h2)*(whatApplyTo[0][1] - whatApplyTo[0][0]) +
                         (q(A1+ 0*h1,B1+ 0*h2) + 2/h1 + 2/h2) * whatApplyTo[0][0];

    // it's (12) equation
    // bottom right
    whatWriteTo[M][0] = 2/(h1*h1)*(whatApplyTo[M][0] - whatApplyTo[M-1][0]) - 
                        2/(h2*h2)*(whatApplyTo[M][1] - whatApplyTo[M][0]) +
                        (q(A1+ M*h1,B1+ 0*h2) + 2/h1 + 2/h2) * whatApplyTo[M][0];
    // it's (13) equation
    // top right
    whatWriteTo[M][N] = 2/(h1*h1)*(whatApplyTo[M][N] - whatApplyTo[M-1][N]) +
                        2/(h2*h2)*(whatApplyTo[M][N] - whatApplyTo[M][N-1]) +
                        (q(A1+ M*h1,B1+ N*h2) + 2/h1 + 2/h2) * whatApplyTo[M][N];
    // it's (14) equation
    // top left
    whatWriteTo[0][N] = -2/(h1*h1)*(whatApplyTo[1][N]- whatApplyTo[0][N])+
                         2/(h2*h2)*(whatApplyTo[0][N]- whatApplyTo[0][N-1])+
                         (q(A1+ 0*h1,B1+ N*h2) + 2/h1 + 2/h2) * whatApplyTo[0][N];
}

void getB(double** whatWriteTo, int M, int N, double h1, double h2, double A1, double A2, double B1, double B2) {
    //with padding, inside "picture"
    for (int i = 1; i < M; ++i) {
        for (int j = 1; j < N; ++j) {
            // here is (7) equation works
            whatWriteTo[i][j] = F(A1+ i*h1,B1+ j*h2);
        }
    }
    for (int i=1; i<M; ++i) {
        // it's (10) equations
        // i=1,M-1
        // top applying
        whatWriteTo[i][N] = psi(A1+ i*h1, B1+ N*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F(A1 + i*h1, B1 + N*h2);
        // bottom applying
        whatWriteTo[i][0] = psi(A1+ i*h1, B1+ 0*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F(A1 + i*h1, B1 + 0*h2);
    }

    for (int j=1; j<N; ++j) {
        // it's (9) equations
        // j=1,N-1
        // right applying
        whatWriteTo[M][j] = psi(A1+ M*h1, B1+ j*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F(A1 + M*h1, B1 + j*h2);
        // left applying
        whatWriteTo[0][j] = psi(A1+ 0*h1, B1+ j*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F(A1 + 0*h1, B1 + j*h2);
    }

    // remaining corner points
    // bottom left
    // it's (11) equation
    whatWriteTo[0][0] = psi(A1+ 0*h1, B1+ 0*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F(A1 + 0*h1, B1 + 0*h2);
    // it's (12) equation
    // bottom right
    whatWriteTo[M][0] = psi(A1+ M*h1, B1+ 0*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F(A1 + M*h1, B1 + 0*h2);
    // it's (13) equation
    // top right
    whatWriteTo[M][N] = psi(A1+ M*h1, B1+ N*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F(A1 + M*h1, B1 + N*h2);
    // it's (14) equation
    // top left
    whatWriteTo[0][N] = psi(A1+ 0*h1, B1+ N*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F(A1 + 0*h1, B1 + N*h2);
}

void minus(double** first, double** second, double** whatWriteTo, double M, double N) {
    for (size_t i = 0; i <= M; ++i) {
            for (size_t j = 0; j <= N; ++j) {
            whatWriteTo[i][j] = first[i][j] - second[i][j];
        }
    }
}

double ro(int index, int M_or_N) {
    if (index == 0 || index == M_or_N) {
        return 0.5;
    } else {
        return 1.0;
    }
}

double scalarProduct(double** first, double** second, double M, double N, double h1, double h2) {
    double sum = 0.0;
    for (size_t i = 0; i <= M; ++i) {
        for (size_t j = 0; j <= N; ++j) {
            sum = sum + h1*h2*ro(i, M)*ro(j, N)*first[i][j] * second[i][j];
        }
    }
    return sum;
}

void multiplyByNum(double** items, double num, double** whatWriteTo, double M, double N) {
    for (size_t i = 0; i <= M; ++i) {
            for (size_t j = 0; j <= N; ++j) {
            whatWriteTo[i][j] = items[i][j]*num;
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
      fprintf(stderr, "Put args!\n");
      return 1;
    }
    const size_t M = atoi(argv[1]);
    const size_t N = atoi(argv[2]);

    clock_t t;
    t = clock();
    
    
    double epsilon = 0.000001;
    // const size_t M = 160;
    // const size_t N = 160;
    
    double A1 = 0.0;
    double A2 = 4.0;
    double B1 = 0.0;
    double B2 = 3.0;
    double h1 = (A2 - A1)/M;
    double h2 = (B2 - B1)/N;
    

    double** omega              = (double**)malloc((M + 1) * sizeof(double*)); for (size_t i = 0; i <= M; ++i) omega[i] = (double*)malloc((N + 1) * sizeof(double));
    double** omega_next         = (double**)malloc((M + 1) * sizeof(double*)); for (size_t i = 0; i <= M; ++i) omega_next[i] = (double*)malloc((N + 1) * sizeof(double));
    double** B                  = (double**)malloc((M + 1) * sizeof(double*)); for (size_t i = 0; i <= M; ++i) B[i] = (double*)malloc((N + 1) * sizeof(double));
    double** A_omega            = (double**)malloc((M + 1) * sizeof(double*)); for (size_t i = 0; i <= M; ++i) A_omega[i] = (double*)malloc((N + 1) * sizeof(double));
    double** r                  = (double**)malloc((M + 1) * sizeof(double*)); for (size_t i = 0; i <= M; ++i) r[i] = (double*)malloc((N + 1) * sizeof(double));
    double** A_r                = (double**)malloc((M + 1) * sizeof(double*)); for (size_t i = 0; i <= M; ++i) A_r[i] = (double*)malloc((N + 1) * sizeof(double));
    double** tau_r              = (double**)malloc((M + 1) * sizeof(double*)); for (size_t i = 0; i <= M; ++i) tau_r[i] = (double*)malloc((N + 1) * sizeof(double));
    double** difference_omega   = (double**)malloc((M + 1) * sizeof(double*)); for (size_t i = 0; i <= M; ++i) difference_omega[i] = (double*)malloc((N + 1) * sizeof(double));

    double tau = 0.0;
    double sq_eps = epsilon * epsilon;
    double squared_difference = sq_eps;
    for (size_t i = 0; i <= M; ++i) {
        for (size_t j = 0; j <= N; ++j) {
            omega[i][j] = 0.0;
            omega_next[i][j] = 0.0;
            B[i][j] = 0.0;
            A_omega[i][j] = 0.0;
            r[i][j] = 0.0;
            A_r[i][j] = 0.0;
            tau_r[i][j] = 0.0;
            difference_omega[i][j] = 0.0;
        }
    }

    getB(B, M, N, h1, h2, A1, A2, B1, B2);

    int count = 0;
    while (squared_difference >= sq_eps)
    {
        for (size_t i = 0; i <= M; ++i) {
            for (size_t j = 0; j <= N; ++j) {
                omega[i][j] = omega_next[i][j];
            }
        }
        
        applyA(omega, A_omega, M, N, h1, h2, A1, B1);
        minus(A_omega, B, r, M, N);
        applyA(r, A_r, M, N, h1, h2, A1, B1);
        tau = scalarProduct(A_r, r, M, N, h1, h2) / scalarProduct(A_r, A_r, M, N, h1, h2);
        multiplyByNum(r, tau, tau_r, M, N);
        minus(omega, tau_r, omega_next, M, N);
        squared_difference = scalarProduct(tau_r, tau_r, M, N, h1, h2);
        if (count % 1000 == 0)
            printf("n:%d, diff:%.10f\n", count,sqrt(squared_difference));
        count++;
    }
    printf("total count: %d\n", count);
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    double max = 0.0;
    for (size_t i = 0; i <= M; ++i) {
        for (size_t j = 0; j <= N; ++j) { 
            double item = fabs(omega_next[i][j] - u(h1*i, h2*j));
            if (item > max) {
                max = item;
            }
        } 
    }
    printf("diff:%.10f\n", max);
    printf("Elapsed time: %.2f seconds\n", time_taken);
    // gcc -o imperative imperative.c && ./imperative
    return 0;
}