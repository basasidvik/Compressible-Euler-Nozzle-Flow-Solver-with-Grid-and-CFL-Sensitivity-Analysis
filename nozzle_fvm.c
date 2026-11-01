/*
 * Time-marching finite-volume solver for the quasi-one-dimensional Euler
 * equations in a convergent-divergent nozzle.
 *
 * Build:
 *   cc -std=c11 -O2 -Wall -Wextra -pedantic nozzle_fvm.c -lm -o nozzle_fvm
 * Run:
 *   ./nozzle_fvm cases/grid_200.cfg
 *
 * The solver writes <output_prefix>.csv and <output_prefix>_residual.csv.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define MAX_LINE 256
#define NVAR 3

typedef struct {
    int nx, max_steps;
    double cfl, residual_tol;
    double gamma, R, p0, T0, p_back;
    double length, throat_x, A_in, A_throat, A_exit;
    char output_prefix[160];
} Config;

typedef struct {
    double rho, u, p;
} Primitive;

static double nozzle_area(const Config *c, double x)
{
    if (x <= c->throat_x) {
        double s = x / c->throat_x;
        return c->A_throat + (c->A_in - c->A_throat) * s * s;
    }
    {
        double s = (x - c->throat_x) / (c->length - c->throat_x);
        return c->A_throat + (c->A_exit - c->A_throat) * s * s;
    }
}

static double area_slope(const Config *c, double x)
{
    if (x < c->throat_x)
        return 2.0 * (c->A_in - c->A_throat) * x
             / (c->throat_x * c->throat_x);
    return 2.0 * (c->A_exit - c->A_throat) * (x - c->throat_x)
         / ((c->length - c->throat_x) * (c->length - c->throat_x));
}

static void primitive_to_U(const Primitive *q, double A, double U[NVAR],
                           double gamma)
{
    double E = q->p / (gamma - 1.0) + 0.5 * q->rho * q->u * q->u;
    U[0] = A * q->rho;
    U[1] = A * q->rho * q->u;
    U[2] = A * E;
}

static Primitive U_to_primitive(const double U[NVAR], double A, double gamma)
{
    Primitive q;
    q.rho = fmax(U[0] / A, 1.0e-10);
    q.u = U[1] / (A * q.rho);
    q.p = fmax((gamma - 1.0) * (U[2] / A - 0.5 * q.rho * q.u * q.u), 1.0e-6);
    return q;
}

static void flux(const Primitive *q, double A, double gamma, double F[NVAR])
{
    double E = q->p / (gamma - 1.0) + 0.5 * q->rho * q->u * q->u;
    F[0] = A * q->rho * q->u;
    F[1] = A * (q->rho * q->u * q->u + q->p);
    F[2] = A * q->u * (E + q->p);
}

static void rusanov_flux(const Primitive *L, const Primitive *R, double A,
                         double gamma, double F[NVAR])
{
    double FL[NVAR], FR[NVAR];
    double aL = sqrt(gamma * L->p / L->rho);
    double aR = sqrt(gamma * R->p / R->rho);
    double s = fmax(fabs(L->u) + aL, fabs(R->u) + aR);
    double UL[NVAR], UR[NVAR];

    primitive_to_U(L, A, UL, gamma);
    primitive_to_U(R, A, UR, gamma);
    flux(L, A, gamma, FL);
    flux(R, A, gamma, FR);
    for (int k = 0; k < NVAR; ++k)
        F[k] = 0.5 * (FL[k] + FR[k]) - 0.5 * s * (UR[k] - UL[k]);
}

static int read_config(const char *filename, Config *c)
{
    FILE *fp = fopen(filename, "r");
    char line[MAX_LINE], key[80], value[160];
    if (!fp) return 0;

    while (fgets(line, sizeof line, fp)) {
        if (line[0] == '#' || sscanf(line, "%79s %159s", key, value) != 2)
            continue;
        if (!strcmp(key, "nx")) c->nx = atoi(value);
        else if (!strcmp(key, "cfl")) c->cfl = atof(value);
        else if (!strcmp(key, "max_steps")) c->max_steps = atoi(value);
        else if (!strcmp(key, "residual_tol")) c->residual_tol = atof(value);
        else if (!strcmp(key, "gamma")) c->gamma = atof(value);
        else if (!strcmp(key, "R")) c->R = atof(value);
        else if (!strcmp(key, "p0")) c->p0 = atof(value);
        else if (!strcmp(key, "T0")) c->T0 = atof(value);
        else if (!strcmp(key, "p_back")) c->p_back = atof(value);
        else if (!strcmp(key, "length")) c->length = atof(value);
        else if (!strcmp(key, "throat_x")) c->throat_x = atof(value);
        else if (!strcmp(key, "A_in")) c->A_in = atof(value);
        else if (!strcmp(key, "A_throat")) c->A_throat = atof(value);
        else if (!strcmp(key, "A_exit")) c->A_exit = atof(value);
        else if (!strcmp(key, "output_prefix"))
            snprintf(c->output_prefix, sizeof c->output_prefix, "%s", value);
    }
    fclose(fp);
    return c->nx > 2 && c->cfl > 0.0 && c->cfl <= 1.0
        && c->gamma > 1.0 && c->R > 0.0 && c->p0 > 0.0 && c->T0 > 0.0
        && c->p_back > 0.0 && c->length > 0.0 && c->throat_x > 0.0
        && c->throat_x < c->length && c->A_throat > 0.0
        && c->A_in >= c->A_throat && c->A_exit >= c->A_throat
        && c->output_prefix[0] != '\0';
}

static void write_solution(const Config *c, const double *U, const double *A,
                           double dx)
{
    char filename[220];
    snprintf(filename, sizeof filename, "%s.csv", c->output_prefix);
    FILE *fp = fopen(filename, "w");
    if (!fp) return;
    fprintf(fp, "x,area,density,velocity,pressure,mach\n");
    for (int i = 0; i < c->nx; ++i) {
        double x = (i + 0.5) * dx;
        Primitive q = U_to_primitive(&U[NVAR * i], A[i], c->gamma);
        double mach = fabs(q.u) / sqrt(c->gamma * q.p / q.rho);
        fprintf(fp, "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                x, A[i], q.rho, q.u, q.p, mach);
    }
    fclose(fp);
}

int main(int argc, char **argv)
{
    Config c = {0};
    if (argc != 2 || !read_config(argv[1], &c)) {
        fprintf(stderr, "Usage: %s cases/file.cfg\n", argv[0]);
        return EXIT_FAILURE;
    }

    const double dx = c.length / c.nx;
    double *A = calloc((size_t)c.nx, sizeof *A);
    double *Af = calloc((size_t)c.nx + 1, sizeof *Af);
    double *U = calloc((size_t)c.nx * NVAR, sizeof *U);
    double *Un = calloc((size_t)c.nx * NVAR, sizeof *Un);
    double *F = calloc((size_t)c.nx + 1, sizeof *F * NVAR);
    double *residuals = calloc((size_t)c.max_steps, sizeof *residuals);
    if (!A || !Af || !U || !Un || !F || !residuals) return EXIT_FAILURE;

    Primitive initial = {c.p0 / (c.R * c.T0), 0.0, c.p0};
    for (int i = 0; i < c.nx; ++i) {
        A[i] = nozzle_area(&c, (i + 0.5) * dx);
        primitive_to_U(&initial, A[i], &U[NVAR * i], c.gamma);
    }
    for (int i = 0; i <= c.nx; ++i) Af[i] = nozzle_area(&c, i * dx);

    char residual_file[220];
    snprintf(residual_file, sizeof residual_file, "%s_residual.csv", c.output_prefix);
    FILE *rf = fopen(residual_file, "w");
    if (rf) fprintf(rf, "step,residual,dt\n");

    int final_step = c.max_steps;
    for (int step = 0; step < c.max_steps; ++step) {
        Primitive *q = calloc((size_t)c.nx, sizeof *q);
        if (!q) break;
        double dt = 1.0e30;
        for (int i = 0; i < c.nx; ++i) {
            q[i] = U_to_primitive(&U[NVAR * i], A[i], c.gamma);
            double wave = fabs(q[i].u) + sqrt(c.gamma * q[i].p / q[i].rho);
            dt = fmin(dt, c.cfl * dx / wave);
        }

        Primitive inlet = {c.p0 / (c.R * c.T0), 0.0, c.p0};
        Primitive outlet = q[c.nx - 1];
        outlet.p = c.p_back;
        rusanov_flux(&inlet, &q[0], Af[0], c.gamma, &F[0]);
        for (int i = 1; i < c.nx; ++i)
            rusanov_flux(&q[i - 1], &q[i], Af[i], c.gamma, &F[NVAR * i]);
        rusanov_flux(&q[c.nx - 1], &outlet, Af[c.nx], c.gamma,
                     &F[NVAR * c.nx]);

        double residual = 0.0;
        for (int i = 0; i < c.nx; ++i) {
            double source[NVAR] = {0.0, q[i].p * area_slope(&c, (i + 0.5) * dx), 0.0};
            for (int k = 0; k < NVAR; ++k) {
                Un[NVAR * i + k] = U[NVAR * i + k]
                    - dt / dx * (F[NVAR * (i + 1) + k] - F[NVAR * i + k])
                    + dt * source[k];
                residual += fabs(Un[NVAR * i + k] - U[NVAR * i + k]);
            }
        }
        residual /= (double)c.nx;
        residuals[step] = residual;
        if (rf) fprintf(rf, "%d,%.12e,%.12e\n", step + 1, residual, dt);
        memcpy(U, Un, (size_t)c.nx * NVAR * sizeof *U);
        free(q);
        if (!isfinite(residual)) {
            final_step = step + 1;
            break;
        }
        if (residual < c.residual_tol) {
            final_step = step + 1;
            break;
        }
    }
    if (rf) fclose(rf);
    write_solution(&c, U, A, dx);
    printf("case=%s steps=%d final_residual=%.6e%s\n", c.output_prefix,
           final_step, residuals[final_step - 1],
           residuals[final_step - 1] < c.residual_tol ? " (converged)"
                                                      : " (not converged)");
    free(A); free(Af); free(U); free(Un); free(F); free(residuals);
    return EXIT_SUCCESS;
}
