/*
 * nozzle_choking.c
 *
 * Locate the choking position in an ideal, quasi-one-dimensional,
 * compressible convergent-divergent nozzle.
 *
 * The steady quasi-1D Euler equations are
 *
 *   d(rho*u*A)/dx = 0
 *   d((rho*u*u + p)*A)/dx = p*dA/dx
 *   d((rho*E + p)*u*A)/dx = 0
 *
 * For an adiabatic, isentropic nozzle these reduce to the critical-area
 * condition A = A* at Mach number M = 1. Therefore the choking position is
 * the minimum-area (throat) position, provided the imposed back pressure is
 * no greater than the critical pressure.
 *
 * Build:
 *   cc -O2 -Wall -Wextra -pedantic nozzle_choking.c -lm -o nozzle_choking
 *
 * Run with defaults:
 *   ./nozzle_choking
 *
 * Or specify:
 *   ./nozzle_choking gamma R p0 T0 p_back L L_convergent A_in A_throat A_exit
 *
 * Units: gamma [-], R [J/(kg K)], pressures [Pa], T0 [K], L and areas [m,m^2].
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

typedef struct {
    double gamma;
    double gas_constant;
    double p0;
    double T0;
    double p_back;
    double length;
    double convergent_length;
    double area_inlet;
    double area_throat;
    double area_exit;
} Nozzle;

static double area(const Nozzle *n, double x)
{
    const double Lc = n->convergent_length;
    const double L = n->length;

    if (x <= 0.0) return n->area_inlet;
    if (x >= L) return n->area_exit;

    if (x <= Lc) {
        /* Smooth cosine contraction: dA/dx is zero at inlet and throat. */
        const double s = x / Lc;
        return n->area_throat
             + 0.5 * (n->area_inlet - n->area_throat)
             * (1.0 + cos(PI * s));
    }

    /* Smooth cosine expansion: dA/dx is zero at throat and exit. */
    {
        const double s = (x - Lc) / (L - Lc);
        return n->area_throat
             + 0.5 * (n->area_exit - n->area_throat)
             * (1.0 - cos(PI * s));
    }
}

/* Golden-section minimization of A(x), useful if the geometry is changed. */
static double locate_minimum_area(const Nozzle *n, double left, double right)
{
    const double golden = 0.6180339887498948482;
    double x1 = right - golden * (right - left);
    double x2 = left + golden * (right - left);
    double f1 = area(n, x1);
    double f2 = area(n, x2);

    for (int iter = 0; iter < 200 && (right - left) > 1.0e-12 * n->length; ++iter) {
        if (f1 < f2) {
            right = x2;
            x2 = x1;
            f2 = f1;
            x1 = right - golden * (right - left);
            f1 = area(n, x1);
        } else {
            left = x1;
            x1 = x2;
            f1 = f2;
            x2 = left + golden * (right - left);
            f2 = area(n, x2);
        }
    }
    return 0.5 * (left + right);
}

static double critical_pressure(const Nozzle *n)
{
    const double g = n->gamma;
    return n->p0 * pow(2.0 / (g + 1.0), g / (g - 1.0));
}

static double critical_temperature(const Nozzle *n)
{
    return n->T0 * 2.0 / (n->gamma + 1.0);
}

static double choked_mass_flow(const Nozzle *n, double throat_area)
{
    const double g = n->gamma;
    const double R = n->gas_constant;
    const double factor = pow(2.0 / (g + 1.0), (g + 1.0) / (2.0 * (g - 1.0)));
    return throat_area * n->p0 / sqrt(R * n->T0) * sqrt(g) * factor;
}

static int read_arguments(int argc, char **argv, Nozzle *n)
{
    if (argc == 1) return 1;
    if (argc != 12) {
        fprintf(stderr,
                "Usage: %s gamma R p0 T0 p_back L L_convergent "
                "A_in A_throat A_exit\n", argv[0]);
        return 0;
    }

    double *values[] = {&n->gamma, &n->gas_constant, &n->p0, &n->T0,
                        &n->p_back, &n->length, &n->convergent_length,
                        &n->area_inlet, &n->area_throat, &n->area_exit};
    for (int i = 0; i < 10; ++i) {
        char *end = NULL;
        *values[i] = strtod(argv[i + 1], &end);
        if (end == argv[i + 1] || *end != '\0') {
            fprintf(stderr, "Invalid numeric argument: %s\n", argv[i + 1]);
            return 0;
        }
    }
    return 1;
}

static int valid_nozzle(const Nozzle *n)
{
    return n->gamma > 1.0 && n->gas_constant > 0.0 && n->p0 > 0.0
        && n->T0 > 0.0 && n->p_back >= 0.0 && n->length > 0.0
        && n->convergent_length > 0.0 && n->convergent_length < n->length
        && n->area_throat > 0.0 && n->area_inlet >= n->area_throat
        && n->area_exit >= n->area_throat;
}

int main(int argc, char **argv)
{
    Nozzle n = {
        .gamma = 1.4,
        .gas_constant = 287.05,
        .p0 = 101325.0,
        .T0 = 300.0,
        .p_back = 20000.0,
        .length = 1.0,
        .convergent_length = 0.4,
        .area_inlet = 0.020,
        .area_throat = 0.010,
        .area_exit = 0.016
    };

    if (!read_arguments(argc, argv, &n) || !valid_nozzle(&n)) {
        fprintf(stderr, "Invalid nozzle or thermodynamic input.\n");
        return EXIT_FAILURE;
    }

    const double x_star = locate_minimum_area(&n, 0.0, n.length);
    const double A_star = area(&n, x_star);
    const double p_star = critical_pressure(&n);
    const double T_star = critical_temperature(&n);
    const int choked = n.p_back <= p_star;

    printf("Quasi-1D Euler convergent-divergent nozzle\n");
    printf("-----------------------------------------\n");
    printf("Throat position x*       = %.9g m\n", x_star);
    printf("Minimum area A*          = %.9g m^2\n", A_star);
    printf("Critical pressure p*     = %.9g Pa\n", p_star);
    printf("Critical temperature T*  = %.9g K\n", T_star);
    printf("Choked mass flow rate    = %.9g kg/s\n", choked_mass_flow(&n, A_star));
    printf("Back pressure            = %.9g Pa\n", n.p_back);
    printf("Flow state               = %s\n", choked ? "CHOKED (M = 1 at x*)"
                                                   : "UNCHOKED (no sonic point)");

    if (!choked) {
        printf("To choke the nozzle, reduce p_back to %.9g Pa or below.\n", p_star);
    }
    if (fabs(A_star - n.area_throat) > 1.0e-8 * n.area_throat) {
        fprintf(stderr, "Warning: the numerical minimum differs from A_throat.\n");
    }
    return EXIT_SUCCESS;
}
