#include <cstdio>
#include <cmath>


void umat(double stress[], double statev[], double ddsdde[][6], double sse,
    double& spd, double scd, double& rpl, double ddsddt[],
    double drplde[], double drpldt, double stran[], double dstran[],
    double time[], double dtime, double temp2[], double dtemp,
    double predef[], double dpred[], const char* cmname, int ndi,
    int nshr, int ntens, int nstatv, double props[], int nprops,
    double coords[], double drot[][3], double pnewdt, double celent,
    double dfgrd0[][3], double dfgrd1[][3], int noel, int npt,
    int layer, int kspt, int jstep, int kinc) {
    const double toler = 1e-6;
    const int newton = 20;

    // Initialization
    int i, j;
    double E, xnu, Sy, xn, eqplas, Sh, Smises, Sf, deqpl, Et, rhs, effg, efflam, effhrd;
    double eplas[6];
    double eelas[6];
    int kewton;



    E = props[0];   // Young's modulus
    xnu = props[1]; // Poisson's ratio
    Sy = props[2];  // Yield stress
    xn = props[3];  // Strain hardening exponent

    // Call rotsig function
    eqplas = statev[2 * ntens];
    double olds[6], oldpl[6], flow[6];

    for (i = 0; i < ntens; ++i) {
        olds[i] = stress[i];
        oldpl[i] = eplas[i];
    }

    // Build elastic stiffness matrix
    double eg = E / (1.0 + xnu) / 2.0;
    double elam = (E / (1.0 - 2.0 * xnu) - 2.0 * eg) / 3.0;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ddsdde[j][i] = elam;
        }
        ddsdde[i][i] = 2.0 * eg + elam;
    }
    for (i = 3; i < ntens; ++i) {
        ddsdde[i][i] = eg;
    }

    // Calculate predictor stress and elastic strain
    for (i = 0; i < ntens; ++i) {
        stress[i] += ddsdde[i][0] * dstran[0] + ddsdde[i][1] * dstran[1] + ddsdde[i][2] * dstran[2];
        eelas[i] += dstran[i];
    }

    // Calculate equivalent von Mises stress
    Smises = pow(stress[0] - stress[1], 2) + pow(stress[1] - stress[2], 2) + pow(stress[2] - stress[0], 2);
    for (i = 3; i < ntens; ++i) {
        Smises += 6.0 * pow(stress[i], 2);
    }
    Smises = sqrt(Smises / 2.0);

    // Get yield stress from the specified hardening curve
    Sf = Sy * pow(1.0 + E * eqplas / Sy, xn);

    // Determine if active yielding
    if (Smises > (1.0 + toler) * Sf) {
        // Calculate the flow direction
        Sh = (stress[0] + stress[1] + stress[2]) / 3.0;
        for (i = 0; i < 3; ++i) {
            flow[i] = (stress[i] - Sh) / Smises;
        }
        for (i = 3; i < ntens; ++i) {
            flow[i] = stress[i] / Smises;
        }

        // Solve for Smises and deqpl using Newton's method
        deqpl = 0.0;
        Et = E * xn * pow(1.0 + E * eqplas / Sy, xn - 1);
        for (kewton = 1; kewton < newton; ++kewton) {
            rhs = Smises - (3.0 * eg) * deqpl - Sf;
            deqpl += rhs / (3.0 * eg + Et);
            Sf = Sy * pow(1.0 + E * (eqplas + deqpl) / Sy, xn);
            Et = E * xn * pow(1.0 + E * (eqplas + deqpl) / Sy, xn - 1);
            if (std::abs(rhs) < toler * Sy) break;
        }
        if (kewton == newton)  std::printf("WARNING: plasticity loop failed\n");

        // Update stresses and strains
        for (i = 0; i < 3; ++i) {
            stress[i] = flow[i] * Sf + Sh;
            eplas[i] += 3.0 / 2.0 * flow[i] * deqpl;
            eelas[i] -= 3.0 / 2.0 * flow[i] * deqpl;
        }
        for (i = 3; i < ntens; ++i) {
            stress[i] = flow[i] * Sf;
            eplas[i] += 3.0 * flow[i] * deqpl;
            eelas[i] -= 3.0 * flow[i] * deqpl;
        }
        eqplas += deqpl;

        // Calculate the plastic strain energy density
        for (i = 0; i < ntens; ++i) {
            spd += (stress[i] + olds[i]) * (eplas[i] - oldpl[i]) / 2.0;
        }

        // Formulate the jacobian (material tangent)
        effg = eg * Sf / Smises;
        efflam = (E / (1.0 - 2.0 * xnu) - 2.0 * effg) / 3.0;
        effhrd = 3.0 * eg * Et / (3.0 * eg + Et) - 3.0 * effg;
        for (i = 0; i < 3; ++i)
        {
            for (j = 0; j < 3; ++j) {
                ddsdde[j][i] = efflam;
            }
            ddsdde[i][i] = 2.0 * effg + efflam;
        }
        for (i = 3; i < ntens; ++i) {
            ddsdde[i][i] = effg;
        }
        for (i = 0; i < ntens; ++i) {
            for (j = 0; j < ntens; ++j) {
                ddsdde[j][i] += effhrd * flow[j] * flow[i];
            }
        }
    }

    // Store strains in state variable array
    for (i = 0; i < ntens; ++i) {
        statev[i] = eelas[i];
        statev[ntens + i] = eplas[i];
    }
    statev[2 * ntens] = eqplas;
}
