#include "raytracingx.h"

CCTK_REAL innerProduct(const CCTK_REAL* U, const CCTK_REAL* V, const Metric* m) { //return g_{\mu\nu} U^\nu V\nu
    printf(("g^tt: " + std::to_string(m->g_tt) + ")\n").c_str());
    return m->g_tt*U[0]*V[0] + m->beta_x*U[1]*V[0] + m->beta_y*U[2]*V[0] + m->beta_z*U[3]*V[0]
         + m->beta_x*U[0]*V[1] + m->g_xx*U[1]*V[1] + m->g_xy*U[2]*V[1] + m->g_xz*U[3]*V[1]
         + m->beta_y*U[0]*V[2] + m->g_xy*U[1]*V[2] + m->g_yy*U[2]*V[2] + m->g_yz*U[3]*V[2]
         + m->beta_z*U[0]*V[3] + m->g_xz*U[1]*V[3] + m->g_yz*U[2]*V[3] + m->g_zz*U[3]*V[3];
}

void generalizedCrossProduct(CCTK_REAL* X, const CCTK_REAL* U, const CCTK_REAL* V, const CCTK_REAL* W, const Metric* m) {//X_\rho = \vareps_{\lambda\mu\nu\rho} U^\lambda V^\mu W^\nu, then raising to get X^\rho
    CCTK_REAL temp[4];
    temp[0] = -U[3]*V[1]*W[2] +  U[2]*V[1]*W[3] +  U[3]*V[2]*W[1] + -U[1]*V[2]*W[3] + -U[2]*V[3]*W[1] +  U[1]*V[3]*W[2];
    temp[1] =  U[3]*V[0]*W[2] + -U[2]*V[0]*W[3] + -U[3]*V[2]*W[0] +  U[0]*V[2]*W[3] +  U[2]*V[3]*W[0] + -U[0]*V[3]*W[2];
    temp[2] = -U[3]*V[0]*W[1] +  U[1]*V[0]*W[3] +  U[3]*V[1]*W[0] + -U[0]*V[1]*W[3] + -U[1]*V[3]*W[0] +  U[0]*V[3]*W[1];
    temp[3] =  U[2]*V[0]*W[1] + -U[1]*V[0]*W[2] + -U[2]*V[1]*W[0] +  U[0]*V[1]*W[2] +  U[1]*V[2]*W[0] + -U[0]*V[2]*W[1];
    oneFormToVector(X, temp, &m);
}

void oneFormToVector(CCTK_REAL* X_vector, const CCTK_REAL* X_oneform, const Metric* m) { //X^\nu = g^{\mu\nu} X_\mu
    X_vector[0] = X_oneform[0]*m->g_upup_tt + X_oneform[1]*m->g_upup_tx + X_oneform[2]*m->g_upup_ty + X_oneform[3]*m->g_upup_tz;
    X_vector[1] = X_oneform[0]*m->g_upup_tx + X_oneform[1]*m->g_upup_xx + X_oneform[2]*m->g_upup_xy + X_oneform[3]*m->g_upup_xz;
    X_vector[2] = X_oneform[0]*m->g_upup_ty + X_oneform[1]*m->g_upup_xy + X_oneform[2]*m->g_upup_yy + X_oneform[3]*m->g_upup_yz;
    X_vector[3] = X_oneform[0]*m->g_upup_tz + X_oneform[1]*m->g_upup_xz + X_oneform[2]*m->g_upup_yz + X_oneform[3]*m->g_upup_zz;
}

void vectorToOneForm(CCTK_REAL* X_oneform, const CCTK_REAL* X_vector, const Metric* m) { //X_\nu = g_{\mu\nu} X^\mu
    X_oneform[0] = X_vector[0]*m->g_tt + X_vector[1]*m->beta_x + X_vector[2]*m->beta_y + X_vector[3]*m->beta_z;
    X_oneform[1] = X_vector[0]*m->beta_x + X_vector[1]*m->g_xx + X_vector[2]*m->g_xy + X_vector[3]*m->g_xz;
    X_oneform[2] = X_vector[0]*m->beta_y + X_vector[1]*m->g_xy + X_vector[2]*m->g_yy + X_vector[3]*m->g_yz;
    X_oneform[3] = X_vector[0]*m->beta_z + X_vector[1]*m->g_xz + X_vector[2]*m->g_yz + X_vector[3]*m->g_zz;
}

void projectUontoV(CCTK_REAL* X, const CCTK_REAL* U, const CCTK_REAL* V, const Metric* m) { //X^\mu = (U \cdot V) / (V \cdot V) V^\mu
    CCTK_REAL r = innerProduct(U, V, &m) / innerProduct(V, V, &m);
    X[0] = r*V[0];
    X[1] = r*V[1];
    X[2] = r*V[2];
    X[3] = r*V[3];
}

void normalize(CCTK_REAL* X_norm, const CCTK_REAL* X, const Metric* m) { //X_{norm}^\mu = X^\mu / (X \cdot X) 
    printf(("X: (" + std::to_string(X[0]) + ", " + std::to_string(X[1]) + ", " + std::to_string(X[2]) + ", " + std::to_string(X[3]) + ")\n").c_str());
    CCTK_REAL mag = sqrt(innerProduct(X, X, &m));
    printf(("X**2: (" + std::to_string(mag) + ")\n").c_str());
    X_norm[0] = X[0] / mag;
    X_norm[1] = X[1] / mag;
    X_norm[2] = X[2] / mag;
    X_norm[3] = X[3] / mag;
}


CCTK_REAL getTimeComponentOf4Velocity(const CCTK_REAL vx, const CCTK_REAL vy, const CCTK_REAL vz, const Metric* m) { //from knowing v^1, v^2, v^3 and having v \cdot v = -1, find v^0
    //g_{\mu\nu} v^\mu v^\nu = -1
    //g_{00} v^0 v^0 + 2v^0 (g_{01} v^1 + g_{02} v^2 + g_{03} v^3) + (g_{ij}) + 1 = 0
    //Ax^2 + Bx + C = 0
    //use quadratic formula

    CCTK_REAL A = m->g_tt;
    CCTK_REAL B = 2*(m->beta_x*vx + m->beta_y*vy + m->beta_z*vz);
    CCTK_REAL C = m->g_xx*vx*vx + 2*m->g_xy*vx*vy + 2*m->g_xz*vx*vz +
                                   m->g_yy*vy*vy + 2*m->g_yz*vy*vz +
                                                    m->g_zz*vz*vz + 1;
    CCTK_REAL v0 = (-B + sqrt(B*B - 4*A*C))/(2*A);
    CCTK_REAL v[4] = {v0, vx, vy, vz};
    if (innerProduct(v, v, &m) + 1 < 0.0000000001) {
        return v0;
    }
    v0 = (-B - sqrt(B*B - 4*A*C))/(2*A);
    v[0] = v0;
    if (innerProduct(v, v, &m) + 1 < 0.0000000001) {
        return v0;
    }
    
    CCTK_VERROR("RaytracingX: Problem encountered with calculating the time component of the camera's 4-velocity from the 3-velocity");
    return -1;
}