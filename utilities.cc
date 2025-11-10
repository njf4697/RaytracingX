#include "raytracingx.h"

CCTK_REAL innerProduct(const CCTK_REAL* U, const CCTK_REAL* V, const Metric m) { //return g_{\mu\nu} U^\nu V\nu
    return m.metric0pr*U[0]*V[0] + m.metric[1]*U[1]*V[0] + m.metric[2]*U[2]*V[0] + m.metric[3]*U[3]*V[0]
         + m.metric[1]*U[0]*V[1] + m.metric[4]*U[1]*V[1] + m.metric[5]*U[2]*V[1] + m.metric[6]*U[3]*V[1]
         + m.metric[2]*U[0]*V[2] + m.metric[5]*U[1]*V[2] + m.metric[7]*U[2]*V[2] + m.metric[8]*U[3]*V[2]
         + m.metric[3]*U[0]*V[3] + m.metric[6]*U[1]*V[3] + m.metric[8]*U[2]*V[3] + m.metric[9]*U[3]*V[3];
}

void generalizedCrossProduct(CCTK_REAL* X, const CCTK_REAL* U, const CCTK_REAL* V, const CCTK_REAL* W, const Metric m) {//X_\rho = \vareps_{\lambda\mu\nu\rho} U^\lambda V^\mu W^\nu, then raising to get X^\rho
    CCTK_REAL* temp;
    temp[0] = -U[3]*V[1]*W[2] +  U[2]*V[1]*W[3] +  U[3]*V[2]*W[1] + -U[1]*V[2]*W[3] + -U[2]*V[3]*W[1] +  U[1]*V[3]*W[2];
    temp[1] =  U[3]*V[0]*W[2] + -U[2]*V[0]*W[3] + -U[3]*V[2]*W[0] +  U[0]*V[2]*W[3] +  U[2]*V[3]*W[0] + -U[0]*V[3]*W[2];
    temp[2] = -U[3]*V[0]*W[1] +  U[1]*V[0]*W[3] +  U[3]*V[1]*W[0] + -U[0]*V[1]*W[3] + -U[1]*V[3]*W[0] +  U[0]*V[3]*W[1];
    temp[3] =  U[2]*V[0]*W[1] + -U[1]*V[0]*W[2] + -U[2]*V[1]*W[0] +  U[0]*V[1]*W[2] +  U[1]*V[2]*W[0] + -U[0]*V[2]*W[1];
    oneFormToVector(X, temp, m);
}

void oneFormToVector(CCTK_REAL* X_vector, const CCTK_REAL* X_oneform, const Metric m) { //X^\nu = g^{\mu\nu} X_\mu
    X_vector[0] = X_oneform[0]*m.metric_inv[0] + X_oneform[1]*m.metric_inv[1] + X_oneform[2]*m.metric_inv[2] + X_oneform[3]*m.metric_inv[3];
    X_vector[1] = X_oneform[0]*m.metric_inv[1] + X_oneform[1]*m.metric_inv[4] + X_oneform[2]*m.metric_inv[5] + X_oneform[3]*m.metric_inv[6];
    X_vector[2] = X_oneform[0]*m.metric_inv[2] + X_oneform[1]*m.metric_inv[5] + X_oneform[2]*m.metric_inv[7] + X_oneform[3]*m.metric_inv[8];
    X_vector[3] = X_oneform[0]*m.metric_inv[3] + X_oneform[1]*m.metric_inv[6] + X_oneform[2]*m.metric_inv[8] + X_oneform[3]*m.metric_inv[9];
}

void vectorToOneForm(CCTK_REAL* X_oneform, const CCTK_REAL* X_vector, const Metric m) { //X_\nu = g_{\mu\nu} X^\mu
    X_oneform[0] = X_vector[0]*m.metric0pr + X_vector[1]*m.metric[1] + X_vector[2]*m.metric[2] + X_vector[3]*m.metric[3];
    X_oneform[1] = X_vector[0]*m.metric[1] + X_vector[1]*m.metric[4] + X_vector[2]*m.metric[5] + X_vector[3]*m.metric[6];
    X_oneform[2] = X_vector[0]*m.metric[2] + X_vector[1]*m.metric[5] + X_vector[2]*m.metric[7] + X_vector[3]*m.metric[8];
    X_oneform[3] = X_vector[0]*m.metric[3] + X_vector[1]*m.metric[6] + X_vector[2]*m.metric[8] + X_vector[3]*m.metric[9];
}

void projectUontoV(CCTK_REAL* X, const CCTK_REAL* U, const CCTK_REAL* V, const Metric m) { //X^\mu = (U \cdot V) / (V \cdot V) V^\mu
    CCTK_REAL* r = innerProduct(U, V, m) / innerProduct(V, V, m);
    X[0] = r*V[0];
    X[1] = r*V[1];
    X[2] = r*V[2];
    X[3] = r*v[3];
}

void normalize(CCTK_REAL* X_norm, const CCTK_REAL* X, const Metric m) { //X_{norm}^\mu = X^\mu / (X \cdot X) 
    CCTK_REAL mag = innerProduct(X, X, m);
    X_norm[0] = X[0] / mag;
    X_norm[1] = X[1] / mag;
    X_norm[2] = X[2] / mag;
    X_norm[3] = X[3] / mag;
}


CCTK_REAL getTimeComponentOf4Velocity(const CCTK_REAL vx, const CCTK_REAL vy, const CCTK_REAL vz, const Metric m) { //from knowing v^1, v^2, v^3 and having v \cdot v = -1, find v^0
    //g_{\mu\nu} v^\mu v^\nu = -1
    //g_{00} v^0 v^0 + 2v^0 (g_{01} v^1 + g_{02} v^2 + g_{03} v^3) + (g_{ij}) + 1 = 0
    //Ax^2 + Bx + C = 0
    //use quadratic formula

    CCTK_REAL A = m.metric0pr;
    CCTK_REAL B = 2*(m.metric[1]*vx + m.metric[2]*vy + m.metric[3]*vz);
    CCTK_REAL C = m.metric[4]*vx*vx + 2*m.metric[5]*vx*vy + 2*m.metric[6]*vx*vz +
                                        m.metric[7]*vy*vy + 2*m.metric[8]*vy*vz +
                                                              m.metric[9]*vz*vz + 1;
    CCTK_REAL v0 = (-B + sqrt(B*B - 4*A*C))/(2*A)
    CCTK_REAL* v[4] = {v0, vx, vy, vz};
    assert(innerProduct(v, m) + 1 < 0.00001);
    
    return v0;
}