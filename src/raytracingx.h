#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <string.h>
#include "util_Table.h"

#ifndef RAYTRACINGX
#define RAYTRACINGX

struct GeodesicInitialConditions { //struct to contain information for the initial conditions for a geodesic
    CCTK_REAL initPos[3]; //initial position of geodesic (X^i), see https://arxiv.org/pdf/1410.7775
    CCTK_REAL initVel[3]; //initial position of geodesic (\Pi_i), see https://arxiv.org/pdf/1410.7775
};

struct Metric { //struct that contains information about the metric interpolated at a point
    CCTK_REAL alpha;
    CCTK_REAL beta_x, beta_y, beta_z; //g_{tx}, g_{ty}, g_{tz}
    CCTK_REAL g_xx, g_xy, g_xz, g_yy, g_yz, g_zz; //g_{ii} <-> \gamma_{ii}
    CCTK_REAL g_tt; //\sqrt{\alpha^2 + \beta_i\beta^i}
    CCTK_REAL beta_xup, beta_yup, beta_zup; //beta^i = \beta_j\gamma_{ij}
    CCTK_REAL g_upup_tt, g_upup_tx, g_upup_ty, g_upup_tz, g_upup_xx, g_upup_xy, g_upup_xz, g_upup_yy, g_upup_yz, g_upup_yz, g_upup_zz; //g^{\mu\nu}

    std::string to_string() {
        return "alpha: " + std::to_string(alpha) + "\n";
               "beta: (" + std::to_string(beta_xup) + ", " + std::to_string(beta_yup) + ", " + std::to_string(beta_up[2]) + ")";
               "metric: " + std::to_string(g_tt) + ", " + std::to_string(beta_x) + ", " + std::to_string(beta_y) + ", " + std::to_string(beta_z) + "\n";
               "        " + std::to_string(beta_x) + ", " + std::to_string(g_xx) + ", " + std::to_string(g_xy) + ", " + std::to_string(g_xz) + "\n";
               "        " + std::to_string(beta_y) + ", " + std::to_string(g_xy) + ", " + std::to_string(g_yy) + ", " + std::to_string(g_yz) + "\n";
               "        " + std::to_string(beta_z) + ", " + std::to_string(g_xz) + ", " + std::to_string(g_yz) + ", " + std::to_string(g_zz);
    }

    void fillInverseMetric() {
        //find inverse via normal matrix inverse (adjoint method)
        CCTK_REAL inv_det_g = 1/(g_xx*g_yy*g_zz+2.*g_xy*g_xz*g_yz-g_xz*g_xz*g_yy-g_yz*g_yz*g_xx-g_xy*g_xy*g_zz);
        CCTK_REAL gamma_upup_xx = (g_yy*g_zz-g_yz*g_yz)*inv_det_g;
        CCTK_REAL gamma_upup_xy = (g_xz*g_yz-g_xy*g_zz)*inv_det_g;
        CCTK_REAL gamma_upup_xz = (g_xy*g_yz-g_yy*g_xz)*inv_det_g;
        CCTK_REAL gamma_upup_yy = (g_xx*g_zz-g_xz*g_xz)*inv_det_g;
        CCTK_REAL gamma_upup_yz = (g_xz*g_xy-g_xx*g_yz)*inv_det_g;
        CCTK_REAL gamma_upup_zz = (g_xx*g_yy-g_xy*g_xy)*inv_det_g;

        g_upup_tt = 1/(pow(alpha,2));
        g_upup_tx = g_upup_tt*beta_xup;
        g_upup_ty = g_upup_tt*beta_yup;
        g_upup_tz = g_upup_tt*beta_up[2];
        g_upup_xx = gamma_upup_xx-g_upup_tx*beta_xup;
        g_upup_xy = gamma_upup_xy-g_upup_tx*beta_yup;
        g_upup_xz = gamma_upup_xz-g_upup_tx*beta_up[2];
        g_upup_yy = gamma_upup_yy-g_upup_ty*beta_yup;
        g_upup_yz = gamma_upup_yz-g_upup_ty*beta_up[2];
        g_upup_zz = gamma_upup_zz-g_upup_tz*beta_up[2];
    }
};

//scheduled
extern "C" void raytraceImage(CCTK_ARGUMENTS); //raytraceImage.cc

//nonscheduled
void gramSchmidtProcess(CCTK_ARGUMENTS, CCTK_REAL* e0, CCTK_REAL* e1, CCTK_REAL* e2, CCTK_REAL* e3, const Metric metric); //raytraceImage.cc
void createGeodesicInitialConditions(CCTK_ARGUMENTS, GeodesicInitialConditions* geodesicArr); //raytraceImage.cc
void interpolateMetricAtPoint(CCTK_ARGUMENTS, const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z, Metric metric_at_point); //interpolateMetric.cc
CCTK_REAL innerProduct(const CCTK_REAL* U, const CCTK_REAL* V, const Metric m); //utilities.cc
void generalizedCrossProduct(CCTK_REAL* X, const CCTK_REAL* U, const CCTK_REAL* V, const CCTK_REAL* W, const Metric m); //utilities.cc
void oneFormToVector(CCTK_REAL* X_vector, const CCTK_REAL* X_oneform, const Metric m); //utilities.cc
void vectorToOneForm(CCTK_REAL* X_oneform, const CCTK_REAL* X_vector, const Metric m); //utilities.cc
void projectUontoV(CCTK_REAL* X, const CCTK_REAL* U, const CCTK_REAL* V, const Metric m); //utilities.cc
void normalize(CCTK_REAL* X_norm, const CCTK_REAL* X, const Metric m); //utilities.cc
CCTK_REAL getTimeComponentOf4Velocity(const CCTK_REAL vx, const CCTK_REAL vy, const CCTK_REAL vz, const Metric m); //utilities.cc

#endif