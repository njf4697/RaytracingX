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
    CCTK_REAL metric[10]; // {\alpha, \beta_0, \beta_1, \beta_2, \gamma_{11}, \gamma_{12}, \gamma_{13}, \gamma_{22}, \gamma_{23}, \gamma_{33}}
    CCTK_REAL metric_inv[10]; // {g^{00}, g^{01}, g^{02}, g^{03}, g^{11}, g^{12}, g^{13}, g^{22}, g^{23}, g^{33}}
    CCTK_REAL metric0pr; // = sqrt(\alpha^2 + \beta_i\beta^i)
    CCTK_REAL beta_up[3]; // \beta^i

    std::string to_string() {
        return "alpha: " + std::to_string(metric[0]) + "\n";
               "beta: (" + std::to_string(beta_up[0]) + ", " + std::to_string(beta_up[1]) + ", " + std::to_string(beta_up[2]) + ")";
               "metric: " + std::to_string(metric0pr) + ", " + std::to_string(metric[1]) + ", " + std::to_string(metric[2]) + ", " + std::to_string(metric[3]) + "\n";
               "        " + std::to_string(metric[1]) + ", " + std::to_string(metric[4]) + ", " + std::to_string(metric[5]) + ", " + std::to_string(metric[6]) + "\n";
               "        " + std::to_string(metric[2]) + ", " + std::to_string(metric[5]) + ", " + std::to_string(metric[7]) + ", " + std::to_string(metric[8]) + "\n";
               "        " + std::to_string(metric[3]) + ", " + std::to_string(metric[6]) + ", " + std::to_string(metric[8]) + ", " + std::to_string(metric[9]);
    }

    void fillInverseMetric() {
        CCTK_REAL inv_det_g = 1/(metric[4]*metric[7]*metric[9]+2.*metric[5]*metric[6]*metric[8]-metric[6]*metric[6]*metric[7]-metric[8]*metric[8]*metric[4]-metric[5]*metric[5]*metric[9]);
        CCTK_REAL inv_spatial_metric[6];
        inv_spatial_metric[0] = (metric[3]*metric[5]-metric[4]*metric[4])*inv_det_g;
        inv_spatial_metric[1] = (metric[4]*metric[2]-metric[1]*metric[5])*inv_det_g;
        inv_spatial_metric[2] = (metric[1]*metric[4]-metric[2]*metric[3])*inv_det_g;
        inv_spatial_metric[3] = (metric[0]*metric[5]-metric[2]*metric[2])*inv_det_g;
        inv_spatial_metric[4] = (metric[2]*metric[1]-metric[0]*metric[4])*inv_det_g;
        inv_spatial_metric[5] = (metric[0]*metric[3]-metric[1]*metric[1])*inv_det_g;

        metric_inv[0] = 1/(pow(metric[0],2));
        metric_inv[1] = metric_inv[0]*beta_up[0];
        metric_inv[2] = metric_inv[0]*beta_up[1];
        metric_inv[3] = metric_inv[0]*beta_up[2];
        metric_inv[4] = inv_spatial_metric[0]-metric_inv[1]*beta_up[0];
        metric_inv[5] = inv_spatial_metric[1]-metric_inv[1]*beta_up[1];
        metric_inv[6] = inv_spatial_metric[2]-metric_inv[1]*beta_up[2];
        metric_inv[7] = inv_spatial_metric[3]-metric_inv[2]*beta_up[1];
        metric_inv[8] = inv_spatial_metric[4]-metric_inv[2]*beta_up[2];
        metric_inv[9] = inv_spatial_metric[5]-metric_inv[3]*beta_up[2];
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