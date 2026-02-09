#include "raytracingx.h"
#include <stdio.h>
#include <AMReX_MFIter.H>
#include "AMReX_ParallelDescriptor.H"

#define DEG2RAD 0.01745329251

void gramSchmidtProcess(CCTK_ARGUMENTS, CCTK_REAL* e0, CCTK_REAL* e1, CCTK_REAL* e2, CCTK_REAL* e3, Metric* metric) { //use Gram-Schmidt Process to turn e0, e1, e2, and e3 into orthonormal vectors for the camera's POV (see https://arxiv.org/pdf/1410.7775)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    //notation borrowed from https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process

    CCTK_REAL u0[4] = {getTimeComponentOf4Velocity(camera_vel[0], camera_vel[1], camera_vel[2], metric), camera_vel[0], camera_vel[1], camera_vel[2]}; // u0=v0 is 4-velocity of camera
    CCTK_REAL u1[4];
    CCTK_REAL u2[4];
    CCTK_REAL u3[4];
    CCTK_REAL v1[4] = {0, camera_point[0], camera_point[1], camera_point[2]}; //v1 is camera pointing direciton
    CCTK_REAL v2[4] = {0, camera_right[0], camera_right[1], camera_right[2]}; //v2 is camera right direction
    

    CCTK_REAL projv1_onto_u0[4];
    projectUontoV(projv1_onto_u0, v1, u0, metric);
    u1[0] = v1[0] - projv1_onto_u0[0];
    u1[1] = v1[1] - projv1_onto_u0[1];
    u1[2] = v1[2] - projv1_onto_u0[2];
    u1[3] = v1[3] - projv1_onto_u0[3];

    CCTK_REAL projv2_onto_u0[4];
    projectUontoV(projv2_onto_u0, v2, u0, metric);
    CCTK_REAL projv2_onto_u1[4];
    projectUontoV(projv2_onto_u1, v2, u1, metric);
    u2[0] = v2[0] - projv2_onto_u0[0] - projv2_onto_u1[0];
    u2[1] = v2[1] - projv2_onto_u0[1] - projv2_onto_u1[1];
    u2[2] = v2[2] - projv2_onto_u0[2] - projv2_onto_u1[2];
    u2[3] = v2[3] - projv2_onto_u0[3] - projv2_onto_u1[3];

    normalize(e0, u0, metric); //e0 should be unit vector in direction of camera 4-vellocity
    normalize(e1, u1, metric); //e1 is (more or less) unit vector in camera facting direction
    normalize(e2, u2, metric); //e2 is (more or less) unit vector in from the center of the camera to the right edge

    generalizedCrossProduct(e3, e0, e1, e2, metric); //e3 is an normalized vector orthogonal to e0, e1, e2, e3
}

void setup_camera_initializer_reals(CCTK_ARGUMENTS, CCTK_REAL* real_params) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    Metric metric;
    interpolateMetricAtPoint(CCTK_PASS_CTOC, camera_pos[0], camera_pos[1], camera_pos[2], &metric); //interpolate metric values and store in Metric struct

    if (CCTK_MyProc(cctkGH) != 0) return;

    CCTK_REAL e0[4];
    CCTK_REAL e1[4];
    CCTK_REAL e2[4];
    CCTK_REAL e3[4];
    gramSchmidtProcess(CCTK_PASS_CTOC, e0, e1, e2, e3, &metric); //create orthonormal basis for camera POV

    CCTK_REAL alpha_h = DEG2RAD * horizontal_fov; //convert FOV to radians
    CCTK_REAL alpha_v = DEG2RAD * vertical_fov;

    real_params[0] = metric.g_tt;
    real_params[1] = metric.beta_x;
    real_params[2] = metric.beta_y;
    real_params[3] = metric.beta_z;
    real_params[4] = metric.g_xx;
    real_params[5] = metric.g_xy;
    real_params[6] = metric.g_xz;
    real_params[7] = metric.g_yy;
    real_params[8] = metric.g_yz;
    real_params[9] = metric.g_zz;
    real_params[10] = e0[0];
    real_params[11] = e0[1];
    real_params[12] = e0[2];
    real_params[13] = e0[3];
    real_params[14] = e1[0];
    real_params[15] = e1[1];
    real_params[16] = e1[2];
    real_params[17] = e1[3];
    real_params[18] = e2[0];
    real_params[19] = e2[1];
    real_params[20] = e2[2];
    real_params[21] = e2[3];
    real_params[22] = e3[0];
    real_params[23] = e3[1];
    real_params[24] = e3[2];
    real_params[25] = e3[3];
    real_params[26] = alpha_h;
    real_params[27] = alpha_v;
    real_params[28] = metric.alpha;
    real_params[29] = camera_pos[0];
    real_params[30] = camera_pos[1];
    real_params[31] = camera_pos[2];
}

void setup_camera_initializer_ints(CCTK_ARGUMENTS, CCTK_INT* int_params) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int_params[0] = num_pixels_width;
    int_params[1] = num_pixels_height;
}