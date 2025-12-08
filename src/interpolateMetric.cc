#include "raytracingx.h"

#define NUM_INTERP_POINTS 1

void inverseSpatialMetric(CCTK_REAL* inv_spatial_metric, const Metric m) { //find \gamma^{ij}
    CCTK_REAL inv_det_g = 1/(m.metric[4]*m.metric[7]*m.metric[9]+2.*m.metric[5]*m.metric[6]*m.metric[8]-m.metric[6]*m.metric[6]*m.metric[7]-m.metric[8]*m.metric[8]*m.metric[4]-m.metric[5]*m.metric[5]*m.metric[9]);
    inv_spatial_metric[0] = (m.metric[3]*m.metric[5]-m.metric[4]*m.metric[4])*inv_det_g;
    inv_spatial_metric[1] = (m.metric[4]*m.metric[2]-m.metric[1]*m.metric[5])*inv_det_g;
    inv_spatial_metric[2] = (m.metric[1]*m.metric[4]-m.metric[2]*m.metric[3])*inv_det_g;
    inv_spatial_metric[3] = (m.metric[0]*m.metric[5]-m.metric[2]*m.metric[2])*inv_det_g;
    inv_spatial_metric[4] = (m.metric[2]*m.metric[1]-m.metric[0]*m.metric[4])*inv_det_g;
    inv_spatial_metric[5] = (m.metric[0]*m.metric[3]-m.metric[1]*m.metric[1])*inv_det_g;
}

void calculateInverseMetric(Metric m) { //find g^{\mu\nu}
    CCTK_REAL h[6];
    inverseSpatialMetric(h, m);  //\gamma^{ij}
    m.metric_inv[0] = 1/(pow(m.metric[0],2));
    m.metric_inv[1] = m.metric_inv[0]*m.beta_up[0];
    m.metric_inv[2] = m.metric_inv[0]*m.beta_up[1];
    m.metric_inv[3] = m.metric_inv[0]*m.beta_up[2];
    m.metric_inv[4] = h[0]-m.metric_inv[1]*m.beta_up[0];
    m.metric_inv[5] = h[1]-m.metric_inv[1]*m.beta_up[1];
    m.metric_inv[6] = h[2]-m.metric_inv[1]*m.beta_up[2];
    m.metric_inv[7] = h[3]-m.metric_inv[2]*m.beta_up[1];
    m.metric_inv[8] = h[4]-m.metric_inv[2]*m.beta_up[2];
    m.metric_inv[9] = h[5]-m.metric_inv[3]*m.beta_up[2];
}

void interpolateMetricAtPoint(CCTK_ARGUMENTS, const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z, Metric metric_at_point) {

    //uses 4th-order generalized polynomial interpolation to find the spacetime quantities at given position
    //see page A126 of https://www.cactuscode.org/documentation/ReferenceManual.pdf

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int NUM_GRID_ARRAYS = 10;

    CCTK_REAL interp_x[NUM_INTERP_POINTS];
    CCTK_REAL interp_y[NUM_INTERP_POINTS];
    CCTK_REAL interp_z[NUM_INTERP_POINTS];
    
    interp_x[0] = x;
    interp_y[0] = y;
    interp_z[0] = z;

    const void* interp_coords[3];
    interp_coords[0] = (const void*) interp_x;
    interp_coords[1] = (const void*) interp_y;
    interp_coords[2] = (const void*) interp_z;

    CCTK_INT variable_indices[NUM_GRID_ARRAYS];
    variable_indices[0] = CCTK_VarIndex("ADMBaseX::alp");
    variable_indices[1] = CCTK_VarIndex("ADMBaseX::betax");
    variable_indices[2] = CCTK_VarIndex("ADMBaseX::betay");
    variable_indices[3] = CCTK_VarIndex("ADMBaseX::betaz");
    variable_indices[4] = CCTK_VarIndex("ADMBaseX::gxx");
    variable_indices[5] = CCTK_VarIndex("ADMBaseX::gxy");
    variable_indices[6] = CCTK_VarIndex("ADMBaseX::gxz");
    variable_indices[7] = CCTK_VarIndex("ADMBaseX::gyy");
    variable_indices[8] = CCTK_VarIndex("ADMBaseX::gyz");
    variable_indices[9] = CCTK_VarIndex("ADMBaseX::gzz");

    static const CCTK_INT output_array_type_codes[NUM_GRID_ARRAYS];
    output_array_type_codes[0] = CCTK_VARIABLE_REAL;
    output_array_type_codes[1] = CCTK_VARIABLE_REAL;
    output_array_type_codes[2] = CCTK_VARIABLE_REAL;
    output_array_type_codes[3] = CCTK_VARIABLE_REAL;
    output_array_type_codes[4] = CCTK_VARIABLE_REAL;
    output_array_type_codes[5] = CCTK_VARIABLE_REAL;
    output_array_type_codes[6] = CCTK_VARIABLE_REAL;
    output_array_type_codes[7] = CCTK_VARIABLE_REAL;
    output_array_type_codes[8] = CCTK_VARIABLE_REAL;
    output_array_type_codes[9] = CCTK_VARIABLE_REAL;

    CCTK_REAL alp_interp[NUM_INTERP_POINTS];
    CCTK_REAL betax_interp[NUM_INTERP_POINTS];
    CCTK_REAL betay_interp[NUM_INTERP_POINTS];
    CCTK_REAL betaz_interp[NUM_INTERP_POINTS];
    CCTK_REAL gxx_interp[NUM_INTERP_POINTS];
    CCTK_REAL gxy_interp[NUM_INTERP_POINTS];
    CCTK_REAL gxz_interp[NUM_INTERP_POINTS];
    CCTK_REAL gyy_interp[NUM_INTERP_POINTS];
    CCTK_REAL gyz_interp[NUM_INTERP_POINTS];
    CCTK_REAL gzz_interp[NUM_INTERP_POINTS];

    void* output_arrays[NUM_GRID_ARRAYS];
    output_arrays[0] = (void*) alp_interp;
    output_arrays[1] = (void*) betax_interp;
    output_arrays[2] = (void*) betay_interp;
    output_arrays[3] = (void*) betaz_interp;
    output_arrays[4] = (void*) gxx_interp;
    output_arrays[5] = (void*) gxy_interp;
    output_arrays[6] = (void*) gxz_interp;
    output_arrays[7] = (void*) gyy_interp;
    output_arrays[8] = (void*) gyz_interp;
    output_arrays[9] = (void*) gzz_interp;

    int operator_handle = CCTK_InterpHandle("generalized polynomial interpolation");
    assert(operator_handle >= 0);

    int coord_system_handle = CCTK_CoordSystemHandle("cart3d");
    assert(coord_system_handle >= 0);

    const cGH *GH;
    int status = CCTK_InterpGridArrays(GH, 3, operator_handle, Util_TableCreateFromString("order=4"), coord_system_handle, NUM_INTERP_POINTS, CCTK_VARIABLE_REAL, interp_coords, 
                                                                                                                            NUM_GRID_ARRAYS, variable_indices,
                                                                                                                            NUM_GRID_ARRAYS, output_array_type_codes, output_arrays);
    assert(status >= 0);

    metric_at_point.metric[0] = output_arrays[0][0]; // \alpha
    metric_at_point.metric[1] = output_arrays[1][0]*output_arrays[4][0] + output_arrays[2][0]*output_arrays[5][0] + output_arrays[3][0]*output_arrays[6][0]; // g_{0j} = \beta_j = \gamma_{ij} \beta^i
    metric_at_point.metric[2] = output_arrays[1][0]*output_arrays[5][0] + output_arrays[2][0]*output_arrays[7][0] + output_arrays[3][0]*output_arrays[8][0];
    metric_at_point.metric[3] = output_arrays[1][0]*output_arrays[6][0] + output_arrays[2][0]*output_arrays[8][0] + output_arrays[3][0]*output_arrays[9][0];
    metric_at_point.metric[4] = output_arrays[4][0]; // g_{xx}
    metric_at_point.metric[5] = output_arrays[5][0]; // g_{xy}
    metric_at_point.metric[6] = output_arrays[6][0]; // g_{xz}
    metric_at_point.metric[7] = output_arrays[7][0]; // g_{yy}
    metric_at_point.metric[8] = output_arrays[8][0]; // g_{yz}
    metric_at_point.metric[9] = output_arrays[9][0]; // g_{zz}

    metric_at_point.beta_up[0] = output_arrays[1][0]; //\beta^i
    metric_at_point.beta_up[1] = output_arrays[2][0];
    metric_at_point.beta_up[2] = output_arrays[3][0];
    metric_at_point.metric0pr = sqrt(pow(metric_at_point.metric[0],2) - metric_at_point.metric[1]*metric_at_point.beta_up[0] - metric_at_point.metric[2]*metric_at_point.beta_up[1] - metric_at_point.metric[3]*metric_at_point.beta_up[2]); //g_{00}

    calculateInverseMetric(metric_at_point); //get g^{\mu\nu}
}