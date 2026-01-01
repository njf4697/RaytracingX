#include "raytracingx.h"
#include <vector>
#include <array>

#define NUM_INTERP_POINTS 1
#define NUM_GRID_ARRAYS 10

void interpolateMetricAtPoint(CCTK_ARGUMENTS, const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z, Metric metric_at_point) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    //if (x < CarpetX::xmin || x > CarpetX::xmax || y < CarpetX::ymin || y > CarpetX::ymax || z < CarpetX::zmin || z > CarpetX::zmax) {
    //  CCTK_VERROR("Camera Location Out of Bounds");
    //}

    // Only Processor 0 interpolates
    const CCTK_INT nPoints = (CCTK_MyProc(cctkGH) == 0) ? 1 : 0;
    printf(std::to_string(CCTK_MyProc(cctkGH)).c_str());

    std::array<std::vector<CCTK_REAL>, 3> location_;
    location_[0].push_back(x);
    location_[1].push_back(y);
    location_[2].push_back(z);

    std::array<std::vector<CCTK_REAL>, 10> metric_;

    // Interpolation coordinates
    const void *interpCoords[3] = {
        location_[0].data(), location_[1].data(), location_[2].data()};

    // Interpolated variables
    const CCTK_INT nInputArrays = 10;
    const CCTK_INT inputArrayIndices[nInputArrays] = {
        CCTK_VarIndex("ADMBaseX::alp"), 
        CCTK_VarIndex("ADMBaseX::betax"), CCTK_VarIndex("ADMBaseX::betay"), CCTK_VarIndex("ADMBaseX::betaz"),
        CCTK_VarIndex("ADMBaseX::gxx"), CCTK_VarIndex("ADMBaseX::gxy"), CCTK_VarIndex("ADMBaseX::gxz"),
        CCTK_VarIndex("ADMBaseX::gyy"), CCTK_VarIndex("ADMBaseX::gyz"), CCTK_VarIndex("ADMBaseX::gzz"),};

    CCTK_POINTER outputArrays[nInputArrays] = {metric_[0].data(), metric_[1].data(), metric_[2].data(), metric_[3].data(), metric_[4].data(),
                                               metric_[5].data(), metric_[6].data(), metric_[7].data(), metric_[8].data(), metric_[9].data()};

    // DriverInterpolate arguments that aren't currently used
    const int coordSystemHandle = 0;
    const CCTK_INT interpCoordsTypeCode = 0;
    const CCTK_INT outputArrayTypes[1] = {0};

    const int interpHandle = CCTK_InterpHandle("CarpetX");
    if (interpHandle < 0) {
      CCTK_WARN(CCTK_WARN_ALERT, "Can't get interpolation handle");
      return;
    }

    // Create parameter table for interpolation
    const int paramTableHandle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
    if (paramTableHandle < 0) {
      CCTK_VERROR("Can't create parameter table: %d", paramTableHandle);
    }

    // Set interpolation order in the parameter table
    int ierr = Util_TableSetInt(paramTableHandle, 1, "order");
    if (ierr < 0) {
      CCTK_VERROR("Can't set order in parameter table: %d", ierr);
    }

    printf(("camera position: x: " + std::to_string(x) + ", y: " + std::to_string(y) + ", z: " + std::to_string(z)).c_str());

    // Perform the interpolation
    ierr = DriverInterpolate(cctkGH, 3, interpHandle, paramTableHandle,
                             coordSystemHandle, nPoints, interpCoordsTypeCode,
                             interpCoords, nInputArrays, inputArrayIndices,
                             nInputArrays, outputArrayTypes, outputArrays);

    if (ierr < 0) {
      CCTK_WARN(CCTK_WARN_ALERT, "Interpolation error");
    }

    metric_at_point.metric[0] = metric_[0].data()[0];
    metric_at_point.beta_up[0] = metric_[1].data()[0];
    metric_at_point.beta_up[1] = metric_[2].data()[0];
    metric_at_point.beta_up[2] = metric_[3].data()[0];
    metric_at_point.metric[4] = metric_[4].data()[0];
    metric_at_point.metric[5] = metric_[5].data()[0];
    metric_at_point.metric[6] = metric_[6].data()[0];
    metric_at_point.metric[7] = metric_[7].data()[0];
    metric_at_point.metric[8] = metric_[8].data()[0];
    metric_at_point.metric[9] = metric_[9].data()[0];

    metric_at_point.metric[1] = metric_at_point.beta_up[0]*metric_at_point.metric[4] + metric_at_point.beta_up[1]*metric_at_point.metric[5] + metric_at_point.beta_up[2]*metric_at_point.metric[6];
    metric_at_point.metric[2] = metric_at_point.beta_up[0]*metric_at_point.metric[5] + metric_at_point.beta_up[1]*metric_at_point.metric[7] + metric_at_point.beta_up[2]*metric_at_point.metric[8];
    metric_at_point.metric[3] = metric_at_point.beta_up[0]*metric_at_point.metric[6] + metric_at_point.beta_up[1]*metric_at_point.metric[8] + metric_at_point.beta_up[2]*metric_at_point.metric[9];

    metric_at_point.metric0pr = sqrt(pow(metric_at_point.metric[0],2) - metric_at_point.metric[1]*metric_at_point.beta_up[0] - metric_at_point.metric[2]*metric_at_point.beta_up[1] - metric_at_point.metric[3]*metric_at_point.beta_up[2]); //g_{00}

    printf(("camera position: x: " + std::to_string(x) + ", y: " + std::to_string(y) + ", z: " + std::to_string(z)).c_str());

    // Destroy the parameter table
    Util_TableDestroy(paramTableHandle);

    metric_at_point.fillInverseMetric(); //get g^{\mu\nu}
} 

/* old version
void interpolateMetricAtPointOld(CCTK_ARGUMENTS, const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z, Metric metric_at_point) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL interp_x[NUM_INTERP_POINTS];
    CCTK_REAL interp_y[NUM_INTERP_POINTS];
    CCTK_REAL interp_z[NUM_INTERP_POINTS];
    
    interp_x[0] = x;
    interp_y[0] = y;
    interp_z[0] = z;

    //assert(x > xmin);
    //assert(x < xmax);
    //assert(y > ymin);
    //assert(y < ymax);
    //assert(z > zmin);
    //assert(z < zmax);

    const void* interp_coords[3] = {(const void*) interp_x, (const void*) interp_y, (const void*) interp_z};

    CCTK_INT variable_indices[NUM_GRID_ARRAYS];
    variable_indices[0] = CCTK_VarIndex("ADMBaseX::alp");
    assert(variable_indices[0] >= 0);
    variable_indices[1] = CCTK_VarIndex("ADMBaseX::betax");
    assert(variable_indices[1] >= 0);
    variable_indices[2] = CCTK_VarIndex("ADMBaseX::betay");
    assert(variable_indices[2] >= 0);
    variable_indices[3] = CCTK_VarIndex("ADMBaseX::betaz");
    assert(variable_indices[3] >= 0);
    variable_indices[4] = CCTK_VarIndex("ADMBaseX::gxx");
    assert(variable_indices[4] >= 0);
    variable_indices[5] = CCTK_VarIndex("ADMBaseX::gxy");
    assert(variable_indices[5] >= 0);
    variable_indices[6] = CCTK_VarIndex("ADMBaseX::gxz");
    assert(variable_indices[6] >= 0);
    variable_indices[7] = CCTK_VarIndex("ADMBaseX::gyy");
    assert(variable_indices[7] >= 0);
    variable_indices[8] = CCTK_VarIndex("ADMBaseX::gyz");
    assert(variable_indices[8] >= 0);
    variable_indices[9] = CCTK_VarIndex("ADMBaseX::gzz");
    assert(variable_indices[9] >= 0);

    //not used by CarpetX:
    //static const CCTK_INT output_array_type_codes[NUM_GRID_ARRAYS] = {CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL};
    static const CCTK_INT output_array_type_codes[NUM_GRID_ARRAYS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

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
    
    int operator_handle = 0; //not used by CarpetX
    //int operator_handle = CCTK_InterpHandle("uniform cartesian");
    //assert(operator_handle >= 0);
    
    int coord_system_handle = 0; //not used by CarpetX
    //int coord_system_handle = CCTK_CoordSystemHandle("cart3d");
    //assert(coord_system_handle >= 0);

    CCTK_INT operand_indices[NUM_GRID_ARRAYS];
    for (int i = 0; i < NUM_GRID_ARRAYS; i++) {
      operand_indices[i] = i;
    }

    CCTK_INT operation_codes[NUM_GRID_ARRAYS];
    for (int i = 0; i < NUM_GRID_ARRAYS; i++) {
      operation_codes[i] = 0;
    }

    //std::string orderstr = "order=" + std::to_string(interpolation_order); //from CarpetX
    std::string orderstr = "order=1";

    int param_table_handle = Util_TableCreateFromString(orderstr.c_str());
    if (param_table_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "bad interpolator parameter(s) \"%s\"!", orderstr.c_str());
    }

    Util_TableSetIntArray(param_table_handle, NUM_GRID_ARRAYS,
                          operand_indices, "operand_indices");

    Util_TableSetIntArray(param_table_handle, NUM_GRID_ARRAYS,
                        operation_codes, "operation_codes");

    printf(("x: " + std::to_string(x) + ", y: " + std::to_string(y) + ", z: " + std::to_string(z)).c_str());

    const cGH *GH;
    int status = CCTK_InterpGridArrays(GH, 3, operator_handle, param_table_handle, coord_system_handle, NUM_INTERP_POINTS, CCTK_VARIABLE_REAL, interp_coords, 
                                                                                                                            NUM_GRID_ARRAYS, variable_indices,
                                                                                                                            NUM_GRID_ARRAYS, output_array_type_codes, output_arrays);
    assert(status >= 0);
    
    CCTK_REAL* alp_interp_filled = (CCTK_REAL*) output_arrays[0];
    CCTK_REAL* betax_interp_filled = (CCTK_REAL*) output_arrays[1];
    CCTK_REAL* betay_interp_filled = (CCTK_REAL*) output_arrays[2];
    CCTK_REAL* betaz_interp_filled = (CCTK_REAL*) output_arrays[3];
    CCTK_REAL* gxx_interp_filled = (CCTK_REAL*) output_arrays[4];
    CCTK_REAL* gxy_interp_filled = (CCTK_REAL*) output_arrays[5];
    CCTK_REAL* gxz_interp_filled = (CCTK_REAL*) output_arrays[6];
    CCTK_REAL* gyy_interp_filled = (CCTK_REAL*) output_arrays[7];
    CCTK_REAL* gyz_interp_filled = (CCTK_REAL*) output_arrays[8];
    CCTK_REAL* gzz_interp_filled = (CCTK_REAL*) output_arrays[9];

    metric_at_point.metric[0] = alp_interp_filled[0]; // \alpha
    metric_at_point.metric[1] = betax_interp_filled[0]*gxx_interp_filled[0] + betay_interp_filled[0]*gxy_interp_filled[0] + betaz_interp_filled[0]*gxz_interp_filled[0]; // g_{0j} = \beta_j = \gamma_{ij} \beta^i
    metric_at_point.metric[2] = betax_interp_filled[0]*gxy_interp_filled[0] + betay_interp_filled[0]*gyy_interp_filled[0] + betaz_interp_filled[0]*gyz_interp_filled[0];
    metric_at_point.metric[3] = betax_interp_filled[0]*gxz_interp_filled[0] + betay_interp_filled[0]*gyz_interp_filled[0] + betaz_interp_filled[0]*gzz_interp_filled[0];
    metric_at_point.metric[4] = gxx_interp_filled[0]; // g_{xx}
    metric_at_point.metric[5] = gxy_interp_filled[0]; // g_{xy}
    metric_at_point.metric[6] = gxz_interp_filled[0]; // g_{xz}
    metric_at_point.metric[7] = gyy_interp_filled[0]; // g_{yy}
    metric_at_point.metric[8] = gyz_interp_filled[0]; // g_{yz}
    metric_at_point.metric[9] = gzz_interp_filled[0]; // g_{zz}

    metric_at_point.beta_up[0] = betax_interp_filled[0]; //\beta^i
    metric_at_point.beta_up[1] = betay_interp_filled[0];
    metric_at_point.beta_up[2] = betaz_interp_filled[0];
    metric_at_point.metric0pr = sqrt(pow(metric_at_point.metric[0],2) - metric_at_point.metric[1]*metric_at_point.beta_up[0] - metric_at_point.metric[2]*metric_at_point.beta_up[1] - metric_at_point.metric[3]*metric_at_point.beta_up[2]); //g_{00}

    calculateInverseMetric(metric_at_point); //get g^{\mu\nu}

    Util_TableDestroy(param_table_handle);
}
*/