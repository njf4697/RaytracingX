#include "raytracingx.h"
#include <vector>
#include <array>

void interpolateMetricAtPoint(CCTK_ARGUMENTS, const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z, Metric* metric_at_point) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    //if (x < CarpetX::xmin || x > CarpetX::xmax || y < CarpetX::ymin || y > CarpetX::ymax || z < CarpetX::zmin || z > CarpetX::zmax) {
    //  CCTK_VERROR("Camera Location Out of Bounds");
    //}

    //only raytrace on one processor
    const CCTK_INT nPoints = 1;//(CCTK_MyProc(cctkGH) == 0) ? 1 : 0;

    std::array<std::vector<CCTK_REAL>, 3> location_;
    location_[0].push_back(x);
    location_[1].push_back(y);
    location_[2].push_back(z);

    std::array<std::vector<CCTK_REAL>, 10> metric_;
    metric_[0].resize(1);
    metric_[1].resize(1);
    metric_[2].resize(1);
    metric_[3].resize(1);
    metric_[4].resize(1);
    metric_[5].resize(1);
    metric_[6].resize(1);
    metric_[7].resize(1);
    metric_[8].resize(1);
    metric_[9].resize(1);
    
    // Interpolation coordinates
    const void *interpCoords[3] = {
        location_[0].data(), location_[1].data(), location_[2].data()};

    // Interpolated variables
    const CCTK_INT nInputArrays = 10;
    const CCTK_INT inputArrayIndices[nInputArrays] = {
        CCTK_VarIndex("ADMBaseX::alp"), 
        CCTK_VarIndex("ADMBaseX::betax"), CCTK_VarIndex("ADMBaseX::betay"), CCTK_VarIndex("ADMBaseX::betaz"),
        CCTK_VarIndex("ADMBaseX::gxx"), CCTK_VarIndex("ADMBaseX::gxy"), CCTK_VarIndex("ADMBaseX::gxz"),
        CCTK_VarIndex("ADMBaseX::gyy"), CCTK_VarIndex("ADMBaseX::gyz"), CCTK_VarIndex("ADMBaseX::gzz")};

    CCTK_POINTER outputArrays[nInputArrays] = {metric_[0].data(), metric_[1].data(), metric_[2].data(), metric_[3].data(), metric_[4].data(),
                                               metric_[5].data(), metric_[6].data(), metric_[7].data(), metric_[8].data(), metric_[9].data()};

    // DriverInterpolate arguments that aren't currently used
    const int coordSystemHandle = 0;
    const CCTK_INT interpCoordsTypeCode = 0;
    const CCTK_INT outputArrayTypes[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

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

    printf(("camera position: x: " + std::to_string(x) + ", y: " + std::to_string(y) + ", z: " + std::to_string(z) + "\n").c_str());

    // Perform the interpolation
    ierr = DriverInterpolate(cctkGH, 3, interpHandle, paramTableHandle,
                             coordSystemHandle, nPoints, interpCoordsTypeCode,
                             interpCoords, nInputArrays, inputArrayIndices,
                             nInputArrays, outputArrayTypes, outputArrays);

    if (ierr < 0) {
      CCTK_WARN(CCTK_WARN_ALERT, "Interpolation error");
    }

    // Destroy the parameter table
    Util_TableDestroy(paramTableHandle);

    printf(std::to_string(metric_[0].data()[0]).c_str());

    metric_at_point->alpha = metric_[0].data()[0];
    metric_at_point->beta_xup = metric_[1].data()[0];
    metric_at_point->beta_yup = metric_[2].data()[0];
    metric_at_point->beta_zup = metric_[3].data()[0];
    metric_at_point->g_xx = metric_[4].data()[0];
    metric_at_point->g_xy = metric_[5].data()[0];
    metric_at_point->g_xz = metric_[6].data()[0];
    metric_at_point->g_yy = metric_[7].data()[0];
    metric_at_point->g_yz = metric_[8].data()[0];
    metric_at_point->g_zz = metric_[9].data()[0];

    //\beta_i = \gamma_ij*\beta^j
    metric_at_point->beta_x = metric_at_point->beta_xup*metric_at_point->g_xx + metric_at_point->beta_yup*metric_at_point->g_xy + metric_at_point->beta_zup*metric_at_point->g_xz;
    metric_at_point->beta_y = metric_at_point->beta_xup*metric_at_point->g_xy + metric_at_point->beta_yup*metric_at_point->g_yy + metric_at_point->beta_zup*metric_at_point->g_yz;
    metric_at_point->beta_z = metric_at_point->beta_xup*metric_at_point->g_xz + metric_at_point->beta_yup*metric_at_point->g_yz + metric_at_point->beta_zup*metric_at_point->g_zz;

    //g_{tt}=\alpha^2-\beta_i\beta^i
    metric_at_point->g_tt = -pow(metric_at_point->alpha,2) + metric_at_point->beta_x*metric_at_point->beta_xup + metric_at_point->beta_y*metric_at_point->beta_yup + metric_at_point->beta_z*metric_at_point->beta_zup; //g_{00}

    metric_at_point->fillInverseMetric(); //get g^{\mu\nu}

    printf(metric_at_point->to_string().c_str());
    printf(metric_at_point->to_string_inv().c_str());
} 