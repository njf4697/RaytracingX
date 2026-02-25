#ifndef RAYTRACINGCONTAINER
#define RAYTRACINGCONTAINER

#include <cctk.h>

#include "Photons.hxx"
#include "PhotonsContainer.hxx"
#include "raytracingx.h"
#include <AMReX_ParallelDescriptor.H>
#include <CParameters.h>

//#define CHECK_OUT_OF_BOUNDS(X,Y) out_of_bounds |= (X0Y > boundarie_hx) || (X0Y < boundarie_lx);\ //old version
//                                 out_of_bounds |= (X1Y > boundarie_hy) || (X1Y < boundarie_ly);\
//                                 out_of_bounds |= (X2Y > boundarie_hz) || (X2Y < boundarie_lz);
#define CHECK_OUT_OF_BOUNDS(X,Y) if (X0Y > boundarie_hx) {\
                                    out_of_bounds = true;      \
                                    tau[i] = -1;               \
                                 }                             \
                                 if (X0Y < boundarie_lx) {\
                                    out_of_bounds = true;      \
                                    tau[i] = -2;               \
                                 }                             \
                                 if (X1Y > boundarie_hy) {\
                                    out_of_bounds = true;      \
                                    tau[i] = -3;               \
                                 }                             \
                                 if (X1Y < boundarie_ly) {\
                                    out_of_bounds = true;      \
                                    tau[i] = -4;               \
                                 }                             \
                                 if (X2Y > boundarie_hz) {\
                                    out_of_bounds = true;      \
                                    tau[i] = -5;               \
                                 }                             \
                                 if (X2Y < boundarie_lz) {\
                                    out_of_bounds = true;      \
                                    tau[i] = -6;               \
                                 }                             

namespace RaytracingPhotons {

struct RaytracingPhotonsData : public Photons::PhotonsData{
    enum {
        vx = 0,       /**< Velocity's lower index on the x direction.*/
        vy,           /**< Velocity's lower index on the y direction.*/
        vz,           /**< Velocity's lower index on the z direction.*/
        ln_E,            /**< Ln Energy value.*/
        tau,             /**< Optical depth value.*/
        index,           /**< Pixel index number.*/
        n_attributes, /**< Total number of attributes*/
    }; // enum
};

}

namespace RaytracingContainers {

template <typename StructType>
class RaytracingPhotonsContainer : public GInX::PhotonsContainer<StructType> {

// ##############################################################################
//                   PhotonsContainer::METHODS DECLARATION
// ##############################################################################

/**
 * \brief Computes the right hand side of the geodesic differential equation.
 *
 * Given differential equation \f[\frac{d}{dt}U = f\left(U, \frac{dU}{dx};
 * t\right)\f] computes
 * \f[f\left(U, \frac{dU}{dx}; t\right)\f]
 *
 * where \f$U\f$ is a vector that contains \f$(x_u, y_u, z_u, vx_d, vy_d, vz_d,
 * E)\f$. That's why each rhs depends on the other components of the vector
 * \f$U\f$. For the position the differential equation  is:
 *
 * \f[\frac{d}{dt} U[i] = \alpha \gamma^{ij} U[3 + j] - \beta^i\f]
 *
 * Where \f$i,j = 0, 1, 2\f$, \f$\gamma\f$ is the induced metric, \f$\alpha\f$
 * is the lapse function and \f$\beta\f$ is the shift vector.
 *
 * For the Velocity_d the differential equation is:
 *
 * \f{eqnarray*}{
 * \frac{d}{dt}U[3 + i] &= -\partial_i\alpha + \left(\gamma^{kj} U[3 + k]
 * \partial_j\alpha - \alpha K_{jk}\gamma^{jl}\gamma^{km}U[3+l]U[3+m]\right) U[3
 * + i]\\ & +
 * \frac{1}{2}\alpha\gamma^{jl}\gamma^{km}U[3+l]U[3+m]\partial_i\gamma{jk} + U[3
 * + j] \partial_i\beta^j
 * \f}
 *
 * and finally, for the \f$ ln E \f$ the differential equation is:
 *
 * \f[ \dfrac{d}{dt} U[6] = \alpha K_{jk}U[3 + l]U[3 + m]\gamma^{lj}\gamma^{mk}
 * - U[3+l]\gamma^{lj}\partial_j\alpha\f]
 *
 *  Where \f$i, j, k, l, m = 0, 1, 2\f$ and we have been using Einstein
 * notation.
 *
 *  @param u A GpuArray of size n_attributes + the coordinates that contains the
 * varaibles needed to evolve.
 *  @param t Current time t.
 *  @param lapse ADM lapse function.
 *  @param shift Shift vector \beta^i
 *  @param metric 3 dimensional ADM metric.
 *  @param curv Extrinsic curvature.
 *  @param dt Timestep.
 *  @param dx Spacestep
 *  @param lev AMR Level of discretization.
 *  @param plo Physical lower bounds of the whole domain.
 *  @return The right hind side of the differential equation.
 */
template <typename StructType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE CCTK_ATTRIBUTE_ALWAYS_INLINE
    amrex::GpuArray<CCTK_REAL, 8>
    PhotonsContainer<StructType>::compute_rhs(
        const amrex::GpuArray<CCTK_REAL, 8> &u, const CCTK_REAL &t,
        amrex::Array4<CCTK_REAL const> const &lapse,
        const amrex::Array4<CCTK_REAL const> &shift,
        const amrex::Array4<CCTK_REAL const> &metric,
        const amrex::Array4<CCTK_REAL const> &curv, 
        const amrex::Array4<CCTK_REAL const> &rho, const CCTK_REAL dt,
        const amrex::GpuArray<double, 3> &dx, const int lev,
        const amrex::GpuArray<double, 3> &plo) {

  amrex::GpuArray<CCTK_REAL, 8> rhs = {0., 0., 0., 0., 0., 0., 0., 0.};

  const long int i0 = amrex::Math::floor((u[0] - plo[0]) / dx[0]);
  const long int j0 = amrex::Math::floor((u[1] - plo[1]) / dx[1]);
  const long int k0 = amrex::Math::floor((u[2] - plo[2]) / dx[2]);

  // Interpolate lapse & partial lapse at \vect{x}
  CCTK_REAL lapse_x;
  amrex::GpuArray<CCTK_REAL, 3> d_lapse_x;
  d_interpolate_array<5>(lapse_x, d_lapse_x, lapse, i0, j0, k0, u[0], u[1],
                         u[2], dx, plo);

  // Interpolate shift & partial shift at \vect{x}
  amrex::GpuArray<CCTK_REAL, 3> shift_x;
  amrex::GpuArray<amrex::GpuArray<CCTK_REAL, 3>, 3> d_shift_x;
  d_interpolate_array<5>(shift_x, d_shift_x, shift, i0, j0, k0, u[0], u[1],
                         u[2], dx, plo);

  // Interpolate metric & partial metric at \vect{x}
  amrex::GpuArray<CCTK_REAL, 6> gamma_x;
  amrex::GpuArray<amrex::GpuArray<CCTK_REAL, 6>, 3> d_gamma_x;
  d_interpolate_array<5>(gamma_x, d_gamma_x, metric, i0, j0, k0, u[0], u[1],
                         u[2], dx, plo);

  // Interpolate Curvature at \vect{x}
  amrex::GpuArray<CCTK_REAL, 6> curv_x;
  interpolate_array<5>(curv_x, curv, i0, j0, k0, u[0], u[1], u[2], dx, plo);

  // Interpolate rho at \vect{x}
  amrex::GpuArray<CCTK_REAL, 1> rho_x;
  interpolate_array<5>(rho_x, rho, i0, j0, k0, u[0], u[1], u[2], dx, plo);

  // Compute the inverse of the metric.
  const CCTK_REAL inv_det_gamma =
      1.0 / (gamma_x[0] * gamma_x[3] * gamma_x[5] +
             2. * gamma_x[1] * gamma_x[2] * gamma_x[4] -
             gamma_x[2] * gamma_x[2] * gamma_x[3] -
             gamma_x[4] * gamma_x[4] * gamma_x[0] -
             gamma_x[1] * gamma_x[1] * gamma_x[5]);

  const amrex::GpuArray<CCTK_REAL, 6> gamma_inv_x = {
      (gamma_x[3] * gamma_x[5] - gamma_x[4] * gamma_x[4]) * inv_det_gamma,
      (gamma_x[4] * gamma_x[2] - gamma_x[1] * gamma_x[5]) * inv_det_gamma,
      (gamma_x[1] * gamma_x[4] - gamma_x[2] * gamma_x[3]) * inv_det_gamma,
      (gamma_x[0] * gamma_x[5] - gamma_x[2] * gamma_x[2]) * inv_det_gamma,
      (gamma_x[2] * gamma_x[1] - gamma_x[0] * gamma_x[4]) * inv_det_gamma,
      (gamma_x[0] * gamma_x[3] - gamma_x[1] * gamma_x[1]) * inv_det_gamma};

  const amrex::GpuArray<CCTK_REAL, 3> V_down = {u[3], u[4], u[5]};

  // Compute the upper index velocity terms.
  const amrex::GpuArray<CCTK_REAL, 3> V_up = {
      gamma_inv_x[0] * u[3] + gamma_inv_x[1] * u[4] + gamma_inv_x[2] * u[5],
      gamma_inv_x[1] * u[3] + gamma_inv_x[3] * u[4] + gamma_inv_x[4] * u[5],
      gamma_inv_x[2] * u[3] + gamma_inv_x[4] * u[4] + gamma_inv_x[5] * u[5]};

  // Compute the rhs for position
  rhs[0] = lapse_x * V_up[0] - shift_x[0];
  rhs[1] = lapse_x * V_up[1] - shift_x[1];
  rhs[2] = lapse_x * V_up[2] - shift_x[2];

  // Compute the rhs for velocity
  for (int i = 0; i < 3; i++) {
    rhs[3 + i] =
        -d_lapse_x[i] +
        (VecVecMul(d_lapse_x, V_up) -
         lapse_x * VecVecMul(SMatVecMul(curv_x, V_up), V_up)) *
            V_down[i] +
        0.5 * lapse_x * VecVecMul(SMatVecMul(d_gamma_x[i], V_up), V_up) +
        VecVecMul(V_down, d_shift_x[i]);
  }

  // Compute the rhs for energy
  rhs[3 + StructType::ln_E] =
      lapse_x * VecVecMul(SMatVecMul(curv_x, V_up), V_up) -
      VecVecMul(V_up, d_lapse_x);

  // Compute the rhs for optical depth
  const CCTK_REAL ds = dx[0] * dx[0] * gamma_inv_x[0] +
                                  dx[1] * dx[1] * gamma_inv_x[3] +
                                  dx[2] * dx[2] * gamma_inv_x[5] +
                                  2.0 * dx[0] * dx[1] * gamma_inv_x[1] +
                                  2.0 * dx[0] * dx[2] * gamma_inv_x[2] +
                                  2.0 * dx[1] * dx[2] * gamma_inv_x[4];
  rhs[3 + StructType::tau] = (0.4 * cgs2cactusOpacity) * (rho[0] * cgs2cactusrho) * (ds / dt);

  return rhs;

} // PhotonsContainer::compute_rhs

/**
 *  \brief Evolving using Runge-Kutta 4.
 *
 *  For the Runge-Kutta 4 we are solving the differential equation
 * \f$\frac{dU}{dt} = f\left(U, \frac{dU}{dx}, t\right)\f$ using:
 *
 *  \f[
 *  U_{n+1} = U_n + \frac{1}{6}\Delta t \left(f_1 + 2f_2 + 2f_3 + f_4\right)
 *  \f]
 *
 *  where:
 *
 *  * \f$f_1 = f(U_n, t),\f$
 *  * \f$f_2 = f\left(U_n + \frac{\Delta t}{2} f_1, t + \frac{\Delta
 * t}{2}\right),\f$
 *  * \f$f_3 = f\left(U_n + \frac{\Delta t}{2} f_2, t + \frac{\Delta
 * t}{2}\right),\f$
 *  * \f$f_4 = f(U_n + \Delta t f_3, t + \Delta t),\f$
 *
 *  And checking if it goes outside of the grid boundaries.
 *
 *  @see compute_rhs()
 *  @param lapse ADM lapse function.
 *  @param shift ADM shift vector.
 *  @param metric ADM 3 dimension metric.
 *  @param curv Extrinsic curvature.
 *  @param dt Timestep.
 *  @param lev AMR Level of discretization.
 */
template <typename StructType>
void PhotonsContainer<StructType>::evolve(const amrex::MultiFab &lapse,
                                          const amrex::MultiFab &shift,
                                          const amrex::MultiFab &metric,
                                          const amrex::MultiFab &curv,
                                          const CCTK_REAL &dt, const int &lev) {

  const auto plo0 = this->Geom(0).ProbLoArray();
  const auto phi0 = this->Geom(0).ProbHiArray();

  const auto dx = this->Geom(lev).CellSizeArray();
  const auto plo = this->Geom(lev).ProbLoArray();

  const CCTK_REAL boundarie_hx = phi0[0] - 0.0 * dx[0];
  const CCTK_REAL boundarie_lx = plo0[0] + 0.0 * dx[0];
  const CCTK_REAL boundarie_hy = phi0[1] - 0.0 * dx[1];
  const CCTK_REAL boundarie_ly = plo0[1] + 0.0 * dx[1];
  const CCTK_REAL boundarie_hz = phi0[2] - 0.0 * dx[2];
  const CCTK_REAL boundarie_lz = plo0[2] + 0.0 * dx[2];

  for (Iterator::ParticleIterator<StructType> pti(*this, lev); pti.isValid();
       ++pti) {

    const int np = pti.numParticles();

    // Get the information relate to the velocities and energy.
    auto &attribs = pti.GetAttributes();
    CCTK_REAL *AMREX_RESTRICT vels_x = attribs[StructType::vx].data();
    CCTK_REAL *AMREX_RESTRICT vels_y = attribs[StructType::vy].data();
    CCTK_REAL *AMREX_RESTRICT vels_z = attribs[StructType::vz].data();
    CCTK_REAL *AMREX_RESTRICT ln_energy = attribs[StructType::ln_E].data();
    CCTK_REAL *AMREX_RESTRICT tau = attribs[StructType::tau].data();
    auto *AMREX_RESTRICT particles = &(pti.GetArrayOfStructs()[0]);

    // Get the array of each parameter.
    auto const lapse_array = lapse.array(pti);
    auto const shift_array = shift.array(pti);
    auto const metric_array = metric.array(pti);
    auto const curv_array = curv.array(pti);
    auto const curv_array = rho.array(pti);

    // Needed for GPU
    auto self = this;

    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int i) noexcept {
      const amrex::GpuArray<CCTK_REAL, 8> U = {
          particles[i].pos(0), particles[i].pos(1), particles[i].pos(2),
          vels_x[i],           vels_y[i],           vels_z[i],
          ln_energy[i], tau[i]};

      bool out_of_bounds = false;

      amrex::GpuArray<CCTK_REAL, 8> U_tmp = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      // f1 = rhs(u , t) for the runge kutta 4 step
      auto k_odd =
          self->compute_rhs(U, 0.0, lapse_array, shift_array, metric_array,
                            curv_array, rho_array, dt, dx, lev, plo0);

      U_tmp[0] = U[0] + 0.5 * dt * k_odd[0];
      U_tmp[1] = U[1] + 0.5 * dt * k_odd[1];
      U_tmp[2] = U[2] + 0.5 * dt * k_odd[2];
      U_tmp[3] = U[3] + 0.5 * dt * k_odd[3];
      U_tmp[4] = U[4] + 0.5 * dt * k_odd[4];
      U_tmp[5] = U[5] + 0.5 * dt * k_odd[5];
      U_tmp[6] = U[6] + 0.5 * dt * k_odd[6];
      U_tmp[7] = U[7] + 0.5 * dt * k_odd[7];
      
      CHECK_OUT_OF_BOUNDS(U_tmp[,])

      if (out_of_bounds) {
        particles[i].id() = -1;
        return;
      }

      // f2 = rhs(u + 0.5 * dt * f1, t) for the runge kutta 4 step
      auto k_even =
          self->compute_rhs(U_tmp, 0.5 * dt, lapse_array, shift_array,
                            metric_array, curv_array, rho_array, dt, dx, lev, plo0);

      // Update particles with the f1 and f2 from RK4
      U_tmp[0] = U[0] + 0.5 * dt * k_even[0];
      U_tmp[1] = U[1] + 0.5 * dt * k_even[1];
      U_tmp[2] = U[2] + 0.5 * dt * k_even[2];
      U_tmp[3] = U[3] + 0.5 * dt * k_even[3];
      U_tmp[4] = U[4] + 0.5 * dt * k_even[4];
      U_tmp[5] = U[5] + 0.5 * dt * k_even[5];
      U_tmp[6] = U[6] + 0.5 * dt * k_even[6];
      U_tmp[7] = U[7] + 0.5 * dt * k_even[7];

      particles[i].pos(0) += (1. / 6.) * dt * (k_odd[0] + 2. * k_even[0]);
      particles[i].pos(1) += (1. / 6.) * dt * (k_odd[1] + 2. * k_even[1]);
      particles[i].pos(2) += (1. / 6.) * dt * (k_odd[2] + 2. * k_even[2]);
      vels_x[i] += (1. / 6.) * dt * (k_odd[3] + 2. * k_even[3]);
      vels_y[i] += (1. / 6.) * dt * (k_odd[4] + 2. * k_even[4]);
      vels_z[i] += (1. / 6.) * dt * (k_odd[5] + 2. * k_even[5]);
      ln_energy[i] += (1. / 6.) * dt * (k_odd[6] + 2. * k_even[6]);
      tau[i] += (1. / 6.) * dt * (k_odd[7] + 2. * k_even[7]);

      CHECK_OUT_OF_BOUNDS(U_tmp[,])

      if (out_of_bounds) {
        particles[i].id() = -1;
        return;
      }

      // f3 = rhs(u + 0.5 * dt * f2, t) for the runge kutta 4 step
      k_odd = self->compute_rhs(U_tmp, 0.5 * dt, lapse_array, shift_array,
                                metric_array, curv_array, dt, dx, lev, plo0);

      U_tmp[0] = U[0] + dt * k_odd[0];
      U_tmp[1] = U[1] + dt * k_odd[1];
      U_tmp[2] = U[2] + dt * k_odd[2];
      U_tmp[3] = U[3] + dt * k_odd[3];
      U_tmp[4] = U[4] + dt * k_odd[4];
      U_tmp[5] = U[5] + dt * k_odd[5];
      U_tmp[6] = U[6] + dt * k_odd[6];

      CHECK_OUT_OF_BOUNDS(U_tmp[,])

      if (out_of_bounds) {
        particles[i].id() = -1;
        return;
      }

      // f4 = rhs(u + dt * f3, t) for the runge kutta 4 step
      k_even = self->compute_rhs(U_tmp, dt, lapse_array, shift_array,
                                 metric_array, curv_array, dt, dx, lev, plo0);

      // Update particles with the f3 and f4 from RK4
      particles[i].pos(0) += (1. / 6.) * dt * (2. * k_odd[0] + k_even[0]);
      particles[i].pos(1) += (1. / 6.) * dt * (2. * k_odd[1] + k_even[1]);
      particles[i].pos(2) += (1. / 6.) * dt * (2. * k_odd[2] + k_even[2]);
      vels_x[i] += (1. / 6.) * dt * (2. * k_odd[3] + k_even[3]);
      vels_y[i] += (1. / 6.) * dt * (2. * k_odd[4] + k_even[4]);
      vels_z[i] += (1. / 6.) * dt * (2. * k_odd[5] + k_even[5]);
      ln_energy[i] += (1. / 6.) * dt * (2. * k_odd[6] + k_even[6]);
      tau[i] += (1. / 6.) * dt * (2. * k_odd[7] + k_even[7]);

      CHECK_OUT_OF_BOUNDS(particles[i].pos(,))
      out_of_bounds |= (tau[i] > 1.);

      if (out_of_bounds) {
        particles[i].id() = -1;
        return;
      }
    });
  }
} // PhotonsContainer::evolve

  /**
   * The check banned zones function check for user defined invalid particles
   * zones.
   *
   * @param level Adaptive Mesh Refinement level
   * @param zones Number of banned zones
   * @param x x-coordinates array for each region
   * @param y y-coordinates array for each region
   * @param z z-coordinates array for each region
   * @param radius Radius array for each region
   */
  void check_banned_zones(const int &level, const CCTK_INT4 &zones,
                          const CCTK_REAL (&x)[10], const CCTK_REAL (&y)[10],
                          const CCTK_REAL (&z)[10],
                          const CCTK_REAL (&radius)[10]) {

    if (!zones) {
      return;
    }

    auto &attribs = pti.GetAttributes();
    CCTK_REAL *AMREX_RESTRICT tau = attribs[StructType::tau].data();

    for (Iterator::ParticleIterator<StructType> pti(*this, level);
         pti.isValid(); ++pti) {
      const int np = pti.numParticles();
      auto *AMREX_RESTRICT particles = &(pti.GetArrayOfStructs()[0]);

      auto self = this;
      amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int i) noexcept {
        for (int check = 0; check < zones; check++) {
          const CCTK_REAL dx = particles[i].pos(0) - x[check];
          const CCTK_REAL dy = particles[i].pos(1) - y[check];
          const CCTK_REAL dz = particles[i].pos(2) - z[check];
          if (dx * dx + dy * dy + dz * dz <= radius[check] * radius[check]) {
            particles[i].id() = -1;
            tau[i] = -check - 7;
          }
        }
      });
    }
  }

}
}

#endif