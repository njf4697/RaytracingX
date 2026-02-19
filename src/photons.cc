#include <cctk.h>

#include "Photons.hxx"
#include "PhotonsContainer.hxx"
#include "raytracingx.h"
#include <AMReX_ParallelDescriptor.H>
#include <CParameters.h>

#include <cstring>
#include <driver.hxx>

#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Types.h>
#include <AMReX_MFIter.H>
#include <cctk_core.h>

#include <iostream>
#include <loop_device.hxx>

using ParticleData = Photons::PhotonsData;
using PC = Containers::PhotonsContainer<ParticleData>;
std::vector<std::unique_ptr<PC>> r_photons;

extern "C" void R_PhotonsContainer_setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int tl = 0;
  const int gi_metric = CCTK_GroupIndex("ADMBaseX::metric");
  assert(gi_metric >= 0 && "Failed to get the metric group index");

  CCTK_INT int_params[2];
  CCTK_REAL real_params[32];
  setup_camera_initializer_reals(CCTK_PASS_CTOC, real_params);
  setup_camera_initializer_ints(CCTK_PASS_CTOC, int_params);

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    if (r_photons.size() < CarpetX::ghext->num_patches()) {
      r_photons.push_back(std::make_unique<PC>(patchdata.amrcore.get()));

      auto &pc = r_photons.at(patch);
      pc->initialize(camera_initializer<ParticleData, PC>,
                     real_params, int_params);
    }
  }

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    auto &pc = r_photons.at(patch);
    auto &pd = CarpetX::ghext->patchdata.at(patch);
    for (int lev = 0; lev < pd.leveldata.size(); ++lev) {

      const auto &ld = pd.leveldata.at(lev);
      const auto &gd_metric = *ld.groupdata.at(gi_metric);
      const amrex::MultiFab &metric = *gd_metric.mfab[tl];

      pc->Redistribute();
      pc->normalize_velocity(metric, lev);
    }
  }
}

extern "C" void R_PhotonsContainer_evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const int tl = 0;
  const int gi_lapse = CCTK_GroupIndex("ADMBaseX::lapse");
  const int gi_shift = CCTK_GroupIndex("ADMBaseX::shift");
  const int gi_metric = CCTK_GroupIndex("ADMBaseX::metric");
  const int gi_curv = CCTK_GroupIndex("ADMBaseX::curv");
  assert(gi_lapse >= 0 && "Failed to get the lapse group index");
  assert(gi_shift >= 0 && "Failed to get the shift group index");
  assert(gi_metric >= 0 && "Failed to get the metric group index");
  assert(gi_curv >= 0 && "Failed to get the curvature group index");

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    auto &pc = r_photons.at(patch);
    auto &pd = CarpetX::ghext->patchdata.at(patch);
    for (int lev = 0; lev < pd.leveldata.size(); ++lev) {
      const auto &ld = pd.leveldata.at(lev);
      const auto &gd_lapse = *ld.groupdata.at(gi_lapse);
      const auto &gd_shift = *ld.groupdata.at(gi_shift);
      const auto &gd_metric = *ld.groupdata.at(gi_metric);
      const auto &gd_curv = *ld.groupdata.at(gi_curv);
      const amrex::MultiFab &lapse = *gd_lapse.mfab[tl];
      const amrex::MultiFab &shift = *gd_shift.mfab[tl];
      const amrex::MultiFab &metric = *gd_metric.mfab[tl];
      const amrex::MultiFab &curv = *gd_curv.mfab[tl];

      pc->evolve(lapse, shift, metric, curv, CCTK_DELTA_TIME, lev);
    }
  }

  // Bounds check
  const CCTK_REAL regions_x[10] = {region_1_position[0], region_2_position[0],
                                   region_3_position[0], region_4_position[0],
                                   region_5_position[0], region_6_position[0],
                                   region_7_position[0], region_8_position[0],
                                   region_9_position[0], region_10_position[0]};
  const CCTK_REAL regions_y[10] = {region_1_position[1], region_2_position[1],
                                   region_3_position[1], region_4_position[1],
                                   region_5_position[1], region_6_position[1],
                                   region_7_position[1], region_8_position[1],
                                   region_9_position[1], region_10_position[1]};
  const CCTK_REAL regions_z[10] = {region_1_position[2], region_2_position[2],
                                   region_3_position[2], region_4_position[2],
                                   region_5_position[2], region_6_position[2],
                                   region_7_position[2], region_8_position[2],
                                   region_9_position[2], region_10_position[2]};
  const CCTK_REAL regions_radius[10] = {
      region_1_radius, region_2_radius, region_3_radius, region_4_radius,
      region_5_radius, region_6_radius, region_7_radius, region_8_radius,
      region_9_radius, region_10_radius};

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    auto &pc = r_photons.at(patch);
    auto &pd = CarpetX::ghext->patchdata.at(patch);
    for (int lev = 0; (lev < pd.leveldata.size()) & banned_regions; ++lev) {
      pc->check_banned_zones(lev, banned_regions, regions_x, regions_y,
                             regions_z, regions_radius);
    }
    pc->Redistribute();
  }
}

extern "C" void R_PhotonsContainer_print(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Printing particles to files");

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    auto &pc = r_photons.at(patch);
    pc->outputParticlesPlot(CCTK_PASS_CTOC, particle_plot_every,
                            std::string(out_dir));
    pc->outputParticlesAscii(CCTK_PASS_CTOC, particle_tsv_every,
                             std::string(out_dir));
  }
}

extern "C" int R_PhotonsContainer_final_cleanup() {
  amrex::Gpu::Device::synchronize();
  r_photons.clear();
  return 0;
}

template <typename StructType, typename ParticleContainerClass>
void camera_initializer(ParticleContainerClass &pc, const CCTK_REAL *real_params, const CCTK_INT *int_params) {
    CCTK_INFO("Initializing particles using the RaytracingX::camera_initializer");
    
    CCTK_REAL e0[4], e1[4], e2[4], e3[4];
    e0[0] = real_params[10];
    e0[1] = real_params[11];
    e0[2] = real_params[12];
    e0[3] = real_params[13];
    e1[0] = real_params[14];
    e1[1] = real_params[15];
    e1[2] = real_params[16];
    e1[3] = real_params[17];
    e2[0] = real_params[18];
    e2[1] = real_params[19];
    e2[2] = real_params[20];
    e2[3] = real_params[21];
    e3[0] = real_params[22];
    e3[1] = real_params[23];
    e3[2] = real_params[24];
    e3[3] = real_params[25];
    CCTK_REAL alpha_h = real_params[26];
    CCTK_REAL alpha_v = real_params[27];
    CCTK_REAL lapse = real_params[28];
    CCTK_INT num_pixels_width = int_params[0];
    CCTK_INT num_pixels_height = int_params[1];
    CCTK_REAL camera_pos[3];
    camera_pos[0] = real_params[29];
    camera_pos[1] = real_params[30];
    camera_pos[2] = real_params[31];

    const CCTK_INT level = 0;
    const CCTK_INT num_pixels = num_pixels_width * num_pixels_height;

    const int n_procs = amrex::ParallelDescriptor::NProcs();
    const int proc_id = amrex::ParallelDescriptor::MyProc();

    const int local_particles_size = num_pixels / n_procs + (proc_id < num_pixels % n_procs);
    const int local_offset = proc_id * num_pixels / n_procs + std::min(proc_id, num_pixels % n_procs);

    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, ("proc: " + std::to_string(proc_id) + " local size: " + std::to_string(local_particles_size) + " local offset: " + std::to_string(local_offset) + " level: " + std::string(level) + " mfi: " + std::string(mfi) + "\n").c_str());

    // Iterating over all the tiles of the particle data structure
    for (amrex::MFIter mfi = pc.MakeMFIter(level); mfi.isValid(); ++mfi) {

        auto &particles = pc.GetParticles(level);
        auto &particle_tile = pc.DefineAndReturnParticleTile(level, mfi);
        assert(particle_tile.GetArrayOfStructs().size() == 0);
        particle_tile.resize(local_particles_size);
        auto arrdata = particle_tile.GetStructOfArrays().realarray();
        auto ptd = particle_tile.getParticleTileData();

        #pragma omp parallel for
        for (int pidx = local_offset; pidx < local_offset + local_particles_size; pidx++) { //create 4-vector \chi parallel to geodesic and fill geodesic initial conditions for each pixel (see https://arxiv.org/pdf/1410.777)
                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, ("proc: " + std::to_string(proc_id) + " pidx: " + std::to_string(pidx) + "\n").c_str());

                int i = pidx / num_pixels_width;
                int j = pidx % num_pixels_width;

                CCTK_REAL a_adj = (2.0 * i / num_pixels_width - 1)*tan(alpha_h / 2.0); // a_{adj} = (2a-1)tan(\alpha_h/2)
                CCTK_REAL b_adj = (2.0 * j / num_pixels_height - 1)*tan(alpha_v / 2.0); // b_{adj} = (2b-1)tan(\alpha_v/2)

                CCTK_REAL C = sqrt(1 + pow(b_adj,2) + pow(a_adj,2));

                CCTK_REAL chi[4];
                chi[0] = C*e0[0] - e1[0] - b_adj*e2[0] - a_adj*e3[0];
                chi[1] = C*e0[1] - e1[1] - b_adj*e2[1] - a_adj*e3[1];
                chi[2] = C*e0[2] - e1[2] - b_adj*e2[2] - a_adj*e3[2];
                chi[3] = C*e0[3] - e1[3] - b_adj*e2[3] - a_adj*e3[3];

                printf("i=%i, j=%i, chi=[%0.2f, %0.2f, %0.2f, %0.2f]\n",i,j,chi[0],chi[1],chi[2],chi[3]);

                CCTK_REAL chi_lower[4];
                vectorToOneFormArr(chi_lower, chi, real_params);

                ptd.id(pidx) = ParticleContainerClass::ParticleType::NextID();
                ptd.cpu(pidx) = amrex::ParallelDescriptor::MyProc();

                ptd.pos(0, pidx) = camera_pos[0];
                ptd.pos(1, pidx) = camera_pos[1];
                ptd.pos(2, pidx) = camera_pos[2];
                CCTK_REAL A = 1 / lapse*chi[0];
                arrdata[StructType::vx][pidx] = chi_lower[0] * A;
                arrdata[StructType::vy][pidx] = chi_lower[1] * A;
                arrdata[StructType::vz][pidx] = chi_lower[2] * A;
                arrdata[StructType::ln_E][pidx] = 0;
        }
    }
    pc.Redistribute();
    pc.SortParticlesByCell();
    CCTK_VINFO("%d particles created", pc.TotalNumberOfParticles()); 
}
