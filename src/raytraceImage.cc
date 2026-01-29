#include "raytracingx.h"
#include <stdio.h>
#include <AMReX_MFIter.H>
#include "AMReX_ParallelDescriptor.H"

#define DEG2RAD 0.01745329251

void gramSchmidtProcess(CCTK_ARGUMENTS, CCTK_REAL* e0, CCTK_REAL* e1, CCTK_REAL* e2, CCTK_REAL* e3, const Metric* metric) { //use Gram-Schmidt Process to turn e0, e1, e2, and e3 into orthonormal vectors for the camera's POV (see https://arxiv.org/pdf/1410.7775)
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

void createGeodesicInitialConditions(CCTK_ARGUMENTS, GeodesicInitialConditions* geodesicArr) { //fills geodesicArr with geodesic initial conditions for each pixel (see https://arxiv.org/pdf/1410.7775)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    Metric metric;
    interpolateMetricAtPoint(CCTK_PASS_CTOC, camera_pos[0], camera_pos[1], camera_pos[2], &metric); //interpolate metric values and store in Metric struct

    //only use values from processer 1
    if (CCTK_MyProc(cctkGH) != 0) return; 

    CCTK_REAL e0[4];
    CCTK_REAL e1[4];
    CCTK_REAL e2[4];
    CCTK_REAL e3[4];
    gramSchmidtProcess(CCTK_PASS_CTOC, e0, e1, e2, e3, &metric); //create orthonormal basis for camera POV

    printf(("e0: (" + std::to_string(e0[0]) + ", " + std::to_string(e0[1]) + ", " + std::to_string(e0[2]) + ", " + std::to_string(e0[3]) + ")\n").c_str());
    printf(("e1: (" + std::to_string(e1[0]) + ", " + std::to_string(e1[1]) + ", " + std::to_string(e1[2]) + ", " + std::to_string(e1[3]) + ")\n").c_str());
    printf(("e2: (" + std::to_string(e2[0]) + ", " + std::to_string(e2[1]) + ", " + std::to_string(e2[2]) + ", " + std::to_string(e2[3]) + ")\n").c_str());
    printf(("e3: (" + std::to_string(e3[0]) + ", " + std::to_string(e3[1]) + ", " + std::to_string(e3[2]) + ", " + std::to_string(e3[3]) + ")\n").c_str());
    
    CCTK_REAL alpha_h = 3.1415926536 / 180 * horizontal_fov; //convert FOV to radians
    CCTK_REAL alpha_v = 3.1415926536 / 180 * vertical_fov;

    //TODO: make parallel for GPU
    #pragma omp parallel for
    for (int i = 0; i < num_pixels_width; i++) {
        for (int j = 0; j < num_pixels_height; j++) { //create 4-vector \chi parallel to geodesic and fill GeodesicInitialConditions struct for each pixel (see https://arxiv.org/pdf/1410.777)
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
            vectorToOneForm(chi_lower, chi, &metric);
            geodesicArr[i*num_pixels_width + j].initPos[0] = camera_pos[0]; 
            geodesicArr[i*num_pixels_width + j].initPos[1] = camera_pos[1]; 
            geodesicArr[i*num_pixels_width + j].initPos[2] = camera_pos[2]; 
            geodesicArr[i*num_pixels_width + j].initVel[0] = chi_lower[0] / (metric.alpha*chi[0]); 
            geodesicArr[i*num_pixels_width + j].initVel[1] = chi_lower[1] / (metric.alpha*chi[0]); 
            geodesicArr[i*num_pixels_width + j].initVel[2] = chi_lower[2] / (metric.alpha*chi[0]);

            printf(("init vel for i=" + std::to_string(i) + " and j=" + std::to_string(j) + ": (" + std::to_string(geodesicArr[i*num_pixels_width + j].initVel[0]) + ", " + std::to_string(geodesicArr[i*num_pixels_width + j].initVel[1]) + ", " + std::to_string(geodesicArr[i*num_pixels_width + j].initVel[2]) + ")\n").c_str());
        }
    }
}

template <typename StructType, typename ParticleContainerClass>
void camera_initializer(ParticleContainerClass &pc, const CCTK_REAL *real_params, const CCTK_INT *int_params) {
    CCTK_INFO("Initializing particles using the RaytracingX::camera_initializer");
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL e0[4];
    CCTK_REAL e1[4];
    CCTK_REAL e2[4];
    CCTK_REAL e3[4];
    gramSchmidtProcess(CCTK_PASS_CTOC, e0, e1, e2, e3, &metric); //create orthonormal basis for camera POV

    CCTK_REAL alpha_h = DEG2RAD * horizontal_fov; //convert FOV to radians
    CCTK_REAL alpha_v = DEG2RAD * vertical_fov;

    const CCTK_INT level = 0;
    const CCTK_INT num_pixels = num_pixels_width * num_pixels_height;

    int iteration = 0;

    // Iterating over all the tiles of the particle data structure
    for (amrex::MFIter mfi = pc.MakeMFIter(level); mfi.isValid(); ++mfi) {
        assert(iteration == 0);

        auto &particles = pc.GetParticles(level);
        auto &particle_tile = pc.DefineAndReturnParticleTile(level, mfi);
        assert(particle_tile.GetArrayOfStructs().size() == 0);
        particle_tile.resize(num_pixels);
        auto arrdata = particle_tile.GetStructOfArrays().realarray();
        auto ptd = particle_tile.GetParticleTileData();

        #pragma omp parallel for
        for (int i = 0; i < num_pixels_width; i++) {
            for (int j = 0; j < num_pixels_height; j++) { //create 4-vector \chi parallel to geodesic and fill GeodesicInitialConditions struct for each pixel (see https://arxiv.org/pdf/1410.777)

                int pidx = i*num_pixels_width + j;

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
                vectorToOneForm(chi_lower, chi, &metric);

                ptd.id(pidx) = ParticleContainerClass::ParticleType::NextID();
                ptd.cpu(pidx) = amrex::ParallelDescriptor::MyProc();

                ptd.pos(0, pidx) = camera_pos[0];
                ptd.pos(1, pidx) = camera_pos[1];
                ptd.pos(2, pidx) = camera_pos[2];
                CCTK_REAL A = 1 / metric.alpha*chi[0];
                arrdata[StructType::vx][pidx] = chi_lower[0] * A;
                arrdata[StructType::vy][pidx] = chi_lower[1] * A;
                arrdata[StructType::vz][pidx] = chi_lower[2] * A;
                arrdata[StructType::ln_E][pidx] = 0;
            }   
        
            pc.Redistribute();
            pc.SortParticlesByCell(); 

            CCTK_VINFO("%d particles created", pc.TotalNumberOfParticles());    
        }

        iteration++;
    }
}