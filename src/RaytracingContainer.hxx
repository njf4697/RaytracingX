#include <cctk.h>

#include "Photons.hxx"
#include "PhotonsContainer.hxx"
#include "raytracingx.h"
#include <AMReX_ParallelDescriptor.H>
#include <CParameters.h>

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
}

}

//namespace RaytracingContainers {
//
//    template <typename StructType>
//    class RaytracingPhotonsContainer : public ConcreteContainer<PhotonsContainer<StructType>, StructType> {
//
//
//
//    }
//}