#ifndef EXOTIC_PROPAGATOR_EXOTICPROPAGATORUTILS_H
#define EXOTIC_PROPAGATOR_EXOTICPROPAGATORUTILS_H

#include <phys-services/I3RandomService.h>
#include <dataclasses/physics/I3MCTreeUtils.h>

namespace I3ExoticPropagatorUtils {

    void PropagateExotic(
            I3FramePtr,
            I3MCTreeConstPtr,
            I3MCTreePtr,
            I3RandomServicePtr,
            // params
            double,
            bool,
            double,
            double,
            bool,
            bool,
            double,
            double,
            bool,
            double,
            double,
            double,
            bool);

    I3ParticlePtr AddExotic(I3Particle, 
                            double,
                            double,
                            bool,
                            double,
                            double,
                            double,
                            bool);

    
    double CalculateNextLength(double, double, double, double, double);

    double EnergyToBeta(double, double);

    void CheckParticle(I3Particle &, double);

    void ExoSpeedProfile(I3MCTreeConstPtr, I3VectorDoublePtr);
}

#endif //EXOTIC_PROPAGATOR_EXOTICPROPAGATORUTILS_H

