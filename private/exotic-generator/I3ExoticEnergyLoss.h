#ifndef I3EXOTICENERGYLOSSUTILS_H
#define I3EXOTICENERGYLOSSUTILS_H

namespace I3ExoticEnergyLoss {
    double CalculateMonopoleEnergyLoss(double, double, double);
    double BetheBlock(double, double,bool);
    double ELEC_LOSS(double, double);
    double NUCL_LOSS(double, double, int, double);
    double SlowIonInIce(double, double);
}

#endif // I3EXOTICENERGYLOSSUTILS_H
