// Must not compile: LogH needs calcGTDriftInv and calcH from the interaction class
// expect: LogHInteraction
#include "test_support.h"
#include "ar_test_system.h"

// force evaluation only
class ForceOnlyInteraction {
public:
    Float gravitational_constant = 1.0;
    template <class Tforce>
    Float calcAccPotAndGTKickInv(Tforce*, Float&, const Particle*, const int, const Particle&, const Perturber&, const Float) { return 1.0; }
    template <class Tforce>
    Float calcInnerAccPotAndGTKickInvTwo(Tforce&, Tforce&, Float&, const Particle&, const Particle&) { return 1.0; }
};

AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, ForceOnlyInteraction, AR::Information<Particle,Particle>, AR::LogH, AR::NoSlowDown> integrator;
