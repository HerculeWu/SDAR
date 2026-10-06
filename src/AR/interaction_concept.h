#pragma once

#include <concepts>
#include "Common/Float.h"
#include "AR/variant.h"
#include "AR/force.h"

namespace AR {

    //! Force evaluation every time transformation requires of the interaction class, with the force type of the transformation
    /*! calcAccPotAndGTKickInv: forces and potential of all members with perturbation, returns the inverse time transformation factor for the kick \n
        calcInnerAccPotAndGTKickInvTwo: the same for one pair without perturbation
     */
    template <class Tmethod, class Tforce, class Tparticle, class Tpcm, class Tpert>
    concept ForceInteraction = requires(Tmethod _method, Tforce* _force, Tforce& _f1, Tforce& _f2, Float& _epot, Tparticle* _particles, Tparticle& _p1, Tparticle& _p2,
                                        const int _n_particle, Tpcm& _particle_cm, Tpert& _perturber, const Float _time) {
        { _method.calcAccPotAndGTKickInv(_force, _epot, _particles, _n_particle, _particle_cm, _perturber, _time) } -> std::convertible_to<Float>;
        { _method.calcInnerAccPotAndGTKickInvTwo(_f1, _f2, _epot, _p1, _p2) } -> std::convertible_to<Float>;
    };

    //! What LogH requires of the interaction class
    /*! Force evaluation with ForceLogH, and the time transformation formulas of the logarithmic Hamiltonian: \n
        calcGTDriftInv: inverse time transformation factor for the drift from (kinetic energy - total energy) \n
        calcH: extended Hamiltonian from (kinetic energy - total energy) and the potential energy
     */
    template <class Tmethod, class Tparticle, class Tpcm, class Tpert>
    concept LogHInteraction = ForceInteraction<Tmethod, ForceLogH, Tparticle, Tpcm, Tpert> &&
        requires(Tmethod _method, const Float _ekin_minus_etot, const Float _epot) {
            { _method.calcGTDriftInv(_ekin_minus_etot) } -> std::convertible_to<Float>;
            { _method.calcH(_ekin_minus_etot, _epot) } -> std::convertible_to<Float>;
        };

    //! What TTL requires of the interaction class
    /*! Force evaluation with ForceTTL, which also fills the gradient of the time transformation function (gtgrad)
     */
    template <class Tmethod, class Tparticle, class Tpcm, class Tpert>
    concept TTLInteraction = ForceInteraction<Tmethod, ForceTTL, Tparticle, Tpcm, Tpert>;

    //! The interaction class provides what the chosen time transformation requires
    template <class Tmethod, class Ttransform, class Tparticle, class Tpcm, class Tpert>
    concept InteractionOf = (std::same_as<Ttransform, LogH> && LogHInteraction<Tmethod, Tparticle, Tpcm, Tpert>) || 
                            (std::same_as<Ttransform, TTL> && TTLInteraction<Tmethod, Tparticle, Tpcm, Tpert>);

}
