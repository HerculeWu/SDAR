#pragma once
// A small few-body system for testing the AR integrator through its public interface,
// with the variant (time transformation and Slow-down scheme) chosen by type.
#include <iomanip>
#include <vector>
#include "Common/binary_tree.h"
#include "AR/symplectic_integrator.h"
#include "AR/information.h"
#include "particle.h"  // sample/AR
#include "perturber.h" // sample/AR

//! Newtonian gravity of isolated point masses, generic over the variant
/*! The force type tells the time transformation: only the TTL force carries the gradient of the time transformation function.
 */
class TestInteraction {
public:
    Float gravitational_constant = 1.0;

    bool checkParams() { return true; }

    void print(std::ostream&) const {}

    template <class Tforce>
    Float calcInnerAccPotAndGTKickInvTwo(Tforce& _f1, Tforce& _f2, Float& _epot, const Particle& _p1, const Particle& _p2) {
        Float dr[3] = {_p2.pos[0]-_p1.pos[0], _p2.pos[1]-_p1.pos[1], _p2.pos[2]-_p1.pos[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float inv_r = 1.0/sqrt(r2);
        Float inv_r3 = inv_r*inv_r*inv_r;
        Float gm1m2 = gravitational_constant*_p1.mass*_p2.mass;
        for (int k=0; k<3; k++) {
            _f1.acc_in[k] =  gravitational_constant*_p2.mass*inv_r3*dr[k];
            _f2.acc_in[k] = -gravitational_constant*_p1.mass*inv_r3*dr[k];
        }
        _f1.pot_in = -gravitational_constant*_p2.mass*inv_r;
        _f2.pot_in = -gravitational_constant*_p1.mass*inv_r;
        if constexpr (requires { _f1.gtgrad; }) {
            for (int k=0; k<3; k++) {
                _f1.gtgrad[k] = gm1m2*inv_r3*dr[k];
                _f2.gtgrad[k] = -_f1.gtgrad[k];
            }
        }
        _epot = -gm1m2*inv_r;
        return gm1m2*inv_r;
    }

    template <class Tforce>
    Float calcInnerAccPotAndGTKickInv(Tforce* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
        _epot = 0.0;
        Float gt_kick_inv = 0.0;
        for (int i=0; i<_n_particle; i++) {
            Float* acci = _force[i].acc_in;
            acci[0] = acci[1] = acci[2] = 0.0;
            if constexpr (requires { _force->gtgrad; }) 
                _force[i].gtgrad[0] = _force[i].gtgrad[1] = _force[i].gtgrad[2] = 0.0;
            Float poti = 0.0;
            for (int j=0; j<_n_particle; j++) {
                if (i==j) continue;
                Float dr[3] = {_particles[j].pos[0]-_particles[i].pos[0], _particles[j].pos[1]-_particles[i].pos[1], _particles[j].pos[2]-_particles[i].pos[2]};
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float inv_r = 1.0/sqrt(r2);
                Float gmor3 = gravitational_constant*_particles[j].mass*inv_r*inv_r*inv_r;
                for (int k=0; k<3; k++) acci[k] += gmor3*dr[k];
                if constexpr (requires { _force->gtgrad; }) 
                    for (int k=0; k<3; k++) _force[i].gtgrad[k] += gravitational_constant*_particles[i].mass*gmor3*dr[k];
                poti -= gravitational_constant*_particles[j].mass*inv_r;
            }
            _epot += poti*_particles[i].mass;
            gt_kick_inv -= poti*_particles[i].mass;
        }
        _epot *= 0.5;
        return 0.5*gt_kick_inv;
    }

    template <class Tforce>
    void calcAccPert(Tforce* _force, const Particle*, const int _n_particle, const Particle&, const Perturber&, const Float) {
        for (int i=0; i<_n_particle; i++) {
            _force[i].acc_pert[0] = _force[i].acc_pert[1] = _force[i].acc_pert[2] = 0.0;
            _force[i].pot_pert = 0.0;
        }
    }

    template <class Tforce>
    Float calcAccPotAndGTKickInv(Tforce* _force, Float& _epot, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const Perturber& _perturber, const Float _time) {
        Float gt_kick_inv;
        if (_n_particle==2) gt_kick_inv = calcInnerAccPotAndGTKickInvTwo(_force[0], _force[1], _epot, _particles[0], _particles[1]);
        else gt_kick_inv = calcInnerAccPotAndGTKickInv(_force, _epot, _particles, _n_particle);
        calcAccPert(_force, _particles, _n_particle, _particle_cm, _perturber, _time);
        return gt_kick_inv;
    }

    static Float calcPertFromBinary(const AR::BinaryTree<Particle>& _bin) {
        Float apo = _bin.semi*(1.0+_bin.ecc);
        return (_bin.m1*_bin.m2)/(apo*apo*apo);
    }

    static Float calcPertFromMR(const Float _r, const Float _mp, const Float _mpert) {
        return (_mp*_mpert)/(_r*_r*_r);
    }

    void calcSlowDownPertOne(Float& _pert_out, Float&, const Particle& _pi, const Particle& _pj) {
        Float dr[3] = {_pj.pos[0]-_pi.pos[0], _pj.pos[1]-_pi.pos[1], _pj.pos[2]-_pi.pos[2]};
        _pert_out += calcPertFromMR(sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]), _pi.mass, _pj.mass);
    }

    void calcSlowDownPert(Float& _pert_out, Float& _t_min_sq, const Float&, const Particle&, const Perturber&) {
        _pert_out = 0.0;
        _t_min_sq = 0.0;
    }

    // LogH only
    Float calcGTDriftInv(Float _ekin_minus_etot) { return _ekin_minus_etot; }

    // LogH only
    Float calcH(Float _ekin_minus_etot, Float _epot) {
        if (_ekin_minus_etot==0.0&&_epot==0.0) return 0;
        else return log(_ekin_minus_etot) - log(-_epot);
    }

    void modifyAndInterruptIter(AR::InterruptBinary<Particle>&, AR::BinaryTree<Particle>&) {}
};

//! one combination of time transformation and Slow-down scheme
template <class Ttransform, class Tslowdown>
struct Variant {
    typedef AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, TestInteraction, AR::Information<Particle,Particle>, Ttransform, Tslowdown> Integrator;
};

//! particle with given mass, position and velocity
inline Particle makeParticle(const Float _mass, const Float _x, const Float _y, const Float _vx, const Float _vy) {
    Particle p;
    p.mass = _mass;
    p.pos[0] = _x; p.pos[1] = _y; p.pos[2] = 0.0;
    p.vel[0] = _vx; p.vel[1] = _vy; p.vel[2] = 0.0;
    return p;
}

//! replace the last body by a binary at apocentre with the same centre of mass (G=1), which makes a hierarchy level
inline void splitLastIntoBinary(std::vector<Particle>& _bodies, const Float _m1, const Float _m2, const Float _semi, const Float _ecc) {
    const Particle cm = _bodies.back();
    _bodies.pop_back();
    const Float mtot = _m1 + _m2;
    const Float r = _semi*(1.0+_ecc);                           // separation along x
    const Float v = sqrt(mtot/_semi*(1.0-_ecc)/(1.0+_ecc));     // relative velocity along y
    _bodies.push_back(makeParticle(_m1, cm.pos[0] - _m2/mtot*r, cm.pos[1], cm.vel[0], cm.vel[1] - _m2/mtot*v));
    _bodies.push_back(makeParticle(_m2, cm.pos[0] + _m1/mtot*r, cm.pos[1], cm.vel[0], cm.vel[1] + _m1/mtot*v));
}

//! AR integrator of one variant with its manager, initialized the way the AR sample does it
template <class Tintegrator>
struct ArSystem {
    AR::TimeTransformedSymplecticManager<TestInteraction> manager;
    Tintegrator ar;

    /*! @param[in] _slowdown_timescale_max: maximum timescale of Slow-down; below the shortest period every Slow-down factor is 1
        @param[in] _slowdown_pert_ratio_ref: Slow-down factor per ratio of inner to outer perturbation
     */
    ArSystem(const std::vector<Particle>& _bodies, const Float _slowdown_timescale_max, const Float _slowdown_pert_ratio_ref=1e-6) {
        manager.interaction.gravitational_constant = 1.0;
        manager.time_step_min = 1e-13;
        manager.time_error_max = 0.25e-13;
        manager.energy_error_relative_max = 1e-10;
        manager.slowdown_pert_ratio_ref = _slowdown_pert_ratio_ref;
        manager.slowdown_timescale_max = _slowdown_timescale_max;
        manager.step_count_max = 1000000;
        manager.ds_scale = 1.0;
        manager.step.initialSymplecticCofficients(-6);
        manager.interrupt_detection_option = 0;

        const int n = int(_bodies.size());
        ar.manager = &manager;
        ar.particles.setMode(COMM::ListMode::local);
        ar.particles.reserveMem(n);
        for (int i=0; i<n; i++) {
            ar.particles.addMember(_bodies[i]);
            ar.particles[i].id = i+1;
        }
        ar.reserveIntegratorMem();
        ar.particles.calcCenterOfMass();
        ar.info.reserveMem(n);
        ar.info.generateBinaryTree(ar.particles, manager.interaction.gravitational_constant);
        ar.info.r_break_crit = 1e-3;
        ar.initialIntegration(0.0);
        ar.info.calcDsAndStepOption(manager.step.getOrder(), manager.interaction.gravitational_constant, manager.ds_scale);
    }

    ~ArSystem() { manager.step.clear(); }

    //! integrate to a Physical time in a number of equal intervals, rebuilding the binary tree after each as the AR sample does
    void integrate(const Float _time_end, const int _n_interval=1) {
        for (int i=1; i<=_n_interval; i++) {
            ar.integrateToTime(_time_end*i/_n_interval);
            ar.info.generateBinaryTree(ar.particles, manager.interaction.gravitational_constant);
        }
    }

    //! largest position difference to another system, particle by particle
    template <class Tother>
    Float distanceTo(ArSystem<Tother>& _other) {
        Float dmax = 0.0;
        for (int i=0; i<ar.particles.getSize(); i++) 
            for (int k=0; k<3; k++) dmax = std::max(dmax, Float(abs(ar.particles[i].pos[k] - _other.ar.particles[i].pos[k])));
        return dmax;
    }
};
