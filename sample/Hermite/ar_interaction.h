#pragma once

#include <cmath>
#include "Common/Float.h"
#include "AR/force.h"
#include "Hermite/neighbor.h"
#include "particle.h"

//! a sample interaction class with newtonian acceleration
class ARInteraction{
    typedef H4::ParticleAR<Particle> ARPtcl;
    typedef H4::ParticleH4<Particle> H4Ptcl;
private:
    //! calculate inner member acceleration, potential and time transformation function gradient and factor for kick (two-body case)
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the time transformation factor (gt_kick) for kick step
    */
    Float calcAccPotAndGTKickTwo(AR::Force* _force, Float& _epot, const ARPtcl* _particles, const int _n_particle) {
        ASSERT(_n_particle==2);

        // acceleration
        const Float mass1 = _particles[0].mass;
        const Float* pos1 = _particles[0].pos;

        const Float mass2 = _particles[1].mass;
        const Float* pos2 = _particles[1].pos;

        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
        ASSERT(r2>0.0);
        Float inv_r = 1.0/sqrt(r2);
        Float inv_r3 = inv_r*inv_r*inv_r;

        Float* acc1 = _force[0].acc_in;
        Float* acc2 = _force[1].acc_in;

        Float mor3_1 = G*mass2*inv_r3;
        acc1[0] = mor3_1 * dr[0];
        acc1[1] = mor3_1 * dr[1];
        acc1[2] = mor3_1 * dr[2];

        Float mor3_2 = G*mass1*inv_r3;
        acc2[0] = - mor3_2 * dr[0];
        acc2[1] = - mor3_2 * dr[1];
        acc2[2] = - mor3_2 * dr[2];

        Float m1m2 = G*mass1*mass2;

#ifdef AR_TTL 
        // trans formation function gradient
        Float m1m2or3 = m1m2*inv_r3;
        Float* gtgrad1 = _force[0].gtgrad;
        Float* gtgrad2 = _force[1].gtgrad;

        gtgrad1[0] = m1m2or3 * dr[0];
        gtgrad1[1] = m1m2or3 * dr[1];
        gtgrad1[2] = m1m2or3 * dr[2];

        gtgrad2[0] = - gtgrad1[0];
        gtgrad2[1] = - gtgrad1[1];
        gtgrad2[2] = - gtgrad1[2];
#endif

        // potential energy
        Float m1m2or = m1m2*inv_r;
        _epot = - m1m2or;

        // transformation factor for kick
        Float gt_kick = 1.0/m1m2or;

        return gt_kick;
    }

    //! calculate inner member acceleration, potential and time transformation function gradient and factor for kick
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the time transformation factor (gt_kick) for kick step
    */
    Float calcAccPotAndGTKick(AR::Force* _force, Float& _epot, const ARPtcl* _particles, const int _n_particle) {
        _epot = Float(0.0);
        Float gt_kick = Float(0.0);
        for (int i=0; i<_n_particle; i++) {
            const Float massi = _particles[i].mass;
            const Float* posi = _particles[i].pos;
            Float* acci = _force[i].acc_in;
            acci[0] = acci[1] = acci[2] = Float(0.0);

#ifdef AR_TTL 
            Float* gtgradi = _force[i].gtgrad;
            gtgradi[0] = gtgradi[1] = gtgradi[2] = Float(0.0);
#endif

            Float poti = Float(0.0);
            Float gtki = Float(0.0);

            for (int j=0; j<_n_particle; j++) {
                if (i==j) continue;
                const Float massj = _particles[j].mass;
                const Float* posj = _particles[j].pos; 
                Float dr[3] = {posj[0] -posi[0],
                               posj[1] -posi[1],
                               posj[2] -posi[2]};
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                Float inv_r = 1.0/sqrt(r2);
                Float inv_r3 = inv_r*inv_r*inv_r;
                Float mor3 = G*massj*inv_r3;
                acci[0] += mor3 * dr[0];
                acci[1] += mor3 * dr[1];
                acci[2] += mor3 * dr[2];

#ifdef AR_TTL                     
                Float mimjor3 = massi*mor3;
                gtgradi[0] += mimjor3 * dr[0];
                gtgradi[1] += mimjor3 * dr[1];
                gtgradi[2] += mimjor3 * dr[2];
#endif

                Float mor = G*massj*inv_r;
                poti -= mor;
                gtki += mor;
                    
            }
            _epot += poti * massi;
            gt_kick += gtki * massi;
        }
        _epot   *= 0.5;
        gt_kick = 2.0/gt_kick;

        return gt_kick;
    }

public:
    Float eps_sq;
    Float G;

    ARInteraction(): eps_sq(Float(-1.0)), G(Float(-1.0)) {}

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(eps_sq>=0.0);
        ASSERT(G>0.0);
        return true;
    }        

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"eps_sq: "<<eps_sq<<std::endl
             <<"G     : "<<G<<std::endl;
    }    


    //! calculate perturbation from c.m. acceleration
    Float calcPertFromAcc(const Float* _acc, const Float _mp, const Float _mpert) {
        Float acc2 = _acc[0]*_acc[0]+_acc[1]*_acc[1]+_acc[2]*_acc[2];
        return acc2/(_mp*_mpert);
    }

    //! calculate perturbation from binary tree
    Float calcPertFromBinary(const COMM::BinaryTree<ARPtcl>& _bin) {
        Float apo = _bin.semi*(1.0+_bin.ecc);
        Float r2 = apo*apo;
        return (_bin.m1*_bin.m2)/(r2*r2);
    }

    //! calculate perturbation from distance to perturber and masses of particle and perturber 
    Float calcPertFromMR(const Float _r, const Float _mp, const Float _mpert) {
        Float r2 = _r*_r;
        return _mp*_mpert/(r2*r2);
    }

#ifdef SLOWDOWN_INTEGRATE

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*! The Force class acc_pert should be updated
      @param[in,out] _slowdown: slowdown class to store perturbation (member: pert: perturbation ([m]/[r^3]); timescale: limited timescale to calculate maximum slowdown factor)
      @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[out] _epot: potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _bin_root: binary tree root
      @param[in] _perturber: pertuber container
      \return perturbation energy to calculate slowdown factor
    */
    Float calcAccEnergyAndSlowDownPert(AR::SlowDown& _slowdown, AR::Force* _force, Float& _epot, const ARPtcl* _particles, const int _n_particle, const H4Ptcl& _particle_cm, const COMM::BinaryTree<ARPtcl>& _bin_root, const H4::Neighbor<Particle>& _perturber) {
        static const Float inv3 = 1.0 / 3.0;

        Float gt_kick;
        if (_n_particle==2) gt_kick = calcAccPotAndGTKickTwo(_force, _epot, _particles, _n_particle);
        else gt_kick = calcAccPotAndGTKick(_force, _epot, _particles, _n_particle);

        // slowdown inner perturbation: m1*m2/apo_in^4
        Float apo_in = _bin_root.semi*(1.0+_bin_root.ecc);
        Float apo_in2 = apo_in*apo_in;
        _slowdown.pert_in = _bin_root.m1*_bin_root.m2/(apo_in2*apo_in2);
        _slowdown.period  = _bin_root.period;

        const Float time = _slowdown.getRealTime();

        const int n_pert = _perturber.neighbor_address.getSize();

        if (n_pert>0) {

            auto* pert_adr = _perturber.neighbor_address.getDataAddress();

            Float xp[n_pert][3], xcm[3], m[n_pert];
            Float vp[n_pert][3], vcm[3];

            for (int j=0; j<n_pert; j++) {
                H4::NBAdr<Particle>::Single* pertj;
                if (pert_adr[j].type==H4::NBType::group) pertj = &(((H4::NBAdr<Particle>::Group*)pert_adr[j].adr)->cm);
                else pertj = (H4::NBAdr<Particle>::Single*)pert_adr[j].adr;
                Float dt = time - pertj->time;
                ASSERT(dt>=0.0);
                xp[j][0] = pertj->pos[0] + dt*(pertj->vel[0] + 0.5*dt*(pertj->acc0[0] + inv3*dt*pertj->acc1[0]));
                xp[j][1] = pertj->pos[1] + dt*(pertj->vel[1] + 0.5*dt*(pertj->acc0[1] + inv3*dt*pertj->acc1[1]));
                xp[j][2] = pertj->pos[2] + dt*(pertj->vel[2] + 0.5*dt*(pertj->acc0[2] + inv3*dt*pertj->acc1[2]));

                vp[j][0] = pertj->vel[0] + dt*(pertj->acc0[0] + 0.5*dt*pertj->acc1[0]);
                vp[j][1] = pertj->vel[1] + dt*(pertj->acc0[1] + 0.5*dt*pertj->acc1[1]);
                vp[j][2] = pertj->vel[2] + dt*(pertj->acc0[2] + 0.5*dt*pertj->acc1[2]);

                m[j] = pertj->mass;
            }

            Float dt = time - _particle_cm.time;
            ASSERT(dt>=0.0);
            xcm[0] = _particle_cm.pos[0] + dt*(_particle_cm.vel[0] + 0.5*dt*(_particle_cm.acc0[0] + inv3*dt*_particle_cm.acc1[0]));
            xcm[1] = _particle_cm.pos[1] + dt*(_particle_cm.vel[1] + 0.5*dt*(_particle_cm.acc0[1] + inv3*dt*_particle_cm.acc1[1]));
            xcm[2] = _particle_cm.pos[2] + dt*(_particle_cm.vel[2] + 0.5*dt*(_particle_cm.acc0[2] + inv3*dt*_particle_cm.acc1[2]));

            vcm[0] = _particle_cm.vel[0] + dt*(_particle_cm.acc0[0] + 0.5*dt*_particle_cm.acc1[0]);
            vcm[1] = _particle_cm.vel[1] + dt*(_particle_cm.acc0[1] + 0.5*dt*_particle_cm.acc1[1]);
            vcm[2] = _particle_cm.vel[2] + dt*(_particle_cm.acc0[2] + 0.5*dt*_particle_cm.acc1[2]);


            Float pert_cm = 0.0, acc_pert_cm[3]={0.0, 0.0, 0.0};
            Float mcm = _particle_cm.mass;
            if (_perturber.need_resolve_flag) {
                // calculate component perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);

                    Float xi[3];
                    xi[0] = pi.pos[0] + xcm[0];
                    xi[1] = pi.pos[1] + xcm[1];
                    xi[2] = pi.pos[2] + xcm[2];

                    Float pert_pot = 0.0;
                    for (int j=0; j<n_pert; j++) {
                        Float dr[3] = {xp[j][0] - xi[0],
                                       xp[j][1] - xi[1],
                                       xp[j][2] - xi[2]};
                        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                        Float r  = sqrt(r2);
                        Float r3 = r*r2;
                        Float mor3 = G*m[j]/r3;
                        pert_pot += mor3/r;

                        acc_pert[0] += mor3 * dr[0];
                        acc_pert[1] += mor3 * dr[1];
                        acc_pert[2] += mor3 * dr[2];
                    }

                    pert_cm += pi.mass * pert_pot;
                    acc_pert_cm[0] += pi.mass *acc_pert[0];
                    acc_pert_cm[1] += pi.mass *acc_pert[1];
                    acc_pert_cm[2] += pi.mass *acc_pert[2];

#ifdef ARC_DEBUG
                    mcm += pi.mass;
#endif
                }
#ifdef ARC_DEBUG
                ASSERT(abs(mcm-_particles_cm.mass)<1e-10);
#endif
                
                // get cm perturbation
                acc_pert_cm[0] /= mcm;
                acc_pert_cm[1] /= mcm;
                acc_pert_cm[2] /= mcm;

                // remove cm. perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    acc_pert[0] -= acc_pert_cm[0]; 
                    acc_pert[1] -= acc_pert_cm[1];        
                    acc_pert[2] -= acc_pert_cm[2]; 
                }

                // get SD time scale
                _slowdown.timescale = NUMERIC_FLOAT_MAX;
                for (int i=0; i<n_pert; i++) {
                    Float dr[3] = {xp[i][0] - xcm[0],
                                   xp[i][1] - xcm[1],
                                   xp[i][2] - xcm[2]};
                    Float dv[3] = {vp[i][0] - vcm[0],
                                   vp[i][1] - vcm[1],
                                   vp[i][2] - vcm[2]};
                    Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                    Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
                    Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
                    Float ti = (drdv*drdv>0.25*v2*r2)? r2*r2/(drdv*drdv) : r2/v2;
                    _slowdown.timescale = std::min(_slowdown.timescale, ti);
                }
                _slowdown.timescale = std::sqrt(_slowdown.timescale);
                _slowdown.timescale = 0.1*std::min(_slowdown.getTimescaleMax(), _slowdown.timescale);
                
                _slowdown.pert_out = pert_cm;
            }
            else {
                // first calculate c.m. acceleration and tidal perturbation
                _slowdown.timescale = NUMERIC_FLOAT_MAX;
                for (int j=0; j<n_pert; j++) {
                    Float dr[3] = {xp[j][0] - xcm[0],
                                   xp[j][1] - xcm[1],
                                   xp[j][2] - xcm[2]};
                    Float dv[3] = {vp[j][0] - vcm[0],
                                   vp[j][1] - vcm[1],
                                   vp[j][2] - vcm[2]};
                    Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                    Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
                    Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
                    Float ti = (drdv*drdv>0.25*v2*r2)? r2*r2/(drdv*drdv) : r2/v2;
                    Float r  = sqrt(r2);
                    Float r3 = r*r2;
                    Float mor3 = G*m[j]/r3;
                    pert_cm += mor3/r;

                    acc_pert_cm[0] += mor3 * dr[0];
                    acc_pert_cm[1] += mor3 * dr[1];
                    acc_pert_cm[2] += mor3 * dr[2];

                    _slowdown.timescale = std::min(_slowdown.timescale, ti);
                }
                _slowdown.timescale = std::sqrt(_slowdown.timescale);
                _slowdown.timescale = 0.1*std::min(_slowdown.getTimescaleMax(), _slowdown.timescale);
            

                // calculate component perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);

                    Float xi[3];
                    xi[0] = pi.pos[0] + xcm[0];
                    xi[1] = pi.pos[1] + xcm[1];
                    xi[2] = pi.pos[2] + xcm[2];

                    // remove c.m. perturbation acc
                    acc_pert[0] = -acc_pert_cm[0]; 
                    acc_pert[1] = -acc_pert_cm[1];        
                    acc_pert[2] = -acc_pert_cm[2]; 

                    for (int j=0; j<n_pert; j++) {
                        Float dr[3] = {xp[j][0] - xi[0],
                                       xp[j][1] - xi[1],
                                       xp[j][2] - xi[2]};
                        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                        Float r  = sqrt(r2);
                        Float r3 = r*r2;
                        Float mor3 = G*m[j]/r3;

                        acc_pert[0] += mor3 * dr[0];
                        acc_pert[1] += mor3 * dr[1];
                        acc_pert[2] += mor3 * dr[2];
                    }
                }   
                _slowdown.pert_out =  _particle_cm.mass*pert_cm;
            }
        }
        else {
            _slowdown.pert_out = 0.0;
            _slowdown.timescale = _slowdown.getTimescaleMax();
        }
        
        return gt_kick;
    }

#else

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*! The Force class acc_pert should be updated
      @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[out] _epot: potential 
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
      \return perturbation energy to calculate slowdown factor
    */
    Float calcAccEnergy(AR::Force* _force, Float& _epot, const ARPtcl* _particles, const int _n_particle, const H4Ptcl& _particle_cm, const H4::Neighbor<Particle>& _perturber, const Float _time) {
        static const Float inv3 = 1.0 / 3.0;

        Float gt_kick;
        if (_n_particle==2) gt_kick = calcAccPotAndGTKickTwo(_force, _epot, _particles, _n_particle);
        else gt_kick = calcAccPotAndGTKick(_force, _epot, _particles, _n_particle);

        const int n_pert = _perturber.neighbor_address.getSize();
        if (n_pert>0) {

            auto* pert_adr = _perturber.neighbor_address.getDataAddress();

            Float xp[n_pert][3], xcm[3], m[n_pert];

            for (int j=0; j<n_pert; j++) {
                H4::NBAdr<Particle>::Single* pertj;
                if (pert_adr[j].type==H4::NBType::group) pertj = &(((H4::NBAdr<Particle>::Group*)pert_adr[j].adr)->cm);
                else pertj = (H4::NBAdr<Particle>::Single*)pert_adr[j].adr;
                Float dt = _time - pertj->time;
                ASSERT(dt>=0.0);
                xp[j][0] = pertj->pos[0] + dt*(pertj->vel[0] + 0.5*dt*(pertj->acc0[0] + inv3*dt*pertj->acc1[0]));
                xp[j][1] = pertj->pos[1] + dt*(pertj->vel[1] + 0.5*dt*(pertj->acc0[1] + inv3*dt*pertj->acc1[1]));
                xp[j][2] = pertj->pos[2] + dt*(pertj->vel[2] + 0.5*dt*(pertj->acc0[2] + inv3*dt*pertj->acc1[2]));

                m[j] = pertj->mass;
            }

            Float dt = _time - _particle_cm.time;
            ASSERT(dt>=0.0);
            xcm[0] = _particle_cm.pos[0] + dt*(_particle_cm.vel[0] + 0.5*dt*(_particle_cm.acc0[0] + inv3*dt*_particle_cm.acc1[0]));
            xcm[1] = _particle_cm.pos[1] + dt*(_particle_cm.vel[1] + 0.5*dt*(_particle_cm.acc0[1] + inv3*dt*_particle_cm.acc1[1]));
            xcm[2] = _particle_cm.pos[2] + dt*(_particle_cm.vel[2] + 0.5*dt*(_particle_cm.acc0[2] + inv3*dt*_particle_cm.acc1[2]));

            Float acc_pert_cm[3]={0.0, 0.0, 0.0};
            Float mcm = _particle_cm.mass;
            if (_perturber.need_resolve_flag) {
                // calculate component perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);

                    Float xi[3];
                    xi[0] = pi.pos[0] + xcm[0];
                    xi[1] = pi.pos[1] + xcm[1];
                    xi[2] = pi.pos[2] + xcm[2];

                    for (int j=0; j<n_pert; j++) {
                        Float dr[3] = {xp[j][0] - xi[0],
                                       xp[j][1] - xi[1],
                                       xp[j][2] - xi[2]};
                        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                        Float r  = sqrt(r2);
                        Float r3 = r*r2;
                        Float mor3 = G*m[j]/r3;

                        acc_pert[0] += mor3 * dr[0];
                        acc_pert[1] += mor3 * dr[1];
                        acc_pert[2] += mor3 * dr[2];
                    }

                    acc_pert_cm[0] += pi.mass *acc_pert[0];
                    acc_pert_cm[1] += pi.mass *acc_pert[1];
                    acc_pert_cm[2] += pi.mass *acc_pert[2];
#ifdef ARC_DEBUG
                    mcm += pi.mass;
#endif
                }
#ifdef ARC_DEBUG
                ASSERT(abs(mcm-_particles_cm.mass)<1e-10);
#endif
                // get cm perturbation
                acc_pert_cm[0] /= mcm;
                acc_pert_cm[1] /= mcm;
                acc_pert_cm[2] /= mcm;

                // remove cm. perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    acc_pert[0] -= acc_pert_cm[0]; 
                    acc_pert[1] -= acc_pert_cm[1];        
                    acc_pert[2] -= acc_pert_cm[2]; 
                }
            }
            else {
                // first calculate c.m. acceleration and tidal perturbation
                for (int j=0; j<n_pert; j++) {
                    Float dr[3] = {xp[j][0] - xcm[0],
                                   xp[j][1] - xcm[1],
                                   xp[j][2] - xcm[2]};
                    Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                    Float r  = sqrt(r2);
                    Float r3 = r*r2;
                    Float mor3 = G*m[j]/r3;

                    acc_pert_cm[0] += mor3 * dr[0];
                    acc_pert_cm[1] += mor3 * dr[1];
                    acc_pert_cm[2] += mor3 * dr[2];
                }

                // calculate component perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);

                    Float xi[3];
                    xi[0] = pi.pos[0] + xcm[0];
                    xi[1] = pi.pos[1] + xcm[1];
                    xi[2] = pi.pos[2] + xcm[2];

                    // remove c.m. perturbation acc
                    acc_pert[0] = -acc_pert_cm[0]; 
                    acc_pert[1] = -acc_pert_cm[1];        
                    acc_pert[2] = -acc_pert_cm[2]; 

                    for (int j=0; j<n_pert; j++) {
                        Float dr[3] = {xp[j][0] - xi[0],
                                       xp[j][1] - xi[1],
                                       xp[j][2] - xi[2]};
                        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                        Float r  = sqrt(r2);
                        Float r3 = r*r2;
                        Float mor3 = G*m[j]/r3;

                        acc_pert[0] += mor3 * dr[0];
                        acc_pert[1] += mor3 * dr[1];
                        acc_pert[2] += mor3 * dr[2];
                    }
                }

            }
        }
        return gt_kick;
        
    }

    //! (Necessary) calculate slowdown perturbation and timescale
    /*!
      @param[in,out] _slowdown: slowdown paramters, (pert_out and timescale should be updated)
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _bin_root: binary tree root
      @param[in] _perturber: pertuber container
    */
    void calcSlowDownPert(AR::SlowDown& _slowdown, const H4Ptcl& _particle_cm, const COMM::BinaryTree<ARPtcl>& _bin_root, const H4::Neighbor<Particle>& _perturber) {
        static const Float inv3 = 1.0 / 3.0;

        // slowdown inner perturbation: m1*m2/apo_in^4
        Float apo_in = _bin_root.semi*(1.0+_bin_root.ecc);
        Float apo_in2 = apo_in*apo_in;
        _slowdown.pert_in = _bin_root.m1*_bin_root.m2/(apo_in2*apo_in2);
        _slowdown.period  = _bin_root.period;

        const Float time = _slowdown.getRealTime();

        const int n_pert = _perturber.neighbor_address.getSize();

        if (n_pert>0) {

            auto* pert_adr = _perturber.neighbor_address.getDataAddress();

            Float xp[n_pert][3], xcm[3], m[n_pert];
            Float vp[n_pert][3], vcm[3];

            Float dt = time - _particle_cm.time;
            //ASSERT(dt>=0.0);
            xcm[0] = _particle_cm.pos[0] + dt*(_particle_cm.vel[0] + 0.5*dt*(_particle_cm.acc0[0] + inv3*dt*_particle_cm.acc1[0]));
            xcm[1] = _particle_cm.pos[1] + dt*(_particle_cm.vel[1] + 0.5*dt*(_particle_cm.acc0[1] + inv3*dt*_particle_cm.acc1[1]));
            xcm[2] = _particle_cm.pos[2] + dt*(_particle_cm.vel[2] + 0.5*dt*(_particle_cm.acc0[2] + inv3*dt*_particle_cm.acc1[2]));

            vcm[0] = _particle_cm.vel[0] + dt*(_particle_cm.acc0[0] + 0.5*dt*_particle_cm.acc1[0]);
            vcm[1] = _particle_cm.vel[1] + dt*(_particle_cm.acc0[1] + 0.5*dt*_particle_cm.acc1[1]);
            vcm[2] = _particle_cm.vel[2] + dt*(_particle_cm.acc0[2] + 0.5*dt*_particle_cm.acc1[2]);

            Float pert_pot = 0.0, t2_min = NUMERIC_FLOAT_MAX;
            Float mcm = _particle_cm.mass;

            for (int j=0; j<n_pert; j++) {
                H4::NBAdr<Particle>::Single* pertj;
                if (pert_adr[j].type==H4::NBType::group) pertj = &(((H4::NBAdr<Particle>::Group*)pert_adr[j].adr)->cm);
                else pertj = (H4::NBAdr<Particle>::Single*)pert_adr[j].adr;
                Float dt = time - pertj->time;
                //ASSERT(dt>=0.0);
                xp[j][0] = pertj->pos[0] + dt*(pertj->vel[0] + 0.5*dt*(pertj->acc0[0] + inv3*dt*pertj->acc1[0]));
                xp[j][1] = pertj->pos[1] + dt*(pertj->vel[1] + 0.5*dt*(pertj->acc0[1] + inv3*dt*pertj->acc1[1]));
                xp[j][2] = pertj->pos[2] + dt*(pertj->vel[2] + 0.5*dt*(pertj->acc0[2] + inv3*dt*pertj->acc1[2]));

                vp[j][0] = pertj->vel[0] + dt*(pertj->acc0[0] + 0.5*dt*pertj->acc1[0]);
                vp[j][1] = pertj->vel[1] + dt*(pertj->acc0[1] + 0.5*dt*pertj->acc1[1]);
                vp[j][2] = pertj->vel[2] + dt*(pertj->acc0[2] + 0.5*dt*pertj->acc1[2]);

                m[j] = pertj->mass;

                Float dr[3] = {xp[j][0] - xcm[0],
                               xp[j][1] - xcm[1],
                               xp[j][2] - xcm[2]};
                Float dv[3] = {vp[j][0] - vcm[0],
                               vp[j][1] - vcm[1],
                               vp[j][2] - vcm[2]};
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                Float r4 = r2*r2;
                pert_pot += m[j]/r4;

                Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
                Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
                Float ti2 = r2*r2/(drdv*drdv+1e-4*v2*r2);

                t2_min = std::min(ti2, t2_min);
            }
            _slowdown.pert_out = mcm* pert_pot;
            _slowdown.timescale = 0.1*std::min(_slowdown.getTimescaleMax(), sqrt(t2_min));
        }
        else{
            _slowdown.pert_out = 0.0;
            _slowdown.timescale = _slowdown.getTimescaleMax();
        }
    }    
        
#endif

#ifndef AR_TTL
    //! (Necessary) calcualte the time transformation factor for drift
    /*! The time transformation factor for drift only depends on (kinetic energy - total energy)
      @param[in] _ekin_minus_etot: ekin - etot
    */
    Float calcGTDrift(Float _ekin_minus_etot) {
        return 1.0/_ekin_minus_etot;
    }
#endif   

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(this, sizeof(*this),1,_fp);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }    
};

