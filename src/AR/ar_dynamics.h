#pragma once

#include "Common/Float.h"
#include "Common/list.h"
#include "Common/particle_group.h"
#include "AR/variant.h"
#include "AR/force.h"
#include "AR/slow_down.h"
#include "AR/profile.h"
#include "AR/information.h"

namespace AR {

    template <class Tmethod> class TimeTransformedSymplecticManager;

    //! State and helpers shared by every AR dynamics
    /*! The integrated and calculated variables that exist for every time transformation and Slow-down scheme,
        and the binary tree helpers that do not depend on them.
     */
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo, class Ttransform>
    class ARDynamicsBase {
    public:
        typedef ForceType<Ttransform> Force; ///< force type, determined by the time transformation

    protected:
        // intergrated variables
        Float time_ = 0.0;   ///< integrated time (not the physical time, it is always 0 initially)
        Float etot_ref_ = 0.0; ///< integrated system energy

        // calculated varaiables
        Float ekin_ = 0.0;   ///< kinetic energy
        Float epot_ = 0.0;   ///< potential

        // cumulative slowdown (inner + outer) energy change
        Float de_change_interrupt_ = 0.0;      // energy change due to interruption
        Float dH_change_interrupt_ = 0.0;      // hamiltonian change due to interruption

        // force array
        COMM::List<Force> force_; ///< acceleration array 

    public:
        TimeTransformedSymplecticManager<Tmethod>* manager = NULL; ///< integration manager
        COMM::ParticleGroup<Tparticle,Tpcm> particles; ///< particle group manager
        Tpert    perturber; ///< perturber class 
        Tinfo    info;   ///< information of the system
        Profile  profile;  ///< profile to measure the performance

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            ASSERT(manager!=NULL);
            ASSERT(manager->checkParams());
            ASSERT(perturber.checkParams());
            ASSERT(info.checkParams());
            return true;
        }

        //! Get current physical time
        /*! \return current physical time
         */
        Float getTime() const {
            return time_;
        }

    protected:
        //! clear the shared data
        void clearBase() {
            time_ = 0.0;
            etot_ref_ =0.0;
            ekin_ = 0.0;
            epot_ = 0.0;
            de_change_interrupt_ = 0.0;
            dH_change_interrupt_ = 0.0;
            force_.clear();
            particles.clear();
            perturber.clear();
            info.clear();
            profile.clear();
        }

        //! copy the shared data (the perturber is not copied)
        void copyBase(const ARDynamicsBase& _sym) {
            time_   = _sym.time_;
            etot_ref_   = _sym.etot_ref_;

            ekin_   = _sym.ekin_;
            epot_   = _sym.epot_;
            de_change_interrupt_= _sym.de_change_interrupt_;
            dH_change_interrupt_= _sym.dH_change_interrupt_;
            force_  = _sym.force_;
            manager = _sym.manager;
            particles = _sym.particles;
            info = _sym.info;
            profile = _sym.profile;
        }

        //! update slowdown velocity iteration function with binary tree
        void updateBinaryVelIter(AR::BinaryTree<Tparticle>& _bin) {
            _bin.vel[0]= _bin.vel[1] = _bin.vel[2] = 0.0;
            Float mcm_inv = 1.0/_bin.mass;
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    updateBinaryVelIter(*bink);
                    _bin.vel[0] += bink->mass*bink->vel[0];
                    _bin.vel[1] += bink->mass*bink->vel[1];
                    _bin.vel[2] += bink->mass*bink->vel[2];
                }
                else {
                    auto* pk = _bin.getMember(k);
                    _bin.vel[0] += pk->mass*pk->vel[0];
                    _bin.vel[1] += pk->mass*pk->vel[1];
                    _bin.vel[2] += pk->mass*pk->vel[2];
                }
            }
            _bin.vel[0] *= mcm_inv;
            _bin.vel[1] *= mcm_inv;
            _bin.vel[2] *= mcm_inv;
        }

        //! update binary cm iteratively
        void updateBinaryCMIter(AR::BinaryTree<Tparticle>& _bin) {
            _bin.pos[0]= _bin.pos[1] = _bin.pos[2] = 0.0;
            _bin.vel[0]= _bin.vel[1] = _bin.vel[2] = 0.0;
            Float mcm_member = 0.0;
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    updateBinaryCMIter(*bink);
                    mcm_member  += bink->mass;
                    _bin.pos[0] += bink->mass*bink->pos[0];
                    _bin.pos[1] += bink->mass*bink->pos[1];
                    _bin.pos[2] += bink->mass*bink->pos[2];
                    _bin.vel[0] += bink->mass*bink->vel[0];
                    _bin.vel[1] += bink->mass*bink->vel[1];
                    _bin.vel[2] += bink->mass*bink->vel[2];
                }
                else {
                    auto* pk = _bin.getMember(k);
                    mcm_member  += pk->mass;
                    _bin.pos[0] += pk->mass*pk->pos[0];
                    _bin.pos[1] += pk->mass*pk->pos[1];
                    _bin.pos[2] += pk->mass*pk->pos[2];
                    _bin.vel[0] += pk->mass*pk->vel[0];
                    _bin.vel[1] += pk->mass*pk->vel[1];
                    _bin.vel[2] += pk->mass*pk->vel[2];
                }
            }
            Float mcm_inv = 1.0/mcm_member;
            _bin.mass = mcm_member;
            _bin.m1 = _bin.getLeftMember()->mass;
            _bin.m2 = _bin.getRightMember()->mass;
            _bin.vel[0] *= mcm_inv;
            _bin.vel[1] *= mcm_inv;
            _bin.vel[2] *= mcm_inv;
            _bin.pos[0] *= mcm_inv;
            _bin.pos[1] *= mcm_inv;
            _bin.pos[2] *= mcm_inv;
        }

        //! set all binary c.m. mass to zero
        void setBinaryCMZeroIter(AR::BinaryTree<Tparticle>& _bin) {
            _bin.mass = 0.0;
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    setBinaryCMZeroIter(*bink);
                }
            }
        }

        //! update binary semi, ecc and period iteratively after 0.25 period for unstable systems
        bool updateBinarySemiEccPeriodIter(AR::BinaryTree<Tparticle>& _bin, const Float& _G, const Float _time, const bool _check=false) {
            bool check = _check;
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    check=updateBinarySemiEccPeriodIter(*bink, _G, _time, _check);
                }
            }
            if ((_time>_bin.stab_check_time&&_bin.stab>1.0&&_bin.m1>0.0&&_bin.m2>0.0)||check) {
                _bin.calcSemiEccPeriod(_G);
                _bin.stab_check_time = _time + _bin.period;
                return true;
            }
            return false;
        }
    };

    //! Helpers shared by the AR dynamics without Slow-down, identical for LogH and TTL
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo, class Ttransform>
    class ARNoSlowDownHelpers: public ARDynamicsBase<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, Ttransform> {
    protected:
        typedef ARDynamicsBase<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, Ttransform> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;

    protected:
        //! Calculate kinetic energy
        inline void calcEKin(){
            ekin_ = Float(0.0);
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            for (int i=0; i<num; i++) {
                const Float *vi=pdat[i].getVel();
                ekin_ += 0.5 * pdat[i].mass * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
            }
        }

        //! kick velocity
        /*! First time step will be calculated, the velocities are kicked
          @param[in] _dt: time size
        */
        inline void kickVel(const Float _dt) {
            const int num = particles.getSize();            
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0; i<num; i++) {
                // kick velocity
                Float* vel = pdat[i].getVel();
                Float* acc = force[i].acc_in;
                Float* pert= force[i].acc_pert;
                // half dv 
                vel[0] += _dt * (acc[0] + pert[0]);
                vel[1] += _dt * (acc[1] + pert[1]);
                vel[2] += _dt * (acc[2] + pert[2]);
            }
        }

        //! drift time and position
        /*! First (real) time is drifted, then positions are drifted
          @param[in] _dt: time step
        */
        inline void driftTimeAndPos(const Float _dt) {
            // drift time 
            time_ += _dt;

            // drift position
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            for (int i=0; i<num; i++) {
                Float* pos = pdat[i].getPos();
                Float* vel = pdat[i].getVel();
                pos[0] += _dt * vel[0];
                pos[1] += _dt * vel[1];
                pos[2] += _dt * vel[2];
            }
        }

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        inline Float calcAccPotAndGTKickInv() {
            Float gt_kick_inv = manager->interaction.calcAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());            

//#ifdef AR_DEBUG
//            // check c.m. force 
//            Force fcm;
//            Float mcm=0.0;
//            for (int i=0; i<particles.getSize(); i++) {
//                for (int k=0; k<3; k++) {
//                    fcm.acc_in[k] += particles[i].mass * force_[i].acc_in[k];
//                }
//                mcm += particles[i].mass;
//            }
//            for (int k=0; k<3; k++) {
//                fcm.acc_in[k] /= mcm;
//                ASSERT(abs(fcm.acc_in[k])<1e-10);
//            }
//#endif
            return gt_kick_inv;
        }

    };

    //! Slow-down energy terms and helpers shared by the AR dynamics with Array or Tree slow-down
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo, class Ttransform>
    class ARSlowDownHelpers: public ARDynamicsBase<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, Ttransform> {
    protected:
        typedef ARDynamicsBase<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, Ttransform> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;

    protected:
        Float ekin_sd_ = 0.0;  ///< slowdown (inner) kinetic energy
        Float epot_sd_ = 0.0;  ///< slowdown (inner) potential energy
        Float etot_sd_ref_ = 0.0;  ///< slowdown (inner) total energy

        Float de_sd_change_cum_ = 0.0;  // slowdown energy change
        Float dH_sd_change_cum_ = 0.0;  // slowdown Hamiltonian change 
        Float de_sd_change_interrupt_ = 0.0;   // slowdown energy change due to interruption
        Float dH_sd_change_interrupt_ = 0.0;   // slowdown energy change due to interruption

        //! clear the Slow-down energy terms
        void clearSlowDownEnergy() {
            ekin_sd_ = 0.0;
            epot_sd_ = 0.0;
            etot_sd_ref_ = 0.0;
            de_sd_change_cum_ = 0.0;
            dH_sd_change_cum_ = 0.0;
            de_sd_change_interrupt_ = 0.0;
            dH_sd_change_interrupt_ = 0.0;
        }

        //! copy the Slow-down energy terms
        void copySlowDownEnergy(const ARSlowDownHelpers& _sym) {
            ekin_sd_= _sym.ekin_sd_;
            epot_sd_= _sym.epot_sd_;
            etot_sd_ref_= _sym.etot_sd_ref_;
            de_sd_change_cum_= _sym.de_sd_change_cum_;
            dH_sd_change_cum_= _sym.dH_sd_change_cum_;
            de_sd_change_interrupt_= _sym.de_sd_change_interrupt_;
            dH_sd_change_interrupt_= _sym.dH_sd_change_interrupt_;
        }

    protected:
        //! iteration function to calculate perturbation and timescale information from binary tree j to binary i
        /*!
          @param[out] _pert_out: perturbation from particle j
          @param[out] _t_min_sq: timescale limit from particle j
          @param[in] _bini: binary i
          @param[in] _binj: binary tree j
         */
        void calcSlowDownPertInnerBinaryIter(Float& _pert_out, Float& _t_min_sq, AR::BinaryTree<Tparticle>& _bini, AR::BinaryTree<Tparticle>& _binj) {
            ASSERT(&_bini != &_binj);
//            ASSERT(_bini.getMemberIndex(0)!=_binj.getMemberIndex(0));

            for (int k=0; k<2; k++) {
                if (_binj.isMemberTree(k)) {
                    auto* bink =  _binj.getMemberAsTree(k);
                    int check_flag = _bini.isSameBranch(*bink);
                    if (check_flag==0) {// no relation
                        if (bink->semi>0.0) manager->interaction.calcSlowDownPertOne(_pert_out, _t_min_sq, _bini, *bink);
                        else calcSlowDownPertInnerBinaryIter(_pert_out, _t_min_sq, _bini, *bink);
                    }
                    else if (check_flag==-2) { // _binj is the upper root tree
                        calcSlowDownPertInnerBinaryIter(_pert_out, _t_min_sq, _bini, *bink);
                    }
                    // other cases (same or sub branch), stop iteration.
                }
                else {
                    auto* pk =  _binj.getMember(k);
                    if (pk->mass>0.0) manager->interaction.calcSlowDownPertOne(_pert_out, _t_min_sq, _bini, *pk);
                }
            }
        }

        //! calculate slowdown factor for inner binary based on other particles and slowdown of system c.m.
        /*!
          @param[in] _bin: binary tree for calculating slowdown
        */
        void calcSlowDownInnerBinary(BinaryTree<Tparticle>& _bin) {
            _bin.slowdown.pert_in = manager->interaction.calcPertFromBinary(_bin);

            Float pert_out = 0.0;
            Float t_min_sq = NUMERIC_FLOAT_MAX;
            auto& bin_root = info.getBinaryTreeRoot();

            calcSlowDownPertInnerBinaryIter(pert_out, t_min_sq, _bin, bin_root);

            _bin.slowdown.pert_out = pert_out + bin_root.slowdown.pert_out;

#ifdef AR_SLOWDOWN_TIMESCALE
            // velocity dependent method
            //Float trv_ave = sdtdat.mtot/sqrt(sdtdat.mvor[0]*sdtdat.mvor[0] + sdtdat.mvor[1]*sdtdat.mvor[1] + sdtdat.mvor[2]*sdtdat.mvor[2]);
            // get min of velocity and force dependent values
            //Float t_min = std::min(trv_ave, sqrt(sdtdat.trf2_min));
            _bin.slowdown.timescale = std::min(_bin.slowdown.getTimescaleMax(), sqrt(t_min_sq));

            // stablility criterion
            // The slowdown factor should not make the system unstable, thus the Qst/Q set the limitation of the increasing of inner semi-major axis.
            if (_bin.stab>0 && _bin.stab != NUMERIC_FLOAT_MAX) {
                Float semi_amplify_max =  std::max(Float(1.0),1.0/_bin.stab);
                Float period_amplify_max = pow(semi_amplify_max,3.0/2.0);
                Float timescale_stab = period_amplify_max*_bin.period;
                _bin.slowdown.timescale = std::min(_bin.slowdown.timescale, timescale_stab);
            }

#else
            _bin.slowdown.timescale = _bin.slowdown.getTimescaleMax();
#endif

            // only set slowdown if semi > 0 and stable
            bool set_sd_flag = true;
            if (_bin.semi>0) {
                for (int k=0; k<2; k++) {
                    if (_bin.isMemberTree(k)) {
                        auto* bink = _bin.getMemberAsTree(k);
                        if (bink->stab>1.0) set_sd_flag = false;
                    }
                }
            }
            else set_sd_flag = false;
            
            if (set_sd_flag) {
                _bin.slowdown.period = _bin.period;
                _bin.slowdown.calcSlowDownFactor();
            }
            else {
                _bin.slowdown.setSlowDownFactor(1.0);
            }
        }

    public:
        //! reset cumulative energy change due to slowdown change
        void resetDESlowDownChangeCum() {
            de_sd_change_cum_ = 0.0;
            dH_sd_change_cum_ = 0.0;
        }

        //! get cumulative energy change due to slowdown change
        Float getDESlowDownChangeCum() const {
            return de_sd_change_cum_;
        }

        //! get cumulative hamiltonian change due to slowdown change
        Float getDHSlowDownChangeCum() const {
            return dH_sd_change_cum_;
        }

        //! reset cumulative energy change due to interruption
        void resetDESlowDownChangeBinaryInterrupt() {
            de_sd_change_interrupt_ = 0.0;
            dH_sd_change_interrupt_ = 0.0;
        }

        //! get cumulative energy change due to interruption
        Float getDESlowDownChangeBinaryInterrupt() const {
            return de_sd_change_interrupt_;
        }

        //! get cumulative hamiltonian change due to interruption
        Float getDHSlowDownChangeBinaryInterrupt() const {
            return dH_sd_change_interrupt_;
        }

        //! Get current kinetic energy with inner slowdown
        /*! \return current kinetic energy with inner slowdown
         */
        Float getEkinSlowDown() const {
            return ekin_sd_;
        }

        //! Get current potential energy with inner slowdown
        /*! \return current potetnial energy with inner slowdown (negative value for bounded systems)
         */
        Float getEpotSlowDown() const {
            return epot_sd_;
        }

        //! Get current total integrated energy with inner slowdown
        /*! \return total integrated energy with inner slowdown
         */
        Float getEtotSlowDownRef() const {
            return etot_sd_ref_;
        }

        //! Get current total energy with inner slowdown from ekin_sdi and epot_sdi
        /*! \return total energy with inner slowdown
         */
        Float getEtotSlowDown() const {
            return ekin_sd_ + epot_sd_;
        }

        //! get energy error with inner slowdown
        /*! \return energy error with inner slowdown
         */
        Float getEnergyErrorSlowDown() const {
            return ekin_sd_ + epot_sd_ - etot_sd_ref_;
        }

        //! get energy error with inner slowdown from backup data 
        Float getEnergyErrorSlowDownFromBackup(Float* _bk) const {
            return -_bk[6] + _bk[7] + _bk[8];
        }

        //! get integrated energy with inner slowdown from backup data
        Float getEtotSlowDownRefFromBackup(Float* _bk) const {
            return _bk[6];
        }

        //! get energy with inner slowdown from backup data (ekin_sdi + epot_sdi)
        Float getEtotSlowDownFromBackup(Float* _bk) const {
            return _bk[7] + _bk[8];
        }

    };

    //! List of slowed binaries and helpers shared by the AR dynamics with Array slow-down, identical for LogH and TTL
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo, class Ttransform>
    class ARArraySlowDownHelpers: public ARSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, Ttransform> {
    protected:
        typedef ARSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, Ttransform> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;
        using Parent::ekin_sd_;
        using Parent::epot_sd_;
        using Parent::etot_sd_ref_;
        using Parent::de_sd_change_cum_;
        using Parent::dH_sd_change_cum_;
        using Parent::de_sd_change_interrupt_;
        using Parent::dH_sd_change_interrupt_;
        using Parent::clearSlowDownEnergy;
        using Parent::copySlowDownEnergy;
        using Parent::calcSlowDownPertInnerBinaryIter;
        using Parent::calcSlowDownInnerBinary;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;
        using Parent::resetDESlowDownChangeCum;
        using Parent::getDESlowDownChangeCum;
        using Parent::getDHSlowDownChangeCum;
        using Parent::resetDESlowDownChangeBinaryInterrupt;
        using Parent::getDESlowDownChangeBinaryInterrupt;
        using Parent::getDHSlowDownChangeBinaryInterrupt;
        using Parent::getEkinSlowDown;
        using Parent::getEpotSlowDown;
        using Parent::getEtotSlowDownRef;
        using Parent::getEtotSlowDown;
        using Parent::getEnergyErrorSlowDown;
        using Parent::getEnergyErrorSlowDownFromBackup;
        using Parent::getEtotSlowDownRefFromBackup;
        using Parent::getEtotSlowDownFromBackup;

    public:
        COMM::List<AR::BinaryTree<Tparticle>*> binary_slowdown; /// binary slowdown, first is root, then others are inner bianries

    protected:
        //! correct postion drift due to inner binary slowdown
        /*! 
          @param[in] _dt: time step
          @param[in] _sd_global_inv: global slowdown factor inverse
         */
        inline void correctPosSlowDownInner(const Float _dt, const Float _sd_global_inv) {
            int n = binary_slowdown.getSize();
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv_m_one = (1.0/kappa - 1.0)*_sd_global_inv;
                Float* velcm = sdi->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi->getMemberIndex(k);
                    ASSERT(j>=0&&j<particles.getSize());
                    Float* pos = particles[j].getPos();
                    Float* vel = particles[j].getVel();

                    // only scale velocity referring to binary c.m.
                    Float vrel[3] = { vel[0] - velcm[0], 
                                      vel[1] - velcm[1], 
                                      vel[2] - velcm[2]}; 
                    pos[0] += _dt * vrel[0] * kappa_inv_m_one;
                    pos[1] += _dt * vrel[1] * kappa_inv_m_one;
                    pos[2] += _dt * vrel[2] * kappa_inv_m_one;
                }
            }
        }

        //! update c.m. for binaries with slowdown inner
        inline void updateCenterOfMassForBinaryWithSlowDownInner() {
            int n = binary_slowdown.getSize();
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                int i1 = sdi->getMemberIndex(0);
                int i2 = sdi->getMemberIndex(1);
                ASSERT(i1>=0&&i1<particles.getSize());
                ASSERT(i2>=0&&i2<particles.getSize());

                Float    m1 = particles[i1].mass;
                Float* pos1 = particles[i1].getPos();
                Float* vel1 = particles[i1].getVel();
                Float    m2 = particles[i2].mass;
                Float* pos2 = particles[i2].getPos();
                Float* vel2 = particles[i2].getVel();
                Float   mcm = m1+m2;

                // first obtain the binary c.m. velocity
                Float mcminv = 1.0/mcm;

                sdi->mass = mcm;
                Float* pos = sdi->getPos();
                pos[0] = (m1*pos1[0] + m2*pos2[0])*mcminv;
                pos[1] = (m1*pos1[1] + m2*pos2[1])*mcminv;
                pos[2] = (m1*pos1[2] + m2*pos2[2])*mcminv;

                Float* vel = sdi->getVel();
                vel[0] = (m1*vel1[0] + m2*vel2[0])*mcminv;
                vel[1] = (m1*vel1[1] + m2*vel2[1])*mcminv;
                vel[2] = (m1*vel1[2] + m2*vel2[2])*mcminv;
            }
        }

        //! calculate kinetic energy with slowdown factor
        /*! @param[in] _ekin: total kinetic energy without slowdown
         */
        inline void calcEkinSlowDownInner(const Float& _ekin) {
            int n = binary_slowdown.getSize();
            if (n==0) return;
            const Float kappa_inv_global = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            Float de = Float(0.0);
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv_m_one = 1.0/kappa - 1.0;
                Float* velcm = sdi->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi->getMemberIndex(k);
                    ASSERT(j>=0&&j<particles.getSize());
                    Float* vel = particles[j].getVel();

                    // only scale velocity referring to binary c.m.
                    Float vrel[3] = { vel[0] - velcm[0], 
                                      vel[1] - velcm[1], 
                                      vel[2] - velcm[2]}; 

                    de += kappa_inv_m_one * particles[j].mass * (vrel[0]*vrel[0] + vrel[1]*vrel[1] + vrel[2]*vrel[2]);
                }
            }
            ekin_sd_ = (_ekin + 0.5*de)*kappa_inv_global;
        }

        //! find inner binaries for slowdown treatment iteration function
        static int findSlowDownInnerBinaryIter(COMM::List<AR::BinaryTree<Tparticle>*>& _binary_slowdown, const int& _c1, const int& _c2, AR::BinaryTree<Tparticle>& _bin) {
            // find leaf binary
            if (_bin.getMemberN()==2 && _bin.semi>0.0) {
                _binary_slowdown.increaseSizeNoInitialize(1);
                auto& sdi_new = _binary_slowdown.getLastMember();
                sdi_new = &_bin;
                return _c1+_c2+1;
            }
            else return _c1+_c2;
        }

        //! find inner binaries for slowdown treatment
        /*! record binary tree address and set update time to _time
          @param[in] _time: next slowdown update time
         */
        void findSlowDownInner(const Float _time) {
            auto& bin_root = info.getBinaryTreeRoot();
            binary_slowdown.resizeNoInitialize(1);
            int ncount[2]={0,0};
            int nsd = bin_root.processTreeIter(binary_slowdown, ncount[0], ncount[1], findSlowDownInnerBinaryIter);
            ASSERT(nsd==binary_slowdown.getSize()-1);

            for (int i=1; i<=nsd; i++) {
#ifdef AR_SLOWDOWN_MASSRATIO
                const Float mass_ratio = manager->slowdown_mass_ref/binary_slowdown[i].bin->mass;
#else 
                const Float mass_ratio = 1.0;
#endif
                binary_slowdown[i]->slowdown.initialSlowDownReference(mass_ratio*manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);
                binary_slowdown[i]->slowdown.setUpdateTime(time_);
            }
        }

        //! Calculate kinetic energy
        inline void calcEKin(){
            ekin_ = Float(0.0);
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            for (int i=0; i<num; i++) {
                const Float *vi=pdat[i].getVel();
                ekin_ += 0.5 * pdat[i].mass * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
            }
            calcEkinSlowDownInner(ekin_);
        }

        //! kick velocity
        /*! First time step will be calculated, the velocities are kicked
          @param[in] _dt: time size
        */
        inline void kickVel(const Float _dt) {
            const int num = particles.getSize();            
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0; i<num; i++) {
                // kick velocity
                Float* vel = pdat[i].getVel();
                Float* acc = force[i].acc_in;
                Float* pert= force[i].acc_pert;
                // half dv 
                vel[0] += _dt * (acc[0] + pert[0]);
                vel[1] += _dt * (acc[1] + pert[1]);
                vel[2] += _dt * (acc[2] + pert[2]);
            }
            // update c.m. of binaries 
            updateCenterOfMassForBinaryWithSlowDownInner();
        }

        //! drift time and position
        /*! First (real) time is drifted, then positions are drifted
          @param[in] _dt: time step
        */
        inline void driftTimeAndPos(const Float _dt) {
            // drift time 
            time_ += _dt;

            // drift position
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            const Float kappa_inv = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            const Float dt_sd = _dt * kappa_inv;

            for (int i=0; i<num; i++) {
                Float* pos = pdat[i].getPos();
                Float* vel = pdat[i].getVel();
                pos[0] += dt_sd * vel[0]; 
                pos[1] += dt_sd * vel[1];
                pos[2] += dt_sd * vel[2];
            }

            // correct postion drift due to inner binary slowdown
            correctPosSlowDownInner(_dt, kappa_inv);
        }

    public:
        //! write back particles with slowdown velocity
        /*! write back particles with slowdown velocity to original address
          @param[in] _particle_cm: center of mass particle to calculate the original frame, different from the particles.cm
         */
        template <class Tptcl>
        void writeBackSlowDownParticles(const Tptcl& _particle_cm) {
            ASSERT(particles.getMode()==COMM::ListMode::copy);
            ASSERT(!particles.isOriginFrame());
            auto* particle_adr = particles.getOriginAddressArray();
            auto* particle_data= particles.getDataAddress();
            const Float kappa_inv = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            for (int i=0; i<particles.getSize(); i++) {
                //ASSERT(particle_adr[i]->mass == particle_data[i].mass);
                particle_adr[i]->mass = particle_data[i].mass;

                particle_adr[i]->pos[0] = particle_data[i].pos[0] + _particle_cm.pos[0];
                particle_adr[i]->pos[1] = particle_data[i].pos[1] + _particle_cm.pos[1];
                particle_adr[i]->pos[2] = particle_data[i].pos[2] + _particle_cm.pos[2];

                particle_adr[i]->vel[0] = particle_data[i].vel[0]*kappa_inv + _particle_cm.vel[0];
                particle_adr[i]->vel[1] = particle_data[i].vel[1]*kappa_inv + _particle_cm.vel[1];
                particle_adr[i]->vel[2] = particle_data[i].vel[2]*kappa_inv + _particle_cm.vel[2];
            }

            // correct inner slowdown velocity
            int nsd= binary_slowdown.getSize();
            for (int i=1; i<nsd; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv_m_one = (1.0/kappa - 1.0)*kappa_inv;
                Float* velcm = sdi->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi->getMemberIndex(k);
                    ASSERT(j>=0&&j<particles.getSize());
                    Float* vel = particle_data[j].getVel();

                    // only scale velocity referring to binary c.m.
                    Float vrel[3] = { vel[0] - velcm[0], 
                                      vel[1] - velcm[1], 
                                      vel[2] - velcm[2]}; 
                    particle_adr[j]->vel[0] += vrel[0] * kappa_inv_m_one;
                    particle_adr[j]->vel[1] += vrel[1] * kappa_inv_m_one;
                    particle_adr[j]->vel[2] += vrel[2] * kappa_inv_m_one;
                }
            }
        }

    };

    //! Helpers shared by the AR dynamics with Tree slow-down, identical for LogH and TTL
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo, class Ttransform>
    class ARTreeSlowDownHelpers: public ARSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, Ttransform> {
    protected:
        typedef ARSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, Ttransform> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;
        using Parent::ekin_sd_;
        using Parent::epot_sd_;
        using Parent::etot_sd_ref_;
        using Parent::de_sd_change_cum_;
        using Parent::dH_sd_change_cum_;
        using Parent::de_sd_change_interrupt_;
        using Parent::dH_sd_change_interrupt_;
        using Parent::clearSlowDownEnergy;
        using Parent::copySlowDownEnergy;
        using Parent::calcSlowDownPertInnerBinaryIter;
        using Parent::calcSlowDownInnerBinary;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;
        using Parent::resetDESlowDownChangeCum;
        using Parent::getDESlowDownChangeCum;
        using Parent::getDHSlowDownChangeCum;
        using Parent::resetDESlowDownChangeBinaryInterrupt;
        using Parent::getDESlowDownChangeBinaryInterrupt;
        using Parent::getDHSlowDownChangeBinaryInterrupt;
        using Parent::getEkinSlowDown;
        using Parent::getEpotSlowDown;
        using Parent::getEtotSlowDownRef;
        using Parent::getEtotSlowDown;
        using Parent::getEnergyErrorSlowDown;
        using Parent::getEnergyErrorSlowDownFromBackup;
        using Parent::getEtotSlowDownRefFromBackup;
        using Parent::getEtotSlowDownFromBackup;

    protected:
        //! Calculate twice (slowdown) kinetic energy iteration function with binary tree
        /*! cumulative ekin_ and ekin_sd_. Notice these two values should be initialized to zero and reduce by two after iteration.
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary tree for kick etot and calc dgt_drift
        */
        void calcTwoEKinIter(const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin){
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
            Float* vel_cm = _bin.getVel();
            for (int k=0; k<2; k++) {
                Float* vk;
                Float  mk;
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    vk = bink->getVel();
                    mk = bink->mass;
                    calcTwoEKinIter(inv_nest_sd, *bink);
                }
                else {
                    auto* pk = _bin.getMember(k);
                    vk = pk->getVel();
                    mk = pk->mass;
                    ekin_ += mk * (vk[0]*vk[0]+vk[1]*vk[1]+vk[2]*vk[2]);
                }
                Float vrel[3] = {vk[0] - vel_cm[0],
                                 vk[1] - vel_cm[1],
                                 vk[2] - vel_cm[2]};
                ekin_sd_ += mk * inv_nest_sd * (vrel[0]*vrel[0] + vrel[1]*vrel[1] + vrel[2]*vrel[2]);
            }
        }

        //! Calculate (slowdown) kinetic energy
        void calcEKin() {
            ekin_ = ekin_sd_ = 0.0;
            auto& bin_root=info.getBinaryTreeRoot();
            Float sd_factor=1.0;
            
            calcTwoEKinIter(sd_factor, bin_root);
            Float* vcm = bin_root.getVel();
            // notice the cm velocity may not be zero after interruption, thus need to be added 
            ekin_sd_ += bin_root.mass*(vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);

            ekin_ *= 0.5;
            ekin_sd_ *= 0.5;
        }

        //! kick velocity
        /*! First time step will be calculated, the velocities are kicked
          @param[in] _dt: time size
        */
        void kickVel(const Float& _dt) {
            const int num = particles.getSize();            
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0; i<num; i++) {
                // kick velocity
                Float* vel = pdat[i].getVel();
                Float* acc = force[i].acc_in;
                Float* pert= force[i].acc_pert;
                // half dv 
                vel[0] += _dt * (acc[0] + pert[0]);
                vel[1] += _dt * (acc[1] + pert[1]);
                vel[2] += _dt * (acc[2] + pert[2]);
            }
            // update binary c.m. velocity interation
            updateBinaryVelIter(info.getBinaryTreeRoot());
        }

        //! drift position with slowdown tree
        /*!
          @param[in] _dt: drift time
          @param[in] _vel_sd_up: upper cm sd vel
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary to drift pos
        */
        void driftPosTreeIter(const Float& _dt, const Float* _vel_sd_up, const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) { 
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();

            Float* vel_cm = _bin.getVel();

            auto driftPos=[&](Float* pos, Float* vel, Float* vel_sd) {
                //scale velocity referring to binary c.m.
                vel_sd[0] = (vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0];
                vel_sd[1] = (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1]; 
                vel_sd[2] = (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2];

                pos[0] += _dt * vel_sd[0];
                pos[1] += _dt * vel_sd[1];
                pos[2] += _dt * vel_sd[2];
            };

            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* pj = _bin.getMemberAsTree(k);
                    Float* pos = pj->getPos();
                    Float* vel = pj->getVel();
                    Float vel_sd[3];
                    driftPos(pos, vel, vel_sd);
                    driftPosTreeIter(_dt, vel_sd, inv_nest_sd, *pj);
//#ifdef AR_DEBUG
//                    auto& pos1= pj->getLeftMember()->pos;
//                    auto& pos2= pj->getRightMember()->pos;
//                    Float m1 = pj->getLeftMember()->mass;
//                    Float m2 = pj->getRightMember()->mass;
//                    Float pos_cm[3] = {(m1*pos1[0]+m2*pos2[0])/(m1+m2),
//                                       (m1*pos1[1]+m2*pos2[1])/(m1+m2),
//                                       (m1*pos1[2]+m2*pos2[2])/(m1+m2)};
//                    Float dpos[3] = {pos[0]-pos_cm[0], 
//                                     pos[1]-pos_cm[1], 
//                                     pos[2]-pos_cm[2]};
//                    Float dr2 = dpos[0]*dpos[0]+dpos[1]*dpos[1]+dpos[2]*dpos[2];
//                    ASSERT(dr2<1e-10);
//#endif
                }
                else {
                    auto* pj = _bin.getMember(k);
                    Float* pos = pj->getPos();
                    Float* vel = pj->getVel();
                    Float vel_sd[3];
                    driftPos(pos, vel, vel_sd);
                }
            }
        }

        //! drift time and position with slowdown tree
        /*! First (real) time is drifted, then positions are drifted
          @param[in] _dt: time step
        */
        void driftTimeAndPos(const Float& _dt) {
            // drift time 
            time_ += _dt;

            // the particle cm velocity is zero (assume in rest-frame)
            ASSERT(!particles.isOriginFrame());
            auto& bin_root=info.getBinaryTreeRoot();
            Float vel_cm[3] = {0.0,0.0,0.0};
            Float sd_factor=1.0;

            driftPosTreeIter(_dt, vel_cm, sd_factor, bin_root);
        }

        /* Reference only: time function from multiplied inverse separations (not compiled).
           This variant of the Tree slow-down TTL time function replaced the kick and drift
           factors by the product of the inverse separations of every binary in the hierarchy.
           It was guarded by AR_TIME_FUNCTION_MULTI_R and never compiled (undefined names
           _bini, calInvR, calcGTDriftInvIter). Its fragments are collected here.

           Helpers:

        //! calc inverse R
        Float calcInvR(Tparticle& _p1, Tparticle& _p2) {
            Float dr[3] = {_p1.pos[0] - _p2.pos[0], 
                           _p1.pos[1] - _p2.pos[1], 
                           _p1.pos[2] - _p2.pos[2]};
            Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            Float invr = 1/sqrt(r2);
            return invr;
        }

        //! calc multiplied inverse R of binary tree
        Float calcMultiInvRIter(AR::BinaryTree<Tparticle>& _bin){
            Float gt_kick_inv= calcInvR(*_bin.getLeftMember(), *_bin.getRightMember());
            for (int k=0; k<2; k++) { 
                if (_bini.isMemberTree(k))  // tree - tree
                    gt_kick_inv *= calcMultiInvRIter(*(_bini.getMemberAsTree(k)));
            }
            return gt_kick_inv;
        }

        //! calc gt_drift_inv for binary tree
        Float calcGtDriftInvIter(AR::BinaryTree<Tparticle>& _bin){
            Float gt_drift_inv= calInvR(*_bin.getLeftMember(), *_bin.getRightMember());
            for (int k=0; k<2; k++) { 
                if (_bini.isMemberTree(k))  // tree - tree
                    gt_kick_inv *= calcMultiInvRIter(*(_bini.getMemberAsTree(k)));
            }
            return gt_kick_inv;
        }

           At the end of calcAccPotAndGTKickInvTreeIter, before returning gt_kick_inv:

            gt_kick_inv = calcMultiInvRIter(_bin);

           In kickEtotAndGTDriftTreeIter, the accumulation
           dgt_drift_inv += vel_sd . gtgrad was skipped.

           In kickEtotAndGTDrift, after the tree iteration:

            dgt_drift_inv = calcGTDriftInvIter(bin_root);
        */

    public:
        //! write back slowdown particle iteration function
        /*!
          @param[in] _particle_cm: center of mass particle to calculate the original frame, different from the particles.cm
          @param[in] _vel_sd_up: upper slowdown cm velocity
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary
         */
        template <class Tptcl>
        void writeBackSlowDownParticlesIter(const Tptcl& _particle_cm, const Float* _vel_sd_up, const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) {
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
            Float* vel_cm = _bin.getVel();
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    Float* vel = bink->getVel();
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]}; 
                    writeBackSlowDownParticlesIter(_particle_cm, vel_sd, inv_nest_sd, *bink);
                }
                else {
                    int i = _bin.getMemberIndex(k);
                    auto& pk = particles[i];
                    auto* pk_adr = particles.getMemberOriginAddress(i);
                    Float* vel = pk.getVel();
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]};
                    pk_adr->mass = pk.mass;

                    pk_adr->pos[0] = pk.pos[0] + _particle_cm.pos[0];
                    pk_adr->pos[1] = pk.pos[1] + _particle_cm.pos[1];
                    pk_adr->pos[2] = pk.pos[2] + _particle_cm.pos[2];

                    pk_adr->vel[0] = vel_sd[0] + _particle_cm.vel[0];
                    pk_adr->vel[1] = vel_sd[1] + _particle_cm.vel[1];
                    pk_adr->vel[2] = vel_sd[2] + _particle_cm.vel[2];
                }
            }
        }

        //! write back particles with slowdown velocity
        /*! write back particles with slowdown velocity to original address
          @param[in] _particle_cm: center of mass particle to calculate the original frame, different from the particles.cm
         */
        template <class Tptcl>
        void writeBackSlowDownParticles(const Tptcl& _particle_cm) {
            //! iteration function using binarytree
            auto& bin_root=info.getBinaryTreeRoot();
            Float vel_cm[3] = {0.0,0.0,0.0};
            Float sd_factor=1.0;
            writeBackSlowDownParticlesIter(_particle_cm, vel_cm, sd_factor, bin_root);
        }

    };

    //! AR dynamics: the equations of motion as AR advances them
    /*! One explicit specialization for each time transformation (LogH or TTL) and Slow-down scheme
        (NoSlowDown, ArraySlowDown or TreeSlowDown). Each holds the state that exists only for its variant,
        kinetic energy evaluation, the velocity kick, the position and time drift, force and potential
        evaluation, the energy kick, the Slow-down update with its energy correction, and the integration of one step.
     */
    template <class Ttransform, class Tslowdown, class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo>
    class ARDynamics;

    //! AR dynamics: LogH without Slow-down
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo>
    class ARDynamics<LogH, NoSlowDown, Tparticle, Tpcm, Tpert, Tmethod, Tinfo>: public ARNoSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, LogH> {
    protected:
        typedef ARNoSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, LogH> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;
        using Parent::calcEKin;
        using Parent::kickVel;
        using Parent::driftTimeAndPos;
        using Parent::calcAccPotAndGTKickInv;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;

    protected:
        //! reserve memory of the data owned by the dynamics
        void reserveDynamicsMem(const int _nmax) {
            (void)_nmax;
        }

        //! clear the data owned by the dynamics
        void clearDynamics() {
        }

        //! copy the data owned by the dynamics
        void copyDynamics(const ARDynamics& _sym) {
            (void)_sym;
        }

        //! kick energy
        /*!
          @param[in] _dt: time step
        */
        inline void kickEtot(const Float _dt) {
            Float de = 0.0;
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                de += mass * (vel[0] * pert[0] +
                              vel[1] * pert[1] +
                              vel[2] * pert[2]);
            }
            etot_ref_ += _dt * de;
        }

    public:
        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time: real physical time to initialize
        */
        void initialIntegration(const Float _time) {
            ASSERT(checkParams());

            // particle number and data address
            const int n_particle = particles.getSize();

            // resize force array
            force_.resizeNoInitialize(n_particle);

            // Initial intgrt value t (avoid confusion of real time when slowdown is used)
            time_ = _time;

            // check particle number
            ASSERT(particles.getSize()>=2);

            // reset particle modification flag
            particles.setModifiedFalse();

            // check the center-of-mass initialization
            if(particles.isOriginFrame()) {
                particles.calcCenterOfMass();
                particles.shiftToCenterOfMassFrame();
                for (int i=0; i<info.binarytree.getSize(); i++) {
                    auto& bin = info.binarytree[i];
                    bin.pos[0] -= particles.cm.pos[0];
                    bin.pos[1] -= particles.cm.pos[1];
                    bin.pos[2] -= particles.cm.pos[2];
                    bin.vel[0] -= particles.cm.vel[0];
                    bin.vel[1] -= particles.cm.vel[1];
                    bin.vel[2] -= particles.cm.vel[2];
                }
            }
            ASSERT(info.getBinaryTreeRoot().pos[0]*info.getBinaryTreeRoot().pos[0]<1e-10);
            ASSERT(info.getBinaryTreeRoot().vel[0]*info.getBinaryTreeRoot().vel[0]<1e-10);

            Tparticle* particle_data = particles.getDataAddress();
            Force* force_data = force_.getDataAddress();

            manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm,  perturber, _time);

            // calculate kinetic energy
            calcEKin();

            // initial total energy
            etot_ref_ = ekin_ + epot_;

        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // inverse time transformation factor for drift
                Float gt_drift_inv = manager->interaction.calcGTDriftInv(ekin_-etot_ref_); // pt = -etot

                // drift
                Float dt_drift = ds_drift/gt_drift_inv;

                // drift time and postion
                driftTimeAndPos(dt_drift);
                _time_table[i] = time_;

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                //! calc force, potential and inverse time transformation factor for kick
                Float gt_kick_inv = calcAccPotAndGTKickInv();

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // kick total energy 
                kickEtot(dt_kick);
                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // calculate kinetic energy
                calcEKin();
            }
        }

        //! integration for two body one step
        /*! For two-body problem the calculation can be much symplified to improve performance. 
          Besides, the slow-down factor calculation is embedded in the Drift (for time) and Kick (for perturbation). 
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().         
        */
        void integrateTwoOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group number
            const int nloop = manager->step.getCDPairSize();
            
            const int n_particle = particles.getSize();
            ASSERT(n_particle==2);


            Tparticle* particle_data = particles.getDataAddress();
            Float mass1 = particle_data[0].mass;
            Float* pos1 = particle_data[0].getPos();
            Float* vel1 = particle_data[0].getVel();

            Float mass2 = particle_data[1].mass;
            Float* pos2 = particle_data[1].getPos();
            Float* vel2 = particle_data[1].getVel();

            Force* force_data = force_.getDataAddress();
            Float* acc1 = force_data[0].acc_in;
            Float* pert1= force_data[0].acc_pert;

            Float* acc2 = force_data[1].acc_in;
            Float* pert2= force_data[1].acc_pert;

#ifdef AR_DEBUG_PRINT_DKD
            std::cout<<"K "<<time_<<" "
                     <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                     <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                     <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds = manager->step.getCK(i)*_ds;
                // inverse time transformation factor for drift
                Float gt_inv = manager->interaction.calcGTDriftInv(ekin_-etot_ref_); // pt = -etot
                // drift
                Float dt = ds/gt_inv;
                ASSERT(!ISNAN(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                _time_table[i] = time_;

                // drift position
                pos1[0] += dt * vel1[0];
                pos1[1] += dt * vel1[1];
                pos1[2] += dt * vel1[2];

                pos2[0] += dt * vel2[0];
                pos2[1] += dt * vel2[1];
                pos2[2] += dt * vel2[2];

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                gt_inv = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm, perturber, _time_table[i]);

                ASSERT(!ISNAN(epot_));

#ifdef AR_DEBUG_PRINT_DKD
                if (i>0)
                    std::cout<<"K "<<time_<<" "
                             <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                             <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                             <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick half step for velocity
                Float dvel1[3], dvel2[3];

                dt = 0.5*ds/gt_inv;

                dvel1[0] = dt * (acc1[0] + pert1[0]);
                dvel1[1] = dt * (acc1[1] + pert1[1]);
                dvel1[2] = dt * (acc1[2] + pert1[2]);

                dvel2[0] = dt * (acc2[0] + pert2[0]);
                dvel2[1] = dt * (acc2[1] + pert2[1]);
                dvel2[2] = dt * (acc2[2] + pert2[2]);

                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];

                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];

#ifdef AR_DEBUG_PRINT_DKD
                std::cout<<"D "<<time_<<" "
                         <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                         <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                         <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick total energy and time transformation factor for drift
                etot_ref_ += 2.0*dt * (mass1* (vel1[0] * pert1[0] + 
                                               vel1[1] * pert1[1] + 
                                               vel1[2] * pert1[2]) +
                                       mass2* (vel2[0] * pert2[0] + 
                                               vel2[1] * pert2[1] + 
                                               vel2[2] * pert2[2]));


                // kick half step for velocity
                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];
                
                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];
                                                                                                                
                // calculate kinetic energy
                ekin_ = 0.5 * (mass1 * (vel1[0]*vel1[0]+vel1[1]*vel1[1]+vel1[2]*vel1[2]) +
                               mass2 * (vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]));

            }

        }

        //! get Hamiltonian
        Float getH() const {
            return manager->interaction.calcH(ekin_ - etot_ref_, epot_);
        }

        //! get Hamiltonian from backup data
        Float getHFromBackup(Float* _bk) const {
            Float& etot_ref =_bk[1];
            Float& ekin = _bk[2];
            Float& epot = _bk[3];
            return manager->interaction.calcH(ekin - etot_ref, epot);
        }

    };

    //! AR dynamics: TTL without Slow-down
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo>
    class ARDynamics<TTL, NoSlowDown, Tparticle, Tpcm, Tpert, Tmethod, Tinfo>: public ARNoSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, TTL> {
    protected:
        typedef ARNoSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, TTL> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;
        using Parent::calcEKin;
        using Parent::kickVel;
        using Parent::driftTimeAndPos;
        using Parent::calcAccPotAndGTKickInv;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;

    protected:
        // transformation factors
        Float gt_drift_inv_ = 0.0;  ///< integrated inverse time transformation factor for drift: dt(drift) = ds/gt_drift_inv_
        Float gt_kick_inv_ = 0.0;   ///< inverse time transformation factor for kick: dt(kick) = ds/gt_kick_inv_

        //! reserve memory of the data owned by the dynamics
        void reserveDynamicsMem(const int _nmax) {
            (void)_nmax;
        }

        //! clear the data owned by the dynamics
        void clearDynamics() {
            gt_drift_inv_ = 0.0;
            gt_kick_inv_ = 0.0;
        }

        //! copy the data owned by the dynamics
        void copyDynamics(const ARDynamics& _sym) {
            gt_drift_inv_ = _sym.gt_drift_inv_;
            gt_kick_inv_ = _sym.gt_kick_inv_;
        }

        //! kick energy and time transformation function for drift
        /*!
          @param[in] _dt: time step
        */

        inline void kickEtotAndGTDrift(const Float _dt) {
            Float de = Float(0.0);
            Float dg = Float(0.0);
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                Float* gtgrad=force[i].gtgrad;
                de += mass * (vel[0] * pert[0] +
                              vel[1] * pert[1] +
                              vel[2] * pert[2]);
                dg +=  (vel[0] * gtgrad[0] +
                        vel[1] * gtgrad[1] +
                        vel[2] * gtgrad[2]);
            }
            etot_ref_ += _dt * de;

            gt_drift_inv_ += _dt * dg;
        }

    public:
        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time: real physical time to initialize
        */
        void initialIntegration(const Float _time) {
            ASSERT(checkParams());

            // particle number and data address
            const int n_particle = particles.getSize();

            // resize force array
            force_.resizeNoInitialize(n_particle);

            // Initial intgrt value t (avoid confusion of real time when slowdown is used)
            time_ = _time;

            // check particle number
            ASSERT(particles.getSize()>=2);

            // reset particle modification flag
            particles.setModifiedFalse();

            // check the center-of-mass initialization
            if(particles.isOriginFrame()) {
                particles.calcCenterOfMass();
                particles.shiftToCenterOfMassFrame();
                for (int i=0; i<info.binarytree.getSize(); i++) {
                    auto& bin = info.binarytree[i];
                    bin.pos[0] -= particles.cm.pos[0];
                    bin.pos[1] -= particles.cm.pos[1];
                    bin.pos[2] -= particles.cm.pos[2];
                    bin.vel[0] -= particles.cm.vel[0];
                    bin.vel[1] -= particles.cm.vel[1];
                    bin.vel[2] -= particles.cm.vel[2];
                }
            }
            ASSERT(info.getBinaryTreeRoot().pos[0]*info.getBinaryTreeRoot().pos[0]<1e-10);
            ASSERT(info.getBinaryTreeRoot().vel[0]*info.getBinaryTreeRoot().vel[0]<1e-10);

            Tparticle* particle_data = particles.getDataAddress();
            Force* force_data = force_.getDataAddress();

            gt_kick_inv_ = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm,  perturber, _time);

            // initially gt_drift 
            gt_drift_inv_ = gt_kick_inv_;


            // calculate kinetic energy
            calcEKin();

            // initial total energy
            etot_ref_ = ekin_ + epot_;

        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // inverse time transformation factor for drift
                Float gt_drift_inv = gt_drift_inv_;

                // drift
                Float dt_drift = ds_drift/gt_drift_inv;

                // drift time and postion
                driftTimeAndPos(dt_drift);
                _time_table[i] = time_;

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                //! calc force, potential and inverse time transformation factor for kick
                Float gt_kick_inv = calcAccPotAndGTKickInv();

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // back up gt_kick 
                gt_kick_inv_ = gt_kick_inv;
                // kick total energy and inverse time transformation factor for drift
                kickEtotAndGTDrift(dt_kick);
                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // calculate kinetic energy
                calcEKin();
            }
        }

        //! integration for two body one step
        /*! For two-body problem the calculation can be much symplified to improve performance. 
          Besides, the slow-down factor calculation is embedded in the Drift (for time) and Kick (for perturbation). 
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().         
        */
        void integrateTwoOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group number
            const int nloop = manager->step.getCDPairSize();
            
            const int n_particle = particles.getSize();
            ASSERT(n_particle==2);


            Tparticle* particle_data = particles.getDataAddress();
            Float mass1 = particle_data[0].mass;
            Float* pos1 = particle_data[0].getPos();
            Float* vel1 = particle_data[0].getVel();

            Float mass2 = particle_data[1].mass;
            Float* pos2 = particle_data[1].getPos();
            Float* vel2 = particle_data[1].getVel();

            Force* force_data = force_.getDataAddress();
            Float* acc1 = force_data[0].acc_in;
            Float* pert1= force_data[0].acc_pert;

            Float* acc2 = force_data[1].acc_in;
            Float* pert2= force_data[1].acc_pert;
            Float* gtgrad1 = force_data[0].gtgrad;
            Float* gtgrad2 = force_data[1].gtgrad;

#ifdef AR_DEBUG_PRINT_DKD
            std::cout<<"K "<<time_<<" "
                     <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                     <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                     <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds = manager->step.getCK(i)*_ds;
                // inverse time transformation factor for drift
                Float gt_inv = gt_drift_inv_;
                // drift
                Float dt = ds/gt_inv;
                ASSERT(!ISNAN(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                _time_table[i] = time_;

                // drift position
                pos1[0] += dt * vel1[0];
                pos1[1] += dt * vel1[1];
                pos1[2] += dt * vel1[2];

                pos2[0] += dt * vel2[0];
                pos2[1] += dt * vel2[1];
                pos2[2] += dt * vel2[2];

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                gt_inv = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm, perturber, _time_table[i]);

                ASSERT(!ISNAN(epot_));

#ifdef AR_DEBUG_PRINT_DKD
                if (i>0)
                    std::cout<<"K "<<time_<<" "
                             <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                             <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                             <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick half step for velocity
                Float dvel1[3], dvel2[3];

                dt = 0.5*ds/gt_inv;

                dvel1[0] = dt * (acc1[0] + pert1[0]);
                dvel1[1] = dt * (acc1[1] + pert1[1]);
                dvel1[2] = dt * (acc1[2] + pert1[2]);

                dvel2[0] = dt * (acc2[0] + pert2[0]);
                dvel2[1] = dt * (acc2[1] + pert2[1]);
                dvel2[2] = dt * (acc2[2] + pert2[2]);

                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];

                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];

#ifdef AR_DEBUG_PRINT_DKD
                std::cout<<"D "<<time_<<" "
                         <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                         <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                         <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick total energy and time transformation factor for drift
                etot_ref_ += 2.0*dt * (mass1* (vel1[0] * pert1[0] + 
                                               vel1[1] * pert1[1] + 
                                               vel1[2] * pert1[2]) +
                                       mass2* (vel2[0] * pert2[0] + 
                                               vel2[1] * pert2[1] + 
                                               vel2[2] * pert2[2]));

                // back up gt_kick_inv
                gt_kick_inv_ = gt_inv;

                // integrate gt_drift_inv
                gt_drift_inv_ +=  2.0*dt* (vel1[0] * gtgrad1[0] +
                                           vel1[1] * gtgrad1[1] +
                                           vel1[2] * gtgrad1[2] +
                                           vel2[0] * gtgrad2[0] +
                                           vel2[1] * gtgrad2[1] +
                                           vel2[2] * gtgrad2[2]);


                // kick half step for velocity
                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];
                
                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];
                                                                                                                
                // calculate kinetic energy
                ekin_ = 0.5 * (mass1 * (vel1[0]*vel1[0]+vel1[1]*vel1[1]+vel1[2]*vel1[2]) +
                               mass2 * (vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]));

            }

        }

        //! get Hamiltonian
        Float getH() const {
            //return (ekin_ - etot_ref_)/gt_drift_inv_ + epot_/gt_kick_inv_;
            return (ekin_ + epot_ - etot_ref_)/gt_kick_inv_;
        }

        //! get Hamiltonian from backup data
        Float getHFromBackup(Float* _bk) const {
            Float& etot_ref =_bk[1];
            Float& ekin = _bk[2];
            Float& epot = _bk[3];
            //Float& gt_drift_inv = _bk[6];
            Float& gt_kick_inv  = _bk[7];
            return (ekin + epot - etot_ref)/gt_kick_inv;
            //return (ekin - etot_ref)/gt_drift_inv + epot/gt_kick_inv;
        }

        //! Get integrated inverse time transformation factor
        /*! In TTF case, it is calculated by integrating \f$ \frac{dg}{dt} = \sum_k \frac{\partial g}{\partial \vec{r_k}} \bullet \vec{v_k} \f$.
          Notice last step is the sub-step in one symplectic loop
          \return inverse time transformation factor for drift
        */
        Float getGTDriftInv() const {
            return gt_drift_inv_;
        }

    };

    //! AR dynamics: LogH with Array slow-down
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo>
    class ARDynamics<LogH, ArraySlowDown, Tparticle, Tpcm, Tpert, Tmethod, Tinfo>: public ARArraySlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, LogH> {
    protected:
        typedef ARArraySlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, LogH> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;
        using Parent::ekin_sd_;
        using Parent::epot_sd_;
        using Parent::etot_sd_ref_;
        using Parent::de_sd_change_cum_;
        using Parent::dH_sd_change_cum_;
        using Parent::de_sd_change_interrupt_;
        using Parent::dH_sd_change_interrupt_;
        using Parent::clearSlowDownEnergy;
        using Parent::copySlowDownEnergy;
        using Parent::calcSlowDownPertInnerBinaryIter;
        using Parent::calcSlowDownInnerBinary;
        using Parent::correctPosSlowDownInner;
        using Parent::updateCenterOfMassForBinaryWithSlowDownInner;
        using Parent::calcEkinSlowDownInner;
        using Parent::findSlowDownInnerBinaryIter;
        using Parent::findSlowDownInner;
        using Parent::calcEKin;
        using Parent::kickVel;
        using Parent::driftTimeAndPos;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;
        using Parent::resetDESlowDownChangeCum;
        using Parent::getDESlowDownChangeCum;
        using Parent::getDHSlowDownChangeCum;
        using Parent::resetDESlowDownChangeBinaryInterrupt;
        using Parent::getDESlowDownChangeBinaryInterrupt;
        using Parent::getDHSlowDownChangeBinaryInterrupt;
        using Parent::getEkinSlowDown;
        using Parent::getEpotSlowDown;
        using Parent::getEtotSlowDownRef;
        using Parent::getEtotSlowDown;
        using Parent::getEnergyErrorSlowDown;
        using Parent::getEnergyErrorSlowDownFromBackup;
        using Parent::getEtotSlowDownRefFromBackup;
        using Parent::getEtotSlowDownFromBackup;
        using Parent::binary_slowdown;
        using Parent::writeBackSlowDownParticles;

    protected:
        //! reserve memory of the data owned by the dynamics
        void reserveDynamicsMem(const int _nmax) {
            binary_slowdown.setMode(COMM::ListMode::local);
            binary_slowdown.reserveMem(_nmax/2+1);
        }

        //! clear the data owned by the dynamics
        void clearDynamics() {
            this->clearSlowDownEnergy();
            binary_slowdown.clear();
        }

        //! copy the data owned by the dynamics
        void copyDynamics(const ARDynamics& _sym) {
            this->copySlowDownEnergy(_sym);
            binary_slowdown = _sym.binary_slowdown;
        }

        //! correct force, potential energy and gt_kick_inv based on slowdown for inner binaries
        /*! 
          @param[in,out] _gt_kick_inv: the inverse time transformation factor for kick step (input), be corrected with slowdown (output)
         */
        inline void correctAccPotGTKickInvSlowDownInner(Float& _gt_kick_inv) {
            int n = binary_slowdown.getSize();
            Float gt_kick_inv_cor = 0.0;
            Float de = 0.0;
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                int i1 = sdi->getMemberIndex(0);
                int i2 = sdi->getMemberIndex(1);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv = 1.0/kappa;
                if (i1>=0) {
                    ASSERT(i2>=0);
                    ASSERT(i1!=i2);
                    ASSERT(i1<particles.getSize());
                    ASSERT(i2<particles.getSize());
                    
                    // calculate pair interaction
                    Force fi[2];
                    Float epoti;
                    Float gt_kick_inv_i = manager->interaction.calcInnerAccPotAndGTKickInvTwo(fi[0], fi[1], epoti, particles[i1], particles[i2]);
                                    // scale binary pair force with slowdown 
                    force_[i1].acc_in[0] += fi[0].acc_in[0]*kappa_inv - fi[0].acc_in[0];
                    force_[i1].acc_in[1] += fi[0].acc_in[1]*kappa_inv - fi[0].acc_in[1];
                    force_[i1].acc_in[2] += fi[0].acc_in[2]*kappa_inv - fi[0].acc_in[2];
                    force_[i2].acc_in[0] += fi[1].acc_in[0]*kappa_inv - fi[1].acc_in[0];
                    force_[i2].acc_in[1] += fi[1].acc_in[1]*kappa_inv - fi[1].acc_in[1];
                    force_[i2].acc_in[2] += fi[1].acc_in[2]*kappa_inv - fi[1].acc_in[2];

                    de += epoti*kappa_inv - epoti;
                    // scale gtgrad with slowdown

                    // gt kick
                    gt_kick_inv_cor += gt_kick_inv_i*(kappa_inv - 1.0);
                }
            }

            // global slowdown
            const Float kappa_inv_global = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            for (int i=0; i<force_.getSize(); i++) {
                force_[i].acc_in[0] *= kappa_inv_global;
                force_[i].acc_in[1] *= kappa_inv_global;
                force_[i].acc_in[2] *= kappa_inv_global;
            }

            epot_sd_ = (epot_ + de)*kappa_inv_global;
            _gt_kick_inv = (_gt_kick_inv + gt_kick_inv_cor)*kappa_inv_global;
        }

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        inline Float calcAccPotAndGTKickInv() {
            Float gt_kick_inv = manager->interaction.calcAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());            

//#ifdef AR_DEBUG
//            // check c.m. force 
//            Force fcm;
//            Float mcm=0.0;
//            for (int i=0; i<particles.getSize(); i++) {
//                for (int k=0; k<3; k++) {
//                    fcm.acc_in[k] += particles[i].mass * force_[i].acc_in[k];
//                }
//                mcm += particles[i].mass;
//            }
//            for (int k=0; k<3; k++) {
//                fcm.acc_in[k] /= mcm;
//                ASSERT(abs(fcm.acc_in[k])<1e-10);
//            }
//#endif
            // slowdown binary acceleration
            correctAccPotGTKickInvSlowDownInner(gt_kick_inv);
//#ifdef AR_DEBUG
//            // check c.m. force 
//            fcm.acc_in[0] = fcm.acc_in[1] = fcm.acc_in[2] = 0.0;
//            for (int i=0; i<particles.getSize(); i++) {
//                for (int k=0; k<3; k++) {
//                    fcm.acc_in[k] += particles[i].mass * force_[i].acc_in[k];
//                }
//            }
//            for (int k=0; k<3; k++) {
//                fcm.acc_in[k] /= mcm;
//                ASSERT(abs(fcm.acc_in[k])<1e-10);
//            }
//#endif
            return gt_kick_inv;
        }

        //! kick energy
        /*!
          @param[in] _dt: time step
        */
        inline void kickEtot(const Float _dt) {
            Float de = 0.0;
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                de += mass * (vel[0] * pert[0] +
                              vel[1] * pert[1] +
                              vel[2] * pert[2]);
            }
            etot_sd_ref_ += _dt * de;
            etot_ref_ += _dt * de;
        }

    public:
        //! update slowdown factor based on perturbation and record slowdown energy change
        /*! Update slowdown inner and global, update gt_inv
            @param [in] _update_energy_flag: Record cumulative slowdown energy change if true;
            @param [in] _stable_check_flag: check whether the binary tree is stable if true;
         */
        void updateSlowDownAndCorrectEnergy(const bool _update_energy_flag, const bool _stable_check_flag) {
            auto& bin_root = *binary_slowdown[0];
            auto& sd_root = bin_root.slowdown;

            // when the maximum inner slowdown is large, the outer should not be slowed down since the system may not be stable.
            //if (inner_sd_change_flag&&sd_org_inner_max<1000.0*manager->slowdown_pert_ratio_ref) sd_root.setSlowDownFactor(1.0);
            //if (time_>=sd_root.getUpdateTime()) {
            sd_root.pert_in = manager->interaction.calcPertFromBinary(bin_root);
            sd_root.pert_out = 0.0;
            Float t_min_sq= NUMERIC_FLOAT_MAX;
            manager->interaction.calcSlowDownPert(sd_root.pert_out, t_min_sq, getTime(), particles.cm, perturber);
            sd_root.timescale = std::min(sd_root.getTimescaleMax(), sqrt(t_min_sq));

            //Float period_amplify_max = NUMERIC_FLOAT_MAX;
            if (_stable_check_flag) {
                // check whether the system is stable for 10000 out period
                Float stab = bin_root.stableCheckIter(bin_root,10000*bin_root.period);
                if (stab<1.0) {
                    sd_root.period = bin_root.period;
                    sd_root.calcSlowDownFactor();
                }
                else sd_root.setSlowDownFactor(1.0);
                // stablility criterion
                // The slowdown factor should not make the system unstable, thus the Qst/Q set the limitation of the increasing of inner semi-major axis.
                //if (stab>0) {
                //    Float semi_amplify_max =  std::max(Float(1.0),1.0/stab);
                //    period_amplify_max = pow(semi_amplify_max,3.0/2.0);
                //}
            }
            else if (bin_root.semi>0) {
                sd_root.period = bin_root.period;
                sd_root.calcSlowDownFactor();
            }
            else sd_root.setSlowDownFactor(1.0);

            //sd_root.increaseUpdateTimeOnePeriod();
            //}


            // inner binary slowdown
            int n = binary_slowdown.getSize();
            bool modified_flag=false;
            for (int i=1; i<n; i++) {
                auto* sdi = binary_slowdown[i];
                //if (time_>=sdi->slowdown.getUpdateTime()) {
                sdi->calcCenterOfMass();
                calcSlowDownInnerBinary(*sdi);
                //sdi->slowdown.increaseUpdateTimeOnePeriod();
                modified_flag=true;
                //}
            }    

            if (_update_energy_flag) {
                Float ekin_sd_bk = ekin_sd_;
                Float epot_sd_bk = epot_sd_;
                Float H_sd_bk = getHSlowDown();

                if (modified_flag) {
                    // initialize the gt_drift_inv_ with new slowdown factor
                    Float gt_kick_inv_sdi = manager->interaction.calcAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());
                    correctAccPotGTKickInvSlowDownInner(gt_kick_inv_sdi);
                    // correct etot_sd_ref_ with new slowdown
                    calcEkinSlowDownInner(ekin_);
                }
                else {
                    Float kappa_inv = 1.0/sd_root.getSlowDownFactor();
                    // only need to correct the total value
                    ekin_sd_ = ekin_*kappa_inv;
                    epot_sd_ = epot_*kappa_inv;
                }

                Float de_sd = (ekin_sd_ - ekin_sd_bk) + (epot_sd_ - epot_sd_bk);
                etot_sd_ref_ += de_sd;

                Float dH_sd = getHSlowDown() - H_sd_bk;

                // add slowdown change to the global slowdown energy
                de_sd_change_cum_ += de_sd;
                dH_sd_change_cum_ += dH_sd;
            }
        }

        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time: real physical time to initialize
        */
        void initialIntegration(const Float _time) {
            ASSERT(checkParams());

            // particle number and data address
            const int n_particle = particles.getSize();

            // resize force array
            force_.resizeNoInitialize(n_particle);

            // Initial intgrt value t (avoid confusion of real time when slowdown is used)
            time_ = _time;

            // check particle number
            ASSERT(particles.getSize()>=2);

            // reset particle modification flag
            particles.setModifiedFalse();

            // check the center-of-mass initialization
            if(particles.isOriginFrame()) {
                particles.calcCenterOfMass();
                particles.shiftToCenterOfMassFrame();
                for (int i=0; i<info.binarytree.getSize(); i++) {
                    auto& bin = info.binarytree[i];
                    bin.pos[0] -= particles.cm.pos[0];
                    bin.pos[1] -= particles.cm.pos[1];
                    bin.pos[2] -= particles.cm.pos[2];
                    bin.vel[0] -= particles.cm.vel[0];
                    bin.vel[1] -= particles.cm.vel[1];
                    bin.vel[2] -= particles.cm.vel[2];
                }
            }
            ASSERT(info.getBinaryTreeRoot().pos[0]*info.getBinaryTreeRoot().pos[0]<1e-10);
            ASSERT(info.getBinaryTreeRoot().vel[0]*info.getBinaryTreeRoot().vel[0]<1e-10);

            binary_slowdown.increaseSizeNoInitialize(1);
            binary_slowdown[0] = &info.getBinaryTreeRoot();

            // set slowdown reference
            SlowDown& slowdown_root = info.getBinaryTreeRoot().slowdown;

            // slowdown for the system
            slowdown_root.initialSlowDownReference(manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);

            if (particles.getSize()>2) {
                findSlowDownInner(time_);
                // update c.m. of binaries 
                //updateCenterOfMassForBinaryWithSlowDownInner();
            }

            updateSlowDownAndCorrectEnergy(false,true);

            calcAccPotAndGTKickInv();

            calcEKin();

            etot_ref_ = ekin_ + epot_;
            etot_sd_ref_ = ekin_sd_ + epot_sd_;

            Float de_sd = etot_sd_ref_ - etot_ref_;

            // add slowdown change to the global slowdown energy
            de_sd_change_cum_ += de_sd;
            dH_sd_change_cum_ = 0.0;

        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // inverse time transformation factor for drift
                Float gt_drift_inv = manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_); // pt = -etot

                // drift
                Float dt_drift = ds_drift/gt_drift_inv;

                // drift time and postion
                driftTimeAndPos(dt_drift);
                _time_table[i] = time_;

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                //! calc force, potential and inverse time transformation factor for kick
                Float gt_kick_inv = calcAccPotAndGTKickInv();

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // kick total energy 
                kickEtot(dt_kick);
                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // calculate kinetic energy
                calcEKin();
            }
        }

        //! integration for two body one step
        /*! For two-body problem the calculation can be much symplified to improve performance. 
          Besides, the slow-down factor calculation is embedded in the Drift (for time) and Kick (for perturbation). 
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().         
        */
        void integrateTwoOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group number
            const int nloop = manager->step.getCDPairSize();
            
            const int n_particle = particles.getSize();
            ASSERT(n_particle==2);

            const Float kappa_inv = 1.0/info.getBinaryTreeRoot().slowdown.getSlowDownFactor();

            Tparticle* particle_data = particles.getDataAddress();
            Float mass1 = particle_data[0].mass;
            Float* pos1 = particle_data[0].getPos();
            Float* vel1 = particle_data[0].getVel();

            Float mass2 = particle_data[1].mass;
            Float* pos2 = particle_data[1].getPos();
            Float* vel2 = particle_data[1].getVel();

            Force* force_data = force_.getDataAddress();
            Float* acc1 = force_data[0].acc_in;
            Float* pert1= force_data[0].acc_pert;

            Float* acc2 = force_data[1].acc_in;
            Float* pert2= force_data[1].acc_pert;

#ifdef AR_DEBUG_PRINT_DKD
            std::cout<<"K "<<time_<<" "
                     <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                     <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                     <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds = manager->step.getCK(i)*_ds;
                // inverse time transformation factor for drift
                Float gt_inv = manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_); // pt = -etot_sd
                // drift
                Float dt = ds/gt_inv;
                ASSERT(!ISNAN(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                _time_table[i] = time_;

                Float dt_sd = dt*kappa_inv;

                // drift position
                pos1[0] += dt_sd * vel1[0];
                pos1[1] += dt_sd * vel1[1];
                pos1[2] += dt_sd * vel1[2];

                pos2[0] += dt_sd * vel2[0];
                pos2[1] += dt_sd * vel2[1];
                pos2[2] += dt_sd * vel2[2];

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                gt_inv = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm, perturber, _time_table[i]);

                ASSERT(!ISNAN(epot_));

#ifdef AR_DEBUG_PRINT_DKD
                if (i>0)
                    std::cout<<"K "<<time_<<" "
                             <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                             <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                             <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick half step for velocity
                Float dvel1[3], dvel2[3];

                // time step for kick
                gt_inv *= kappa_inv;

                dt = 0.5*ds/gt_inv;

                dvel1[0] = dt * (acc1[0]*kappa_inv + pert1[0]);
                dvel1[1] = dt * (acc1[1]*kappa_inv + pert1[1]);
                dvel1[2] = dt * (acc1[2]*kappa_inv + pert1[2]);

                dvel2[0] = dt * (acc2[0]*kappa_inv + pert2[0]);
                dvel2[1] = dt * (acc2[1]*kappa_inv + pert2[1]);
                dvel2[2] = dt * (acc2[2]*kappa_inv + pert2[2]);

                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];

                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];

#ifdef AR_DEBUG_PRINT_DKD
                std::cout<<"D "<<time_<<" "
                         <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                         <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                         <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick total energy and time transformation factor for drift
                etot_ref_ += 2.0*dt * (mass1* (vel1[0] * pert1[0] + 
                                               vel1[1] * pert1[1] + 
                                               vel1[2] * pert1[2]) +
                                       mass2* (vel2[0] * pert2[0] + 
                                               vel2[1] * pert2[1] + 
                                               vel2[2] * pert2[2]));


                // kick half step for velocity
                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];
                
                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];
                                                                                                                
                // calculate kinetic energy
                ekin_ = 0.5 * (mass1 * (vel1[0]*vel1[0]+vel1[1]*vel1[1]+vel1[2]*vel1[2]) +
                               mass2 * (vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]));

            }

            // make consistent slowdown inner energy 
            etot_sd_ref_ = etot_ref_*kappa_inv;
            ekin_sd_ = ekin_*kappa_inv;
            epot_sd_ = epot_*kappa_inv;
        }

        //! get Hamiltonian
        Float getH() const {
            return manager->interaction.calcH(ekin_ - etot_ref_, epot_);
        }

        //! get Hamiltonian from backup data
        Float getHFromBackup(Float* _bk) const {
            Float& etot_ref =_bk[1];
            Float& ekin = _bk[2];
            Float& epot = _bk[3];
            return manager->interaction.calcH(ekin - etot_ref, epot);
        }

        //! get slowdown Hamiltonian
        Float getHSlowDown() const {
            return manager->interaction.calcH(ekin_sd_ - etot_sd_ref_, epot_sd_);
        }

        //! get slowdown Hamiltonian from backup data
        Float getHSlowDownFromBackup(Float* _bk) const {
            Float& etot_sd_ref =_bk[6];
            Float& ekin_sd = _bk[7];
            Float& epot_sd = _bk[8];
            return manager->interaction.calcH(ekin_sd - etot_sd_ref, epot_sd);            
        }

    };

    //! AR dynamics: TTL with Array slow-down
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo>
    class ARDynamics<TTL, ArraySlowDown, Tparticle, Tpcm, Tpert, Tmethod, Tinfo>: public ARArraySlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, TTL> {
    protected:
        typedef ARArraySlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, TTL> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;
        using Parent::ekin_sd_;
        using Parent::epot_sd_;
        using Parent::etot_sd_ref_;
        using Parent::de_sd_change_cum_;
        using Parent::dH_sd_change_cum_;
        using Parent::de_sd_change_interrupt_;
        using Parent::dH_sd_change_interrupt_;
        using Parent::clearSlowDownEnergy;
        using Parent::copySlowDownEnergy;
        using Parent::calcSlowDownPertInnerBinaryIter;
        using Parent::calcSlowDownInnerBinary;
        using Parent::correctPosSlowDownInner;
        using Parent::updateCenterOfMassForBinaryWithSlowDownInner;
        using Parent::calcEkinSlowDownInner;
        using Parent::findSlowDownInnerBinaryIter;
        using Parent::findSlowDownInner;
        using Parent::calcEKin;
        using Parent::kickVel;
        using Parent::driftTimeAndPos;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;
        using Parent::resetDESlowDownChangeCum;
        using Parent::getDESlowDownChangeCum;
        using Parent::getDHSlowDownChangeCum;
        using Parent::resetDESlowDownChangeBinaryInterrupt;
        using Parent::getDESlowDownChangeBinaryInterrupt;
        using Parent::getDHSlowDownChangeBinaryInterrupt;
        using Parent::getEkinSlowDown;
        using Parent::getEpotSlowDown;
        using Parent::getEtotSlowDownRef;
        using Parent::getEtotSlowDown;
        using Parent::getEnergyErrorSlowDown;
        using Parent::getEnergyErrorSlowDownFromBackup;
        using Parent::getEtotSlowDownRefFromBackup;
        using Parent::getEtotSlowDownFromBackup;
        using Parent::binary_slowdown;
        using Parent::writeBackSlowDownParticles;

    protected:
        // transformation factors
        Float gt_drift_inv_ = 0.0;  ///< integrated inverse time transformation factor for drift: dt(drift) = ds/gt_drift_inv_
        Float gt_kick_inv_ = 0.0;   ///< inverse time transformation factor for kick: dt(kick) = ds/gt_kick_inv_

        //! reserve memory of the data owned by the dynamics
        void reserveDynamicsMem(const int _nmax) {
            binary_slowdown.setMode(COMM::ListMode::local);
            binary_slowdown.reserveMem(_nmax/2+1);
        }

        //! clear the data owned by the dynamics
        void clearDynamics() {
            this->clearSlowDownEnergy();
            gt_drift_inv_ = 0.0;
            gt_kick_inv_ = 0.0;
            binary_slowdown.clear();
        }

        //! copy the data owned by the dynamics
        void copyDynamics(const ARDynamics& _sym) {
            this->copySlowDownEnergy(_sym);
            gt_drift_inv_ = _sym.gt_drift_inv_;
            gt_kick_inv_ = _sym.gt_kick_inv_;
            binary_slowdown = _sym.binary_slowdown;
        }

        //! correct force, potential energy and gt_kick_inv based on slowdown for inner binaries
        /*! 
          @param[in,out] _gt_kick_inv: the inverse time transformation factor for kick step (input), be corrected with slowdown (output)
         */
        inline void correctAccPotGTKickInvSlowDownInner(Float& _gt_kick_inv) {
            int n = binary_slowdown.getSize();
            Float gt_kick_inv_cor = 0.0;
            Float de = 0.0;
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                int i1 = sdi->getMemberIndex(0);
                int i2 = sdi->getMemberIndex(1);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv = 1.0/kappa;
                if (i1>=0) {
                    ASSERT(i2>=0);
                    ASSERT(i1!=i2);
                    ASSERT(i1<particles.getSize());
                    ASSERT(i2<particles.getSize());
                    
                    // calculate pair interaction
                    Force fi[2];
                    Float epoti;
                    Float gt_kick_inv_i = manager->interaction.calcInnerAccPotAndGTKickInvTwo(fi[0], fi[1], epoti, particles[i1], particles[i2]);
                                    // scale binary pair force with slowdown 
                    force_[i1].acc_in[0] += fi[0].acc_in[0]*kappa_inv - fi[0].acc_in[0];
                    force_[i1].acc_in[1] += fi[0].acc_in[1]*kappa_inv - fi[0].acc_in[1];
                    force_[i1].acc_in[2] += fi[0].acc_in[2]*kappa_inv - fi[0].acc_in[2];
                    force_[i2].acc_in[0] += fi[1].acc_in[0]*kappa_inv - fi[1].acc_in[0];
                    force_[i2].acc_in[1] += fi[1].acc_in[1]*kappa_inv - fi[1].acc_in[1];
                    force_[i2].acc_in[2] += fi[1].acc_in[2]*kappa_inv - fi[1].acc_in[2];

                    de += epoti*kappa_inv - epoti;
                    // scale gtgrad with slowdown
                    force_[i1].gtgrad[0] += fi[0].gtgrad[0]*kappa_inv - fi[0].gtgrad[0];
                    force_[i1].gtgrad[1] += fi[0].gtgrad[1]*kappa_inv - fi[0].gtgrad[1];
                    force_[i1].gtgrad[2] += fi[0].gtgrad[2]*kappa_inv - fi[0].gtgrad[2];
                    force_[i2].gtgrad[0] += fi[1].gtgrad[0]*kappa_inv - fi[1].gtgrad[0];
                    force_[i2].gtgrad[1] += fi[1].gtgrad[1]*kappa_inv - fi[1].gtgrad[1];
                    force_[i2].gtgrad[2] += fi[1].gtgrad[2]*kappa_inv - fi[1].gtgrad[2];

                    // gt kick
                    gt_kick_inv_cor += gt_kick_inv_i*(kappa_inv - 1.0);
                }
            }

            // global slowdown
            const Float kappa_inv_global = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            for (int i=0; i<force_.getSize(); i++) {
                force_[i].acc_in[0] *= kappa_inv_global;
                force_[i].acc_in[1] *= kappa_inv_global;
                force_[i].acc_in[2] *= kappa_inv_global;
                force_[i].gtgrad[0] *= kappa_inv_global;
                force_[i].gtgrad[1] *= kappa_inv_global;
                force_[i].gtgrad[2] *= kappa_inv_global;
            }

            epot_sd_ = (epot_ + de)*kappa_inv_global;
            _gt_kick_inv = (_gt_kick_inv + gt_kick_inv_cor)*kappa_inv_global;
        }

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        inline Float calcAccPotAndGTKickInv() {
            Float gt_kick_inv = manager->interaction.calcAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());            

//#ifdef AR_DEBUG
//            // check c.m. force 
//            Force fcm;
//            Float mcm=0.0;
//            for (int i=0; i<particles.getSize(); i++) {
//                for (int k=0; k<3; k++) {
//                    fcm.acc_in[k] += particles[i].mass * force_[i].acc_in[k];
//                }
//                mcm += particles[i].mass;
//            }
//            for (int k=0; k<3; k++) {
//                fcm.acc_in[k] /= mcm;
//                ASSERT(abs(fcm.acc_in[k])<1e-10);
//            }
//#endif
            // slowdown binary acceleration
            correctAccPotGTKickInvSlowDownInner(gt_kick_inv);
//#ifdef AR_DEBUG
//            // check c.m. force 
//            fcm.acc_in[0] = fcm.acc_in[1] = fcm.acc_in[2] = 0.0;
//            for (int i=0; i<particles.getSize(); i++) {
//                for (int k=0; k<3; k++) {
//                    fcm.acc_in[k] += particles[i].mass * force_[i].acc_in[k];
//                }
//            }
//            for (int k=0; k<3; k++) {
//                fcm.acc_in[k] /= mcm;
//                ASSERT(abs(fcm.acc_in[k])<1e-10);
//            }
//#endif
            return gt_kick_inv;
        }

        //! kick energy and time transformation function for drift
        /*!
          @param[in] _dt: time step
        */

        inline void kickEtotAndGTDrift(const Float _dt) {
            Float de = Float(0.0);
            Float dg = Float(0.0);
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                Float* gtgrad=force[i].gtgrad;
                de += mass * (vel[0] * pert[0] +
                              vel[1] * pert[1] +
                              vel[2] * pert[2]);
                dg +=  (vel[0] * gtgrad[0] +
                        vel[1] * gtgrad[1] +
                        vel[2] * gtgrad[2]);
            }
            etot_ref_ += _dt * de;

            etot_sd_ref_ += _dt * de;

            // correct gt_drift_inv
            const Float kappa_inv_global = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            int n = binary_slowdown.getSize();
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv = 1.0/kappa;
                Float* velcm = sdi->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi->getMemberIndex(k);
                    if (j>=0) {
                        ASSERT(j<particles.getSize());
                        Float* gtgrad=force_[j].gtgrad;
                        Float* vel = particles[j].getVel();
                        Float vrel[3] = { vel[0] - velcm[0], 
                                          vel[1] - velcm[1], 
                                          vel[2] - velcm[2]}; 
                        dg +=  (vrel[0] * (kappa_inv-1)* gtgrad[0] +
                                vrel[1] * (kappa_inv-1)* gtgrad[1] +
                                vrel[2] * (kappa_inv-1)* gtgrad[2]);
                    }
                }
            }
            gt_drift_inv_ += _dt * dg *kappa_inv_global;
        }

    public:
        //! update slowdown factor based on perturbation and record slowdown energy change
        /*! Update slowdown inner and global, update gt_inv
            @param [in] _update_energy_flag: Record cumulative slowdown energy change if true;
            @param [in] _stable_check_flag: check whether the binary tree is stable if true;
         */
        void updateSlowDownAndCorrectEnergy(const bool _update_energy_flag, const bool _stable_check_flag) {
            auto& bin_root = *binary_slowdown[0];
            auto& sd_root = bin_root.slowdown;
            Float sd_backup = sd_root.getSlowDownFactor();

            // when the maximum inner slowdown is large, the outer should not be slowed down since the system may not be stable.
            //if (inner_sd_change_flag&&sd_org_inner_max<1000.0*manager->slowdown_pert_ratio_ref) sd_root.setSlowDownFactor(1.0);
            //if (time_>=sd_root.getUpdateTime()) {
            sd_root.pert_in = manager->interaction.calcPertFromBinary(bin_root);
            sd_root.pert_out = 0.0;
            Float t_min_sq= NUMERIC_FLOAT_MAX;
            manager->interaction.calcSlowDownPert(sd_root.pert_out, t_min_sq, getTime(), particles.cm, perturber);
            sd_root.timescale = std::min(sd_root.getTimescaleMax(), sqrt(t_min_sq));

            //Float period_amplify_max = NUMERIC_FLOAT_MAX;
            if (_stable_check_flag) {
                // check whether the system is stable for 10000 out period
                Float stab = bin_root.stableCheckIter(bin_root,10000*bin_root.period);
                if (stab<1.0) {
                    sd_root.period = bin_root.period;
                    sd_root.calcSlowDownFactor();
                }
                else sd_root.setSlowDownFactor(1.0);
                // stablility criterion
                // The slowdown factor should not make the system unstable, thus the Qst/Q set the limitation of the increasing of inner semi-major axis.
                //if (stab>0) {
                //    Float semi_amplify_max =  std::max(Float(1.0),1.0/stab);
                //    period_amplify_max = pow(semi_amplify_max,3.0/2.0);
                //}
            }
            else if (bin_root.semi>0) {
                sd_root.period = bin_root.period;
                sd_root.calcSlowDownFactor();
            }
            else sd_root.setSlowDownFactor(1.0);

            //sd_root.increaseUpdateTimeOnePeriod();
            //}


            // inner binary slowdown
            int n = binary_slowdown.getSize();
            bool modified_flag=false;
            for (int i=1; i<n; i++) {
                auto* sdi = binary_slowdown[i];
                //if (time_>=sdi->slowdown.getUpdateTime()) {
                sdi->calcCenterOfMass();
                calcSlowDownInnerBinary(*sdi);
                //sdi->slowdown.increaseUpdateTimeOnePeriod();
                modified_flag=true;
                //}
            }    

            if (_update_energy_flag) {
                Float ekin_sd_bk = ekin_sd_;
                Float epot_sd_bk = epot_sd_;
                Float H_sd_bk = getHSlowDown();

                if (modified_flag) {
                    // initialize the gt_drift_inv_ with new slowdown factor
                    Float gt_kick_inv_sdi = manager->interaction.calcAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());
                    correctAccPotGTKickInvSlowDownInner(gt_kick_inv_sdi);
                    gt_drift_inv_ += gt_kick_inv_sdi - gt_kick_inv_;
                    gt_kick_inv_ = gt_kick_inv_sdi;
                    // correct etot_sd_ref_ with new slowdown
                    calcEkinSlowDownInner(ekin_);
                }
                else {
                    Float kappa_inv = 1.0/sd_root.getSlowDownFactor();
                    Float gt_kick_inv_sdi = gt_kick_inv_*sd_backup*kappa_inv;
                    gt_drift_inv_ += gt_kick_inv_sdi - gt_kick_inv_;
                    gt_kick_inv_ = gt_kick_inv_sdi;
                    // only need to correct the total value
                    ekin_sd_ = ekin_*kappa_inv;
                    epot_sd_ = epot_*kappa_inv;
                }

                Float de_sd = (ekin_sd_ - ekin_sd_bk) + (epot_sd_ - epot_sd_bk);
                etot_sd_ref_ += de_sd;

                Float dH_sd = getHSlowDown() - H_sd_bk;

                // add slowdown change to the global slowdown energy
                de_sd_change_cum_ += de_sd;
                dH_sd_change_cum_ += dH_sd;
            }
        }

        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time: real physical time to initialize
        */
        void initialIntegration(const Float _time) {
            ASSERT(checkParams());

            // particle number and data address
            const int n_particle = particles.getSize();

            // resize force array
            force_.resizeNoInitialize(n_particle);

            // Initial intgrt value t (avoid confusion of real time when slowdown is used)
            time_ = _time;

            // check particle number
            ASSERT(particles.getSize()>=2);

            // reset particle modification flag
            particles.setModifiedFalse();

            // check the center-of-mass initialization
            if(particles.isOriginFrame()) {
                particles.calcCenterOfMass();
                particles.shiftToCenterOfMassFrame();
                for (int i=0; i<info.binarytree.getSize(); i++) {
                    auto& bin = info.binarytree[i];
                    bin.pos[0] -= particles.cm.pos[0];
                    bin.pos[1] -= particles.cm.pos[1];
                    bin.pos[2] -= particles.cm.pos[2];
                    bin.vel[0] -= particles.cm.vel[0];
                    bin.vel[1] -= particles.cm.vel[1];
                    bin.vel[2] -= particles.cm.vel[2];
                }
            }
            ASSERT(info.getBinaryTreeRoot().pos[0]*info.getBinaryTreeRoot().pos[0]<1e-10);
            ASSERT(info.getBinaryTreeRoot().vel[0]*info.getBinaryTreeRoot().vel[0]<1e-10);

            binary_slowdown.increaseSizeNoInitialize(1);
            binary_slowdown[0] = &info.getBinaryTreeRoot();

            // set slowdown reference
            SlowDown& slowdown_root = info.getBinaryTreeRoot().slowdown;

            // slowdown for the system
            slowdown_root.initialSlowDownReference(manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);

            if (particles.getSize()>2) {
                findSlowDownInner(time_);
                // update c.m. of binaries 
                //updateCenterOfMassForBinaryWithSlowDownInner();
            }

            updateSlowDownAndCorrectEnergy(false,true);

            gt_kick_inv_ = calcAccPotAndGTKickInv();

            // initially gt_drift 
            gt_drift_inv_ = gt_kick_inv_;


            calcEKin();

            etot_ref_ = ekin_ + epot_;
            etot_sd_ref_ = ekin_sd_ + epot_sd_;

            Float de_sd = etot_sd_ref_ - etot_ref_;

            // add slowdown change to the global slowdown energy
            de_sd_change_cum_ += de_sd;
            dH_sd_change_cum_ = 0.0;

        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // inverse time transformation factor for drift
                Float gt_drift_inv = gt_drift_inv_;

                // drift
                Float dt_drift = ds_drift/gt_drift_inv;

                // drift time and postion
                driftTimeAndPos(dt_drift);
                _time_table[i] = time_;

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                //! calc force, potential and inverse time transformation factor for kick
                Float gt_kick_inv = calcAccPotAndGTKickInv();

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // back up gt_kick 
                gt_kick_inv_ = gt_kick_inv;
                // kick total energy and inverse time transformation factor for drift
                kickEtotAndGTDrift(dt_kick);
                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // calculate kinetic energy
                calcEKin();
            }
        }

        //! integration for two body one step
        /*! For two-body problem the calculation can be much symplified to improve performance. 
          Besides, the slow-down factor calculation is embedded in the Drift (for time) and Kick (for perturbation). 
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().         
        */
        void integrateTwoOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group number
            const int nloop = manager->step.getCDPairSize();
            
            const int n_particle = particles.getSize();
            ASSERT(n_particle==2);

            const Float kappa_inv = 1.0/info.getBinaryTreeRoot().slowdown.getSlowDownFactor();

            Tparticle* particle_data = particles.getDataAddress();
            Float mass1 = particle_data[0].mass;
            Float* pos1 = particle_data[0].getPos();
            Float* vel1 = particle_data[0].getVel();

            Float mass2 = particle_data[1].mass;
            Float* pos2 = particle_data[1].getPos();
            Float* vel2 = particle_data[1].getVel();

            Force* force_data = force_.getDataAddress();
            Float* acc1 = force_data[0].acc_in;
            Float* pert1= force_data[0].acc_pert;

            Float* acc2 = force_data[1].acc_in;
            Float* pert2= force_data[1].acc_pert;
            Float* gtgrad1 = force_data[0].gtgrad;
            Float* gtgrad2 = force_data[1].gtgrad;

#ifdef AR_DEBUG_PRINT_DKD
            std::cout<<"K "<<time_<<" "
                     <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                     <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                     <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds = manager->step.getCK(i)*_ds;
                // inverse time transformation factor for drift
                Float gt_inv = gt_drift_inv_;
                // drift
                Float dt = ds/gt_inv;
                ASSERT(!ISNAN(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                _time_table[i] = time_;

                Float dt_sd = dt*kappa_inv;

                // drift position
                pos1[0] += dt_sd * vel1[0];
                pos1[1] += dt_sd * vel1[1];
                pos1[2] += dt_sd * vel1[2];

                pos2[0] += dt_sd * vel2[0];
                pos2[1] += dt_sd * vel2[1];
                pos2[2] += dt_sd * vel2[2];

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                gt_inv = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm, perturber, _time_table[i]);

                ASSERT(!ISNAN(epot_));

#ifdef AR_DEBUG_PRINT_DKD
                if (i>0)
                    std::cout<<"K "<<time_<<" "
                             <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                             <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                             <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick half step for velocity
                Float dvel1[3], dvel2[3];

                // time step for kick
                gt_inv *= kappa_inv;

                dt = 0.5*ds/gt_inv;

                dvel1[0] = dt * (acc1[0]*kappa_inv + pert1[0]);
                dvel1[1] = dt * (acc1[1]*kappa_inv + pert1[1]);
                dvel1[2] = dt * (acc1[2]*kappa_inv + pert1[2]);

                dvel2[0] = dt * (acc2[0]*kappa_inv + pert2[0]);
                dvel2[1] = dt * (acc2[1]*kappa_inv + pert2[1]);
                dvel2[2] = dt * (acc2[2]*kappa_inv + pert2[2]);

                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];

                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];

#ifdef AR_DEBUG_PRINT_DKD
                std::cout<<"D "<<time_<<" "
                         <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                         <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                         <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick total energy and time transformation factor for drift
                etot_ref_ += 2.0*dt * (mass1* (vel1[0] * pert1[0] + 
                                               vel1[1] * pert1[1] + 
                                               vel1[2] * pert1[2]) +
                                       mass2* (vel2[0] * pert2[0] + 
                                               vel2[1] * pert2[1] + 
                                               vel2[2] * pert2[2]));

                // back up gt_kick_inv
                gt_kick_inv_ = gt_inv;

                // integrate gt_drift_inv
                gt_drift_inv_ +=  2.0*dt*kappa_inv*kappa_inv* (vel1[0] * gtgrad1[0] +
                                                               vel1[1] * gtgrad1[1] +
                                                               vel1[2] * gtgrad1[2] +
                                                               vel2[0] * gtgrad2[0] +
                                                               vel2[1] * gtgrad2[1] +
                                                               vel2[2] * gtgrad2[2]);


                // kick half step for velocity
                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];
                
                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];
                                                                                                                
                // calculate kinetic energy
                ekin_ = 0.5 * (mass1 * (vel1[0]*vel1[0]+vel1[1]*vel1[1]+vel1[2]*vel1[2]) +
                               mass2 * (vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]));

            }

            // make consistent slowdown inner energy 
            etot_sd_ref_ = etot_ref_*kappa_inv;
            ekin_sd_ = ekin_*kappa_inv;
            epot_sd_ = epot_*kappa_inv;
        }

        //! get Hamiltonian
        Float getH() const {
            //return (ekin_ - etot_ref_)/gt_drift_inv_ + epot_/gt_kick_inv_;
            return (ekin_ + epot_ - etot_ref_)/gt_kick_inv_;
        }

        //! get Hamiltonian from backup data
        Float getHFromBackup(Float* _bk) const {
            Float& etot_ref =_bk[1];
            Float& ekin = _bk[2];
            Float& epot = _bk[3];
            //Float& gt_drift_inv = _bk[13];
            Float& gt_kick_inv  = _bk[14];
            return (ekin + epot - etot_ref)/gt_kick_inv;
            //return (ekin - etot_ref)/gt_drift_inv + epot/gt_kick_inv;
        }

        //! get slowdown Hamiltonian
        Float getHSlowDown() const {
            //return (ekin_sd_ - etot_sd_ref_)/gt_drift_inv_ + epot_sd_/gt_kick_inv_;
            return (ekin_sd_ + epot_sd_ - etot_sd_ref_)/gt_kick_inv_;
        }

        //! get slowdown Hamiltonian from backup data
        Float getHSlowDownFromBackup(Float* _bk) const {
            Float& etot_sd_ref =_bk[6];
            Float& ekin_sd = _bk[7];
            Float& epot_sd = _bk[8];
            //Float& gt_drift_inv = _bk[13];
            Float& gt_kick_inv  = _bk[14];
            //return (ekin_sd - etot_sd_ref)/gt_drift_inv + epot_sd/gt_kick_inv;
            return (ekin_sd + epot_sd - etot_sd_ref)/gt_kick_inv;
        }

        //! Get integrated inverse time transformation factor
        /*! In TTF case, it is calculated by integrating \f$ \frac{dg}{dt} = \sum_k \frac{\partial g}{\partial \vec{r_k}} \bullet \vec{v_k} \f$.
          Notice last step is the sub-step in one symplectic loop
          \return inverse time transformation factor for drift
        */
        Float getGTDriftInv() const {
            return gt_drift_inv_;
        }

    };

    //! AR dynamics: LogH with Tree slow-down
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo>
    class ARDynamics<LogH, TreeSlowDown, Tparticle, Tpcm, Tpert, Tmethod, Tinfo>: public ARTreeSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, LogH> {
    protected:
        typedef ARTreeSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, LogH> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;
        using Parent::ekin_sd_;
        using Parent::epot_sd_;
        using Parent::etot_sd_ref_;
        using Parent::de_sd_change_cum_;
        using Parent::dH_sd_change_cum_;
        using Parent::de_sd_change_interrupt_;
        using Parent::dH_sd_change_interrupt_;
        using Parent::clearSlowDownEnergy;
        using Parent::copySlowDownEnergy;
        using Parent::calcSlowDownPertInnerBinaryIter;
        using Parent::calcSlowDownInnerBinary;
        using Parent::calcTwoEKinIter;
        using Parent::calcEKin;
        using Parent::kickVel;
        using Parent::driftPosTreeIter;
        using Parent::driftTimeAndPos;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;
        using Parent::resetDESlowDownChangeCum;
        using Parent::getDESlowDownChangeCum;
        using Parent::getDHSlowDownChangeCum;
        using Parent::resetDESlowDownChangeBinaryInterrupt;
        using Parent::getDESlowDownChangeBinaryInterrupt;
        using Parent::getDHSlowDownChangeBinaryInterrupt;
        using Parent::getEkinSlowDown;
        using Parent::getEpotSlowDown;
        using Parent::getEtotSlowDownRef;
        using Parent::getEtotSlowDown;
        using Parent::getEnergyErrorSlowDown;
        using Parent::getEnergyErrorSlowDownFromBackup;
        using Parent::getEtotSlowDownRefFromBackup;
        using Parent::getEtotSlowDownFromBackup;
        using Parent::writeBackSlowDownParticlesIter;
        using Parent::writeBackSlowDownParticles;

    protected:
        //! reserve memory of the data owned by the dynamics
        void reserveDynamicsMem(const int _nmax) {
            (void)_nmax;
        }

        //! clear the data owned by the dynamics
        void clearDynamics() {
            this->clearSlowDownEnergy();
        }

        //! copy the data owned by the dynamics
        void copyDynamics(const ARDynamics& _sym) {
            this->copySlowDownEnergy(_sym);
        }

        //! calc force, potential and inverse time transformation factor for one pair of particles
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _i: particle i index
          @param[in] _j: particle j index
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvTwo(const Float& _inv_nest_sd, const int _i, const int _j) {
            ASSERT(_i>=0&&_i<particles.getSize());
            ASSERT(_j>=0&&_j<particles.getSize());
            
            // calculate pair interaction
            Force fij[2];
            Float epotij;
            Float gt_kick_inv = manager->interaction.calcInnerAccPotAndGTKickInvTwo(fij[0], fij[1], epotij, particles[_i], particles[_j]);

            // scale binary pair force with slowdown 
            force_[_i].acc_in[0] += fij[0].acc_in[0]*_inv_nest_sd;
            force_[_i].acc_in[1] += fij[0].acc_in[1]*_inv_nest_sd;
            force_[_i].acc_in[2] += fij[0].acc_in[2]*_inv_nest_sd;
            force_[_j].acc_in[0] += fij[1].acc_in[0]*_inv_nest_sd;
            force_[_j].acc_in[1] += fij[1].acc_in[1]*_inv_nest_sd;
            force_[_j].acc_in[2] += fij[1].acc_in[2]*_inv_nest_sd;

            epot_    += epotij;
            epot_sd_ += epotij*_inv_nest_sd;

            
            return gt_kick_inv * _inv_nest_sd;
        }

        //! calc force, potential and inverse time transformation factor for one particle by walking binary tree
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _i: particle index
          @param[in] _bin: binary tree for walking
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvOneTreeIter(const Float& _inv_nest_sd, const int _i, AR::BinaryTree<Tparticle>& _bin) {
            Float gt_kick_inv=0.0;
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) // particle - tree
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(_inv_nest_sd, _i, *(_bin.getMemberAsTree(k)));
                else  // particle - particle
                    gt_kick_inv += calcAccPotAndGTKickInvTwo(_inv_nest_sd, _i, _bin.getMemberIndex(k));
            }
            return gt_kick_inv;
        }

        //! calc crossing force, potential and inverse time transformation factor between binary tree i and binary tree j
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _bini: binary tree i for walking
          @param[in] _binj: binary tree j for walking
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvCrossTreeIter(const Float& _inv_nest_sd, AR::BinaryTree<Tparticle>& _bini, AR::BinaryTree<Tparticle>& _binj) {
            ASSERT(&_bini!=&_binj);
            Float gt_kick_inv=0.0;
            for (int k=0; k<2; k++) { 
                if (_bini.isMemberTree(k))  // tree - tree
                    gt_kick_inv += calcAccPotAndGTKickInvCrossTreeIter(_inv_nest_sd, *(_bini.getMemberAsTree(k)), _binj);
                else  // particle - tree
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(_inv_nest_sd, _bini.getMemberIndex(k), _binj);
            }
            return gt_kick_inv;
        }

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary to drift pos
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvTreeIter(const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) {
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
            Float gt_kick_inv = 0.0;

            // check left 
            if (_bin.isMemberTree(0)) { // left is tree
                auto* bin_left = _bin.getMemberAsTree(0);
                // inner interaction of left tree
                gt_kick_inv += calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_left);

                if (_bin.isMemberTree(1)) { // right is tree
                    auto* bin_right = _bin.getMemberAsTree(1);

                    // inner interaction of right tree
                    gt_kick_inv += calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_right);
                    
                    // cross interaction
                    gt_kick_inv += calcAccPotAndGTKickInvCrossTreeIter(inv_nest_sd, *bin_left, *bin_right);
                }
                else { // right is particle
                    // cross interaction from particle j to tree left
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(inv_nest_sd,  _bin.getMemberIndex(1), *bin_left);
                }
            }
            else { // left is particle
                if (_bin.isMemberTree(1)) { // right is tree
                    auto* bin_right = _bin.getMemberAsTree(1);
                    // inner interaction of right tree
                    gt_kick_inv += calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_right);

                    // cross interaction from particle i to tree right
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(inv_nest_sd, _bin.getMemberIndex(0), *bin_right);
                }
                else { // right is particle
                    // particle - particle interaction
                    gt_kick_inv += calcAccPotAndGTKickInvTwo(inv_nest_sd, _bin.getMemberIndex(0), _bin.getMemberIndex(1));
                }
            }


            return gt_kick_inv;
        }

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        inline Float calcAccPotAndGTKickInv() {
            epot_ = 0.0;
            epot_sd_ = 0.0;
            for (int i=0; i<force_.getSize(); i++) force_[i].clear();

            Float gt_kick_inv = calcAccPotAndGTKickInvTreeIter(1.0, info.getBinaryTreeRoot());
            // pertuber force
            manager->interaction.calcAccPert(force_.getDataAddress(), particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());
            return gt_kick_inv;
        }

        //! kick energy
        /*!
          @param[in] _dt: time step
        */
        inline void kickEtot(const Float _dt) {
            Float de = 0.0;
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                de += mass * (vel[0] * pert[0] +
                              vel[1] * pert[1] +
                              vel[2] * pert[2]);
            }
            etot_sd_ref_ += _dt * de;
            etot_ref_ += _dt * de;
        }

    public:
        //! update slowdown factor based on perturbation and record slowdown energy change
        /*! Update slowdown inner and global.
            @param [in] _update_energy_flag: Record cumulative slowdown energy change if true;
            @param [in] _stable_check_flag: check whether the binary tree is stable if true;
         */
        void updateSlowDownAndCorrectEnergy(const bool _update_energy_flag, const bool _stable_check_flag) {
            auto& bin_root = info.getBinaryTreeRoot();
            auto& sd_root = bin_root.slowdown;


            // when the maximum inner slowdown is large, the outer should not be slowed down since the system may not be stable.
            //if (inner_sd_change_flag&&sd_org_inner_max<1000.0*manager->slowdown_pert_ratio_ref) sd_root.setSlowDownFactor(1.0);
            //if (time_>=sd_root.getUpdateTime()) {
            sd_root.pert_in = manager->interaction.calcPertFromBinary(bin_root);
            sd_root.pert_out = 0.0;
            Float t_min_sq= NUMERIC_FLOAT_MAX;
            manager->interaction.calcSlowDownPert(sd_root.pert_out, t_min_sq, getTime(), particles.cm, perturber);
            sd_root.timescale = std::min(sd_root.getTimescaleMax(), sqrt(t_min_sq));

            //Float period_amplify_max = NUMERIC_FLOAT_MAX;
            if (_stable_check_flag) {
                // check whether the system is stable for 10000 out period and the apo-center is below break criterion
                Float stab = bin_root.stableCheckIter(bin_root,10000*bin_root.period);
                Float apo = bin_root.semi*(1+bin_root.ecc);
                if (stab<1.0 && apo<info.r_break_crit) {
                    sd_root.period = bin_root.period;
                    sd_root.calcSlowDownFactor();
                }
                else sd_root.setSlowDownFactor(1.0);

                // stablility criterion
                // The slowdown factor should not make the system unstable, thus the Qst/Q set the limitation of the increasing of inner semi-major axis.
                //if (stab>0 && stab != NUMERIC_FLOAT_MAX) {
                //    Float semi_amplify_max =  std::max(Float(1.0),1.0/stab);
                //    period_amplify_max = pow(semi_amplify_max,3.0/2.0);
                //}
            }
            else if (bin_root.semi>0) {
                sd_root.period = bin_root.period;
                sd_root.calcSlowDownFactor();
            }
            else sd_root.setSlowDownFactor(1.0);

            //sd_root.increaseUpdateTimeOnePeriod();
            //}

            // inner binary slowdown
            Float sd_org_inner_max = 0.0;
            bool inner_sd_change_flag=false;
            int n_bin = info.binarytree.getSize();
            for (int i=0; i<n_bin-1; i++) {
                auto& bini = info.binarytree[i];
                //if (time_>=bini.slowdown.getUpdateTime()) {
                bini.calcCenterOfMass();
                calcSlowDownInnerBinary(bini);

                //sdi->slowdown.increaseUpdateTimeOnePeriod();
                sd_org_inner_max = std::max(bini.slowdown.getSlowDownFactorOrigin(),sd_org_inner_max);
                inner_sd_change_flag=true;
                //}
            }


            if (_update_energy_flag) {
                Float ekin_sd_bk = ekin_sd_;
                Float epot_sd_bk = epot_sd_;
                Float H_sd_bk = getHSlowDown();
                if(inner_sd_change_flag) {
                    calcAccPotAndGTKickInv();
                    calcEKin();
                }
                else {
                    Float kappa_inv = 1.0/sd_root.getSlowDownFactor();
                    ekin_sd_ = ekin_*kappa_inv;
                    epot_sd_ = epot_*kappa_inv;
                }
                Float de_sd = (ekin_sd_ - ekin_sd_bk) + (epot_sd_ - epot_sd_bk);
                etot_sd_ref_ += de_sd;

                Float dH_sd = getHSlowDown() - H_sd_bk;

                // add slowdown change to the global slowdown energy
                de_sd_change_cum_ += de_sd;
                dH_sd_change_cum_ += dH_sd;
            }
        }

        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time: real physical time to initialize
        */
        void initialIntegration(const Float _time) {
            ASSERT(checkParams());

            // particle number and data address
            const int n_particle = particles.getSize();

            // resize force array
            force_.resizeNoInitialize(n_particle);

            // Initial intgrt value t (avoid confusion of real time when slowdown is used)
            time_ = _time;

            // check particle number
            ASSERT(particles.getSize()>=2);

            // reset particle modification flag
            particles.setModifiedFalse();

            // check the center-of-mass initialization
            if(particles.isOriginFrame()) {
                particles.calcCenterOfMass();
                particles.shiftToCenterOfMassFrame();
                for (int i=0; i<info.binarytree.getSize(); i++) {
                    auto& bin = info.binarytree[i];
                    bin.pos[0] -= particles.cm.pos[0];
                    bin.pos[1] -= particles.cm.pos[1];
                    bin.pos[2] -= particles.cm.pos[2];
                    bin.vel[0] -= particles.cm.vel[0];
                    bin.vel[1] -= particles.cm.vel[1];
                    bin.vel[2] -= particles.cm.vel[2];
                }
            }
            ASSERT(info.getBinaryTreeRoot().pos[0]*info.getBinaryTreeRoot().pos[0]<1e-10);
            ASSERT(info.getBinaryTreeRoot().vel[0]*info.getBinaryTreeRoot().vel[0]<1e-10);

            for (int i=0; i<info.binarytree.getSize(); i++) 
                info.binarytree[i].slowdown.initialSlowDownReference(manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);

            updateSlowDownAndCorrectEnergy(false,true);

            calcAccPotAndGTKickInv();

            calcEKin();

            etot_ref_ = ekin_ + epot_;
            etot_sd_ref_ = ekin_sd_ + epot_sd_;

            Float de_sd = etot_sd_ref_ - etot_ref_;

            // add slowdown change to the global slowdown energy
            de_sd_change_cum_ += de_sd;
            dH_sd_change_cum_ = 0.0;

        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // inverse time transformation factor for drift
                Float gt_drift_inv = manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_); // pt = -etot

                // drift
                Float dt_drift = ds_drift/gt_drift_inv;

                // drift time and postion
                driftTimeAndPos(dt_drift);
                _time_table[i] = time_;

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                //! calc force, potential and inverse time transformation factor for kick
                Float gt_kick_inv = calcAccPotAndGTKickInv();

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // kick total energy 
                kickEtot(dt_kick);
                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // calculate kinetic energy
                calcEKin();
            }
        }

        //! integration for two body one step
        /*! For two-body problem the calculation can be much symplified to improve performance. 
          Besides, the slow-down factor calculation is embedded in the Drift (for time) and Kick (for perturbation). 
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().         
        */
        void integrateTwoOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group number
            const int nloop = manager->step.getCDPairSize();
            
            const int n_particle = particles.getSize();
            ASSERT(n_particle==2);

            const Float kappa_inv = 1.0/info.getBinaryTreeRoot().slowdown.getSlowDownFactor();

            Tparticle* particle_data = particles.getDataAddress();
            Float mass1 = particle_data[0].mass;
            Float* pos1 = particle_data[0].getPos();
            Float* vel1 = particle_data[0].getVel();

            Float mass2 = particle_data[1].mass;
            Float* pos2 = particle_data[1].getPos();
            Float* vel2 = particle_data[1].getVel();

            Force* force_data = force_.getDataAddress();
            Float* acc1 = force_data[0].acc_in;
            Float* pert1= force_data[0].acc_pert;

            Float* acc2 = force_data[1].acc_in;
            Float* pert2= force_data[1].acc_pert;

#ifdef AR_DEBUG_PRINT_DKD
            std::cout<<"K "<<time_<<" "
                     <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                     <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                     <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds = manager->step.getCK(i)*_ds;
                // inverse time transformation factor for drift
                Float gt_inv = manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_); // pt = -etot_sd
                // drift
                Float dt = ds/gt_inv;
                ASSERT(!ISNAN(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                _time_table[i] = time_;

                Float dt_sd = dt*kappa_inv;

                // drift position
                pos1[0] += dt_sd * vel1[0];
                pos1[1] += dt_sd * vel1[1];
                pos1[2] += dt_sd * vel1[2];

                pos2[0] += dt_sd * vel2[0];
                pos2[1] += dt_sd * vel2[1];
                pos2[2] += dt_sd * vel2[2];

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                gt_inv = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm, perturber, _time_table[i]);

                ASSERT(!ISNAN(epot_));

#ifdef AR_DEBUG_PRINT_DKD
                if (i>0)
                    std::cout<<"K "<<time_<<" "
                             <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                             <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                             <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick half step for velocity
                Float dvel1[3], dvel2[3];

                // time step for kick
                gt_inv *= kappa_inv;

                dt = 0.5*ds/gt_inv;

                dvel1[0] = dt * (acc1[0]*kappa_inv + pert1[0]);
                dvel1[1] = dt * (acc1[1]*kappa_inv + pert1[1]);
                dvel1[2] = dt * (acc1[2]*kappa_inv + pert1[2]);

                dvel2[0] = dt * (acc2[0]*kappa_inv + pert2[0]);
                dvel2[1] = dt * (acc2[1]*kappa_inv + pert2[1]);
                dvel2[2] = dt * (acc2[2]*kappa_inv + pert2[2]);

                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];

                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];

#ifdef AR_DEBUG_PRINT_DKD
                std::cout<<"D "<<time_<<" "
                         <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                         <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                         <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick total energy and time transformation factor for drift
                etot_ref_ += 2.0*dt * (mass1* (vel1[0] * pert1[0] + 
                                               vel1[1] * pert1[1] + 
                                               vel1[2] * pert1[2]) +
                                       mass2* (vel2[0] * pert2[0] + 
                                               vel2[1] * pert2[1] + 
                                               vel2[2] * pert2[2]));


                // kick half step for velocity
                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];
                
                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];
                                                                                                                
                // calculate kinetic energy
                ekin_ = 0.5 * (mass1 * (vel1[0]*vel1[0]+vel1[1]*vel1[1]+vel1[2]*vel1[2]) +
                               mass2 * (vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]));

            }

            // make consistent slowdown inner energy 
            etot_sd_ref_ = etot_ref_*kappa_inv;
            ekin_sd_ = ekin_*kappa_inv;
            epot_sd_ = epot_*kappa_inv;
        }

        //! get Hamiltonian
        Float getH() const {
            return manager->interaction.calcH(ekin_ - etot_ref_, epot_);
        }

        //! get Hamiltonian from backup data
        Float getHFromBackup(Float* _bk) const {
            Float& etot_ref =_bk[1];
            Float& ekin = _bk[2];
            Float& epot = _bk[3];
            return manager->interaction.calcH(ekin - etot_ref, epot);
        }

        //! get slowdown Hamiltonian
        Float getHSlowDown() const {
            return manager->interaction.calcH(ekin_sd_ - etot_sd_ref_, epot_sd_);
        }

        //! get slowdown Hamiltonian from backup data
        Float getHSlowDownFromBackup(Float* _bk) const {
            Float& etot_sd_ref =_bk[6];
            Float& ekin_sd = _bk[7];
            Float& epot_sd = _bk[8];
            return manager->interaction.calcH(ekin_sd - etot_sd_ref, epot_sd);            
        }

    };

    //! AR dynamics: TTL with Tree slow-down
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo>
    class ARDynamics<TTL, TreeSlowDown, Tparticle, Tpcm, Tpert, Tmethod, Tinfo>: public ARTreeSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, TTL> {
    protected:
        typedef ARTreeSlowDownHelpers<Tparticle, Tpcm, Tpert, Tmethod, Tinfo, TTL> Parent;
        using Parent::time_;
        using Parent::etot_ref_;
        using Parent::ekin_;
        using Parent::epot_;
        using Parent::de_change_interrupt_;
        using Parent::dH_change_interrupt_;
        using Parent::force_;
        using Parent::updateBinaryVelIter;
        using Parent::updateBinaryCMIter;
        using Parent::setBinaryCMZeroIter;
        using Parent::updateBinarySemiEccPeriodIter;
        using Parent::ekin_sd_;
        using Parent::epot_sd_;
        using Parent::etot_sd_ref_;
        using Parent::de_sd_change_cum_;
        using Parent::dH_sd_change_cum_;
        using Parent::de_sd_change_interrupt_;
        using Parent::dH_sd_change_interrupt_;
        using Parent::clearSlowDownEnergy;
        using Parent::copySlowDownEnergy;
        using Parent::calcSlowDownPertInnerBinaryIter;
        using Parent::calcSlowDownInnerBinary;
        using Parent::calcTwoEKinIter;
        using Parent::calcEKin;
        using Parent::kickVel;
        using Parent::driftPosTreeIter;
        using Parent::driftTimeAndPos;

    public:
        typedef typename Parent::Force Force;
        using Parent::manager;
        using Parent::particles;
        using Parent::perturber;
        using Parent::info;
        using Parent::profile;
        using Parent::checkParams;
        using Parent::getTime;
        using Parent::resetDESlowDownChangeCum;
        using Parent::getDESlowDownChangeCum;
        using Parent::getDHSlowDownChangeCum;
        using Parent::resetDESlowDownChangeBinaryInterrupt;
        using Parent::getDESlowDownChangeBinaryInterrupt;
        using Parent::getDHSlowDownChangeBinaryInterrupt;
        using Parent::getEkinSlowDown;
        using Parent::getEpotSlowDown;
        using Parent::getEtotSlowDownRef;
        using Parent::getEtotSlowDown;
        using Parent::getEnergyErrorSlowDown;
        using Parent::getEnergyErrorSlowDownFromBackup;
        using Parent::getEtotSlowDownRefFromBackup;
        using Parent::getEtotSlowDownFromBackup;
        using Parent::writeBackSlowDownParticlesIter;
        using Parent::writeBackSlowDownParticles;

    protected:
        // transformation factors
        Float gt_drift_inv_ = 0.0;  ///< integrated inverse time transformation factor for drift: dt(drift) = ds/gt_drift_inv_
        Float gt_kick_inv_ = 0.0;   ///< inverse time transformation factor for kick: dt(kick) = ds/gt_kick_inv_

        //! reserve memory of the data owned by the dynamics
        void reserveDynamicsMem(const int _nmax) {
            (void)_nmax;
        }

        //! clear the data owned by the dynamics
        void clearDynamics() {
            this->clearSlowDownEnergy();
            gt_drift_inv_ = 0.0;
            gt_kick_inv_ = 0.0;
        }

        //! copy the data owned by the dynamics
        void copyDynamics(const ARDynamics& _sym) {
            this->copySlowDownEnergy(_sym);
            gt_drift_inv_ = _sym.gt_drift_inv_;
            gt_kick_inv_ = _sym.gt_kick_inv_;
        }

        //! calc force, potential and inverse time transformation factor for one pair of particles
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _i: particle i index
          @param[in] _j: particle j index
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvTwo(const Float& _inv_nest_sd, const int _i, const int _j) {
            ASSERT(_i>=0&&_i<particles.getSize());
            ASSERT(_j>=0&&_j<particles.getSize());
            
            // calculate pair interaction
            Force fij[2];
            Float epotij;
            Float gt_kick_inv = manager->interaction.calcInnerAccPotAndGTKickInvTwo(fij[0], fij[1], epotij, particles[_i], particles[_j]);

            // scale binary pair force with slowdown 
            force_[_i].acc_in[0] += fij[0].acc_in[0]*_inv_nest_sd;
            force_[_i].acc_in[1] += fij[0].acc_in[1]*_inv_nest_sd;
            force_[_i].acc_in[2] += fij[0].acc_in[2]*_inv_nest_sd;
            force_[_j].acc_in[0] += fij[1].acc_in[0]*_inv_nest_sd;
            force_[_j].acc_in[1] += fij[1].acc_in[1]*_inv_nest_sd;
            force_[_j].acc_in[2] += fij[1].acc_in[2]*_inv_nest_sd;

            epot_    += epotij;
            epot_sd_ += epotij*_inv_nest_sd;

            // scale gtgrad with slowdown
            force_[_i].gtgrad[0] += fij[0].gtgrad[0]*_inv_nest_sd;
            force_[_i].gtgrad[1] += fij[0].gtgrad[1]*_inv_nest_sd;
            force_[_i].gtgrad[2] += fij[0].gtgrad[2]*_inv_nest_sd;
            force_[_j].gtgrad[0] += fij[1].gtgrad[0]*_inv_nest_sd;
            force_[_j].gtgrad[1] += fij[1].gtgrad[1]*_inv_nest_sd;
            force_[_j].gtgrad[2] += fij[1].gtgrad[2]*_inv_nest_sd;
            
            return gt_kick_inv * _inv_nest_sd;
        }

        //! calc force, potential and inverse time transformation factor for one particle by walking binary tree
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _i: particle index
          @param[in] _bin: binary tree for walking
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvOneTreeIter(const Float& _inv_nest_sd, const int _i, AR::BinaryTree<Tparticle>& _bin) {
            Float gt_kick_inv=0.0;
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) // particle - tree
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(_inv_nest_sd, _i, *(_bin.getMemberAsTree(k)));
                else  // particle - particle
                    gt_kick_inv += calcAccPotAndGTKickInvTwo(_inv_nest_sd, _i, _bin.getMemberIndex(k));
            }
            return gt_kick_inv;
        }

        //! calc crossing force, potential and inverse time transformation factor between binary tree i and binary tree j
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _bini: binary tree i for walking
          @param[in] _binj: binary tree j for walking
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvCrossTreeIter(const Float& _inv_nest_sd, AR::BinaryTree<Tparticle>& _bini, AR::BinaryTree<Tparticle>& _binj) {
            ASSERT(&_bini!=&_binj);
            Float gt_kick_inv=0.0;
            for (int k=0; k<2; k++) { 
                if (_bini.isMemberTree(k))  // tree - tree
                    gt_kick_inv += calcAccPotAndGTKickInvCrossTreeIter(_inv_nest_sd, *(_bini.getMemberAsTree(k)), _binj);
                else  // particle - tree
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(_inv_nest_sd, _bini.getMemberIndex(k), _binj);
            }
            return gt_kick_inv;
        }

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary to drift pos
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvTreeIter(const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) {
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
            Float gt_kick_inv = 0.0;

            // check left 
            if (_bin.isMemberTree(0)) { // left is tree
                auto* bin_left = _bin.getMemberAsTree(0);
                // inner interaction of left tree
                gt_kick_inv += calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_left);

                if (_bin.isMemberTree(1)) { // right is tree
                    auto* bin_right = _bin.getMemberAsTree(1);

                    // inner interaction of right tree
                    gt_kick_inv += calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_right);
                    
                    // cross interaction
                    gt_kick_inv += calcAccPotAndGTKickInvCrossTreeIter(inv_nest_sd, *bin_left, *bin_right);
                }
                else { // right is particle
                    // cross interaction from particle j to tree left
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(inv_nest_sd,  _bin.getMemberIndex(1), *bin_left);
                }
            }
            else { // left is particle
                if (_bin.isMemberTree(1)) { // right is tree
                    auto* bin_right = _bin.getMemberAsTree(1);
                    // inner interaction of right tree
                    gt_kick_inv += calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_right);

                    // cross interaction from particle i to tree right
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(inv_nest_sd, _bin.getMemberIndex(0), *bin_right);
                }
                else { // right is particle
                    // particle - particle interaction
                    gt_kick_inv += calcAccPotAndGTKickInvTwo(inv_nest_sd, _bin.getMemberIndex(0), _bin.getMemberIndex(1));
                }
            }


            return gt_kick_inv;
        }

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        inline Float calcAccPotAndGTKickInv() {
            epot_ = 0.0;
            epot_sd_ = 0.0;
            for (int i=0; i<force_.getSize(); i++) force_[i].clear();

            Float gt_kick_inv = calcAccPotAndGTKickInvTreeIter(1.0, info.getBinaryTreeRoot());
            // pertuber force
            manager->interaction.calcAccPert(force_.getDataAddress(), particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());
            return gt_kick_inv;
        }

        //! kick energy and time transformation function for drift of binary tree 
        /*!
          @param[in] _dt: time step
          @param[in] _vel_sd_up: upper cm sd vel
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary tree for kick etot and calc dgt_drift
          \return gt_drift_inv change 
        */
        Float kickEtotAndGTDriftTreeIter(const Float& _dt, const Float* _vel_sd_up, const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) {
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
            Float dgt_drift_inv = 0.0;
            Float de = 0.0;

            Float* vel_cm = _bin.getVel();

            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    Float* vel = bink->getVel();
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]}; 
                    dgt_drift_inv += kickEtotAndGTDriftTreeIter(_dt, vel_sd, inv_nest_sd, *bink);
                }
                else {
                    int i = _bin.getMemberIndex(k);
                    ASSERT(i>=0&&i<particles.getSize());
                    ASSERT(&particles[i]==_bin.getMember(k));
                    
                    Float* gtgrad = force_[i].gtgrad;
                    Float* pert   = force_[i].acc_pert;
                    Float* vel = particles[i].getVel();
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]};

                    de += particles[i].mass * (vel[0] * pert[0] +
                                               vel[1] * pert[1] +
                                               vel[2] * pert[2]);
                    dgt_drift_inv +=  (vel_sd[0] * gtgrad[0] +
                                       vel_sd[1] * gtgrad[1] +
                                       vel_sd[2] * gtgrad[2]);
                }
            }
            etot_ref_     += _dt * de;
            etot_sd_ref_ += _dt * de;

            return dgt_drift_inv;
        }

        //! kick energy and time transformation function for drift
        /*!
          @param[in] _dt: time step
        */
        void kickEtotAndGTDrift(const Float _dt) {
            // the particle cm velocity is zero (assume in rest-frame)
            ASSERT(!particles.isOriginFrame());
            auto& bin_root=info.getBinaryTreeRoot();
            Float vel_cm[3] = {0.0,0.0,0.0};
            Float sd_factor=1.0;

            Float dgt_drift_inv = kickEtotAndGTDriftTreeIter(_dt, vel_cm, sd_factor, bin_root);
            gt_drift_inv_ += dgt_drift_inv*_dt;
        }

    public:
        //! update slowdown factor based on perturbation and record slowdown energy change
        /*! Update slowdown inner and global.
            @param [in] _update_energy_flag: Record cumulative slowdown energy change if true;
            @param [in] _stable_check_flag: check whether the binary tree is stable if true;
         */
        void updateSlowDownAndCorrectEnergy(const bool _update_energy_flag, const bool _stable_check_flag) {
            auto& bin_root = info.getBinaryTreeRoot();
            auto& sd_root = bin_root.slowdown;

            Float sd_backup = sd_root.getSlowDownFactor();

            // when the maximum inner slowdown is large, the outer should not be slowed down since the system may not be stable.
            //if (inner_sd_change_flag&&sd_org_inner_max<1000.0*manager->slowdown_pert_ratio_ref) sd_root.setSlowDownFactor(1.0);
            //if (time_>=sd_root.getUpdateTime()) {
            sd_root.pert_in = manager->interaction.calcPertFromBinary(bin_root);
            sd_root.pert_out = 0.0;
            Float t_min_sq= NUMERIC_FLOAT_MAX;
            manager->interaction.calcSlowDownPert(sd_root.pert_out, t_min_sq, getTime(), particles.cm, perturber);
            sd_root.timescale = std::min(sd_root.getTimescaleMax(), sqrt(t_min_sq));

            //Float period_amplify_max = NUMERIC_FLOAT_MAX;
            if (_stable_check_flag) {
                // check whether the system is stable for 10000 out period and the apo-center is below break criterion
                Float stab = bin_root.stableCheckIter(bin_root,10000*bin_root.period);
                Float apo = bin_root.semi*(1+bin_root.ecc);
                if (stab<1.0 && apo<info.r_break_crit) {
                    sd_root.period = bin_root.period;
                    sd_root.calcSlowDownFactor();
                }
                else sd_root.setSlowDownFactor(1.0);

                // stablility criterion
                // The slowdown factor should not make the system unstable, thus the Qst/Q set the limitation of the increasing of inner semi-major axis.
                //if (stab>0 && stab != NUMERIC_FLOAT_MAX) {
                //    Float semi_amplify_max =  std::max(Float(1.0),1.0/stab);
                //    period_amplify_max = pow(semi_amplify_max,3.0/2.0);
                //}
            }
            else if (bin_root.semi>0) {
                sd_root.period = bin_root.period;
                sd_root.calcSlowDownFactor();
            }
            else sd_root.setSlowDownFactor(1.0);

            //sd_root.increaseUpdateTimeOnePeriod();
            //}

            // inner binary slowdown
            Float sd_org_inner_max = 0.0;
            bool inner_sd_change_flag=false;
            int n_bin = info.binarytree.getSize();
            for (int i=0; i<n_bin-1; i++) {
                auto& bini = info.binarytree[i];
                //if (time_>=bini.slowdown.getUpdateTime()) {
                bini.calcCenterOfMass();
                calcSlowDownInnerBinary(bini);

                //sdi->slowdown.increaseUpdateTimeOnePeriod();
                sd_org_inner_max = std::max(bini.slowdown.getSlowDownFactorOrigin(),sd_org_inner_max);
                inner_sd_change_flag=true;
                //}
            }


            if (_update_energy_flag) {
                Float ekin_sd_bk = ekin_sd_;
                Float epot_sd_bk = epot_sd_;
                Float H_sd_bk = getHSlowDown();
                if(inner_sd_change_flag) {
                    Float gt_kick_inv_new = calcAccPotAndGTKickInv();
                    gt_drift_inv_ += gt_kick_inv_new - gt_kick_inv_;
                    gt_kick_inv_ = gt_kick_inv_new;
                    calcEKin();
                }
                else {
                    Float kappa_inv = 1.0/sd_root.getSlowDownFactor();
                    Float gt_kick_inv_new = gt_kick_inv_*sd_backup*kappa_inv;
                    gt_drift_inv_ += gt_kick_inv_new - gt_kick_inv_;
                    gt_kick_inv_ = gt_kick_inv_new;
                    ekin_sd_ = ekin_*kappa_inv;
                    epot_sd_ = epot_*kappa_inv;
                }
                Float de_sd = (ekin_sd_ - ekin_sd_bk) + (epot_sd_ - epot_sd_bk);
                etot_sd_ref_ += de_sd;

                Float dH_sd = getHSlowDown() - H_sd_bk;

                // add slowdown change to the global slowdown energy
                de_sd_change_cum_ += de_sd;
                dH_sd_change_cum_ += dH_sd;
            }
        }

        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time: real physical time to initialize
        */
        void initialIntegration(const Float _time) {
            ASSERT(checkParams());

            // particle number and data address
            const int n_particle = particles.getSize();

            // resize force array
            force_.resizeNoInitialize(n_particle);

            // Initial intgrt value t (avoid confusion of real time when slowdown is used)
            time_ = _time;

            // check particle number
            ASSERT(particles.getSize()>=2);

            // reset particle modification flag
            particles.setModifiedFalse();

            // check the center-of-mass initialization
            if(particles.isOriginFrame()) {
                particles.calcCenterOfMass();
                particles.shiftToCenterOfMassFrame();
                for (int i=0; i<info.binarytree.getSize(); i++) {
                    auto& bin = info.binarytree[i];
                    bin.pos[0] -= particles.cm.pos[0];
                    bin.pos[1] -= particles.cm.pos[1];
                    bin.pos[2] -= particles.cm.pos[2];
                    bin.vel[0] -= particles.cm.vel[0];
                    bin.vel[1] -= particles.cm.vel[1];
                    bin.vel[2] -= particles.cm.vel[2];
                }
            }
            ASSERT(info.getBinaryTreeRoot().pos[0]*info.getBinaryTreeRoot().pos[0]<1e-10);
            ASSERT(info.getBinaryTreeRoot().vel[0]*info.getBinaryTreeRoot().vel[0]<1e-10);

            for (int i=0; i<info.binarytree.getSize(); i++) 
                info.binarytree[i].slowdown.initialSlowDownReference(manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);

            updateSlowDownAndCorrectEnergy(false,true);

            gt_kick_inv_ = calcAccPotAndGTKickInv();

            // initially gt_drift 
            gt_drift_inv_ = gt_kick_inv_;


            calcEKin();

            etot_ref_ = ekin_ + epot_;
            etot_sd_ref_ = ekin_sd_ + epot_sd_;

            Float de_sd = etot_sd_ref_ - etot_ref_;

            // add slowdown change to the global slowdown energy
            de_sd_change_cum_ += de_sd;
            dH_sd_change_cum_ = 0.0;

        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // inverse time transformation factor for drift
                Float gt_drift_inv = gt_drift_inv_;

                // drift
                Float dt_drift = ds_drift/gt_drift_inv;

                // drift time and postion
                driftTimeAndPos(dt_drift);
                _time_table[i] = time_;

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                //! calc force, potential and inverse time transformation factor for kick
                Float gt_kick_inv = calcAccPotAndGTKickInv();

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // back up gt_kick 
                gt_kick_inv_ = gt_kick_inv;
                // kick total energy and inverse time transformation factor for drift
                kickEtotAndGTDrift(dt_kick);
                // kick half step for velocity
                kickVel(0.5*dt_kick);

                // calculate kinetic energy
                calcEKin();
            }
        }

        //! integration for two body one step
        /*! For two-body problem the calculation can be much symplified to improve performance. 
          Besides, the slow-down factor calculation is embedded in the Drift (for time) and Kick (for perturbation). 
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().         
        */
        void integrateTwoOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group number
            const int nloop = manager->step.getCDPairSize();
            
            const int n_particle = particles.getSize();
            ASSERT(n_particle==2);

            const Float kappa_inv = 1.0/info.getBinaryTreeRoot().slowdown.getSlowDownFactor();

            Tparticle* particle_data = particles.getDataAddress();
            Float mass1 = particle_data[0].mass;
            Float* pos1 = particle_data[0].getPos();
            Float* vel1 = particle_data[0].getVel();

            Float mass2 = particle_data[1].mass;
            Float* pos2 = particle_data[1].getPos();
            Float* vel2 = particle_data[1].getVel();

            Force* force_data = force_.getDataAddress();
            Float* acc1 = force_data[0].acc_in;
            Float* pert1= force_data[0].acc_pert;

            Float* acc2 = force_data[1].acc_in;
            Float* pert2= force_data[1].acc_pert;
            Float* gtgrad1 = force_data[0].gtgrad;
            Float* gtgrad2 = force_data[1].gtgrad;

#ifdef AR_DEBUG_PRINT_DKD
            std::cout<<"K "<<time_<<" "
                     <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                     <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                     <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds = manager->step.getCK(i)*_ds;
                // inverse time transformation factor for drift
                Float gt_inv = gt_drift_inv_;
                // drift
                Float dt = ds/gt_inv;
                ASSERT(!ISNAN(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                _time_table[i] = time_;

                Float dt_sd = dt*kappa_inv;

                // drift position
                pos1[0] += dt_sd * vel1[0];
                pos1[1] += dt_sd * vel1[1];
                pos1[2] += dt_sd * vel1[2];

                pos2[0] += dt_sd * vel2[0];
                pos2[1] += dt_sd * vel2[1];
                pos2[2] += dt_sd * vel2[2];

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                gt_inv = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm, perturber, _time_table[i]);

                ASSERT(!ISNAN(epot_));

#ifdef AR_DEBUG_PRINT_DKD
                if (i>0)
                    std::cout<<"K "<<time_<<" "
                             <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                             <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                             <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick half step for velocity
                Float dvel1[3], dvel2[3];

                // time step for kick
                gt_inv *= kappa_inv;

                dt = 0.5*ds/gt_inv;

                dvel1[0] = dt * (acc1[0]*kappa_inv + pert1[0]);
                dvel1[1] = dt * (acc1[1]*kappa_inv + pert1[1]);
                dvel1[2] = dt * (acc1[2]*kappa_inv + pert1[2]);

                dvel2[0] = dt * (acc2[0]*kappa_inv + pert2[0]);
                dvel2[1] = dt * (acc2[1]*kappa_inv + pert2[1]);
                dvel2[2] = dt * (acc2[2]*kappa_inv + pert2[2]);

                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];

                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];

#ifdef AR_DEBUG_PRINT_DKD
                std::cout<<"D "<<time_<<" "
                         <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                         <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                         <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick total energy and time transformation factor for drift
                etot_ref_ += 2.0*dt * (mass1* (vel1[0] * pert1[0] + 
                                               vel1[1] * pert1[1] + 
                                               vel1[2] * pert1[2]) +
                                       mass2* (vel2[0] * pert2[0] + 
                                               vel2[1] * pert2[1] + 
                                               vel2[2] * pert2[2]));

                // back up gt_kick_inv
                gt_kick_inv_ = gt_inv;

                // integrate gt_drift_inv
                gt_drift_inv_ +=  2.0*dt*kappa_inv*kappa_inv* (vel1[0] * gtgrad1[0] +
                                                               vel1[1] * gtgrad1[1] +
                                                               vel1[2] * gtgrad1[2] +
                                                               vel2[0] * gtgrad2[0] +
                                                               vel2[1] * gtgrad2[1] +
                                                               vel2[2] * gtgrad2[2]);


                // kick half step for velocity
                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];
                
                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];
                                                                                                                
                // calculate kinetic energy
                ekin_ = 0.5 * (mass1 * (vel1[0]*vel1[0]+vel1[1]*vel1[1]+vel1[2]*vel1[2]) +
                               mass2 * (vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]));

            }

            // make consistent slowdown inner energy 
            etot_sd_ref_ = etot_ref_*kappa_inv;
            ekin_sd_ = ekin_*kappa_inv;
            epot_sd_ = epot_*kappa_inv;
        }

        //! get Hamiltonian
        Float getH() const {
            //return (ekin_ - etot_ref_)/gt_drift_inv_ + epot_/gt_kick_inv_;
            return (ekin_ + epot_ - etot_ref_)/gt_kick_inv_;
        }

        //! get Hamiltonian from backup data
        Float getHFromBackup(Float* _bk) const {
            Float& etot_ref =_bk[1];
            Float& ekin = _bk[2];
            Float& epot = _bk[3];
            //Float& gt_drift_inv = _bk[13];
            Float& gt_kick_inv  = _bk[14];
            return (ekin + epot - etot_ref)/gt_kick_inv;
            //return (ekin - etot_ref)/gt_drift_inv + epot/gt_kick_inv;
        }

        //! get slowdown Hamiltonian
        Float getHSlowDown() const {
            //return (ekin_sd_ - etot_sd_ref_)/gt_drift_inv_ + epot_sd_/gt_kick_inv_;
            return (ekin_sd_ + epot_sd_ - etot_sd_ref_)/gt_kick_inv_;
        }

        //! get slowdown Hamiltonian from backup data
        Float getHSlowDownFromBackup(Float* _bk) const {
            Float& etot_sd_ref =_bk[6];
            Float& ekin_sd = _bk[7];
            Float& epot_sd = _bk[8];
            //Float& gt_drift_inv = _bk[13];
            Float& gt_kick_inv  = _bk[14];
            //return (ekin_sd - etot_sd_ref)/gt_drift_inv + epot_sd/gt_kick_inv;
            return (ekin_sd + epot_sd - etot_sd_ref)/gt_kick_inv;
        }

        //! Get integrated inverse time transformation factor
        /*! In TTF case, it is calculated by integrating \f$ \frac{dg}{dt} = \sum_k \frac{\partial g}{\partial \vec{r_k}} \bullet \vec{v_k} \f$.
          Notice last step is the sub-step in one symplectic loop
          \return inverse time transformation factor for drift
        */
        Float getGTDriftInv() const {
            return gt_drift_inv_;
        }

    };

}
