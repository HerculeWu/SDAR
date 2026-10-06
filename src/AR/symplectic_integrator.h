#pragma once

#include <functional>
#include "Common/list.h"
#include "Common/particle_group.h"
#include "AR/symplectic_step.h"
#include "AR/regularization_step_control.h"
#include "AR/force.h"
#include "AR/slow_down.h"
#include "AR/profile.h"
#include "AR/information.h"
#include "AR/interrupt.h"
#include "AR/variant.h"
#include "AR/ar_dynamics.h"
#include "AR/interaction_concept.h"

//! Algorithmic regularization (time transformed explicit symplectic integrator) namespace
/*!
  All major AR classes and related acceleration functions (typedef) are defined
*/
namespace AR {

    //! print features
    void printFeatures(std::ostream & fout) {
        if (std::same_as<DefaultTimeTransformation, TTL>) fout<<"Use AR TTL method\n";
        else fout<<"Use AR LogH method\n";
        if (std::same_as<DefaultSlowDownScheme, TreeSlowDown>) fout<<"Use slowdown Tree method\n";
        if (std::same_as<DefaultSlowDownScheme, ArraySlowDown>) fout<<"Use slowdown array method\n";
#ifdef AR_SLOWDOWN_TIMESCALE
        fout<<"Use slowdown timescale criterion\n";         
#endif
#ifdef AR_SLOWDOWN_MASSRATIO
        fout<<"Use slowdown mass ratio criterion\n";
#endif
    }

    //! print debug features
    void printDebugFeatures(std::ostream & fout) {
#ifdef AR_DEBUG
        fout<<"Debug mode: AR\n";
#endif        
    }

    //! print reference to cite
    void printReference(std::ostream & fout, const int offset=4) {
        for (int i=0; i<offset; i++) fout<<" ";
        fout<<"SDAR: Wang L., Nitadori K., Makino J., 2020, MNRAS, 493, 3398"
            <<std::endl;
    }

    //! Time Transformed Symplectic integrator manager
    /*! Tmethod is the class contain the interaction function, see sample of interaction.h:\n
     */
    template <class Tmethod>
    class TimeTransformedSymplecticManager {
    public:
        Float time_error_max; ///> maximum time error (absolute), should be positive and larger than round-off error 
        Float energy_error_relative_max; ///> maximum energy error requirement 
        Float time_step_min;        ///> minimum real time step allown
        Float ds_scale;            ///> scaling factor to determine ds
        Float slowdown_pert_ratio_ref;   ///> slowdown perturbation /inner ratio reference factor
#ifdef AR_SLOWDOWN_MASSRATIO
        Float slowdown_mass_ref;         ///> slowdown mass factor reference
#endif
        Float slowdown_timescale_max;       ///> slowdown maximum timescale to calculate maximum slowdown factor
        long long unsigned int step_count_max; ///> maximum step counts
        int interrupt_detection_option;    ///> 1: detect interruption; 0: no detection
        
        Tmethod interaction; ///> class contain interaction function
        SymplecticStep step;  ///> class to manager kick drift step

        //! constructor
        TimeTransformedSymplecticManager(): time_error_max(Float(-1.0)), energy_error_relative_max(Float(-1.0)), time_step_min(Float(-1.0)), ds_scale(1.0), slowdown_pert_ratio_ref(Float(-1.0)), 
#ifdef AR_SLOWDOWN_MASSRATIO
                                            slowdown_mass_ref(Float(-1.0)), 
#endif
                                            slowdown_timescale_max(0.0),
                                            step_count_max(0), interrupt_detection_option(0), interaction(), step() {}

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            //ASSERT(time_error_max>ROUND_OFF_ERROR_LIMIT);
            ASSERT(time_error_max>0.0);
            ASSERT(energy_error_relative_max>ROUND_OFF_ERROR_LIMIT);
            //ASSERT(time_step_min>ROUND_OFF_ERROR_LIMIT);
            ASSERT(time_step_min>0.0);
            ASSERT(ds_scale>0.0);
            ASSERT(slowdown_pert_ratio_ref>0.0);
#ifdef AR_SLOWDOWN_MASSRATIO
            ASSERT(slowdown_mass_ref>0.0);
#endif
            ASSERT(slowdown_timescale_max>0);
            ASSERT(step_count_max>0);
            ASSERT(step.getOrder()>0);
            ASSERT(interaction.checkParams());
            return true;
        }

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBinary(FILE *_fout) {
            size_t size = sizeof(*this) - sizeof(interaction) - sizeof(step);
            fwrite(this, size, 1,_fout);
            interaction.writeBinary(_fout);
            step.writeBinary(_fout);
        }

        //! read class data with BINARY format and initial the array
        /*! @param[in] _fin: file IO for read
          @param[in] _version: version for reading. 0: default; 1: missing ds_scale
         */
        void readBinary(FILE *_fin, int _version=0) {
            if (_version==0) {
                size_t size = sizeof(*this) - sizeof(interaction) - sizeof(step);
                size_t rcount = fread(this, size, 1, _fin);
                if (rcount<1) {
                    std::cerr<<"Error: TimeTransformedSymplecticManager parameter reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                    abort();
                }
            }
            else if (_version==1) {
                size_t rcount = fread(this, sizeof(Float), 3, _fin);
                if (rcount<3) {
                    std::cerr<<"Error: TimeTransformedSymplecticManager parameter data reading fails! requiring data number is 3, only obtain "<<rcount<<".\n";
                    abort();
                }
                ds_scale=1.0;
                size_t size = sizeof(*this) - sizeof(interaction) - sizeof(step) - 4*sizeof(Float);
                rcount = fread(&slowdown_pert_ratio_ref, size, 1, _fin);
                if (rcount<1) {
                    std::cerr<<"Error: TimeTransformedSymplecticManager parameter data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                    abort();
                }
            }
            else {
                std::cerr<<"Error: TimeTransformedSymplecticManager.readBinary unknown version "<<_version<<", should be 0 or 1."<<std::endl;
                abort();
            }
            interaction.readBinary(_fin);
            step.readBinary(_fin);
        }

        //! print parameters
        void print(std::ostream & _fout) const{
            _fout<<"time_error_max            : "<<time_error_max<<std::endl
                 <<"energy_error_relative_max : "<<energy_error_relative_max<<std::endl 
                 <<"time_step_min             : "<<time_step_min<<std::endl
                 <<"slowdown_pert_ratio_ref   : "<<slowdown_pert_ratio_ref<<std::endl
#ifdef AR_SLOWDOWN_MASSRATIO
                 <<"slowdown_mass_ref         : "<<slowdown_mass_ref<<std::endl
#endif
                 <<"slowdown_timescale_max    : "<<slowdown_timescale_max<<std::endl
                 <<"step_count_max            : "<<step_count_max<<std::endl
                 <<"ds_scale                  : "<<ds_scale<<std::endl;
            interaction.print(_fout);
            step.print(_fout);
        }
    };

    //! Time Transformed Symplectic integrator class for a group of particles
    /*! The basic steps to use the integrator \n
      1. Add particles (particles.addParticle/particles.linkParticleList)  \n
      2. Initial system (initial) \n
      3. Integration (integrateOneStep/integrateToTime) \n
      Requirement for Tparticle class, public memebers: pos[3], vel[3], mass\n
      Template dependence: Tparticle: particle type; Tpcm: particle cm type  Tpert: perturber class type, Tmethod: interaction class;
    */
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo, TimeTransformation Ttransform = DefaultTimeTransformation, SlowDownScheme Tslowdown = DefaultSlowDownScheme>
    class TimeTransformedSymplecticIntegrator: public ARDynamics<Ttransform, Tslowdown, Tparticle, Tpcm, Tpert, Tmethod, Tinfo> {
    private:
        typedef ARDynamics<Ttransform, Tslowdown, Tparticle, Tpcm, Tpert, Tmethod, Tinfo> Dynamics;

        static constexpr bool is_ttl = std::same_as<Ttransform, TTL>;
        static constexpr bool is_slowdown_array = std::same_as<Tslowdown, ArraySlowDown>;
        static constexpr bool is_slowdown_tree = std::same_as<Tslowdown, TreeSlowDown>;
        static constexpr bool has_slowdown = is_slowdown_array || is_slowdown_tree;

        // the time transformation states what it requires of the interaction class
        static_assert(InteractionOf<Tmethod, Ttransform, Tparticle, Tpcm, Tpert>, 
                      "the interaction class does not provide what the time transformation requires: see AR::LogHInteraction and AR::TTLInteraction");

    protected:
        using Dynamics::time_;
        using Dynamics::etot_ref_;
        using Dynamics::ekin_;
        using Dynamics::epot_;
        using Dynamics::de_change_interrupt_;
        using Dynamics::dH_change_interrupt_;
        using Dynamics::force_;

    public:
        typedef typename Dynamics::Force Force;
        using Dynamics::manager;
        using Dynamics::particles;
        using Dynamics::perturber;
        using Dynamics::info;
        using Dynamics::profile;
        
        //! Constructor
        TimeTransformedSymplecticIntegrator() {}

        //! reserve memory for force
        /*! The size of force depends on the particle data size.Thus particles should be added first before call this function
        */
        void reserveIntegratorMem() {
            // force array always allocated local memory
            int nmax = particles.getSizeMax();
            ASSERT(nmax>0);
            force_.setMode(COMM::ListMode::local);
            force_.reserveMem(nmax);
            this->reserveDynamicsMem(nmax);
        }

        //! Clear function
        /*! Free dynamical memory space allocated
         */
        void clear() {
            this->clearBase();
            this->clearDynamics();
        }

        //! destructor
        ~TimeTransformedSymplecticIntegrator() {
            clear();
        }

        //! operator = 
        /*! Copy function will remove the local data and also copy the particle data or the link
         */
        TimeTransformedSymplecticIntegrator& operator = (const TimeTransformedSymplecticIntegrator& _sym) {
            clear();
            this->copyBase(_sym);
            this->copyDynamics(_sym);

            return *this;
        }

        //! whether the time transformation is TTL (otherwise LogH)
        static constexpr bool isTTL() { return is_ttl; }

        //! whether a Slow-down scheme is used
        static constexpr bool hasSlowDown() { return has_slowdown; }

        //! whether the Slow-down scheme is Array slow-down
        static constexpr bool isSlowDownArray() { return is_slowdown_array; }

        //! whether the Slow-down scheme is Tree slow-down
        static constexpr bool isSlowDownTree() { return is_slowdown_tree; }

        //! number of inner slowed binaries, whatever the Slow-down scheme
        /*! Array slow-down: the list of slowed binaries; Tree slow-down: every binary of the tree; no Slow-down: zero
         */
        int getSlowDownInnerNumber() const {
            if constexpr (is_slowdown_array) return this->binary_slowdown.getSize();
            else if constexpr (is_slowdown_tree) return info.binarytree.getSize();
            else return 0;
        }

        //! Slow-down of the i-th inner slowed binary, whatever the Slow-down scheme (see getSlowDownInnerNumber)
        SlowDown& getSlowDownInner(const int _i) requires has_slowdown {
            if constexpr (is_slowdown_array) return this->binary_slowdown[_i]->slowdown;
            else return info.binarytree[_i].slowdown;
        }

        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time: real physical time to initialize
        */
        void initialIntegration(const Float _time) {
            ASSERT(this->checkParams());

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

            if constexpr (has_slowdown) {
                if constexpr (is_slowdown_tree) {
                    for (int i=0; i<info.binarytree.getSize(); i++) 
                        info.binarytree[i].slowdown.initialSlowDownReference(manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);
                }
                else {
                    this->binary_slowdown.increaseSizeNoInitialize(1);
                    this->binary_slowdown[0] = &info.getBinaryTreeRoot();

                    // set slowdown reference
                    SlowDown& slowdown_root = info.getBinaryTreeRoot().slowdown;

                    // slowdown for the system
                    slowdown_root.initialSlowDownReference(manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);

                    if (particles.getSize()>2) {
                        this->findSlowDownInner(time_);
                        // update c.m. of binaries 
                        //updateCenterOfMassForBinaryWithSlowDownInner();
                    }
                }

                this->updateSlowDownAndCorrectEnergy(false,true);

                if constexpr (is_ttl) {
                    this->gt_kick_inv_ = this->calcAccPotAndGTKickInv();

                    // initially gt_drift 
                    this->gt_drift_inv_ = this->gt_kick_inv_;

                }
                else {
                    this->calcAccPotAndGTKickInv();
                }

                this->calcEKin();

                etot_ref_ = ekin_ + epot_;
                this->etot_sd_ref_ = this->ekin_sd_ + this->epot_sd_;

                Float de_sd = this->etot_sd_ref_ - etot_ref_;

                // add slowdown change to the global slowdown energy
                this->de_sd_change_cum_ += de_sd;
                this->dH_sd_change_cum_ = 0.0;

            }
            else {
                Tparticle* particle_data = particles.getDataAddress();
                Force* force_data = force_.getDataAddress();

                if constexpr (is_ttl) {
                    this->gt_kick_inv_ = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm,  perturber, _time);

                    // initially gt_drift 
                    this->gt_drift_inv_ = this->gt_kick_inv_;

                }
                else {
                    manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm,  perturber, _time);
                }

                // calculate kinetic energy
                this->calcEKin();

                // initial total energy
                etot_ref_ = ekin_ + epot_;

            }
        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(this->checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // inverse time transformation factor for drift
                Float gt_drift_inv;
                if constexpr (is_ttl) {
                    gt_drift_inv = this->gt_drift_inv_;
                }
                else {
                    if constexpr (has_slowdown) {
                        gt_drift_inv = manager->interaction.calcGTDriftInv(this->ekin_sd_-this->etot_sd_ref_); // pt = -etot
                    }
                    else {
                        gt_drift_inv = manager->interaction.calcGTDriftInv(ekin_-etot_ref_); // pt = -etot
                    }
                }

                // drift
                Float dt_drift = ds_drift/gt_drift_inv;

                // drift time and postion
                this->driftTimeAndPos(dt_drift);
                _time_table[i] = time_;

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                //! calc force, potential and inverse time transformation factor for kick
                Float gt_kick_inv = this->calcAccPotAndGTKickInv();

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                this->kickVel(0.5*dt_kick);

                if constexpr (is_ttl) {
                    // back up gt_kick 
                    this->gt_kick_inv_ = gt_kick_inv;
                    // kick total energy and inverse time transformation factor for drift
                    this->kickEtotAndGTDrift(dt_kick);
                }
                else {
                    // kick total energy 
                    this->kickEtot(dt_kick);
                }
                // kick half step for velocity
                this->kickVel(0.5*dt_kick);

                // calculate kinetic energy
                this->calcEKin();
            }
        }


        //! integration for two body one step
        /*! For two-body problem the calculation can be much symplified to improve performance. 
          Besides, the slow-down factor calculation is embedded in the Drift (for time) and Kick (for perturbation). 
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().         
        */
        void integrateTwoOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(this->checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group number
            const int nloop = manager->step.getCDPairSize();
            
            const int n_particle = particles.getSize();
            ASSERT(n_particle==2);

            [[maybe_unused]] Float kappa_inv = 1.0;
            if constexpr (has_slowdown) {
                kappa_inv = 1.0/info.getBinaryTreeRoot().slowdown.getSlowDownFactor();
            }

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
            [[maybe_unused]] Float* gtgrad1 = NULL;
            [[maybe_unused]] Float* gtgrad2 = NULL;
            if constexpr (is_ttl) {
                gtgrad1 = force_data[0].gtgrad;
                gtgrad2 = force_data[1].gtgrad;
            }

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
                Float gt_inv;
                if constexpr (is_ttl) {
                    gt_inv = this->gt_drift_inv_;
                }
                else {
                    if constexpr (has_slowdown) {
                        gt_inv = manager->interaction.calcGTDriftInv(this->ekin_sd_-this->etot_sd_ref_); // pt = -etot_sd
                    }
                    else {
                        gt_inv = manager->interaction.calcGTDriftInv(ekin_-etot_ref_); // pt = -etot
                    }
                }
                // drift
                Float dt = ds/gt_inv;
                ASSERT(!ISNAN(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                _time_table[i] = time_;

                if constexpr (has_slowdown) {
                    Float dt_sd = dt*kappa_inv;

                    // drift position
                    pos1[0] += dt_sd * vel1[0];
                    pos1[1] += dt_sd * vel1[1];
                    pos1[2] += dt_sd * vel1[2];

                    pos2[0] += dt_sd * vel2[0];
                    pos2[1] += dt_sd * vel2[1];
                    pos2[2] += dt_sd * vel2[2];
                }
                else {
                    // drift position
                    pos1[0] += dt * vel1[0];
                    pos1[1] += dt * vel1[1];
                    pos1[2] += dt * vel1[2];

                    pos2[0] += dt * vel2[0];
                    pos2[1] += dt * vel2[1];
                    pos2[2] += dt * vel2[2];
                }

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

                if constexpr (has_slowdown) {
                    // time step for kick
                    gt_inv *= kappa_inv;

                    dt = 0.5*ds/gt_inv;

                    dvel1[0] = dt * (acc1[0]*kappa_inv + pert1[0]);
                    dvel1[1] = dt * (acc1[1]*kappa_inv + pert1[1]);
                    dvel1[2] = dt * (acc1[2]*kappa_inv + pert1[2]);

                    dvel2[0] = dt * (acc2[0]*kappa_inv + pert2[0]);
                    dvel2[1] = dt * (acc2[1]*kappa_inv + pert2[1]);
                    dvel2[2] = dt * (acc2[2]*kappa_inv + pert2[2]);
                }
                else {
                    dt = 0.5*ds/gt_inv;

                    dvel1[0] = dt * (acc1[0] + pert1[0]);
                    dvel1[1] = dt * (acc1[1] + pert1[1]);
                    dvel1[2] = dt * (acc1[2] + pert1[2]);

                    dvel2[0] = dt * (acc2[0] + pert2[0]);
                    dvel2[1] = dt * (acc2[1] + pert2[1]);
                    dvel2[2] = dt * (acc2[2] + pert2[2]);
                }

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

                if constexpr (is_ttl) {
                    // back up gt_kick_inv
                    this->gt_kick_inv_ = gt_inv;

                    if constexpr (has_slowdown) {
                        // integrate gt_drift_inv
                        this->gt_drift_inv_ +=  2.0*dt*kappa_inv*kappa_inv* (vel1[0] * gtgrad1[0] +
                                                                       vel1[1] * gtgrad1[1] +
                                                                       vel1[2] * gtgrad1[2] +
                                                                       vel2[0] * gtgrad2[0] +
                                                                       vel2[1] * gtgrad2[1] +
                                                                       vel2[2] * gtgrad2[2]);
                    }
                    else {
                        // integrate gt_drift_inv
                        this->gt_drift_inv_ +=  2.0*dt* (vel1[0] * gtgrad1[0] +
                                                   vel1[1] * gtgrad1[1] +
                                                   vel1[2] * gtgrad1[2] +
                                                   vel2[0] * gtgrad2[0] +
                                                   vel2[1] * gtgrad2[1] +
                                                   vel2[2] * gtgrad2[2]);
                    }

                }

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

            if constexpr (has_slowdown) {
                // make consistent slowdown inner energy 
                this->etot_sd_ref_ = etot_ref_*kappa_inv;
                this->ekin_sd_ = ekin_*kappa_inv;
                this->epot_sd_ = epot_*kappa_inv;
            }
        }

    private:
        //! Working data of one integrateToTime call, shared by the stepper operations
        struct IntegrateToTimeData {
            Float time_end;        // target Physical time without offset
            Float dt_full;         // full Physical time interval of the call
            Float* backup_data;    // for backup chain data
            int bk_data_size;
            Float* time_table;     // for storing sub-integrated time 
            int n_particle;
            InterruptBinary<Tparticle> bin_interrupt;
            InterruptBinary<Tparticle> bin_interrupt_return;
            bool warning_print_once; // warning print flag
            // last step, for messages
            Float energy_error;
            Float energy_error_bk;
            Float etot_ref_bk;
            Float H;

            Float getEnergyErrorRelAbs() const { return abs((energy_error - energy_error_bk)/etot_ref_bk); }
            // H should be zero initially
            Float getIntegrationErrorCumAbs() const { return abs(H); }
        };

        //! Stepper over the integrator, driven by Regularization step control in integrateToTime
        struct ToTimeStepper {
            TimeTransformedSymplecticIntegrator& integrator;
            IntegrateToTimeData& data;

            Float getTime() const { return integrator.time_; }

            CheckpointAction checkpoint(StepCheckpoint& _checkpoint) { return integrator.checkpointToTime(_checkpoint, data); }

            void save() {
                int bk_return_size = integrator.backupIntData(data.backup_data);
                ASSERT(bk_return_size == data.bk_data_size);
                (void)bk_return_size;
            }

            void restore() {
                int bk_return_size = integrator.restoreIntData(data.backup_data);
                ASSERT(bk_return_size == data.bk_data_size);
                (void)bk_return_size;
                // binary c.m. is not backup, thus recalculate to get correct c.m. velocity for position drift correction due to slowdown inner (the first drift in integrateonestep assume c.m. vel is up to date)
                integrator.updateBinaryCMIter(integrator.info.getBinaryTreeRoot());
            }

            StepOutcome step(const Float _ds) { return integrator.stepToTime(_ds, data); }

            void observe(const StepReport& _report) { integrator.observeToTime(_report, data); }
        };

        //! request a restart of the Regularization step with a new estimate at a checkpoint
        void restartRegularizationStep(StepCheckpoint& _checkpoint, const RestartCause _cause, const Float _ds_estimate) {
#ifdef AR_DEBUG_PRINT
            Float ds_init = _checkpoint.getStepInit();
            Float ds_now = _checkpoint.getStepNow();
#endif
            bool accepted = _checkpoint.requestRestart(_cause, _ds_estimate);
            info.ds = _checkpoint.getPersistentStep();
#ifdef AR_DEBUG_PRINT
            if (accepted) 
                std::cerr<<(_cause==RestartCause::interrupt ? "Change ds after interruption" : "Change ds after update binary orbit")
                         <<": ds(init): "<<ds_init<<" ds(new): "<<_ds_estimate<<" ds(now): "<<ds_now<<std::endl;
#endif
            (void)accepted;
        }

        //! checkpoint of integrateToTime: Interrupt detection and binary update, before the data are backuped
        /*! 
          \return stop when an Interrupt is handed back to the caller, abandon when nothing is left to integrate, otherwise proceed
         */
        CheckpointAction checkpointToTime(StepCheckpoint& _checkpoint, IntegrateToTimeData& _data) {
            const Float _time_end = _data.time_end;
            const int n_particle = _data.n_particle;
            auto& bin_interrupt = _data.bin_interrupt;
            auto& bin_interrupt_return = _data.bin_interrupt_return;
            (void)n_particle;

            bool binary_update_flag=false;
            auto& bin_root = info.getBinaryTreeRoot();
            auto& G = manager->interaction.gravitational_constant;

            // check interrupt condiction, ensure that time end not reach
            if (manager->interrupt_detection_option>0 && !_checkpoint.isSynchronizing()) {
                bin_interrupt.time_now = time_ + info.time_offset;
                bin_interrupt.time_end = _time_end + info.time_offset;
                // calc perturbation energy
                //Float epert=0.0;
                //for (int i=0; i<n_particle; i++) {
                //    epert += force_[i].pot_pert*particles[i].mass;
                //}
                manager->interaction.modifyAndInterruptIter(bin_interrupt, bin_root);
                //InterruptBinary<Tparticle>* bin_intr_ptr = &bin_interrupt;
                //bin_intr_ptr = bin_root.processRootIter(bin_intr_ptr, Tmethod::modifyAndInterruptIter);
                ASSERT(bin_interrupt.checkParams());
                if (bin_interrupt.status!=InterruptStatus::none) {
                    // the mode return back to the root scope
                    if (manager->interrupt_detection_option==2) {
                        return CheckpointAction::stop;
                    }
                    else {

                        // check whether destroy appears (all masses becomes zero)
                        if (bin_interrupt.status==InterruptStatus::destroy) {
                            // all particles become zero masses
#ifdef AR_DEBUG
                            for (int j=0; j<n_particle; j++) {
                                ASSERT(particles[j].mass==0.0);
                            }
#endif
                            de_change_interrupt_ -= etot_ref_;
                            dH_change_interrupt_ -= this->getH();
                            ekin_ = epot_ = etot_ref_ = 0.0;
                            if constexpr (has_slowdown) {
                                this->de_sd_change_cum_ -= this->etot_sd_ref_;
                                this->dH_sd_change_interrupt_ -= this->getHSlowDown();
                                this->ekin_sd_ = this->epot_sd_ = this->etot_sd_ref_ = 0.0;
                            }
#ifdef AR_DEBUG_PRINT
                            std::cerr<<"Interrupt condition triggered! Destroy";
                            std::cerr<<" Time: "<<time_;
                            bin_interrupt.adr->printColumnTitle(std::cerr);
                            std::cerr<<std::endl;
                            bin_interrupt.adr->printColumn(std::cerr);
                            std::cerr<<std::endl;
                            Tparticle::printColumnTitle(std::cerr);
                            std::cerr<<std::endl;
                            for (int j=0; j<2; j++) {
                                bin_interrupt.adr->getMember(j)->printColumn(std::cerr);
                                std::cerr<<std::endl;
                            }
#endif

                            // set binary tree mass to zero
                            this->setBinaryCMZeroIter(bin_root);

                            Float dt = _time_end - time_;
                            time_ += dt;

                            return CheckpointAction::abandon;
                        }

                        Float ekin_bk = ekin_;
                        Float epot_bk = epot_;
                        Float H_bk = this->getH();

                        [[maybe_unused]] Float ekin_sd_bk = 0.0, epot_sd_bk = 0.0, H_sd_bk = 0.0;
                        if constexpr (has_slowdown) {
                            ekin_sd_bk = this->ekin_sd_;
                            epot_sd_bk = this->epot_sd_;
                            H_sd_bk = this->getHSlowDown();
                        }
                        
                        // update binary tree mass
                        info.generateBinaryTree(particles, G);
                        //this->updateBinaryCMIter(bin_root);
                        //this->updateBinarySemiEccPeriodIter(bin_root, G, time_, true);
                        binary_update_flag = true;
                        //bool stable_check=
                        //if (stable_check) bin_root.stableCheckIter(bin_root, 10000*bin_root.period);
                        
                        // should do later, original mass still needed
                        //particles.cm.mass += bin_interrupt.dm;

                        if constexpr (is_ttl) {
                            Float gt_kick_inv_new = this->calcAccPotAndGTKickInv();
                            Float d_gt_kick_inv = gt_kick_inv_new - this->gt_kick_inv_;
                            // when the change is large, initialize this->gt_drift_inv_ to avoid large error
                            if (fabs(d_gt_kick_inv)/std::max(fabs(this->gt_kick_inv_),fabs(gt_kick_inv_new)) >1e-3) 
                                this->gt_drift_inv_ = gt_kick_inv_new;
                            else 
                                this->gt_drift_inv_ += d_gt_kick_inv;
                            this->gt_kick_inv_ = gt_kick_inv_new;
                        }
                        else {
                            this->calcAccPotAndGTKickInv();
                        }
                        // calculate kinetic energy
                        this->calcEKin();

                        // Notice initially etot_ref_ does not include epert. The perturbation effect is accumulated in the integration. Here instance change of mass does not create any work. So no need to add de_pert
                        // get perturbation energy change due to mass change
                        //Float epert_new = 0.0;
                        //for (int i=0; i<n_particle; i++) {
                        //    epert_new += force_[i].pot_pert*particles[i].mass;
                        //}
                        //Float de_pert = epert_new - epert; // notice this is double perturbation potential

                        // get energy change
                        Float de = (ekin_ - ekin_bk) + (epot_ - epot_bk); //+ de_pert;
                        etot_ref_ += de;
                        de_change_interrupt_ += de;
                        dH_change_interrupt_ += this->getH() - H_bk;

                        [[maybe_unused]] Float de_sd = 0.0, dH_sd = 0.0;
                        if constexpr (has_slowdown) {
                            de_sd = (this->ekin_sd_ - ekin_sd_bk) + (this->epot_sd_ - epot_sd_bk);// + de_pert;
                            this->etot_sd_ref_ += de_sd;

                            dH_sd = this->getHSlowDown() - H_sd_bk;

                            // add slowdown change to the global slowdown energy
                            this->de_sd_change_interrupt_ += de_sd;
                            this->dH_sd_change_interrupt_ += dH_sd;
                            this->de_sd_change_cum_ += de_sd;
                            this->dH_sd_change_cum_ += dH_sd;
                        }

#ifdef AR_DEBUG_PRINT
                        std::cerr<<"Interrupt condition triggered!";
                        std::cerr<<" Time: "<<time_;
                        if constexpr (has_slowdown) {
                            std::cerr<<" Energy change: dE_SD: "<<de_sd<<" dH_SD: "<<dH_sd;
                            std::cerr<<" Slowdown: "<<bin_root.slowdown.getSlowDownFactor()<<std::endl;
                        }
                        bin_interrupt.adr->printColumnTitle(std::cerr);
                        std::cerr<<std::endl;
                        bin_interrupt.adr->printColumn(std::cerr);
                        std::cerr<<std::endl;
                        Tparticle::printColumnTitle(std::cerr);
                        std::cerr<<std::endl;
                        for (int j=0; j<2; j++) {
                            bin_interrupt.adr->getMember(j)->printColumn(std::cerr);
                            std::cerr<<std::endl;
                        }
#endif

                        // change fix step option to make safety if energy change is large
                        //info.fix_step_option=FixStepOption::none;
                        
                        // if time_end flag set, reset it to be safety
                        //time_end_flag = false;

                        // check merger case
                        if (bin_interrupt.status==InterruptStatus::merge) {
                            // count particle having mass
                            int count_mass=0;
                            int index_mass_last=-1;
                            for (int j=0; j<n_particle; j++) {
                                if (particles[j].mass>0.0) {
                                    count_mass++;
                                    index_mass_last=j;
                                }
                            }
                            // only one particle has mass, drift directly
                            if (count_mass==1) {
                                ASSERT(index_mass_last<n_particle&&index_mass_last>=0);
                                auto& p = particles[index_mass_last];
                                Float dt = _time_end - time_;
                                p.pos[0] += dt * p.vel[0];
                                p.pos[1] += dt * p.vel[1];
                                p.pos[2] += dt * p.vel[2];

                                time_ += dt;

                                return CheckpointAction::abandon;
                            }
                            // if only two particles have mass, switch off auto ds adjustment
                            if (count_mass==2) {
                                info.fix_step_option=FixStepOption::later;
                                _checkpoint.setStepOption(info.fix_step_option);
                            }
                            //else {
                            //    info.generateBinaryTree(particles, G);
                            //}
                        }

                        if constexpr (has_slowdown) {
                            this->updateSlowDownAndCorrectEnergy(true, true);
                        }

                        restartRegularizationStep(_checkpoint, RestartCause::interrupt, info.calcDsKeplerBinaryTree(*bin_interrupt.adr, manager->step.getOrder(), G, manager->ds_scale));

                        // return one should be the top root
                        if (bin_interrupt_return.status!=InterruptStatus::none) {
                            if (bin_interrupt_return.adr!= bin_interrupt.adr) {
                                // give root address if interrupted binaries are different from previous one
                                bin_interrupt_return.adr = &(info.getBinaryTreeRoot());
                            }
                            if (bin_interrupt.status==InterruptStatus::merge) 
                                bin_interrupt_return.status = InterruptStatus::merge;
                        }
                        else bin_interrupt_return = bin_interrupt;
                    }
                    bin_interrupt.clear();
                }
            }


            // update binary orbit and ds if unstable
            if (!_checkpoint.isSynchronizing()&&!binary_update_flag) {
                bool update_flag=this->updateBinarySemiEccPeriodIter(bin_root, G, time_);

                if constexpr (has_slowdown) {
                    this->updateSlowDownAndCorrectEnergy(true, true);
                }

                if (update_flag) {
            // update slowdown and correct slowdown energy and gt_inv

#ifdef AR_DEBUG_PRINT
                    std::cerr<<"Update binary tree orbits, time= "<<time_<<"\n";
#endif
                    restartRegularizationStep(_checkpoint, RestartCause::binary_update, info.calcDsKeplerBinaryTree(bin_root, manager->step.getOrder(), G, manager->ds_scale));
                }
            }

            // in case the step option is changed above
            _checkpoint.setStepOption(info.fix_step_option);

            return CheckpointAction::proceed;
        }

        //! one step of integrateToTime
        /*! 
          \return the integration error (change of the extended Hamiltonian) and the sub-step time table
         */
        StepOutcome stepToTime(const Float _ds, IntegrateToTimeData& _data) {
            // integrate one step
            if(_data.n_particle==2) integrateTwoOneStep(_ds, _data.time_table);
            else integrateOneStep(_ds, _data.time_table);

            // energy check
            Float energy_error_bk, etot_ref_bk, energy_error, H_bk, H;
            if constexpr (has_slowdown) {
                energy_error_bk = this->getEnergyErrorSlowDownFromBackup(_data.backup_data);
                etot_ref_bk = this->getEtotSlowDownRefFromBackup(_data.backup_data);
                energy_error = this->getEnergyErrorSlowDown();
                H_bk = this->getHSlowDownFromBackup(_data.backup_data);
                H = this->getHSlowDown();
            }
            else {
                energy_error_bk = getEnergyErrorFromBackup(_data.backup_data);
                etot_ref_bk = getEtotRefFromBackup(_data.backup_data);
                energy_error = getEnergyError();
                H_bk = this->getHFromBackup(_data.backup_data);
                H = this->getH();
            }
            _data.energy_error = energy_error;
            _data.energy_error_bk = energy_error_bk;
            _data.etot_ref_bk = etot_ref_bk;
            _data.H = H;

            // get integration error for extended Hamiltonian
            Float integration_error_rel_abs = abs(H-H_bk);

            return StepOutcome{integration_error_rel_abs, _data.time_table};
        }

        //! error message print of integrateToTime for the last step
        void printMessageToTime(const char* message, const StepReport& report, const IntegrateToTimeData& _data) {
            std::cerr<<message<<std::endl;
            std::cerr<<"  T: "<<time_
                     <<"  dT_err/T: "<<(_data.time_end - time_)/_data.dt_full
                     <<"  ds: "<<report.ds
                     <<"  ds_init: "<<report.ds_init
                     <<"  |Int_err/E|: "<<report.integration_error
                     <<"  |Int_err_cum/E|: "<<_data.getIntegrationErrorCumAbs()
                     <<"  |dE/E|: "<<_data.getEnergyErrorRelAbs()
                     <<"  dE_cum: "<<_data.energy_error
                     <<"  Etot_sd: "<<_data.etot_ref_bk
                     <<"  T_end_flag: "<<report.synchronizing
                     <<"  Step_count: "<<report.step_count;
            switch (report.fix_step_option) {
            case FixStepOption::always:
                std::cerr<<"  Fix:  always"<<std::endl;
                break;
            case FixStepOption::later:
                std::cerr<<"  Fix:  later"<<std::endl;
                break;
            case FixStepOption::none:
                std::cerr<<"  Fix:  none"<<std::endl;
                break;
            default:
                break;
            }
        }

#ifdef AR_COLLECT_DS_MODIFY_INFO
        void collectDsModifyInfo(const char* error_message, const StepReport& _report, const IntegrateToTimeData& _data) {
            std::cerr<<error_message<<": "
                     <<"time "<<time_<<" " 
                     <<"ds_new "<<_report.ds_next<<" "
                     <<"ds_init "<<_report.ds_init<<" "
                     <<"modify "<<_report.step_modify_factor<<" "
                     <<"steps "<<_report.step_count<<" "
                     <<"n_mods "<<_report.reduce_count<<" "
                     <<"err "<<_report.integration_error<<" "
                     <<"err/max "<<1.0/_report.integration_error_ratio<<" "
                     <<"errcum/E "<<_data.getIntegrationErrorCumAbs()<<" "
                     <<"dt "<<_report.dt<<" "
                     <<"n_ptcl "<<_data.n_particle<<" ";
            for (int i=0; i<info.binarytree.getSize(); i++) {
                auto& bini = info.binarytree[i];
                std::cerr<<"semi "<<bini.semi<<" "
                         <<"ecc "<<bini.ecc<<" "
                         <<"period "<<bini.period<<" "
                         <<"m1 "<<bini.m1<<" "
                         <<"m2 "<<bini.m2<<" "
                         <<"stab "<<bini.stab<<" "
                         <<"sd "<<bini.slowdown.getSlowDownFactor()<<" "
                         <<"sd_org "<<bini.slowdown.getSlowDownFactorOrigin()<<" "
                         <<"pert_in "<<bini.slowdown.getPertIn()<<" "
                         <<"pert_out "<<bini.slowdown.getPertOut()<<" ";
            }
            std::cerr<<std::endl;
        }
#endif 

        //! per-step report of integrateToTime: all warning and debug printing
        void observeToTime(const StepReport& _report, IntegrateToTimeData& _data) {
            switch (_report.event) {
            case StepEvent::step_taken:
//#ifdef AR_WARN
                // warning for large number of steps
                if(_data.warning_print_once&&_report.step_count>=manager->step_count_max) {
                    if(_report.step_count%manager->step_count_max==0) {
                        printMessageToTime("Warning: step count is signficiant large", _report, _data);
                        for (int i=0; i<info.binarytree.getSize(); i++){
                            auto& bin = info.binarytree[i];
                            std::cerr<<"  Binary["<<i<<"]: "
                                     <<"  i1="<<bin.getMemberIndex(0)
                                     <<"  i2="<<bin.getMemberIndex(1)
                                     <<"  m1="<<bin.m1
                                     <<"  m2="<<bin.m2
                                     <<"  semi= "<<bin.semi
                                     <<"  ecc= "<<bin.ecc
                                     <<"  period= "<<bin.period
                                     <<"  stab= "<<bin.stab
                                     <<"  SD= "<<bin.slowdown.getSlowDownFactor()
                                     <<"  SD_org= "<<bin.slowdown.getSlowDownFactorOrigin()
                                     <<"  Tscale= "<<bin.slowdown.timescale
                                     <<"  pert_in= "<<bin.slowdown.pert_in
                                     <<"  pert_out= "<<bin.slowdown.pert_out;
                            std::cerr<<std::endl;
                            _data.warning_print_once = false;
                        }
#ifdef AR_DEBUG_DUMP
                        if (!info.dump_flag) {
                            DATADUMP("dump_large_step");
                            info.dump_flag=true;
                        }
#endif
                    }
                }
//#endif

#ifdef AR_DEEP_DEBUG
                printMessageToTime("", _report, _data);
                std::cerr<<"Timetable: ";
                for (int i=0; i<manager->step.getCDPairSize(); i++) std::cerr<<" "<<_data.time_table[manager->step.getSortCumSumCKIndex(i)];
                std::cerr<<std::endl;
#endif
                break;
            case StepEvent::large_energy_error:
#ifdef AR_COLLECT_DS_MODIFY_INFO
                collectDsModifyInfo("Large_energy_error", _report, _data);
#endif
                break;
            case StepEvent::negative_step:
#ifdef AR_COLLECT_DS_MODIFY_INFO
                collectDsModifyInfo("Negative_step", _report, _data);
#endif
                break;
            case StepEvent::reuse_backup_step:
#ifdef AR_COLLECT_DS_MODIFY_INFO
                collectDsModifyInfo("Reuse_backup_ds", _report, _data);
#endif
                break;
            case StepEvent::increase_step:
                // the persistent Regularization step follows the growth
                info.ds = _report.ds_next;
#ifdef AR_DEBUG_PRINT
                std::cerr<<"Energy error is small enough for increase step, integration_error_rel_abs="<<_report.integration_error
                         <<" energy_error_rel_max="<<manager->energy_error_relative_max<<" step_modify_factor="<<_report.step_modify_factor<<" new ds="<<_report.ds_next<<std::endl;
#endif
                break;
            case StepEvent::sync_enlarge_step:
#ifdef AR_DEEP_DEBUG
                std::cerr<<"Time step dt(real) "<<_report.dt<<" <0.3*(time_end-time)(real) "<<_data.time_end - time_<<" enlarge step factor: "<<_report.step_modify_factor<<" new ds: "<<_report.ds_next<<std::endl;
#endif
                break;
            case StepEvent::sync_overshoot_first:
#ifdef AR_DEEP_DEBUG
                std::cerr<<"Time_end reach, time[k]= "<<_report.time_next<<" time= "<<time_<<" time_end/time[k]="<<_data.time_end/_report.time_next<<" CumSum_CK="<<_report.cck<<" ds(next) = "<<_report.ds<<" ds(next_next) = "<<_report.ds_next<<"\n";
#endif
                break;
            case StepEvent::sync_overshoot_between:
#ifdef AR_DEEP_DEBUG
                std::cerr<<"Time_end reach, time_prev= "<<_report.time_prev<<" time[k]= "<<_report.time_next<<" time= "<<time_<<" (time_end-time_prev)/dt="<<(_data.time_end-_report.time_prev)/_report.dt<<" CumSum_CK="<<_report.cck<<" CumSum_CK(prev)="<<_report.cck_prev<<" ds(next) = "<<_report.ds<<" ds(next_next) = "<<_report.ds_next<<" \n";
#endif
                break;
            case StepEvent::synchronization_failed:
                // When time sychronization steps too large, the integration is aborted
                printMessageToTime("Error! step count after time synchronization is too large", _report, _data);
                printColumnTitle(std::cerr,20,info.binarytree.getSize());
                std::cerr<<std::endl;
                printColumn(std::cerr,20,info.binarytree.getSize());
                std::cerr<<std::endl;
                break;
            case StepEvent::finish:
#ifdef AR_DEEP_DEBUG
                std::cerr<<"Finish, time_diff_rel = "<<(_data.time_end - time_)/_data.dt_full<<" integration_error_rel_abs = "<<_report.integration_error<<std::endl;
#endif
                break;
            }
        }

    public:
        // Integrate the system to a given time
        /*! The Regularization step of every step is chosen by RegularizationStepControl, which drives the integration through ToTimeStepper
          @param[in] _time_end: the expected finishing time without offset
          \return binary tree of the pair which triggers interruption condition
         */
        InterruptBinary<Tparticle> integrateToTime(const Float _time_end) {
            ASSERT(this->checkParams());

            // backup data size
            const int bk_data_size = getBackupDataSize();
            
            Float backup_data[bk_data_size]; // for backup chain data
#ifdef AR_DEBUG_DUMP
            Float backup_data_init[bk_data_size]; // for backup initial data
#endif

            // time table
            const int cd_pair_size = manager->step.getCDPairSize();
            Float time_table[cd_pair_size]; // for storing sub-integrated time 

            // Regularization step control, starting from the persistent step
            const StepControlLimits step_limits = {manager->time_error_max, manager->energy_error_relative_max, manager->step_count_max};
            RegularizationStepControl step_control(manager->step, step_limits, info.ds, info.fix_step_option);

            IntegrateToTimeData data;
            data.time_end = _time_end;
            data.dt_full = _time_end - time_; // real full time step
            data.backup_data = backup_data;
            data.bk_data_size = bk_data_size;
            data.time_table = time_table;
            data.bin_interrupt.time_now=time_ + info.time_offset;
            data.bin_interrupt.time_end=_time_end + info.time_offset;
            data.bin_interrupt_return = data.bin_interrupt;
            data.warning_print_once = true;
            
            // particle data
            const int n_particle = particles.getSize();
            data.n_particle = n_particle;

/* This must suppress since after findslowdowninner, slowdown inner is reset to 1.0, recalculate ekin_sdi give completely wrong value for energy correction for slowdown change later
#ifdef AR_DEBUG
            Float ekin_check = ekin_;
            this->calcEKin();
            ASSERT(abs(ekin_check-ekin_)<1e-10);
            ekin_ = ekin_check;
#endif
*/
#ifdef AR_DEBUG_DUMP
            // back up initial data
            backupIntData(backup_data_init);
#endif
      
            if constexpr (is_slowdown_array) {
                // find new inner slowdown binaries, the binary tree data may be modified, thus it is safer to recheck slowdown inner binary at beginning to avoid memory issue (bin is pointer).
                if constexpr (is_ttl) {
                    if (n_particle >2) {
                        int nold = this->binary_slowdown.getSize();
                        this->findSlowDownInner(time_);
                        int nnew = this->binary_slowdown.getSize();
                        // in case slowdown is disabled in the next step, this->gt_drift_inv_ should be re-initialized
                        if (nold>0&&nnew==0) {
                            this->gt_kick_inv_ = manager->interaction.calcAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, time_);
                            this->gt_drift_inv_ = this->gt_kick_inv_;
                        }
                    }
                }
                else {
                    if (n_particle >2) this->findSlowDownInner(time_);
                }
            }

            if constexpr (is_slowdown_tree) {
                // update slowdown and correct slowdown energy and gt_inv
                this->updateSlowDownAndCorrectEnergy(true, true);
            }


            // reset binary stab_check_time
            for (int i=0; i<info.binarytree.getSize(); i++)
                info.binarytree[i].stab_check_time = time_;


            ToTimeStepper stepper = {*this, data};
            StepControlStatus status = step_control.advanceTo(stepper, _time_end);

            // the persistent Regularization step and the step option return to the group information on every exit
            info.ds = step_control.getPersistentStep();
            info.fix_step_option = step_control.getStepOption();

            // When time sychronization steps too large, abort
            if (status==StepControlStatus::failed_to_synchronize) {
#ifdef AR_DEBUG_DUMP
                if (!info.dump_flag) {
                    DATADUMP("dump_large_step");
                    info.dump_flag=true;
                }
#endif
                abort();
            }

            // cumulative step count 
            profile.step_count = step_control.getStepCount();
            profile.step_count_tsyn = step_control.getStepCountSync();
            profile.step_count_sum += step_control.getStepCount();
            profile.step_count_tsyn_sum += step_control.getStepCountSync();

            // an Interrupt that ends the integration early is returned as it is
            if (status!=StepControlStatus::reached) return data.bin_interrupt;

            return data.bin_interrupt_return;
        }


    public:
        //! correct CM drift
        /*! calculate c.m. and correct the member data to the c.m. frame.
          This is used after the perturbation, in case the c.m. drift when members are in c.m. frame
         */
        void correctCenterOfMassDrift() {
            ASSERT(!particles.isOriginFrame());
            Float mcm=0.0, pos_cm[3]={0.0,0.0,0.0}, vel_cm[3]={0.0,0.0,0.0};
            auto* particle_data= particles.getDataAddress();
            for (int i=0; i<particles.getSize(); i++) {
                const Float *ri = particle_data[i].pos;
                const Float *vi = particle_data[i].getVel();
                const Float mi  = particle_data[i].mass;

                pos_cm[0] += ri[0] * mi;
                pos_cm[1] += ri[1] * mi;
                pos_cm[2] += ri[2] * mi;

                vel_cm[0] += vi[0] * mi;
                vel_cm[1] += vi[1] * mi;
                vel_cm[2] += vi[2] * mi;
                mcm += mi;
            }
            pos_cm[0] /= mcm; 
            pos_cm[1] /= mcm; 
            pos_cm[2] /= mcm; 
            vel_cm[0] /= mcm; 
            vel_cm[1] /= mcm; 
            vel_cm[2] /= mcm;

            for (int i=0; i<particles.getSize(); i++) {
                Float *ri = particle_data[i].pos;
                Float *vi = particle_data[i].getVel();

                ri[0] -= pos_cm[0]; 
                ri[1] -= pos_cm[1]; 
                ri[2] -= pos_cm[2]; 
                vi[0] -= vel_cm[0]; 
                vi[1] -= vel_cm[1]; 
                vi[2] -= vel_cm[2]; 
            }
        }

        //! write back particles to original address
        /*! If particles are in center-off-mass frame, write back the particle in original frame but not modify local copies to avoid roundoff error
         */
        template <class Tptcl>
        void writeBackParticlesOriginFrame() {
            ASSERT(particles.getMode()==COMM::ListMode::copy);
            auto* particle_adr = particles.getOriginAddressArray();
            auto* particle_data= particles.getDataAddress();
            if (particles.isOriginFrame()) {
                for (int i=0; i<particles.getSize(); i++) {
                    *(Tptcl*)particle_adr[i] = particle_data[i];
                }
            }
            else {
                for (int i=0; i<particles.getSize(); i++) {
                    Tptcl pc = particle_data[i];

                    pc.pos[0] = particle_data[i].pos[0] + particles.cm.pos[0];
                    pc.pos[1] = particle_data[i].pos[1] + particles.cm.pos[1];
                    pc.pos[2] = particle_data[i].pos[2] + particles.cm.pos[2];

                    pc.vel[0] = particle_data[i].vel[0] + particles.cm.vel[0];
                    pc.vel[1] = particle_data[i].vel[1] + particles.cm.vel[1];
                    pc.vel[2] = particle_data[i].vel[2] + particles.cm.vel[2];
                    
                    *(Tptcl*)particle_adr[i] = pc;
                }
            }
        }



        //! Get current kinetic energy
        /*! \return current kinetic energy
         */
        Float getEkin() const {
            return ekin_;
        }

        //! Get current potential energy
        /*! \return current potetnial energy (negative value for bounded systems)
         */
        Float getEpot() const {
            return epot_;
        }

        //! Get current total integrated energy 
        /*! \return total integrated energy 
         */
        Float getEtotRef() const {
            return etot_ref_;
        }

        //! Get current total energy from ekin and epot
        /*! \return total integrated energy 
         */
        Float getEtot() const {
            return ekin_ + epot_;
        }

        //! get perturbation potential energy
        /*! \return perturbation potential energy
         */
        Float getEpert() const {
            Float epert=0.0;
            int n_particle = particles.getSize();
            for (int i=0; i<n_particle; i++) {
                epert += force_[i].pot_pert*particles[i].mass;
            }
            return epert;
        }

        //! get energy error 
        /*! \return energy error
         */
        Float getEnergyError() const {
            return ekin_ + epot_ - etot_ref_;
        }

        //! get energy error from backup data
        Float getEnergyErrorFromBackup(Float* _bk) const {
            return -_bk[1] + _bk[2] + _bk[3];
        }

        //! get integrated energy from backup data
        Float getEtotRefFromBackup(Float* _bk) const {
            return _bk[1];
        }

        //! get total energy from backup data (ekin+epot)
        Float getEtotFromBackup(Float* _bk) const {
            return _bk[2] + _bk[3];
        }

        //! reset cumulative energy/hamiltonian change due to interruption
        void resetDEChangeBinaryInterrupt() {
            de_change_interrupt_ = 0.0;
            dH_change_interrupt_ = 0.0;
        }

        //! get cumulative energy change due to interruption
        Float getDEChangeBinaryInterrupt() const {
            return de_change_interrupt_;
        }

        //! get cumulative hamiltonian change due to interruption
        Float getDHChangeBinaryInterrupt() const {
            return dH_change_interrupt_;
        }

        //! get backup data size
        int getBackupDataSize() const {
            int bk_size = 6;
            if constexpr (has_slowdown) {
                bk_size += 7; 
                //bk_size += SlowDown::getBackupDataSize();
            }
            if constexpr (is_ttl) {
                bk_size += 2;
            }
            bk_size += particles.getBackupDataSize();
            return bk_size;
        }


        //! Backup integration data 
        /*! Backup $time_, $etot_, $ekin_, $epot_, $gt_drift_, $this->gt_kick_inv_, #particles, $slowdown to one Float data array
          \return backup array size
        */
        int backupIntData(Float* _bk) {
            int bk_size=0;
            _bk[bk_size++] = time_;       //0
            _bk[bk_size++] = etot_ref_;   //1
            _bk[bk_size++] = ekin_;       //2
            _bk[bk_size++] = epot_;       //3 
            _bk[bk_size++] = de_change_interrupt_;     //4
            _bk[bk_size++] = dH_change_interrupt_;     //5
            if constexpr (has_slowdown) {
                _bk[bk_size++] = this->etot_sd_ref_; //6
                _bk[bk_size++] = this->ekin_sd_;     //7
                _bk[bk_size++] = this->epot_sd_;     //8
                _bk[bk_size++] = this->de_sd_change_cum_; //9
                _bk[bk_size++] = this->dH_sd_change_cum_; //10
                _bk[bk_size++] = this->de_sd_change_interrupt_; //11
                _bk[bk_size++] = this->dH_sd_change_interrupt_; //12
            }

            if constexpr (is_ttl) {
                _bk[bk_size++] = this->gt_drift_inv_;  //13 / 6
                _bk[bk_size++] = this->gt_kick_inv_;   //14 / 7
            }

            bk_size += particles.backupParticlePosVel(&_bk[bk_size]); 
//            bk_size += info.getBinaryTreeRoot().slowdown.backup(&_bk[bk_size]); // slowdownfactor
            return bk_size;
        }

        //! Restore integration data
        /*! restore $time_, $etot_, $ekin_, $epot_, $gt_drift_, $this->gt_kick_inv_, #particles, $slowdown from one Float data array
          \return backup array size
        */
        int restoreIntData(Float* _bk) {
            int bk_size = 0;
            time_     = _bk[bk_size++];
            etot_ref_ = _bk[bk_size++];
            ekin_     = _bk[bk_size++];
            epot_     = _bk[bk_size++];
            de_change_interrupt_= _bk[bk_size++];
            dH_change_interrupt_= _bk[bk_size++];

            if constexpr (has_slowdown) {
                this->etot_sd_ref_ = _bk[bk_size++];
                this->ekin_sd_     = _bk[bk_size++];
                this->epot_sd_     = _bk[bk_size++];

                this->de_sd_change_cum_= _bk[bk_size++];
                this->dH_sd_change_cum_= _bk[bk_size++];
                this->de_sd_change_interrupt_= _bk[bk_size++];
                this->dH_sd_change_interrupt_= _bk[bk_size++];
            }
            if constexpr (is_ttl) {
                this->gt_drift_inv_  = _bk[bk_size++];
                this->gt_kick_inv_   = _bk[bk_size++];
            }
            bk_size += particles.restoreParticlePosVel(&_bk[bk_size]);
//            bk_size += info.getBinaryTreeRoot().slowdown.restore(&_bk[bk_size]);
            return bk_size;
        }


        //! print group information 
        /*! Message, Number of members, time, binary tree printing interation
          @param[in] _type: 0: new group (if pair id is same, no printing); 1: end group (always print and reset pair id)
          @param[in] _fout: FILE IO
          @param[in] _width: print width
          @param[in] _pcm: center of mass particle to calculate origin position and velocity, if NULL, assume cm pos and vel are zero
        */
        template<class Tptcl>
        void printGroupInfo(const int _type, std::ostream& _fout, const int _width, const Tptcl* _pcm=NULL) {
            auto& bin_root = info.getBinaryTreeRoot();
            //auto* p1 = bin_root.getLeftMember();
            //auto* p2 = bin_root.getRightMember();
            
            bool reset_flag = (_type==1 && bin_root.semi<0 && bin_root.ecca>0);

            if (info.checkAndSetBinaryPairIDIter(bin_root, reset_flag)) {
                if (_type==0) return; // if it is new but already existed binary, do not print
                else if (!reset_flag) return; // in the end case, if the system is still bound, do not print 
            }

            Float pos_cm[3], vel_cm[3];
            auto& pcm_loc = particles.cm;
            if (_pcm!=NULL) {
                pos_cm[0] = pcm_loc.pos[0] + _pcm->pos[0];
                pos_cm[1] = pcm_loc.pos[1] + _pcm->pos[1];
                pos_cm[2] = pcm_loc.pos[2] + _pcm->pos[2];
                vel_cm[0] = pcm_loc.vel[0] + _pcm->vel[0];
                vel_cm[1] = pcm_loc.vel[1] + _pcm->vel[1];
                vel_cm[2] = pcm_loc.vel[2] + _pcm->vel[2];
            }
            else {
                pos_cm[0] = pcm_loc.pos[0]; 
                pos_cm[1] = pcm_loc.pos[1]; 
                pos_cm[2] = pcm_loc.pos[2]; 
                vel_cm[0] = pcm_loc.vel[0]; 
                vel_cm[1] = pcm_loc.vel[1]; 
                vel_cm[2] = pcm_loc.vel[2]; 
            }
#pragma omp critical
            {
                _fout<<std::setw(_width)<<_type
                     <<std::setw(_width)<<bin_root.getMemberN()
                     <<std::setw(_width)<<time_ + info.time_offset;
                _fout<<std::setw(_width)<<pos_cm[0]
                     <<std::setw(_width)<<pos_cm[1]
                     <<std::setw(_width)<<pos_cm[2]
                     <<std::setw(_width)<<vel_cm[0]
                     <<std::setw(_width)<<vel_cm[1]
                     <<std::setw(_width)<<vel_cm[2];
                bin_root.printBinaryTreeIter(_fout, _width);
                _fout<<std::endl;
            }
            //if (_type==0) { // register pair id to avoid repeating printing
            //    p1->setBinaryPairID(p2->id);
            //    p2->setBinaryPairID(p1->id);
            //}
            //else { // break case reset pair id
            //    p1->setBinaryPairID(0);
            //    p2->setBinaryPairID(0);
            //}
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width 
          @param[in] _n_sd: slowdown inner group
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20, const int _n_sd=0) {
            _fout<<std::setw(_width)<<"Time"
                 <<std::setw(_width)<<"dE"
                 <<std::setw(_width)<<"Etot"
                 <<std::setw(_width)<<"Ekin"
                 <<std::setw(_width)<<"Epot"
                 <<std::setw(_width)<<"Gt_drift"
                 <<std::setw(_width)<<"H"
                 <<std::setw(_width)<<"dE_intr"
                 <<std::setw(_width)<<"dH_intr";
            perturber.printColumnTitle(_fout, _width);
            info.printColumnTitle(_fout, _width);
            profile.printColumnTitle(_fout, _width);
            if constexpr (has_slowdown) {
                _fout<<std::setw(_width)<<"dE_SD" 
                     <<std::setw(_width)<<"Etot_SD" 
                     <<std::setw(_width)<<"Ekin_SD" 
                     <<std::setw(_width)<<"Epot_SD"
                     <<std::setw(_width)<<"dE_SDC_cum" 
                     <<std::setw(_width)<<"dH_SDC_cum"
                     <<std::setw(_width)<<"dE_SDC_intr" 
                     <<std::setw(_width)<<"dH_SDC_intr"; 
                _fout<<std::setw(_width)<<"N_SD";
                for (int i=0; i<_n_sd; i++) {
                    _fout<<std::setw(_width)<<"I1"
                         <<std::setw(_width)<<"I2";
                    SlowDown::printColumnTitle(_fout, _width);
                }
            }
            particles.printColumnTitle(_fout, _width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width 
          @param[in] _n_sd: slowdown inner group
        */
        void printColumn(std::ostream & _fout, const int _width=20, const int _n_sd=0){
            // inverse time transformation factor for drift and Hamiltonian, with Slow-down if it is used
            Float gt_drift_inv, H;
            if constexpr (is_ttl) {
                gt_drift_inv = this->gt_drift_inv_;
            }
            else {
                if constexpr (has_slowdown) {
                    gt_drift_inv = manager->interaction.calcGTDriftInv(this->ekin_sd_-this->etot_sd_ref_);
                }
                else {
                    gt_drift_inv = manager->interaction.calcGTDriftInv(ekin_-etot_ref_);
                }
            }
            if constexpr (has_slowdown) {
                H = this->getHSlowDown();
            }
            else {
                H = this->getH();
            }

            _fout<<std::setw(_width)<<this->getTime()
                 <<std::setw(_width)<<getEnergyError()
                 <<std::setw(_width)<<etot_ref_
                 <<std::setw(_width)<<ekin_
                 <<std::setw(_width)<<epot_
                 <<std::setw(_width)<<1.0/gt_drift_inv
                 <<std::setw(_width)<<H
                 <<std::setw(_width)<<de_change_interrupt_
                 <<std::setw(_width)<<dH_change_interrupt_;
            perturber.printColumn(_fout, _width);
            info.printColumn(_fout, _width);
            profile.printColumn(_fout, _width);
            if constexpr (has_slowdown) {
                _fout<<std::setw(_width)<<this->getEnergyErrorSlowDown()
                     <<std::setw(_width)<<this->etot_sd_ref_ 
                     <<std::setw(_width)<<this->ekin_sd_ 
                     <<std::setw(_width)<<this->epot_sd_
                     <<std::setw(_width)<<this->de_sd_change_cum_
                     <<std::setw(_width)<<this->dH_sd_change_cum_
                     <<std::setw(_width)<<this->de_sd_change_interrupt_
                     <<std::setw(_width)<<this->dH_sd_change_interrupt_;
                SlowDown sd_empty;
                if constexpr (is_slowdown_array) {
                    int n_sd_now = this->binary_slowdown.getSize();
                    _fout<<std::setw(_width)<<n_sd_now;
                    for (int i=0; i<_n_sd; i++) {
                        if (i<n_sd_now) {
                            _fout<<std::setw(_width)<<this->binary_slowdown[i]->getMemberIndex(0)
                                 <<std::setw(_width)<<this->binary_slowdown[i]->getMemberIndex(1);
                            this->binary_slowdown[i]->slowdown.printColumn(_fout, _width);
                        }
                        else {
                            _fout<<std::setw(_width)<<-1
                                 <<std::setw(_width)<<-1;
                            sd_empty.printColumn(_fout, _width);
                        }
                    }
                }
                else {
                    int n_sd_now = info.binarytree.getSize();
                    _fout<<std::setw(_width)<<n_sd_now;
                    for (int i=0; i<_n_sd; i++) {
                        if (i<n_sd_now) {
                            _fout<<std::setw(_width)<<info.binarytree[i].getMemberIndex(0)
                                 <<std::setw(_width)<<info.binarytree[i].getMemberIndex(1);
                            info.binarytree[i].slowdown.printColumn(_fout, _width);
                        }
                        else {
                            _fout<<std::setw(_width)<<-1
                                 <<std::setw(_width)<<-1;
                            sd_empty.printColumn(_fout, _width);
                        }
                    }
                }
            }
            particles.printColumn(_fout, _width);
        }

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBinary(FILE *_fout) {
            fwrite(&time_, sizeof(Float), 1, _fout);
            fwrite(&etot_ref_, sizeof(Float), 1, _fout);
            fwrite(&ekin_, sizeof(Float), 1, _fout);
            fwrite(&epot_, sizeof(Float), 1, _fout);
            if constexpr (is_ttl) {
                fwrite(&this->gt_drift_inv_, sizeof(Float), 1, _fout);
            }
            int size = force_.getSize();
            fwrite(&size, sizeof(int), 1, _fout);
            for (int i=0; i<size; i++) force_[i].writeBinary(_fout);
            
            particles.writeBinary(_fout);
            perturber.writeBinary(_fout);
            info.writeBinary(_fout);
            profile.writeBinary(_fout);
        }

        //! read class data with BINARY format and initial the array
        /*! @param[in] _fin: file IO for read
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(&time_, sizeof(Float), 1, _fin);
            rcount += fread(&etot_ref_, sizeof(Float), 1, _fin);
            rcount += fread(&ekin_, sizeof(Float), 1, _fin);
            rcount += fread(&epot_, sizeof(Float), 1, _fin);
            if (rcount<4) {
                std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
                abort();
            }
            if constexpr (is_ttl) {
                rcount = fread(&this->gt_drift_inv_, sizeof(Float), 1, _fin);
                if (rcount<1) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                    abort();
                }
            }
            int size;
            rcount = fread(&size, sizeof(int),1, _fin);
            if(rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            if(size<0) {
                std::cerr<<"Error: array size <0 "<<size<<"<=0!\n";
                abort();
            }
            if (size>0) {
                force_.setMode(COMM::ListMode::local);
                force_.reserveMem(size);
                force_.resizeNoInitialize(size);
                for (int i=0; i<size; i++) force_[i].readBinary(_fin);
            }
            
            particles.setMode(COMM::ListMode::local);
            particles.readBinary(_fin);
            perturber.readBinary(_fin);
            info.readBinary(_fin);
            profile.readBinary(_fin);
        }

    };
}
