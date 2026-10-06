#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include "Common/Float.h"
#include "AR/fix_step_option.h"
#include "AR/symplectic_step.h"

namespace AR {

    //! Numeric limits read by Regularization step control
    struct StepControlLimits {
        Float time_error_max;                  ///> time synchronization (absolute) error limit
        Float energy_error_relative_max;       ///> integration error limit of one step
        long long unsigned int step_count_max; ///> maximum step counts
    };

    //! How an advance to a target Physical time ended
    /*! reached: the target Physical time is reached within the time error limit \n
        stopped: the stepper asked to stop at a checkpoint before the target is reached \n
        abandoned: the stepper reported at a checkpoint that nothing is left to integrate \n
        failed_to_synchronize: the time synchronization steps exceed the maximum step count
     */
    enum class StepControlStatus {reached, stopped, abandoned, failed_to_synchronize};

    //! What the stepper wants after a checkpoint
    enum class CheckpointAction {proceed, stop, abandon};

    //! Why a restart of the Regularization step is requested at a checkpoint
    enum class RestartCause {interrupt, binary_update};

    //! Outcome of one step reported by the stepper
    struct StepOutcome {
        Float integration_error; ///> integration error of the step, compared with the energy error limit
        const Float* time_table; ///> Physical time after each sub-step (storage order of the symplectic step table), read only when the step overshoots the target
    };

    //! Reason code of a per-step report
    /*! step_taken: a step is finished, before any decision \n
        large_energy_error: the step is rejected due to its integration error; a reduced step is retried \n
        negative_step: the step is rejected since Physical time moved backwards; a reduced step is retried \n
        reuse_backup_step: a previously reduced Regularization step is recovered \n
        increase_step: the Regularization step grows, the persistent one follows \n
        sync_enlarge_step: during time synchronization the step is too small to reach the target and is enlarged \n
        sync_overshoot_first: the target is passed inside the first sub-step; the step is rejected \n
        sync_overshoot_between: the target is passed between two sub-steps; the step is rejected \n
        finish: the target is reached
     */
    enum class StepEvent {step_taken, large_energy_error, negative_step, reuse_backup_step, increase_step, sync_enlarge_step, sync_overshoot_first, sync_overshoot_between, finish};

    //! Per-step report handed to the stepper, which does all printing
    struct StepReport {
        StepEvent event;
        Float ds;                 ///> Regularization step of the step just taken (after modification when the step is rejected)
        Float ds_next;            ///> Regularization step in the other slot, used after the next one
        Float ds_init;            ///> reference Regularization step of this advance
        Float dt;                 ///> Physical time advance of the step just taken
        Float integration_error;
        Float integration_error_ratio; ///> energy error limit / integration error
        Float step_modify_factor; 
        int reduce_count;         ///> number of reductions directly after a recovery
        bool synchronizing;       ///> whether the target was passed once (time synchronization phase)
        long long unsigned int step_count;
        long long unsigned int step_count_tsyn;
        FixStepOption fix_step_option;
        Float time_prev;          ///> sync_overshoot_between: sub-step time before the target
        Float time_next;          ///> sync_overshoot_*: first sub-step time after the target
        Float cck_prev;           ///> sync_overshoot_between: cumulative coefficient of the sub-step before the target
        Float cck;                ///> sync_overshoot_*: cumulative coefficient of the first sub-step after the target
    };

    class RegularizationStepControl;

    //! Handle given to the stepper during a checkpoint, the only place where a restart is possible
    class StepCheckpoint {
    private:
        RegularizationStepControl& control_;
        friend class RegularizationStepControl;
        StepCheckpoint(RegularizationStepControl& _control): control_(_control) {}
    public:
        //! whether the target was passed once, i.e. time synchronization is going on
        bool isSynchronizing() const;

        //! request a restart with a new estimate of the Regularization step, return true if the estimate is accepted
        bool requestRestart(const RestartCause _cause, const Float _ds_estimate);

        //! change the step option 
        void setStepOption(const FixStepOption _option);

        //! reference Regularization step of this advance
        Float getStepInit() const;

        //! Regularization step in the first slot (for messages)
        Float getStepNow() const;
    };

    //! What Regularization step control needs from the integration it drives
    /*! getTime: current Physical time \n
        checkpoint: called before state is saved; may request a restart or a new step option through the handle, and may ask to stop \n
        save / restore: backup the state, or return to the last backup \n
        step: advance one step of the given Regularization step \n
        observe: receive a per-step report
     */
    template <class Tstepper>
    concept RegularizationStepper = requires(Tstepper _stepper, const Float _ds, const StepReport& _report, StepCheckpoint& _checkpoint) {
        { _stepper.getTime() } -> std::convertible_to<Float>;
        { _stepper.checkpoint(_checkpoint) } -> std::same_as<CheckpointAction>;
        { _stepper.save() };
        { _stepper.restore() };
        { _stepper.step(_ds) } -> std::same_as<StepOutcome>;
        { _stepper.observe(_report) };
    };

    //! Regularization step control
    /*! Chooses the Regularization step of every AR step while advancing to a target Physical time:
        error-driven reduction and recovery, growth, and time synchronization at the target.
        It drives the integration through a stepper (see RegularizationStepper) and performs no IO.
        The behaviour is the one inherited from the integration loop of TimeTransformedSymplecticIntegrator; inherited quirks are marked.
     */
    class RegularizationStepControl {
    private:
        //! reduce level number limit
        static const int n_reduce_level_max=10;

        //! Store of reduced Regularization steps waiting for recovery
        struct ReduceLevelStore {
            Float ds_backup[n_reduce_level_max+1];
            int n_step_wait_recover_ds[n_reduce_level_max+1];
            int n_reduce_level; 

            ReduceLevelStore(const Float _ds): ds_backup{}, n_step_wait_recover_ds{} { initial(_ds), n_reduce_level = -1; }

            // Inherited: only ten of the eleven slots are initialized and the level counter is not reset
            void initial(const Float _ds) {
                for (int i=0; i<n_reduce_level_max; i++) {
                    ds_backup[i] = _ds;
                    n_step_wait_recover_ds[i] = 0;
                }
            }

            // record ds to backup
            void backup(const Float _ds, const Float _modify_factor) {
                if (n_reduce_level==n_reduce_level_max) 
                    n_step_wait_recover_ds[n_reduce_level] += 2*to_int(1.0/_modify_factor);
                else {
                    n_reduce_level++;
                    ds_backup[n_reduce_level] = _ds;
                    n_step_wait_recover_ds[n_reduce_level] = 2*to_int(1.0/_modify_factor);
                }
            }

            // count step (return false) and recover ds if necessary (return true)
            bool countAndRecover(Float &_ds, Float &_modify_factor, const bool _recover_flag) {
                if (n_reduce_level>=0) {
                    if (n_step_wait_recover_ds[n_reduce_level] ==0) {
                        if (_recover_flag) {
                            _modify_factor = ds_backup[n_reduce_level]/_ds;
                            _ds = ds_backup[n_reduce_level];
                            n_step_wait_recover_ds[n_reduce_level] = -1;
                            n_reduce_level--;
                            return true;
                        }
                    }
                    else {
                        n_step_wait_recover_ds[n_reduce_level]--;
                        return false;
                    }
                }
                return false;
            }
        };

        const SymplecticStep& step_;
        StepControlLimits limits_;
        Float ds_persistent_;           // persistent Regularization step, seen by callers between advances
        FixStepOption fix_step_option_;

        // state of one advance
        Float ds_[2];                   // two switch step control: step with a buffer
        Float ds_init_;                 // reference step
        int   ds_switch_;               // 0 or 1
        ReduceLevelStore ds_backup_;
        bool time_end_flag_;            // indicate whether time reach the end
        long long unsigned int step_count_;      // integration step 
        long long unsigned int step_count_tsyn_; // time synchronization step

        //! regular block time step modification factor
        static Float regularStepFactor(const Float _fac) {
            Float fac = 1.0;
            if (_fac<1) while (fac>_fac) fac *= 0.5;
            else {
                while (fac<=_fac) fac *= 2.0;
                fac *= 0.5;
            }
            return fac;
        }

        //! use the new estimate as reference step
        void acceptRestart() {
            ASSERT(ds_persistent_>0);
            ds_[0] = std::min(ds_[0], ds_persistent_);
            ds_[1] = std::min(ds_[1], ds_persistent_);
            ds_backup_.initial(ds_persistent_);
            ds_init_ = ds_persistent_;
        }

        //! restart with a new estimate of the Regularization step, return true if the estimate is accepted
        bool restart(const RestartCause _cause, const Float _ds_estimate) {
            ds_persistent_ = _ds_estimate;
            if (_cause==RestartCause::interrupt) {
                Float ds_max = step_.calcStepModifyFactorFromErrorRatio(2.0)*ds_init_;
                Float ds_min = step_.calcStepModifyFactorFromErrorRatio(0.5)*ds_init_;
                if (ds_persistent_>ds_max || ds_persistent_<ds_min) {
                    acceptRestart();
                    return true;
                }
                // Inherited: a rejected estimate after an Interrupt reverts the persistent step
                else ds_persistent_ = ds_init_;
            }
            else {
                // Inherited: a rejected estimate after a binary update stays in the persistent step
                if (abs(ds_init_-ds_persistent_)/ds_init_>0.1) {
                    acceptRestart();
                    return true;
                }
            }
            return false;
        }

        friend class StepCheckpoint;

    public:
        //! constructor
        /*!
          @param[in] _step: symplectic step table
          @param[in] _limits: numeric limits
          @param[in] _ds: persistent Regularization step 
          @param[in] _fix_step_option: step option
         */
        RegularizationStepControl(const SymplecticStep& _step, const StepControlLimits& _limits, const Float _ds, const FixStepOption _fix_step_option): 
            step_(_step), limits_(_limits), ds_persistent_(_ds), fix_step_option_(_fix_step_option), 
            ds_{_ds,_ds}, ds_init_(_ds), ds_switch_(0), ds_backup_(_ds), time_end_flag_(false), step_count_(0), step_count_tsyn_(0) {}

        //! persistent Regularization step: it changes only on growth and on restart, it is not the last working step
        Float getPersistentStep() const { return ds_persistent_; }

        //! step option
        FixStepOption getStepOption() const { return fix_step_option_; }

        //! number of steps of the last advance
        long long unsigned int getStepCount() const { return step_count_; }

        //! number of time synchronization steps of the last advance
        long long unsigned int getStepCountSync() const { return step_count_tsyn_; }

        //! advance the stepper to a target Physical time
        /*!
          @param[in] _stepper: the integration to drive
          @param[in] _time_end: target Physical time
          \return how the advance ended
         */
        template <RegularizationStepper Tstepper>
        StepControlStatus advanceTo(Tstepper& _stepper, const Float _time_end) {
            // real full time step
            const Float dt_full = _time_end - _stepper.getTime();

            // time error 
            const Float time_error = limits_.time_error_max;

            // energy error limit
            const Float energy_error_rel_max = limits_.energy_error_relative_max;

            bool backup_flag=true; // flag for backup or restore

            const int cd_pair_size = step_.getCDPairSize();

            // start from the persistent step
            ds_[0] = ds_[1] = ds_persistent_;
            ds_init_   = ds_persistent_;
            ds_switch_ = 0;
            ds_backup_ = ReduceLevelStore(ds_persistent_);

            int reduce_ds_count=0; // number of reduce ds (ignore first few steps)
            Float step_modify_factor=1.0; // step modify factor 
            Float previous_step_modify_factor=1.0; // step modify factor 
            Float previous_error_ratio=-1; // previous error ratio when ds is reduced
            bool previous_is_restore=false; // previous step is reduced or not

            // time end control
            int n_step_end=0;  // number of steps integrated to reach the time end for one during the time sychronization sub steps
            time_end_flag_=false;

            // step count
            step_count_=0;
            step_count_tsyn_=0;

            Float*  ds = ds_;
            int&    ds_switch = ds_switch_;
            bool&   time_end_flag = time_end_flag_;

            StepReport report;
            report.time_prev = report.time_next = report.cck_prev = report.cck = 0.0;
            report.dt = report.integration_error = report.integration_error_ratio = 0.0;

            auto observe = [&](const StepEvent _event, const Float _step_modify_factor) {
                report.event = _event;
                report.ds = ds[ds_switch];
                report.ds_next = ds[1-ds_switch];
                report.ds_init = ds_init_;
                report.step_modify_factor = _step_modify_factor;
                report.reduce_count = reduce_ds_count;
                report.synchronizing = time_end_flag;
                report.step_count = step_count_;
                report.step_count_tsyn = step_count_tsyn_;
                report.fix_step_option = fix_step_option_;
                _stepper.observe(report);
            };

            // integration loop
            while(true) {
                if(backup_flag) {
                    StepCheckpoint checkpoint(*this);
                    CheckpointAction action = _stepper.checkpoint(checkpoint);
                    if (action==CheckpointAction::stop) return StepControlStatus::stopped;
                    if (action==CheckpointAction::abandon) return StepControlStatus::abandoned;
                    _stepper.save();
                }
                else _stepper.restore();

                // get real time 
                Float dt = _stepper.getTime();

                // integrate one step
                ASSERT(!ISINF(ds[ds_switch]));
                StepOutcome outcome = _stepper.step(ds[ds_switch]);
                const Float time = _stepper.getTime();
                const Float* time_table = outcome.time_table;

                // real step size
                dt =  time - dt;
                
                step_count_++;

                // get integration error for extended Hamiltonian
                Float integration_error_rel_abs = outcome.integration_error;

                Float integration_error_ratio = energy_error_rel_max/integration_error_rel_abs;
      
                Float error_increase_ratio_regular = step_.calcErrorRatioFromStepModifyFactor(2.0);

                report.dt = dt;
                report.integration_error = integration_error_rel_abs;
                report.integration_error_ratio = integration_error_ratio;
                observe(StepEvent::step_taken, step_modify_factor);

                // When time sychronization steps too large, give up
                if(step_count_tsyn_>limits_.step_count_max) return StepControlStatus::failed_to_synchronize;

                ASSERT(!ISNAN(integration_error_rel_abs));

                // modify step if energy error is large
                if(integration_error_rel_abs>energy_error_rel_max && fix_step_option_!=FixStepOption::always) {

                    bool check_flag = true;

                    // check whether already modified
                    if (previous_step_modify_factor!=1.0) {
                        ASSERT(previous_error_ratio>0.0);

                        // Inherited: if the error did not grow to twice the one before the previous reduction, the over-tolerance step is accepted
                        // if error does not reduce much, do not modify step anymore
                        if (integration_error_ratio>0.5*previous_error_ratio) check_flag=false;
                    }

                    if (check_flag) {
                        // for initial steps, reduce step permanently 
                        if(step_count_<5) {

                            // estimate the modification factor based on the symplectic order
                            // limit step_modify_factor to 0.125
                            step_modify_factor = std::max(regularStepFactor(step_.calcStepModifyFactorFromErrorRatio(integration_error_ratio)), Float(0.125));
                            ASSERT(step_modify_factor>0.0);

                            previous_step_modify_factor = step_modify_factor;
                            previous_error_ratio = integration_error_ratio;

                            ds[ds_switch] *= step_modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            // Inherited: the persistent step is not reduced
                            ds_backup_.initial(ds_persistent_);

                            backup_flag = false;
                            observe(StepEvent::large_energy_error, step_modify_factor);
                            continue;
                        }
                        // for big energy error, reduce step temparely
                        else if (fix_step_option_==FixStepOption::none) {

                            // estimate the modification factor based on the symplectic order
                            // limit step_modify_factor to 0.125
                            step_modify_factor = std::max(regularStepFactor(step_.calcStepModifyFactorFromErrorRatio(integration_error_ratio)), Float(0.125));
                            ASSERT(step_modify_factor>0.0);

                            previous_step_modify_factor = step_modify_factor;
                            previous_error_ratio = integration_error_ratio;

                            ds_backup_.backup(ds[ds_switch], step_modify_factor);
                            if(previous_is_restore) reduce_ds_count++;

                            ds[ds_switch] *= step_modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            ASSERT(!ISINF(ds[ds_switch]));

                            backup_flag = false;
                            observe(StepEvent::large_energy_error, step_modify_factor);
                            continue;
                        }
                    }
                }

                // if negative step, reduce step size
                if(!time_end_flag&&dt<0) {
                    // limit step_modify_factor to 0.125
                    step_modify_factor = std::min(std::max(regularStepFactor(step_.calcStepModifyFactorFromErrorRatio(abs(_time_end/dt))), Float(0.0625)),Float(0.5)); 
                    ASSERT(step_modify_factor>0.0);
                    previous_step_modify_factor = step_modify_factor;
                    previous_error_ratio = integration_error_ratio;

                    ds[ds_switch] *= step_modify_factor;
                    ds[1-ds_switch] = ds[ds_switch];
                    ASSERT(!ISINF(ds[ds_switch]));

                    // for initial steps, reduce step permanently
                    if (step_count_<5) {
                        ds_backup_.initial(ds_persistent_);
                    }
                    else { // reduce step temparely
                        ds_backup_.backup(ds[ds_switch], step_modify_factor);
                    }

                    backup_flag = false;

                    observe(StepEvent::negative_step, step_modify_factor);
                    continue;
                }

                // if no modification, reset previous values
                previous_step_modify_factor = 1.0;
                previous_error_ratio = -1.0;

                // check integration time
                if(time < _time_end - time_error){
                    // step increase depend on n_step_wait_recover_ds
                    if(fix_step_option_==FixStepOption::none && !time_end_flag) {
                        // waiting step count reach
                        previous_is_restore=ds_backup_.countAndRecover(ds[1-ds_switch], step_modify_factor, integration_error_ratio>error_increase_ratio_regular);
                        if (previous_is_restore) {
                            observe(StepEvent::reuse_backup_step, step_modify_factor);
                        }
                        // increase step size if energy error is small
                        else if(integration_error_rel_abs<0.5*energy_error_rel_max && dt>0.0 && dt_full/dt>std::max(100.0,0.02*limits_.step_count_max)) {
                            Float integration_error_ratio = energy_error_rel_max/integration_error_rel_abs;
                            Float step_modify_factor = std::min(Float(100.0),step_.calcStepModifyFactorFromErrorRatio(integration_error_ratio));
                            ASSERT(step_modify_factor>0.0);
                            ds[1-ds_switch] *= step_modify_factor;
                            ds_persistent_ = ds[1-ds_switch];
                            ASSERT(!ISINF(ds[1-ds_switch]));
                            observe(StepEvent::increase_step, step_modify_factor);
                        }
                    }

                    // time sychronization on case, when step size too small to reach time end, increase step size
                    if(time_end_flag && ds[ds_switch]==ds[1-ds_switch]) {
                        step_count_tsyn_++;

                        Float dt_end = _time_end - time;
                        if (dt<0) {
                            // limit step_modify_factor to 0.125
                            step_modify_factor = std::min(std::max(regularStepFactor(step_.calcStepModifyFactorFromErrorRatio(abs(_time_end/dt))), Float(0.0625)),Float(0.5)); 
                            ASSERT(step_modify_factor>0.0);

                            ds[ds_switch] *= step_modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            ASSERT(!ISINF(ds[ds_switch]));
                        }
                        else if (n_step_end>1 && dt<0.3*dt_end) {
                            // dt should be >0.0
                            ds[1-ds_switch] = ds[ds_switch] * dt_end/dt;
                            ASSERT(!ISINF(ds[1-ds_switch]));
                            observe(StepEvent::sync_enlarge_step, dt_end/dt);
                        }
                        else n_step_end++;
                    }

                    // when used once, update to the new step
                    ds[ds_switch] = ds[1-ds_switch]; 
                    ASSERT(!ISINF(ds[ds_switch]));
                    ds_switch = 1-ds_switch;

                    // Inherited: a step that does not advance Physical time is retried, not accepted
                    if (dt>0) backup_flag = true;
                    else backup_flag = false;
                }
                else if(time > _time_end + time_error) {
                    time_end_flag = true;
                    backup_flag = false;

                    step_count_tsyn_++;
                    n_step_end=0;

                    // check timetable
                    int i=-1,k=0; // i indicate the increasing time index, k is the corresponding index in time_table
                    for(i=0; i<cd_pair_size; i++) {
                        k = step_.getSortCumSumCKIndex(i);
                        if(_time_end<time_table[k]) break;
                    }
                    if (i==0) { // first step case
                        ASSERT(time_table[k]>0.0);
                        ds[ds_switch] *= step_.getSortCumSumCK(i)*_time_end/time_table[k];
                        ds[1-ds_switch] = ds[ds_switch];
                        ASSERT(!ISINF(ds[ds_switch]));
                        report.time_next = time_table[k];
                        report.cck = step_.getSortCumSumCK(i);
                        observe(StepEvent::sync_overshoot_first, step_modify_factor);
                    }
                    else { // not first step case, get the interval time 
                        // previous integrated sub time in time table
                        Float time_prev = time_table[step_.getSortCumSumCKIndex(i-1)];
                        Float dt_k = time_table[k] - time_prev;
                        Float ds_tmp = ds[ds_switch];
                        // get cumsum CK factor for two steps near the time_end
                        Float cck_prev = step_.getSortCumSumCK(i-1);
                        Float cck = step_.getSortCumSumCK(i);
                        // in case the time is between two sub step, first scale the next step with the previous step CumSum CK cck(i-1)
                        ASSERT(!ISINF(cck_prev));
                        ds[ds_switch] *= cck_prev;  
                        ASSERT(!ISINF(ds[ds_switch]));
                        // then next next step, scale with the CumSum CK between two step: cck(i) - cck(i-1) 
                        ASSERT(dt_k>0.0);
                        ds[1-ds_switch] = ds_tmp*(cck-cck_prev)*std::min(Float(1.0),(_time_end-time_prev+time_error)/dt_k); 
                        ASSERT(!ISINF(ds[1-ds_switch]));

                        report.time_prev = time_prev;
                        report.time_next = time_table[k];
                        report.cck_prev = cck_prev;
                        report.cck = cck;
                        observe(StepEvent::sync_overshoot_between, step_modify_factor);
                    }
                }
                else {
                    observe(StepEvent::finish, step_modify_factor);
                    break;
                }
            }

            return StepControlStatus::reached;
        }
    };

    inline bool StepCheckpoint::isSynchronizing() const { return control_.time_end_flag_; }

    inline bool StepCheckpoint::requestRestart(const RestartCause _cause, const Float _ds_estimate) { return control_.restart(_cause, _ds_estimate); }

    inline void StepCheckpoint::setStepOption(const FixStepOption _option) { control_.fix_step_option_ = _option; }

    inline Float StepCheckpoint::getStepInit() const { return control_.ds_init_; }

    inline Float StepCheckpoint::getStepNow() const { return control_.ds_[0]; }

}
