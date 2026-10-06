// Regularization step control driven by a scripted stepper
#include "doctest.h"
#include "test_support.h"
#include "scripted_stepper.h"

using AR::CheckpointAction;
using AR::FixStepOption;
using AR::RestartCause;
using AR::StepControlStatus;
using AR::StepEvent;

namespace {

    const Float E_MAX = 1e-10;       // energy error limit
    const Float T_ERR = 1e-9;        // time error limit
    const Float SMALL = E_MAX/100.0; // error small enough to recover a reduced step
    const Float FAIR  = 0.6*E_MAX;   // error under the limit, too large to recover or grow
    const Float T_TOL = 1.001*T_ERR; // time error limit with room for round-off

    struct Fixture {
        AR::SymplecticStep table;
        AR::StepControlLimits limits = {T_ERR, E_MAX, 1000000};
        ScriptedStepper stepper;

        Fixture(): stepper((table.initialSymplecticCofficients(-6), table)) {}
        ~Fixture() { table.clear(); }

        AR::RegularizationStepControl control(const Float _ds, const FixStepOption _option) {
            return AR::RegularizationStepControl(table, limits, _ds, _option);
        }
    };

    //! one unit of Physical time per step whatever the Regularization step, so that the step never grows
    ScriptedStep unitStep(const Float _error) { return ScriptedStep{1.0, _error}; }
}

TEST_CASE_FIXTURE(Fixture, "a step with fine error reaching the target ends the advance") {
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 3.0)==StepControlStatus::reached);
    CHECK(stepper.calls=="CSPCSPCSP");
    CHECK(stepper.ds==std::vector<Float>{1.0, 1.0, 1.0});
    CHECK(c.getStepCount()==3);
    CHECK(c.getStepCountSync()==0);
    CHECK(c.getPersistentStep()==1.0);
}

TEST_CASE_FIXTURE(Fixture, "reduction: an over-tolerance step among the first steps is retried with a smaller step") {
    stepper.on_step = [](int _n, Float) { return unitStep(_n==1 ? 100.0*E_MAX : FAIR); };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 4.0)==StepControlStatus::reached);
    // error ratio 0.01 scales the step by 0.01^(1/6)=0.46, rounded down to a power of two
    CHECK(stepper.ds==std::vector<Float>{1.0, 1.0, 0.25, 0.25, 0.25});
    CHECK(stepper.calls=="CSPCSPTPCSPCSP");
    CHECK(stepper.events.front()==StepEvent::large_energy_error);
    // the persistent step is not the working step
    CHECK(c.getPersistentStep()==1.0);
}

TEST_CASE_FIXTURE(Fixture, "reduction: one reduction is limited to a factor 0.125") {
    stepper.on_step = [](int _n, Float) { return unitStep(_n==0 ? 1e12*E_MAX : FAIR); };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 2.0)==StepControlStatus::reached);
    CHECK(stepper.ds==std::vector<Float>{1.0, 0.125, 0.125});
}

TEST_CASE_FIXTURE(Fixture, "reduction: after the first steps the step is reduced only when the step option is none") {
    stepper.on_step = [](int _n, Float) { return unitStep(_n==6 ? 100.0*E_MAX : FAIR); };

    SUBCASE("none: reduced") {
        auto c = control(1.0, FixStepOption::none);
        CHECK(c.advanceTo(stepper, 8.0)==StepControlStatus::reached);
        CHECK(stepper.ds[6]==1.0);
        CHECK(stepper.ds[7]==0.25);
        CHECK(stepper.calls.substr(3*6)=="CSPTPCSP");
    }
    SUBCASE("later: accepted as it is") {
        auto c = control(1.0, FixStepOption::later);
        CHECK(c.advanceTo(stepper, 8.0)==StepControlStatus::reached);
        CHECK(stepper.ds==std::vector<Float>(8, 1.0));
        CHECK(stepper.events==std::vector<StepEvent>{StepEvent::finish});
    }
}

TEST_CASE_FIXTURE(Fixture, "reduction: never with step option always") {
    stepper.on_step = [](int, Float) { return unitStep(100.0*E_MAX); };
    auto c = control(1.0, FixStepOption::always);
    CHECK(c.advanceTo(stepper, 8.0)==StepControlStatus::reached);
    CHECK(stepper.ds==std::vector<Float>(8, 1.0));
}

TEST_CASE_FIXTURE(Fixture, "recovery: a reduced step returns after its wait once the error is small") {
    // reduced by 0.25 at step 6: the wait is 2/0.25 = 8 accepted steps
    stepper.on_step = [](int _n, Float) { return unitStep(_n==6 ? 100.0*E_MAX : SMALL); };
    auto c = control(1.0, FixStepOption::none);
    CHECK(c.advanceTo(stepper, 20.0)==StepControlStatus::reached);
    for (int i=7; i<=15; i++) CHECK(stepper.ds[i]==0.25);
    CHECK(stepper.ds[16]==1.0);
    CHECK(stepper.saw(StepEvent::reuse_backup_step));
    CHECK(c.getPersistentStep()==1.0);
}

TEST_CASE_FIXTURE(Fixture, "recovery: waits further while the error is not small") {
    stepper.on_step = [](int _n, Float) { return unitStep(_n==6 ? 100.0*E_MAX : (_n<20 ? FAIR : SMALL)); };
    auto c = control(1.0, FixStepOption::none);
    CHECK(c.advanceTo(stepper, 25.0)==StepControlStatus::reached);
    CHECK(stepper.ds[20]==0.25);
    CHECK(stepper.ds[21]==1.0);
}

TEST_CASE_FIXTURE(Fixture, "growth: a small error far from the target grows the step, and the persistent step with it") {
    // error ratio 1e6 scales the step by 1e6^(1/6) = 10
    stepper.on_step = [](int _n, Float _ds) { return ScriptedStep{_ds, _n==0 ? Float(1e-6*E_MAX) : FAIR}; };
    auto c = control(1.0, FixStepOption::none);
    const Float time_end = 1e6;
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint&) { return _n<3 ? CheckpointAction::proceed : CheckpointAction::stop; };
    CHECK(c.advanceTo(stepper, time_end)==StepControlStatus::stopped);
    CHECK(stepper.ds[0]==1.0);
    CHECK(stepper.ds[1]==doctest::Approx(10.0));
    CHECK(stepper.events==std::vector<StepEvent>{StepEvent::increase_step});
    CHECK(c.getPersistentStep()==stepper.ds[1]);
}

TEST_CASE_FIXTURE(Fixture, "growth: limited to a factor 100 per step") {
    stepper.on_step = [](int _n, Float _ds) { return ScriptedStep{_ds, _n==0 ? Float(0.0) : FAIR}; };
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint&) { return _n<2 ? CheckpointAction::proceed : CheckpointAction::stop; };
    auto c = control(1.0, FixStepOption::none);
    CHECK(c.advanceTo(stepper, 1e9)==StepControlStatus::stopped);
    CHECK(stepper.ds[1]==100.0);
    CHECK(c.getPersistentStep()==100.0);
}

TEST_CASE_FIXTURE(Fixture, "growth: none when few steps are left, or when the step option is not none") {
    stepper.on_step = [](int, Float _ds) { return ScriptedStep{_ds, 0.0}; };

    SUBCASE("the target is within max(100, 2% of the step count limit) steps") {
        auto c = control(1.0, FixStepOption::none);
        CHECK(c.advanceTo(stepper, 20000.0)==StepControlStatus::reached);
        CHECK(c.getPersistentStep()==1.0);
        CHECK(c.getStepCount()==20000);
    }
    SUBCASE("step option later") {
        stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint&) { return _n<3 ? CheckpointAction::proceed : CheckpointAction::stop; };
        auto c = control(1.0, FixStepOption::later);
        CHECK(c.advanceTo(stepper, 1e9)==StepControlStatus::stopped);
        CHECK(stepper.ds==std::vector<Float>(3, 1.0));
    }
}

TEST_CASE_FIXTURE(Fixture, "synchronization: an overshooting step is rejected and split at the sub-step before the target") {
    auto c = control(1.0, FixStepOption::later);
    const Float time_end = 2.7; // between the sub-steps at 2.61 and 2.90 of the third step
    CHECK(c.advanceTo(stepper, time_end)==StepControlStatus::reached);
    CHECK(abs(stepper.time-time_end)<=T_TOL);

    // the third step passes the target: restored, then retried in two parts
    CHECK(stepper.calls.substr(0, 11)=="CSPCSPCSPTP");
    REQUIRE(stepper.events.front()==StepEvent::sync_overshoot_between);
    const AR::StepReport& report = stepper.reports.front();
    CHECK(report.time_prev<time_end);
    CHECK(report.time_next>time_end);
    // first part: up to the last sub-step before the target
    CHECK(stepper.ds[3]==1.0*report.cck_prev);
    // second part: the fraction of the sub-step interval that is left
    CHECK(stepper.ds[4]==1.0*(report.cck-report.cck_prev)*((time_end-report.time_prev+T_ERR)/(report.time_next-report.time_prev)));
    CHECK(c.getStepCountSync()>=1);
    // synchronization does not touch the persistent step
    CHECK(c.getPersistentStep()==1.0);
}

TEST_CASE_FIXTURE(Fixture, "synchronization: overshoot inside the first sub-step scales the whole step") {
    auto c = control(1.0, FixStepOption::later);
    // smallest cumulative coefficient of the table: the first sub-step in time
    const Float cck0 = table.getSortCumSumCK(0);
    REQUIRE(cck0>0.0);
    const Float time_end = 0.5*cck0;
    CHECK(c.advanceTo(stepper, time_end)==StepControlStatus::reached);
    REQUIRE(stepper.events.front()==StepEvent::sync_overshoot_first);
    CHECK(stepper.ds[1]==1.0*(cck0*time_end/stepper.reports.front().time_next));
    CHECK(abs(stepper.time-time_end)<=T_TOL);
}

TEST_CASE_FIXTURE(Fixture, "synchronization: a step that keeps falling short of the target is enlarged") {
    // after the overshoot, Physical time advances ten times slower than expected
    bool slow = false;
    stepper.on_step = [&](int, Float _ds) { return ScriptedStep{slow ? 0.1*_ds : _ds, FAIR}; };
    auto c = control(1.0, FixStepOption::later);
    const Float time_end = 2.7;
    stepper.on_checkpoint = [&](int, AR::StepCheckpoint& _checkpoint) { 
        slow = _checkpoint.isSynchronizing(); 
        return CheckpointAction::proceed; 
    };
    CHECK(c.advanceTo(stepper, time_end)==StepControlStatus::reached);
    REQUIRE(stepper.saw(StepEvent::sync_enlarge_step));
    CHECK(abs(stepper.time-time_end)<=T_TOL);
    // the same step is tried three times before it is enlarged
    const int n = stepper.countSteps();
    CHECK(stepper.ds[n-2]==stepper.ds[n-3]);
    CHECK(stepper.ds[n-3]==stepper.ds[n-4]);
    CHECK(stepper.ds[n-1]>stepper.ds[n-2]);
}

TEST_CASE_FIXTURE(Fixture, "synchronization: failure is returned as a status") {
    limits.step_count_max = 5;
    // every step lands far beyond the target
    stepper.on_step = [&](int, Float) { return ScriptedStep{10.0, FAIR}; };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 2.5)==StepControlStatus::failed_to_synchronize);
    CHECK(c.getStepCountSync()==6);
    CHECK(c.getStepCount()==7);
    // the last report tells the stepper, which does the printing
    CHECK(stepper.events.back()==StepEvent::synchronization_failed);
}

TEST_CASE_FIXTURE(Fixture, "negative step: Physical time moving backwards before the target is retried with a smaller step") {
    stepper.on_step = [](int _n, Float) { return _n==1 ? ScriptedStep{-1.0, FAIR} : unitStep(FAIR); };
    auto c = control(1.0, FixStepOption::later);
    // |time_end/dt| = 4 scales the step by 4^(1/6) = 1.26, which is limited to 0.5
    CHECK(c.advanceTo(stepper, 4.0)==StepControlStatus::reached);
    CHECK(stepper.calls.substr(0, 8)=="CSPCSPTP");
    CHECK(stepper.ds[2]==0.5);
    CHECK(stepper.events.front()==StepEvent::negative_step);
}

TEST_CASE_FIXTURE(Fixture, "restart after an Interrupt: an estimate outside the band replaces the step") {
    // band: [0.5^(1/6), 2^(1/6)] = [0.89, 1.12] times the reference step
    bool accepted = false;
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint& _checkpoint) {
        if (_n==2) accepted = _checkpoint.requestRestart(RestartCause::interrupt, 0.5);
        return CheckpointAction::proceed;
    };
    stepper.on_step = [](int, Float) { return unitStep(FAIR); };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 4.0)==StepControlStatus::reached);
    CHECK(accepted);
    CHECK(stepper.ds==std::vector<Float>{1.0, 1.0, 0.5, 0.5});
    CHECK(c.getPersistentStep()==0.5);
}

TEST_CASE_FIXTURE(Fixture, "restart after an Interrupt: a larger estimate becomes the persistent step but does not enlarge the working step") {
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint& _checkpoint) {
        if (_n==2) CHECK(_checkpoint.requestRestart(RestartCause::interrupt, 3.0));
        return CheckpointAction::proceed;
    };
    stepper.on_step = [](int, Float) { return unitStep(FAIR); };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 4.0)==StepControlStatus::reached);
    CHECK(stepper.ds==std::vector<Float>(4, 1.0));
    CHECK(c.getPersistentStep()==3.0);
}

TEST_CASE_FIXTURE(Fixture, "restart after an Interrupt, inherited: an estimate inside the band is dropped and the persistent step reverts") {
    bool accepted = true;
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint& _checkpoint) {
        if (_n==2) accepted = _checkpoint.requestRestart(RestartCause::interrupt, 1.1);
        return CheckpointAction::proceed;
    };
    stepper.on_step = [](int, Float) { return unitStep(FAIR); };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 4.0)==StepControlStatus::reached);
    CHECK(!accepted);
    CHECK(stepper.ds==std::vector<Float>(4, 1.0));
    CHECK(c.getPersistentStep()==1.0);
}

TEST_CASE_FIXTURE(Fixture, "restart after a binary update: an estimate more than 10% away replaces the step") {
    bool accepted = false;
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint& _checkpoint) {
        if (_n==2) accepted = _checkpoint.requestRestart(RestartCause::binary_update, 0.89);
        return CheckpointAction::proceed;
    };
    stepper.on_step = [](int, Float) { return unitStep(FAIR); };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 4.0)==StepControlStatus::reached);
    // 0.89 is inside the band of an Interrupt restart, but outside the one of a binary update
    CHECK(accepted);
    CHECK(stepper.ds==std::vector<Float>{1.0, 1.0, 0.89, 0.89});
    CHECK(c.getPersistentStep()==0.89);
}

TEST_CASE_FIXTURE(Fixture, "restart after a binary update, inherited: an estimate within 10% is not used but stays in the persistent step") {
    bool accepted = true;
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint& _checkpoint) {
        if (_n==2) accepted = _checkpoint.requestRestart(RestartCause::binary_update, 0.95);
        return CheckpointAction::proceed;
    };
    stepper.on_step = [](int, Float) { return unitStep(FAIR); };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 4.0)==StepControlStatus::reached);
    CHECK(!accepted);
    CHECK(stepper.ds==std::vector<Float>(4, 1.0));
    CHECK(c.getPersistentStep()==0.95);
}

TEST_CASE_FIXTURE(Fixture, "a step-option change at a checkpoint takes effect for the following steps") {
    stepper.on_step = [](int _n, Float) { return unitStep((_n==6||_n==10) ? 100.0*E_MAX : FAIR); };
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint& _checkpoint) {
        if (_n==9) _checkpoint.setStepOption(FixStepOption::later);
        return CheckpointAction::proceed;
    };
    auto c = control(1.0, FixStepOption::none);
    CHECK(c.advanceTo(stepper, 11.0)==StepControlStatus::reached);
    // reduced while the option is none ...
    CHECK(stepper.ds[7]==0.25);
    // ... and accepted over tolerance once it is later
    CHECK(stepper.calls.substr(stepper.calls.size()-6)=="CSPCSP");
    CHECK(c.getStepOption()==FixStepOption::later);
}

TEST_CASE_FIXTURE(Fixture, "a checkpoint can stop or abandon the advance") {
    auto c = control(1.0, FixStepOption::later);
    SUBCASE("stop") {
        stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint&) { return _n==2 ? CheckpointAction::stop : CheckpointAction::proceed; };
        CHECK(c.advanceTo(stepper, 5.0)==StepControlStatus::stopped);
    }
    SUBCASE("abandon") {
        stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint&) { return _n==2 ? CheckpointAction::abandon : CheckpointAction::proceed; };
        CHECK(c.advanceTo(stepper, 5.0)==StepControlStatus::abandoned);
    }
    // nothing is saved or stepped after the checkpoint
    CHECK(stepper.calls=="CSPCSPC");
    CHECK(c.getStepCount()==2);
}

TEST_CASE_FIXTURE(Fixture, "each advance starts again from the persistent step") {
    stepper.on_step = [](int _n, Float) { return unitStep(_n==0 ? 100.0*E_MAX : FAIR); };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 2.0)==StepControlStatus::reached);
    CHECK(stepper.ds.back()==0.25);
    CHECK(c.advanceTo(stepper, 4.0)==StepControlStatus::reached);
    CHECK(stepper.ds.back()==1.0);
    CHECK(c.getStepCount()==2);
}

// ---- inherited quirks: kept on purpose, see the spec "Regularization step control and AR dynamics as typed variants"

TEST_CASE_FIXTURE(Fixture, "inherited: an over-tolerance retry is accepted unless its error doubled") {
    // step 6 is rejected with 100 times the limit; what happens to the retry (step 7) depends on its error
    Float retry_error = 0.0;
    stepper.on_step = [&](int _n, Float) { return unitStep(_n==6 ? 100.0*E_MAX : (_n==7 ? retry_error : FAIR)); };
    auto c = control(1.0, FixStepOption::none);

    SUBCASE("1.5 times the previous error: accepted although over the limit") {
        retry_error = 150.0*E_MAX;
        CHECK(c.advanceTo(stepper, 8.0)==StepControlStatus::reached);
        CHECK(stepper.calls.substr(3*6)=="CSPTPCSP");
        CHECK(stepper.ds[8]==0.25);
    }
    SUBCASE("3 times the previous error: reduced again") {
        retry_error = 300.0*E_MAX;
        CHECK(c.advanceTo(stepper, 8.0)==StepControlStatus::reached);
        CHECK(stepper.calls.substr(3*6)=="CSPTPTPCSP");
        CHECK(stepper.ds[8]==0.0625);
    }
}

TEST_CASE_FIXTURE(Fixture, "inherited: a step that does not advance Physical time is retried, not accepted") {
    stepper.on_step = [](int _n, Float) { return _n==1 ? ScriptedStep{0.0, FAIR} : unitStep(FAIR); };
    auto c = control(1.0, FixStepOption::later);
    CHECK(c.advanceTo(stepper, 3.0)==StepControlStatus::reached);
    // restored and repeated with the same step
    CHECK(stepper.calls=="CSPCSPTPCSP");
    CHECK(stepper.ds==std::vector<Float>(4, 1.0));
}

TEST_CASE_FIXTURE(Fixture, "inherited: a step moving Physical time backwards during synchronization is reduced and retried") {
    // the first step after the checkpoint that follows the overshoot goes backwards
    bool synchronizing = false;
    int n_backward = -1;
    stepper.on_checkpoint = [&](int, AR::StepCheckpoint& _checkpoint) { 
        synchronizing = _checkpoint.isSynchronizing(); 
        return CheckpointAction::proceed; 
    };
    stepper.on_step = [&](int _n, Float _ds) { 
        if (synchronizing && n_backward<0) {
            n_backward = _n;
            return ScriptedStep{-0.01, FAIR};
        }
        return ScriptedStep{_ds, FAIR}; 
    };
    auto c = control(1.0, FixStepOption::later);
    const Float time_end = 2.7;
    CHECK(c.advanceTo(stepper, time_end)==StepControlStatus::reached);
    REQUIRE(n_backward>0);
    // not accepted: state is restored before the next step
    int n_step = 0;
    for (std::size_t i=0; i<stepper.calls.size(); i++) {
        if (stepper.calls[i]=='P' && n_step++==n_backward) {
            CHECK(stepper.calls[i+1]=='T');
            break;
        }
    }
    // |time_end/dt| = 270 scales the step by 270^(1/6) = 2.5, which is limited to 0.5
    CHECK(stepper.ds[n_backward+1]==0.5*stepper.ds[n_backward]);
    CHECK(abs(stepper.time-time_end)<=T_TOL);
}

TEST_CASE_FIXTURE(Fixture, "inherited: a restart keeps the reduce level, so a reduced step recovers to the estimate without waiting") {
    // reduced by 0.25 at step 6 (wait 8); restart with estimate 0.5 before step 8
    stepper.on_step = [](int _n, Float) { return unitStep(_n==6 ? 100.0*E_MAX : SMALL); };
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint& _checkpoint) {
        if (_n==8) CHECK(_checkpoint.requestRestart(RestartCause::binary_update, 0.5));
        return CheckpointAction::proceed;
    };
    auto c = control(1.0, FixStepOption::none);
    CHECK(c.advanceTo(stepper, 12.0)==StepControlStatus::reached);
    CHECK(stepper.ds[8]==0.25);
    CHECK(stepper.ds[9]==0.5);
    CHECK(stepper.ds.back()==0.5);
}

TEST_CASE_FIXTURE(Fixture, "inherited: a restart initializes ten of the eleven reduce levels, the eleventh keeps its step and wait") {
    // eleven reductions by 0.25 (steps 6, 8, ..., 26), each followed by an accepted retry, fill all levels
    const int n_restart = 28;
    stepper.on_step = [&](int _n, Float) { 
        if (_n>=6 && _n<n_restart) return unitStep(_n%2==0 ? 100.0*E_MAX : FAIR);
        return unitStep(SMALL);
    };
    stepper.on_checkpoint = [&](int _n, AR::StepCheckpoint& _checkpoint) {
        if (_n==n_restart) CHECK(_checkpoint.requestRestart(RestartCause::binary_update, 1e-3));
        return CheckpointAction::proceed;
    };
    auto c = control(1.0, FixStepOption::none);
    CHECK(c.advanceTo(stepper, 50.0)==StepControlStatus::reached);

    const Float ds_level10 = std::pow(0.25, 10); // step stored in the eleventh level
    const Float ds_reduced = std::pow(0.25, 11); // working step at the restart
    REQUIRE(stepper.ds[n_restart]==ds_reduced);
    // the eleventh level still waits: 8 steps minus the accepted retry
    for (int i=0; i<=7; i++) CHECK(stepper.ds[n_restart+i]==ds_reduced);
    // then it recovers its stale step, not the estimate
    CHECK(stepper.ds[n_restart+8]==ds_level10);
    // the ten initialized levels recover to the estimate at once
    CHECK(stepper.ds[n_restart+9]==1e-3);
    CHECK(stepper.ds.back()==1e-3);
}
