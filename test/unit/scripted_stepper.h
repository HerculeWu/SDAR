#pragma once
#include <functional>
#include <string>
#include <vector>
#include "AR/regularization_step_control.h"

//! What a scripted step does: advance of Physical time and integration error
struct ScriptedStep {
    Float dt;
    Float error;
};

//! Fake stepper driven by a script, recording what Regularization step control asked for
/*! Physical time advances by the scripted dt; within a step it is taken to advance linearly with the
    cumulative coefficients of the symplectic step table, which gives the sub-step time table.
 */
struct ScriptedStepper {
    const AR::SymplecticStep& table;
    Float time = 0.0;
    Float time_saved = 0.0;

    //! script of the n-th step (counted from 0) with the requested Regularization step; default: dt = ds, tiny error
    std::function<ScriptedStep(int, Float)> on_step = [](int, Float _ds) { return ScriptedStep{_ds, 0.0}; };
    //! script of the checkpoint before the n-th step
    std::function<AR::CheckpointAction(int, AR::StepCheckpoint&)> on_checkpoint = 
        [](int, AR::StepCheckpoint&) { return AR::CheckpointAction::proceed; };

    std::vector<Float> ds;               // Regularization step of every step
    std::string calls;                   // C: checkpoint, S: save, T: restore, P: step
    std::vector<AR::StepEvent> events;   // every reported event except step_taken
    std::vector<AR::StepReport> reports; // the reports of those events
    std::vector<Float> time_table;

    ScriptedStepper(const AR::SymplecticStep& _table): table(_table), time_table(_table.getCDPairSize()) {}

    int countSteps() const { return int(ds.size()); }

    bool saw(const AR::StepEvent _event) const {
        return std::find(events.begin(), events.end(), _event)!=events.end();
    }

    Float getTime() const { return time; }

    AR::CheckpointAction checkpoint(AR::StepCheckpoint& _checkpoint) {
        calls += 'C';
        return on_checkpoint(countSteps(), _checkpoint);
    }

    void save() {
        calls += 'S';
        time_saved = time;
    }

    void restore() {
        calls += 'T';
        time = time_saved;
    }

    AR::StepOutcome step(const Float _ds) {
        ScriptedStep script = on_step(countSteps(), _ds);
        calls += 'P';
        ds.push_back(_ds);
        for (int i=0; i<table.getCDPairSize(); i++) 
            time_table[table.getSortCumSumCKIndex(i)] = time + script.dt*table.getSortCumSumCK(i);
        time += script.dt;
        return AR::StepOutcome{script.error, time_table.data()};
    }

    void observe(const AR::StepReport& _report) {
        if (_report.event!=AR::StepEvent::step_taken) {
            events.push_back(_report.event);
            reports.push_back(_report);
        }
    }
};
