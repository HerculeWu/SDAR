// Replay of the trace fixtures recorded from the integration loop before
// Regularization step control was extracted (format: test/fixtures/step_control/README.md)
#include "doctest.h"
#include "test_support.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "AR/regularization_step_control.h"

namespace {

    using Line = std::vector<std::string>;

    Float hex(const std::string& _s) { return std::strtod(_s.c_str(), NULL); }

    //! Stepper that replays one recorded integrateToTime call and checks every call of step control against it
    struct ReplayStepper {
        const std::vector<Line>& lines;
        std::size_t next;
        const AR::RegularizationStepControl* control;
        Float time, time_saved, ds_persistent;
        std::vector<Float> time_table;
        std::string error; // first mismatch

        ReplayStepper(const std::vector<Line>& _lines, const int _cd_pair_size): 
            lines(_lines), next(1), control(NULL), time(hex(_lines[0][3])), time_saved(time), ds_persistent(hex(_lines[0][1])), time_table(_cd_pair_size) {}

        void fail(const std::string& _what) {
            if (error.empty()) error = _what + " at event " + std::to_string(next);
        }

        // "D" lines say what the persistent Regularization step is from here on
        void checkPersistentStep() {
            while (next<lines.size() && lines[next][0]=="D") ds_persistent = hex(lines[next++][1]);
            if (control->getPersistentStep()!=ds_persistent) fail("persistent Regularization step differs");
        }

        bool expect(const char* _key) {
            if (next<lines.size() && lines[next][0]==_key) return true;
            fail(std::string("expected ") + (next<lines.size()? lines[next][0] : "nothing") + ", step control asked for " + _key);
            return false;
        }

        Float getTime() const { return time; }

        AR::CheckpointAction checkpoint(AR::StepCheckpoint& _checkpoint) {
            checkPersistentStep();
            if (!expect("C")) return AR::CheckpointAction::stop;
            next++;
            while (next<lines.size()) {
                const Line& line = lines[next];
                if (line[0]=="O") _checkpoint.setStepOption(AR::FixStepOption(std::stoi(line[1])));
                else if (line[0]=="R") _checkpoint.requestRestart(line[1]=="I" ? AR::RestartCause::interrupt : AR::RestartCause::binary_update, hex(line[2]));
                else if (line[0]=="X") { next++; return AR::CheckpointAction::stop; }
                else break;
                next++;
            }
            return AR::CheckpointAction::proceed;
        }

        void save() {
            checkPersistentStep();
            if (expect("S")) next++;
            time_saved = time;
        }

        void restore() {
            checkPersistentStep();
            if (expect("T")) next++;
            time = time_saved;
        }

        AR::StepOutcome step(const Float _ds) {
            AR::StepOutcome outcome = {0.0, time_table.data()};
            if (!expect("P")) throw std::runtime_error(error);
            const Line& line = lines[next++];
            if (_ds!=hex(line[1])) fail("Regularization step differs");
            time = hex(line[2]);
            outcome.integration_error = hex(line[3]);
            for (std::size_t i=0; i<time_table.size(); i++) 
                time_table[i] = line.size()>4 ? hex(line[4+i]) : std::nan("");
            return outcome;
        }

        void observe(const AR::StepReport&) {}
    };

    //! replay every record of one fixture file, return the number of records
    int replay(const std::string& _name) {
        std::ifstream fin(std::string(FIXTURE_DIR) + "/step_control/" + _name);
        REQUIRE(fin.is_open());
        std::vector<Line> lines;
        std::string text;
        int n_record = 0;
        while (std::getline(fin, text)) {
            std::istringstream words(text);
            Line line;
            for (std::string word; words>>word;) line.push_back(word);
            lines.push_back(line);
            if (line[0]!="END") continue;

            const Line& begin = lines[0];
            REQUIRE(begin[0]=="BEGIN");
            const int order = std::stoi(begin[8]), cd_pair_size = std::stoi(begin[9]);
            AR::SymplecticStep table;
            // order 6 with 8 pairs and order 8 with 16 pairs are the second Yoshida solution
            table.initialSymplecticCofficients(((order==6&&cd_pair_size==8)||(order==8&&cd_pair_size==16)) ? -order : order);
            REQUIRE(table.getCDPairSize()==cd_pair_size);
            AR::StepControlLimits limits = {hex(begin[5]), hex(begin[6]), std::stoull(begin[7])};
            AR::RegularizationStepControl control(table, limits, hex(begin[1]), AR::FixStepOption(std::stoi(begin[2])));
            ReplayStepper stepper(lines, cd_pair_size);
            stepper.control = &control;

            AR::StepControlStatus status = control.advanceTo(stepper, hex(begin[4]));
            stepper.checkPersistentStep();

            INFO(_name << " record " << n_record);
            CHECK(stepper.error=="");
            CHECK(stepper.next==lines.size()-1); // every recorded event was consumed
            CHECK((status==AR::StepControlStatus::stopped) == (line[1]=="stop"));
            CHECK(control.getPersistentStep()==hex(line[2]));
            CHECK(int(control.getStepOption())==std::stoi(line[3]));
            CHECK(control.getStepCount()==std::stoull(line[4]));
            CHECK(control.getStepCountSync()==std::stoull(line[5]));

            table.clear();
            lines.clear();
            n_record++;
        }
        return n_record;
    }
}

TEST_CASE("replay: restart after a binary update, reduction, recovery and overshoot (Tree slow-down)") {
    CHECK(replay("binary-update.trace")==4);
}

TEST_CASE("replay: the same run under Array slow-down") {
    CHECK(replay("array-slow-down.trace")==2);
}

TEST_CASE("replay: growth of the Regularization step") {
    CHECK(replay("growth.trace")==1);
}

TEST_CASE("replay: restart after an Interrupt and step-option change") {
    CHECK(replay("interrupt.trace")==4);
}

TEST_CASE("replay: early stop by an Interrupt") {
    CHECK(replay("merge-stop.trace")==1);
}

TEST_CASE("replay: short calls driven by Hermite") {
    CHECK(replay("hermite.trace")==8);
}
