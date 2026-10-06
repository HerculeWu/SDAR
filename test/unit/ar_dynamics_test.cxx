// The AR integrator through its public interface, for every combination of
// time transformation and Slow-down scheme, all compiled into this one program
#include "doctest.h"
#include "test_support.h"
#include "ar_test_system.h"

typedef Variant<AR::LogH, AR::NoSlowDown>    LogHNone;
typedef Variant<AR::TTL,  AR::NoSlowDown>    TTLNone;
typedef Variant<AR::LogH, AR::ArraySlowDown> LogHArray;
typedef Variant<AR::TTL,  AR::ArraySlowDown> TTLArray;
typedef Variant<AR::LogH, AR::TreeSlowDown>  LogHTree;
typedef Variant<AR::TTL,  AR::TreeSlowDown>  TTLTree;

TYPE_TO_STRING_AS("LogH", LogHNone);
TYPE_TO_STRING_AS("TTL", TTLNone);
TYPE_TO_STRING_AS("LogH + Array slow-down", LogHArray);
TYPE_TO_STRING_AS("TTL + Array slow-down", TTLArray);
TYPE_TO_STRING_AS("LogH + Tree slow-down", LogHTree);
TYPE_TO_STRING_AS("TTL + Tree slow-down", TTLTree);

#define ALL_VARIANTS LogHNone, TTLNone, LogHArray, TTLArray, LogHTree, TTLTree

namespace {
    const Float PI = 3.14159265358979323846;
    const Float NO_SLOWDOWN = 1e-12; // maximum Slow-down timescale below every period: all Slow-down factors are 1

    //! equal-mass binary, semi-major axis 1, eccentricity 0.5, at apocentre; the period is 2 pi
    std::vector<Particle> keplerBinary() {
        std::vector<Particle> bodies = {makeParticle(1.0, 0.0, 0.0, 0.0, 0.0)};
        splitLastIntoBinary(bodies, 0.5, 0.5, 1.0, 0.5);
        return bodies;
    }

    //! hierarchical triple: inner binary (semi 0.01, period 0.0063) orbiting a third body (semi 1, period 4.4)
    std::vector<Particle> hierarchicalTriple() {
        std::vector<Particle> bodies = {makeParticle(2.0, 0.0, 0.0, 0.0, 0.0)};
        splitLastIntoBinary(bodies, 1.0, 1.0, 1.0, 0.2);
        splitLastIntoBinary(bodies, 0.5, 0.5, 0.01, 0.3);
        return bodies;
    }

    //! two-level hierarchy below the root: ((binary, single), single) with separations 0.001, 0.1 and 10
    std::vector<Particle> hierarchicalQuadruple() {
        std::vector<Particle> bodies = {makeParticle(3.0, 0.0, 0.0, 0.0, 0.0)};
        splitLastIntoBinary(bodies, 1.0, 2.0, 10.0, 0.1);
        splitLastIntoBinary(bodies, 1.0, 1.0, 0.1, 0.1);
        splitLastIntoBinary(bodies, 0.5, 0.5, 0.001, 0.1);
        return bodies;
    }

    const Float TRIPLE_TIME = 0.05; // eight inner periods
}

TEST_CASE_TEMPLATE("Kepler orbit returns to its start after one period and conserves energy", V, LogHNone, TTLNone, TTLArray, TTLTree) {
    const std::vector<Particle> start = keplerBinary();
    ArSystem<typename V::Integrator> system(start, NO_SLOWDOWN);
    const Float etot = system.ar.getEtotRef();
    CHECK(etot==doctest::Approx(-0.125));

    system.integrate(2.0*PI);

    CHECK(system.ar.getTime()==doctest::Approx(2.0*PI).epsilon(1e-12));
    for (int i=0; i<2; i++) {
        for (int k=0; k<2; k++) {
            CHECK(abs(system.ar.particles[i].pos[k]-start[i].pos[k])<1e-7);
            CHECK(abs(system.ar.particles[i].vel[k]-start[i].vel[k])<1e-7);
        }
    }
    CHECK(abs(system.ar.getEnergyError()/etot)<1e-9);
}

// Inherited: the two-body step refreshes the Slow-down energies only after its last sub-step, while the LogH
// time transformation reads them at every sub-step. A Kepler orbit under LogH with a Slow-down scheme therefore
// carries an energy error of order 1e-3 even when its Slow-down factor is 1. Found while writing these tests and kept
// as it is, because results must stay bitwise-identical to the previous version.
TEST_CASE_TEMPLATE("inherited: Kepler orbit under LogH with a Slow-down scheme is only approximate", V, LogHArray, LogHTree) {
    const std::vector<Particle> start = keplerBinary();
    ArSystem<typename V::Integrator> system(start, NO_SLOWDOWN);
    const Float etot = system.ar.getEtotRef();

    system.integrate(2.0*PI);

    CHECK(system.ar.getTime()==doctest::Approx(2.0*PI).epsilon(1e-12));
    for (int i=0; i<2; i++) 
        for (int k=0; k<2; k++) CHECK(abs(system.ar.particles[i].pos[k]-start[i].pos[k])<1e-2);
    CHECK(abs(system.ar.getEnergyError()/etot)<1e-2);
    // a fix of the two-body step shows up here
    CHECK(abs(system.ar.getEnergyError()/etot)>1e-6);
}

TEST_CASE_TEMPLATE("hierarchical triple conserves energy", V, ALL_VARIANTS) {
    ArSystem<typename V::Integrator> system(hierarchicalTriple(), NO_SLOWDOWN);
    const Float etot = system.ar.getEtotRef();

    system.integrate(TRIPLE_TIME, 5);

    CHECK(system.ar.getTime()==doctest::Approx(TRIPLE_TIME).epsilon(1e-12));
    CHECK(abs(system.ar.getEnergyError()/etot)<1e-8);
    // it did move: the inner binary went around eight times
    CHECK(system.ar.profile.step_count_sum>100);
}

TEST_CASE("the three Slow-down schemes agree when every Slow-down factor is 1") {
    ArSystem<LogHNone::Integrator>  none (hierarchicalTriple(), NO_SLOWDOWN);
    ArSystem<LogHArray::Integrator> array(hierarchicalTriple(), NO_SLOWDOWN);
    ArSystem<LogHTree::Integrator>  tree (hierarchicalTriple(), NO_SLOWDOWN);
    none.integrate(TRIPLE_TIME, 5);
    array.integrate(TRIPLE_TIME, 5);
    tree.integrate(TRIPLE_TIME, 5);

    // the schemes are different algorithms: agreement is within round-off amplified over eight inner orbits
    CHECK(none.distanceTo(array)<1e-8);
    CHECK(none.distanceTo(tree)<1e-8);

    SUBCASE("also with TTL") {
        ArSystem<TTLNone::Integrator>  none_ttl (hierarchicalTriple(), NO_SLOWDOWN);
        ArSystem<TTLArray::Integrator> array_ttl(hierarchicalTriple(), NO_SLOWDOWN);
        ArSystem<TTLTree::Integrator>  tree_ttl (hierarchicalTriple(), NO_SLOWDOWN);
        none_ttl.integrate(TRIPLE_TIME, 5);
        array_ttl.integrate(TRIPLE_TIME, 5);
        tree_ttl.integrate(TRIPLE_TIME, 5);
        CHECK(none_ttl.distanceTo(array_ttl)<1e-8);
        CHECK(none_ttl.distanceTo(tree_ttl)<1e-8);
    }
}

TEST_CASE("LogH and TTL agree on the same problem under each Slow-down scheme") {
    // two time transformations take different steps: agreement is within the integration error
    SUBCASE("no Slow-down") {
        ArSystem<LogHNone::Integrator> logh(hierarchicalTriple(), NO_SLOWDOWN);
        ArSystem<TTLNone::Integrator>  ttl (hierarchicalTriple(), NO_SLOWDOWN);
        logh.integrate(TRIPLE_TIME, 5);
        ttl.integrate(TRIPLE_TIME, 5);
        CHECK(logh.distanceTo(ttl)<1e-7);
    }
    SUBCASE("Array slow-down") {
        ArSystem<LogHArray::Integrator> logh(hierarchicalTriple(), NO_SLOWDOWN);
        ArSystem<TTLArray::Integrator>  ttl (hierarchicalTriple(), NO_SLOWDOWN);
        logh.integrate(TRIPLE_TIME, 5);
        ttl.integrate(TRIPLE_TIME, 5);
        CHECK(logh.distanceTo(ttl)<1e-7);
    }
    SUBCASE("Tree slow-down") {
        ArSystem<LogHTree::Integrator> logh(hierarchicalTriple(), NO_SLOWDOWN);
        ArSystem<TTLTree::Integrator>  ttl (hierarchicalTriple(), NO_SLOWDOWN);
        logh.integrate(TRIPLE_TIME, 5);
        ttl.integrate(TRIPLE_TIME, 5);
        CHECK(logh.distanceTo(ttl)<1e-7);
    }
    SUBCASE("Kepler orbit") {
        ArSystem<LogHNone::Integrator> logh(keplerBinary(), NO_SLOWDOWN);
        ArSystem<TTLNone::Integrator>  ttl (keplerBinary(), NO_SLOWDOWN);
        logh.integrate(1.0);
        ttl.integrate(1.0);
        CHECK(logh.distanceTo(ttl)<1e-8);
    }
}

TEST_CASE_TEMPLATE("Array slow-down: a weakly perturbed inner binary is slowed and the Slow-down energy stays within bound", V, LogHArray, TTLArray) {
    // the Slow-down timescale may reach the outer period
    ArSystem<typename V::Integrator> system(hierarchicalTriple(), 4.0, 1e-3);
    const Float etot_sd = system.ar.getEtotSlowDownRef();

    // the list holds the root and the inner binary
    REQUIRE(system.ar.getSlowDownInnerNumber()==2);
    REQUIRE(system.ar.binary_slowdown.getSize()==2);
    CHECK(system.ar.binary_slowdown[1]->slowdown.getSlowDownFactor()>1.0);

    system.integrate(TRIPLE_TIME, 5);

    CHECK(system.ar.binary_slowdown[1]->slowdown.getSlowDownFactor()>1.0);
    CHECK(abs(system.ar.getEnergyErrorSlowDown()/etot_sd)<1e-8);
}

TEST_CASE_TEMPLATE("Tree slow-down: the inner and the outer binary of a two-level hierarchy can both be slowed", V, LogHTree, TTLTree) {
    ArSystem<typename V::Integrator> system(hierarchicalQuadruple(), 1.0, 1e-3);
    const Float etot_sd = system.ar.getEtotSlowDownRef();

    // three binaries: the root (separation 10), the middle (0.1) and the inner (0.001)
    auto& tree = system.ar.info.binarytree;
    REQUIRE(system.ar.getSlowDownInnerNumber()==3);
    int n_slowed = 0;
    for (int i=0; i<tree.getSize(); i++) {
        if (tree[i].semi<1.0) {
            CHECK(tree[i].slowdown.getSlowDownFactor()>1.0);
            n_slowed++;
        }
    }
    CHECK(n_slowed==2);

    system.integrate(0.01, 2);
    CHECK(abs(system.ar.getEnergyErrorSlowDown()/etot_sd)<1e-8);
}

TEST_CASE_TEMPLATE("the number of inner slowed binaries is zero without Slow-down", V, LogHNone, TTLNone) {
    ArSystem<typename V::Integrator> system(hierarchicalTriple(), NO_SLOWDOWN);
    CHECK(system.ar.getSlowDownInnerNumber()==0);
}
