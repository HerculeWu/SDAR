#include <iostream>
#include <cmath>
#include <cassert>
#include "doctest.h"

#define ASSERT(x) assert(x)

#include "Common/Float.h"
#include "AR/symplectic_step.h"

TEST_CASE("step-modification factor and error ratio are inverse of each other") {
    AR::SymplecticStep step;
    step.initialSymplecticCofficients(-6);
    const AR::SymplecticStep& table = step;

    CHECK(table.getOrder() == 6);
    CHECK(table.calcErrorRatioFromStepModifyFactor(0.5) == doctest::Approx(std::pow(0.5, 6)));
    Float factor = table.calcStepModifyFactorFromErrorRatio(table.calcErrorRatioFromStepModifyFactor(0.5));
    CHECK(factor == doctest::Approx(0.5));
}
