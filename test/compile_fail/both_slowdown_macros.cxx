// Must not compile: the two Slow-down macros select different schemes
// expect: define at most one of them
#define AR_SLOWDOWN_ARRAY
#define AR_SLOWDOWN_TREE
#include "AR/variant.h"
