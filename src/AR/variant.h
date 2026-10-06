#pragma once

#include <concepts>

namespace AR {

    //! Time transformation: logarithmic Hamiltonian, the transformation during the drift comes from the kinetic and binding energy
    struct LogH {};

    //! Time transformation: TTL (Mikkola & Aarseth 2002), the transformation during the drift comes from an integrated auxiliary variable
    struct TTL {};

    //! Slow-down scheme: no Slow-down
    struct NoSlowDown {};

    //! Slow-down scheme: Array slow-down, inner binaries listed in a flat array
    struct ArraySlowDown {};

    //! Slow-down scheme: Tree slow-down, every binary of the tree with nested Slow-down factors
    struct TreeSlowDown {};

    //! A time transformation of AR
    template <class T>
    concept TimeTransformation = std::same_as<T, LogH> || std::same_as<T, TTL>;

    //! A Slow-down scheme of AR
    template <class T>
    concept SlowDownScheme = std::same_as<T, NoSlowDown> || std::same_as<T, ArraySlowDown> || std::same_as<T, TreeSlowDown>;

    // ---- Compatibility shim: the variant macros select the default variant types.
    // This is the only place of the library where these macros are read.

#ifdef AR_TTL
    typedef TTL DefaultTimeTransformation;
#else
    typedef LogH DefaultTimeTransformation;
#endif

#if (defined AR_SLOWDOWN_ARRAY) && (defined AR_SLOWDOWN_TREE)
#error "AR_SLOWDOWN_ARRAY and AR_SLOWDOWN_TREE select different Slow-down schemes; define at most one of them"
#endif

#if defined AR_SLOWDOWN_ARRAY
    typedef ArraySlowDown DefaultSlowDownScheme;
#elif defined AR_SLOWDOWN_TREE
    typedef TreeSlowDown DefaultSlowDownScheme;
#else
    typedef NoSlowDown DefaultSlowDownScheme;
#endif

}
