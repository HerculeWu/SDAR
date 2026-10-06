# C++20 is the minimum language standard

The library was written in C++11 and selected its AR variants (LogH or TTL, and no, Array or Tree slow-down) with preprocessor macros, which forked the integrator's code and state so that no two variants could be compiled or tested together. We raise the minimum standard to C++20 so those variants can be expressed as compile-time types checked by concepts, accepting that downstream codes must build with a C++20 compiler.

## Consequences

The variant macros (`AR_TTL`, `AR_SLOWDOWN_ARRAY`, `AR_SLOWDOWN_TREE`) keep working as a thin compatibility shim that picks the default variant types, so existing downstream builds and interaction classes compile unchanged.
