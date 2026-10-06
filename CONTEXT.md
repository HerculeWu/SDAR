# SDAR

A library for integrating few-body gravitational systems, combining a Hermite integrator for the global system with algorithmic regularization (AR) for compact subsystems.

## Language

### Methods

**AR**:
Algorithmic regularization: integration of a few-body system with an explicit symplectic integrator applied to time-transformed equations of motion, so that close encounters and highly eccentric orbits need no special treatment.
_Avoid_: Symplectic integrator, TSI

**SDAR**:
The combination of **AR** with **Slow-down**; also the name of this library.

**LogH**:
The logarithmic Hamiltonian variant of **AR**, in which the time transformation during the drift is computed from the system's kinetic and binding energy.

**TTL**:
The variant of **AR** after Mikkola & Aarseth (2002), in which the time transformation during the drift is computed from an auxiliary variable integrated alongside the system instead of from the kinetic energy.

**AR dynamics**:
The equations of motion as **AR** advances them for one choice of time transformation (**LogH** or **TTL**) and one choice of **Slow-down** scheme (none, **Array slow-down** or **Tree slow-down**).

### Particles and groups

**Group**:
A set of particles integrated together by one AR integrator and represented to the Hermite system by its centre of mass.
_Avoid_: Subsystem, cluster

**Member**:
A particle inside a **Group**.

**Single**:
A particle integrated directly by Hermite, not belonging to any **Group**.

**Resolved group**:
A **Group** in the state where Hermite sees its **Members** individually instead of its centre of mass. It is a state of a **Group**, not a separate kind of thing.

### Time

**Physical time**:
The time that appears in the equations of motion, in whatever unit system the user has chosen.
_Avoid_: Real time

**Integrated time**:
**Physical time** measured from a shifted origin; the two differ only by a constant offset.
_Avoid_: Integration time

**Regularization time**:
The fictitious independent variable that replaces **Physical time** in the time-transformed equations of motion; AR advances in steps of it.
_Avoid_: Step, step size

**Regularization step**:
The increment of **Regularization time** advanced by one AR step. It names the increment, never the time variable itself.
_Avoid_: Step size, ds

**Slow-down time**:
The time by which a slowed-down binary's internal orbit has actually advanced, which is shorter than the elapsed **Physical time** by the **Slow-down factor**.
_Avoid_: Slowdown time

### Interrupts

**Interrupt**:
A physical phenomenon outside the gravitational equations of motion that intervenes on a binary during integration. It names the physics, not a pause in the integration: an interrupt may or may not stop the integrator.
_Avoid_: Interruption

**Interrupt cause**:
The physical process behind an **Interrupt**, such as a collision or a common envelope. The set is open-ended and defined by the user of the library.

**Interrupt outcome**:
What an **Interrupt** did to the binary: _change_ (members altered, binary survives), _merge_ (the two members become one) or _destroy_ (the binary no longer exists).

### Slow-down

**Slow-down**:
The method of artificially slowing the internal orbital motion of a weakly perturbed binary, so that fewer orbits need to be integrated while the long-term effect of the perturbation is kept.
_Avoid_: Slowdown, slow down

**Slow-down factor**:
How many times slower a binary's internal motion runs under **Slow-down**; 1 means no slow-down. Each binary has exactly one.
_Avoid_: Kappa, original slow-down factor

**Array slow-down**:
**Slow-down** applied only to the innermost binaries of a hierarchy: for particles a, b, c where a-b is a binary and its centre of mass ab forms a binary with c, only a-b is slowed.

**Tree slow-down**:
**Slow-down** applied at every level of a hierarchy: in the same system both a-b and ab-c are slowed.
