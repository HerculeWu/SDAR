# Step-control trace fixtures

Each `.trace` file records how the integration loop of `AR::TimeTransformedSymplecticIntegrator::integrateToTime`
chose the Regularization step, before that logic was moved into the Regularization step control module.
Replaying a record through the module must reproduce every Regularization step, every save/restore decision
and every value of the persistent Regularization step exactly.

## Format

Plain text, one event per line, fields separated by one space. Floating-point values are C99 hexadecimal
(`%a`), which `strtod` reads back exactly. A file holds one or more records; a record is one `integrateToTime` call.

| Line | Meaning |
|---|---|
| `BEGIN ds opt t0 t_end t_err e_err n_max order cd` | Start of a call. `ds`: persistent Regularization step on entry. `opt`: step option (0 always, 1 later, 2 none). `t0`: Physical time on entry. `t_end`: target Physical time. `t_err`: time synchronization error limit. `e_err`: relative energy error limit. `n_max`: maximum step count. `order`: symplectic order. `cd`: size of the sub-step time table (order 6 with 8 entries and order 8 with 16 entries are the second Yoshida solution, i.e. `initialSymplecticCofficients(-order)`). |
| `C` | A checkpoint begins: the loop is about to save state. The lines up to the next `S`, `X` are what happened inside it. |
| `O opt` | Inside a checkpoint: the step option changed to `opt`. |
| `R I ds_est` | Inside a checkpoint: restart requested after an Interrupt, with the Kepler-based estimate `ds_est`. |
| `R B ds_est` | Inside a checkpoint: restart requested after a binary update, with estimate `ds_est`. |
| `X` | Inside a checkpoint: the integrator stopped (early return caused by an Interrupt). `END` follows. |
| `S` | State saved (end of a checkpoint). |
| `T` | State restored instead of a checkpoint: the previous step is retried. Physical time returns to its value at the last `S`. |
| `P ds t err [tt_0 ... tt_{cd-1}]` | One step taken with Regularization step `ds`. `t`: Physical time after the step. `err`: integration error, the absolute change of the extended Hamiltonian over the step. The time advance is `t` minus the Physical time before the step. The sub-step time table `tt` (in storage order, not sorted) is present only when `t > t_end + t_err`, the one case in which it is read. |
| `D ds` | The persistent Regularization step (the value callers see between calls) became `ds`. Written at the point of the loop where the change is visible: before the next `C`/`T`, before `S`, or before `END`. |
| `END status ds opt n_step n_step_tsyn` | End of the call. `status`: `reach` (target time reached) or `stop`. `ds`, `opt`: persistent Regularization step and step option on exit. `n_step`, `n_step_tsyn`: steps taken and time-synchronization steps. |

Between two `P` lines there is exactly one of: `C ... S` (the step was accepted) or `T` (it was rejected).

## Files

| File | Source run | Contains |
|---|---|---|
| `binary-update.trace` | `ar.ttl.sd.t -S -t 2e-6 -n 20 triple.unstable`, calls 0-3 | restart after a binary update, error-driven reduction in the first steps and later, recovery of a reduced step, overshoot |
| `array-slow-down.trace` | `ar.logh.sd.a -S -t 2e-6 -n 20 triple.unstable`, calls 0-1 | the same under Array slow-down |
| `growth.trace` | `ar.ttl -S --fix-step-option none -e 1e-7 --n-step-max 2000 -t 2e-6 -n 4 triple.unstable`, call 0 | growth of the Regularization step, reduction, recovery |
| `interrupt.trace` | `ar.logh -S -i 1 -t 0.0002 -n 4 ../../golden/input/triple.collide`, all calls | restart after an Interrupt, step-option change |
| `merge-stop.trace` | `ar.logh.sd.a -S -i 2 -t 0.002 -n 20 ../../golden/input/binary.collide` | early stop by an Interrupt |
| `hermite.trace` | `hermite --dt-max-power 6 --dt-min-power 34 -o 3 -t 0.5 ../../../sample/input/triple.stable.lowm3`, calls 0-7 | short calls driven by Hermite, reduction in the first steps, overshoot |

`triple.unstable` is the input of the first three runs: the inner binary of the sample triple with a
close, massive third body.

The traces were recorded once, with Apple clang 17, by a temporary hook compiled with `-D AR_STEP_TRACE`
that appended to the file named by the environment variable `AR_STEP_TRACE_FILE`.
