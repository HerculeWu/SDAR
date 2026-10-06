# Golden-output comparison

One command proves that a change to the library altered no numerical result:

    test/golden/golden.py

It builds every AR sample variant (`sample/AR/ar.*`) and the Hermite sample at the
sample Makefiles' flags, runs each on fixed inputs, and compares the column output,
the binary dumps (`.last`) and the parameter files (`.par`) byte-for-byte with the
files in `expected/`. It exits non-zero and names the variant and file when any byte
differs. Text on the error stream is not compared.

## Cases

Per AR variant (`expected/ar.<variant>/`):

| Case | What it exercises |
|---|---|
| `step` | fixed Regularization step on the hierarchical triple |
| `restart` | reading the `step` dump back and continuing (not for Array slow-down, whose sample aborts on load today) |
| `sync` | time-end synchronization over many inner orbits; Slow-down in the `sd` variants; carries the run time |
| `interrupt-modify` | an Interrupt that modifies the inner binary without stopping |
| `interrupt-merge` | an Interrupt that merges a two-body Group and returns early |

Hermite sample (`expected/hermite/`, built TTL with Tree slow-down):

| Case | What it exercises |
|---|---|
| `triple` | one hierarchical Group driven by Hermite, with Slow-down |
| `groups` | a small break radius, so Groups are broken and re-formed thousands of times |

## Run time

The run time per variant is printed beside the baseline in `expected/TIMING`, which is
the time on the machine that recorded the golden output. Compare ratios on that machine
only; use `--repeat 3` to damp noise.

## Recording

    test/golden/golden.py --record

overwrites `expected/`. Do this only on a tree whose results are known to be right.
The golden files are the bitwise results of one compiler, named in `expected/COMPILER`;
with another compiler the script warns and differences are expected. Set `CXX` to
choose the compiler.
