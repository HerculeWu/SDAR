#!/usr/bin/env python3
"""Golden-output comparison for the sample programs.

Builds every AR sample variant and the Hermite sample, runs each on fixed
inputs and compares column output (stdout), binary dumps (.last) and the
parameter files (.par) byte-for-byte with the recorded files in expected/.
Text on the error stream is not compared.

    test/golden/golden.py            build, run, compare (exit 1 on any difference)
    test/golden/golden.py --record   build, run, overwrite expected/ and the timing baseline

Golden files are bitwise results of one compiler at the sample Makefiles'
flags; expected/COMPILER names it. See README.md beside this script.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
INPUT = os.path.join(HERE, "input")
EXPECTED = os.path.join(HERE, "expected")
AR_DIR = os.path.join(ROOT, "sample", "AR")
HERMITE_DIR = os.path.join(ROOT, "sample", "Hermite")

TRIPLE = os.path.join(ROOT, "sample", "input", "triple.stable.lowm3")

# AR sample binaries: ar.<variant>
AR_VARIANTS = ["logh", "logh.sd.a", "logh.sd.t", "ttl", "ttl.sd.a", "ttl.sd.t"]

# Regularization step of the "step" case, needed to continue from its dump
# because the sample does not restore it on load.
RESTART_DS = "2.3651431997916e-07"


class Case:
    """One run of one binary: copy `source` to `<name>` in the work directory, run on it."""

    def __init__(self, name, args, source, dump=True, source_is_dump_of=None):
        self.name = name
        self.args = args
        self.source = source
        self.dump = dump
        self.source_is_dump_of = source_is_dump_of

    def outputs(self):
        files = [self.name + ".out", self.name + ".par"]
        if self.dump:
            files.append(self.name + ".last")
        return files


def ar_cases(variant):
    cases = [
        # fixed Regularization step, hierarchical triple, no synchronization
        Case("step", ["-n", "200"], TRIPLE),
        # time-end synchronization and Slow-down over many inner orbits
        Case("sync", ["-S", "-t", "0.01", "-n", "100"], TRIPLE),
        # Interrupt that modifies the inner binary without stopping the integration
        Case("interrupt-modify", ["-S", "-i", "1", "-t", "0.002", "-n", "20"],
             os.path.join(INPUT, "triple.collide")),
        # Interrupt that merges a two-body Group and returns early
        Case("interrupt-merge", ["-S", "-i", "2", "-t", "0.002", "-n", "20"],
             os.path.join(INPUT, "binary.collide")),
    ]
    # Reading a dump back. The Array slow-down sample aborts on load today
    # (its list of slowed binaries is not rebuilt), so it has no restart case.
    if not variant.endswith("sd.a"):
        cases.insert(1, Case("restart", ["-l", "-s", RESTART_DS, "-n", "100"], None,
                             source_is_dump_of="step"))
    return cases


def hermite_cases():
    return [
        # all three particles in one Group: Hermite drives a hierarchical AR Group with Slow-down
        Case("triple", ["--dt-max-power", "6", "--dt-min-power", "34", "-o", "3", "-t", "32"],
             TRIPLE, dump=False),
        # small break radius: Groups are broken and re-formed thousands of times
        Case("groups", ["--dt-max-power", "10", "--dt-min-power", "30", "-o", "10", "-t", "0.02",
                        "-r", "1e-4"],
             os.path.join(INPUT, "triple.inner-group"), dump=False),
    ]


def targets():
    """(label, binary path, cases) for everything compared."""
    result = [("ar." + v, os.path.join(AR_DIR, "ar." + v), ar_cases(v)) for v in AR_VARIANTS]
    result.append(("hermite", os.path.join(HERMITE_DIR, "hermite"), hermite_cases()))
    return result


def compiler_id(cxx):
    out = subprocess.run([cxx, "--version"], capture_output=True, text=True).stdout
    return out.splitlines()[0].strip() if out else "unknown"


def build(cxx, jobs):
    make_args = ["-j", str(jobs)] + (["CXX=" + cxx] if cxx != "g++" else [])
    for directory, make_targets in ((AR_DIR, ["ar." + v for v in AR_VARIANTS]),
                                    (HERMITE_DIR, ["hermite"])):
        proc = subprocess.run(["make", "-C", directory] + make_args + make_targets,
                              capture_output=True, text=True)
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout + proc.stderr)
            sys.exit("golden: build failed in " + os.path.relpath(directory, ROOT))


def run_case(binary, case, workdir):
    """Run one case; return (seconds, error message or None). Outputs land in workdir."""
    data = os.path.join(workdir, case.name)
    if case.source_is_dump_of:
        shutil.copyfile(os.path.join(workdir, case.source_is_dump_of + ".last"), data)
    else:
        shutil.copyfile(case.source, data)
    with open(data + ".out", "wb") as out, open(data + ".err", "wb") as err:
        start = time.perf_counter()
        proc = subprocess.run([binary] + case.args + [case.name], cwd=workdir,
                              stdout=out, stderr=err)
        seconds = time.perf_counter() - start
    if proc.returncode != 0:
        with open(data + ".err", "rb") as err:
            tail = err.read().decode(errors="replace").splitlines()[-3:]
        return seconds, "exit status %d: %s" % (proc.returncode, " | ".join(tail))
    return seconds, None


def first_difference(a, b):
    with open(a, "rb") as fa, open(b, "rb") as fb:
        da, db = fa.read(), fb.read()
    if da == db:
        return None
    n = min(len(da), len(db))
    offset = next((i for i in range(n) if da[i] != db[i]), n)
    return "first difference at byte %d (sizes %d and %d)" % (offset, len(db), len(da))


def read_timing(path):
    timing = {}
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    label, seconds = line.split()
                    timing[label] = float(seconds)
    return timing


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--record", action="store_true",
                        help="overwrite the golden files and the timing baseline")
    parser.add_argument("--repeat", type=int, default=1,
                        help="run every case N times and report the fastest (default 1)")
    parser.add_argument("--no-build", action="store_true", help="use the binaries as they are")
    parser.add_argument("--keep", action="store_true", help="keep the work directory")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()

    cxx = os.environ.get("CXX", "g++")
    if not args.no_build:
        build(cxx, args.jobs)

    compiler = compiler_id(cxx)
    compiler_file = os.path.join(EXPECTED, "COMPILER")
    timing_file = os.path.join(EXPECTED, "TIMING")
    if args.record:
        os.makedirs(EXPECTED, exist_ok=True)
        with open(compiler_file, "w") as f:
            f.write(compiler + "\n")
    elif not os.path.exists(compiler_file):
        sys.exit("golden: no golden output recorded; run with --record on a known-good tree")
    else:
        with open(compiler_file) as f:
            recorded = f.read().strip()
        if recorded != compiler:
            print("golden: WARNING compiler differs from the one that recorded the golden output\n"
                  "  recorded: %s\n  current : %s\n"
                  "  differences below may come from the compiler, not the code" % (recorded, compiler))

    baseline = read_timing(timing_file)
    work_root = tempfile.mkdtemp(prefix="sdar-golden-")
    failures = []
    timing = {}
    n_files = 0
    for label, binary, cases in targets():
        workdir = os.path.join(work_root, label)
        expected_dir = os.path.join(EXPECTED, label)
        os.makedirs(workdir)
        if args.record:
            shutil.rmtree(expected_dir, ignore_errors=True)
            os.makedirs(expected_dir)
        total = 0.0
        for case in cases:
            best = None
            error = None
            for _ in range(max(1, args.repeat)):
                seconds, error = run_case(binary, case, workdir)
                best = seconds if best is None else min(best, seconds)
                if error:
                    break
            total += best
            if error:
                failures.append("%s %s: %s" % (label, case.name, error))
                continue
            for name in case.outputs():
                produced = os.path.join(workdir, name)
                golden = os.path.join(expected_dir, name)
                n_files += 1
                if args.record:
                    shutil.copyfile(produced, golden)
                elif not os.path.exists(golden):
                    failures.append("%s %s: no golden file" % (label, name))
                else:
                    difference = first_difference(produced, golden)
                    if difference:
                        failures.append("%s %s: %s" % (label, name, difference))
        timing[label] = total

    print("%-14s %10s %10s %8s" % ("variant", "time[s]", "baseline", "ratio"))
    for label, seconds in timing.items():
        if label in baseline and baseline[label] > 0:
            print("%-14s %10.2f %10.2f %8.3f" % (label, seconds, baseline[label],
                                                 seconds / baseline[label]))
        else:
            print("%-14s %10.2f %10s %8s" % (label, seconds, "-", "-"))
    if baseline and not args.record:
        print("%-14s %10.2f %10.2f %8.3f" % ("total", sum(timing.values()), sum(baseline.values()),
                                             sum(timing.values()) / sum(baseline.values())))

    if args.record and not failures:
        with open(timing_file, "w") as f:
            f.write("# seconds per variant when the golden output was recorded (machine-specific)\n")
            for label, seconds in timing.items():
                f.write("%s %.3f\n" % (label, seconds))

    if args.keep or failures:
        print("golden: work directory kept at " + work_root)
    else:
        shutil.rmtree(work_root)

    if failures:
        print("golden: FAILED")
        for failure in failures:
            print("  " + failure)
        return 1
    print("golden: %s %d files" % ("recorded" if args.record else "OK,", n_files)
          + ("" if args.record else " identical"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
