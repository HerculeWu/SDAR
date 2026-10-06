## Unit tests and golden-output comparison
#  make          build and run the unit tests (fails when a test fails), and check that
#                the sources in compile_fail/ are rejected by the compiler
#  make golden   run the golden-output comparison of the sample programs
#  make all      both

TARGET=unit_test

VPATH=../src

SOURCE=${shell ls ../src/*/*.h} ${shell ls unit/*.h 2>/dev/null} ../sample/AR/particle.h ../sample/AR/perturber.h
TEST_SOURCE=${shell ls unit/*.cxx}
TEST_OBJECT=$(TEST_SOURCE:.cxx=.o)

## compiler-------------------------------------------#
CXX= g++
CXXFLAGS += -O2 -Wall

## Flag ----------------------------------------------#
CXXFLAGS += -I../src -I./third_party -I./unit -I../sample/AR
CXXFLAGS += -std=c++20
CXXFLAGS += -D FIXTURE_DIR='"$(CURDIR)/fixtures"'

## Target --------------------------------------------#

test: $(TARGET) compile_fail
	./$(TARGET)

# sources that must be rejected by the compiler, with the message named in their "// expect:" line
compile_fail: ${SOURCE}
	@for f in compile_fail/*.cxx; do \
		expect=`sed -n 's|^// expect: ||p' $$f`; \
		if $(CXX) $(CXXFLAGS) -fsyntax-only $$f 2>compile_fail.log; then echo "compile_fail: $$f compiles, but must not"; exit 1; fi; \
		if ! grep -q "$$expect" compile_fail.log; then echo "compile_fail: $$f fails without the message '$$expect'"; cat compile_fail.log; exit 1; fi; \
		echo "compile_fail: $$f rejected as expected ($$expect)"; \
	done; rm -f compile_fail.log

all: test golden

golden:
	./golden/golden.py

$(TARGET): main.o $(TEST_OBJECT)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(CXXLIBS)

main.o: main.cxx third_party/doctest.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

unit/%.o: unit/%.cxx third_party/doctest.h ${SOURCE}
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) main.o unit/*.o

.PHONY: test all golden clean compile_fail
