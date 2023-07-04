COMPILER = g++-10
FLAGS = -O3 -std=c++17 -Wall -Wfatal-errors -fopenmp

CPP := main.cpp model.cpp solving_algorithms.cpp tw_preprocessor.cpp local_search.cpp
REQ := $(CPP) arg_parser.hpp bitsets.hpp constants.h cost_models.h instances.h local_search.h model.h msa.hpp scc.hpp\
 solving_algorithms.h static_array.hpp tw_preprocessor.h util.hpp

CPP := $(addprefix src/,$(CPP))
REQ := $(addprefix src/,$(REQ))

sizes = 128 192 256 320 384 448 512

TDTSPTW = tdtsptw $(foreach size,$(sizes),tdtsptw$(size))

ifneq ($(strip $(MARCH)),)
FLAGS := $(FLAGS) -march=$(MARCH)
BIN_DIR := bin-$(MARCH)/
TDTSPTW := $(addprefix $(BIN_DIR),$(TDTSPTW))
$(shell mkdir -p $(BIN_DIR))
RM_BIN := rm -rf $(BIN_DIR)
endif

small: $(firstword $(TDTSPTW))
large: $(lastword $(TDTSPTW))
all: $(TDTSPTW)

.PHONY: small large all clean

$(TDTSPTW): $(BIN_DIR)tdtsptw%: $(REQ)
	$(COMPILER) $(CPP) $(FLAGS) -o $@ $(if $*,-DMAX_INSTANCE_SIZE=$*,) -DNDEBUG

clean:
	rm -f $(TDTSPTW)
	$(RM_BIN)
