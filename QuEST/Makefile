ROOT = .
SRC_DIR = $(ROOT)/src/
MISC_DIR = $(ROOT)/misc/
BIN_DIR = $(ROOT)/bin/
#../misc_files/
CC = g++
#CCFLAGS = -Wall -ansi -pedantic -D DEBUG -O4 -funroll-loops -pthread
INCLUDE = $(MISC_DIR)
COPY = cp

include config.mk

LOCAL_TARGETS =	generate_profile peak_caller quick_window_scan calibrate_peak_shift generate_profile_views metrics profile_2_wig align_2_QuEST collapse_reads bin_align_2_bedgraph report_bad_control_regions

LOCAL_EXECS = $(LOCAL_TARGETS)
LOCAL_OBJECTS = $(LOCAL_TARGETS:=.o)
OTHER_DEPENDENCIES = params string_utils seq_contig sequence_utility_functions
OTHER_OBJECTS = $(OTHER_DEPENDENCIES:=.o)

EXPANDED_OTHER_OBJECTS = $(addprefix $(MISC_DIR), $(OTHER_OBJECTS))

default: all

.PHONY: default all clean


all: other_dependencies local_objects local_execs

place:
	mv $(LOCAL_TARGETS) $(BIN_DIR)

local_execs: $(LOCAL_EXECS)

local_objects:
	cd $(SRC_DIR) && $(MAKE) $(LOCAL_OBJECTS)

$(LOCAL_EXECS) : %: $(SRC_DIR)/%.o
	echo $(INCLUDE)
	$(CC) -I$(INCLUDE) $(CCFLAGS) $(SRC_DIR)$@.o $(EXPANDED_OTHER_OBJECTS) -o $@

other_dependencies: other_objects

other_objects:
	cd $(MISC_DIR) && $(MAKE) $(OTHER_OBJECTS)

clean:
	rm -rf $(SRC_DIR)*.o
	cd $(MISC_DIR) && $(MAKE) clean	

all_clean:
	rm -rf $(SRC_DIR)*.o
	cd $(MISC_DIR) && $(MAKE) clean
