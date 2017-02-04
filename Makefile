# GCC = g++
GCC = $(CXX)

BIN = bin/minimizerindex

CC_FLAGS_RELEASE = -O3 -Wall -c -fmessage-length=0 -pthread -std=c++11 -lpthread
CC_FLAGS_DEBUG = -O0 -g -Wall -c -fmessage-length=0 -pthread -std=c++11 -lpthread
CC_LIBS = 
INCLUDE = -Isrc/

LD_FLAGS = -lpthread -pthread -Wl,--no-as-needed
LIB_DIRS = 
LD_LIBS =

SOURCE_DIR = src
OBJ_DEBUG_DIR = obj_debug
OBJ_RELEASE_DIR = obj_release
CC_FILES :=  $(wildcard $(SOURCE_DIR)/*.cc)
H_FILES := $(wildcard $(SOURCE_DIR)/*.h) $(wildcard $(SOURCE_DIR)/*.hpp)
OBJ_FILES_DEBUG := $(addprefix $(OBJ_DEBUG_DIR)/,$(CC_FILES:.cc=.o))
OBJ_FILES_RELEASE := $(addprefix $(OBJ_RELEASE_DIR)/,$(CC_FILES:.cc=.o))

all: debug

clean: cleanbuild cleantest


debug: $(OBJ_FILES_DEBUG) 
	@echo [LD DEBUG] $<
	@mkdir -p $(dir $(BIN))
	@$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN) $(OBJ_FILES_DEBUG) $(LD_LIBS)

$(OBJ_DEBUG_DIR)/%.o: %.cc $(H_FILES)
	@echo [CP DEBUG] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<

release: $(OBJ_FILES_RELEASE)
	@echo [LD RELEASE] $<
	@mkdir -p $(dir $(BIN))
	@$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN) $(OBJ_FILES_RELEASE) $(LD_LIBS)

$(OBJ_RELEASE_DIR)/%.o: %.cc $(H_FILES)
	@echo [CP RELEASE] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<

cleanbuild:
	rm -rf $(OBJ_RELEASE_DIR) $(OBJ_DEBUG_DIR) $(BIN)

#######################
### Compiling tests ###
#######################
BIN_TEST = $(BIN)_test

GTEST_DIR = lib/googletest
INCLUDE_GTEST = -I$(GTEST_DIR)/include -Itest/

TEST_SOURCE_DIR = test_src
OBJ_TEST_DIR = obj_test
# Remove the main source file from the SOURCE_DIR
TEST_CC_FILES :=  $(wildcard $(TEST_SOURCE_DIR)/*.cc) $(filter-out $(SOURCE_DIR)/main_test.cc, $(CC_FILES))
TEST_H_FILES := $(wildcard $(TEST_SOURCE_DIR)/*.h) $(wildcard $(TEST_SOURCE_DIR)/*.hpp) $(H_FILES)
OBJ_FILES_TEST := $(addprefix $(OBJ_TEST_DIR)/,$(TEST_CC_FILES:.cc=.o)) $(GTEST_DIR)/build/gtest.a
TEST_MACROS = -DRUN_ALL_TESTS_

test: $(GTEST_DIR)/build/gtest.a $(OBJ_FILES_TEST)
	@echo [LD TESTS] $<
	@mkdir -p $(dir $(BIN_TEST))
	@$(GCC) $(LD_FLAGS) $(LIB_DIRS) $(TEST_MACROS) -o $(BIN_TEST) $(OBJ_FILES_TEST) $(LD_LIBS)

$(OBJ_TEST_DIR)/%.o: %.cc $(TEST_H_FILES)
	@echo [CP TESTS] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(INCLUDE_GTEST) $(CC_FLAGS_DEBUG) $(TEST_MACROS) -o $@ $<

$(GTEST_DIR)/build/gtest.a:
	cd $(GTEST_DIR); make

cleantest:
	rm -rf $(OBJ_TEST_DIR) $(BIN_TEST)

#######################
