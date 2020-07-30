TARGET = countK
BUILD_DIR = ./build
SRC_DIR = ./src

SRCS := $(shell find $(SRC_DIR) -type f -name *.cpp)
OBJS := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SRCS:.cpp=.o))

# Compiler
CXX = g++
CC = $(CXX)

# Set this to include SeqAn libraries, either system wide
# or download into current folder and set to .
SEQAN_LIB = ./src
CXXFLAGS += -I$(SEQAN_LIB) #-DSEQAN_HAS_ZLIB=1 -DSEQAN_DISABLE_VERSION_CHECK -std=c++14

# # Date and version number from git
# DATE := on $(shell git log --pretty=format:"%cd" --date=iso | cut -f 1,2 -d " " | head -n 1)
# VERSION := 0.0.1-alpha-$(shell git log --pretty=format:"%h" --date=iso | head -n 1)
# CXXFLAGS += -DDATE=\""$(DATE)"\" -DVERSION=\""$(VERSION)"\"

# Enable warnings
# CXXFLAGS += -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result

# DEBUG build
#CXXFLAGS += -g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
# CXXFLAGS += -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0

LDLIBS = -lz -lpthread -lrt

all: $(TARGET)

# $(TARGET): countKmers.cpp $(BUILD_DIR)/functions.o #$(SRC_DIR)/functions.h
# 	$(CC) countKmers.cpp $(BUILD_DIR)/functions.o -o $@ $(LDLIBS)
#
# $(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h
# 	$(CC) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDLIBS)

# $(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
