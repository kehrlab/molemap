TARGET = molemap
BUILD_DIR = ./build
SRC_DIR = ./src

SRCS := $(shell find $(SRC_DIR) -type f -name "*.cpp")
OBJS := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SRCS:.cpp=.o))

SEQAN_LIB =.

# Compiler
CXX = g++ -std=c++14
CC = $(CXX)

CXXFLAGS += -O3 -I $(SEQAN_LIB) -DSEQAN_HAS_ZLIB=1 #-DSEQAN_DISABLE_VERSION_CHECK -std=c++14

# # Date and version number from git
DATE := on $(shell git log --pretty=format:"%cd" --date=iso | cut -f 1,2 -d " " | head -n 1)
VERSION := 1.3.0-$(shell git log --pretty=format:"%h" --date=iso | head -n 1)
CXXFLAGS += -DDATE=\""$(DATE)"\" -DVERSION=\""$(VERSION)"\"

# Enable warnings
# CXXFLAGS += -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result

# DEBUG build
#CXXFLAGS += -g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
# CXXFLAGS += -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0

LDLIBS = -lrt -lpthread -fopenmp -fPIE -lz

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDLIBS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -f $(OBJS) $(TARGET)
