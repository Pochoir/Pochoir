CC := icpc
CXX := icpc
CFLAGS = -Wall -m64 
CXXFLAGS := -Wall -m64
LDFLAGS = -DBUILD_64 -lrt
CFLAGS_RELEASE = -O3 -DNDEBUG
CXXFLAGS_RELEASE = -O3 -DNDEBUG
CFLAGS_DEBUG = -g -DDEBUG -O0
CXXFLAGS_DEBUG = -g -DDEBUG -O0

ifeq ($(DEBUG),1)
CFLAGS := $(CFLAGS_DEBUG) $(CFLAGS)
CXXFLAGS := $(CXXFLAGS_DEBUG) $(CXXFLAGS)
else
CFLAGS := $(CFLAGS_RELEASE) $(CFLAGS)
CXXFLAGS := $(CXXFLAGS_RELEASE) $(CXXFLAGS)
endif

TARGET = rbq_simple
COMMON_HEADERS := ktiming.h rmq_2d_simple.hpp

OBJS := \
	ktiming.o \
    test_rbq_2d.o \

.PHONY: all
all: $(TARGET)

rbq_simple: $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

# pattern rule for building objects
%.o: %.cxx %.h $(COMMON_HEADERS) .buildmode Makefile
	$(CXX) $(CXXFLAGS) -c $< -o $@
%.o: %.c %.h $(COMMON_HEADERS) .buildmode Makefile
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJS)
