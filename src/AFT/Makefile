address = ../ADT
CC = g++
CFLAGS = -c -std=c++11 -fpic -O3
LDFLAGS = -shared
#LDFLAGS += -L$(address) -Wl,-rpath=$(address)
SOURCES=AFT.cpp\
	Edge.cpp\
	Preset.cpp\
	Trimming.cpp\
	Triangle.cpp\
	Point.cpp\
	ReBlanking.cpp\
	Front.cpp\
	NewGrid.cpp\
	FinalGrid.cpp	
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=libAFT.so

all: $(TARGET)
	
$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(TARGET)
