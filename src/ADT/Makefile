CC = g++ -g
CFLAGS = -c -std=c++11 -fpic -O3
LDFLAGS = -shared
SOURCES=ADT.cpp	
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=libADT.so

all: $(TARGET)
	
$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(TARGET)
