# the compiler: gcc for C program
CC = gcc

#compiler flags:
# -g adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS = -w
DBGCFLAGS = -g
#  -lz links the gz library
LIBS = -lz
#sources
SOURCES = remove_redundant.c options_filter.c 
HASHMAP = hashmap.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = remove_redundant

all: $(TARGET)
$(TARGET): $(TARGET).c
	$(CC) -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)

clean:
	$(RM) $(TARGET)
