# the compiler: gcc for C program
CC = gcc

#compiler flags:
# -g adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS = -w
DBGCFLAGS = -g
# -lm links the math library
#  -lz links the gz library
LIBS = -lm -lz
#sources
SOURCES = sarscov2_imputation.c math.c options.c 
HASHMAP = hashmap.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = sarscov2_imputation

all: $(TARGET)
$(TARGET): $(TARGET).c
	$(CC) -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)

clean:
	$(RM) $(TARGET)
