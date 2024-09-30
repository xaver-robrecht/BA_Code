# Directories
SRC_DIR = src
INCLUDE_DIR = include
BIN_DIR = bin
BUILD_DIR = build
DATA_DIR = data

# Compiler and flags
CC := gcc
CFLAGS := -I $(INCLUDE_DIR) -O3 -march=native -flto -ftree-parallelize-loops=8
LDFLAGS := -lgsl -lgslcblas -lm

# Source files
SRC_FILES = $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC_FILES))

all: clean build
# Executable
EXEC = $(BIN_DIR)/main.exe

# Build target: Compiles the main executable from all object files
build: $(EXEC)

# Link all object files to create the main executable
$(EXEC): $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile each .c file in the src directory into an object file in the build directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

# Clean target: Removes all executables and object files
clean:
	rm -f $(BIN_DIR)/*.exe
	rm -f $(BUILD_DIR)/*.o
	rm -f $(TEST_DIR)/*.exe
	rm -f $(DATA_DIR)/*.txt

# Debug target: Adds debugging information and builds the main executable
debug: CFLAGS += -g
debug: build

run: clean build $(EXEC)
	for r in 4 6 8; do \
		for lambda in 1e-6 0.01 0.02; do \
			./bin/main.exe $$r 0.01 0 $$lambda; \
		done \
	done
	./bin/main.exe 6 0.01 0 0.1;
	./bin/main.exe 8 0.01 0 0.1;
# Declare phony targets to avoid conflicts with files of the same name
.PHONY: build clean test debug ass run
