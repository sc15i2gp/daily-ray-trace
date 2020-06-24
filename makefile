CXX = g++
OPTS = -lgdi32
EXE_DIR = bin
SRC_DIR = src
OUT_DIR = output
SRC = $(wildcard $(SRC_DIR)/*.cpp)
EXE = $(EXE_DIR)\raytrace.exe
DEBUG_STR = debug_raytrace
DEBUG_EXE = $(EXE_DIR)/$(DEBUG_STR).exe

$(info [${SRC}])


.PHONY: clean build debug all
.DEFAULT: build

build: $(EXE)

debug: $(DEBUG_EXE)

all: build debug

clean:
	del $(EXE_DIR)\*.exe $(DEBUG_STR).* *.pdb *.obj

$(EXE):
	g++ $(SRC) $(OPTS) -o $(EXE)
	
$(DEBUG_EXE):
	cl /Zi $(SRC) /Fe:$(DEBUG_EXE) user32.lib gdi32.lib 
