CXX = g++
OPTS = -lgdi32
SRC = $(wildcard *.cpp)
EXE = raytrace.exe
DEBUG_STR = debug_raytrace
DEBUG_EXE = $(DEBUG_STR).exe

.PHONY: clean build debug all
.DEFAULT: build

build: $(EXE)

debug: $(DEBUG_EXE)

all: build debug

clean:
	del $(EXE) $(DEBUG_STR).* *.pdb *.obj

$(EXE):
	g++ $(SRC) $(OPTS) -o $(EXE)
	
$(DEBUG_EXE):
	cl /Zi $(SRC) user32.lib gdi32.lib /Fe:$(DEBUG_EXE)
