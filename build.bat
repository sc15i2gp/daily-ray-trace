@echo off
IF NOT EXIST src (mkdir src)
IF NOT EXIST output (mkdir output)
IF NOT EXIST bin (mkdir bin)
del *.obj bin\*.ilk bin\*.exe bin\*.pdb
g++ src\main.cpp -static-libgcc -static-libstdc++ -lgdi32 -march=native -O3 -o bin\raytrace.exe
cl /DDEBUG_BUILD /Zi /MD src\main.cpp /Fe:bin\debug_raytrace.exe user32.lib gdi32.lib /link /debug:full
cl /Zi /MD /O2it src\main.cpp /Fe:bin\profile_raytrace.exe user32.lib gdi32.lib /link /debug:full
