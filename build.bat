@echo off
IF NOT EXIST src (mkdir src)
IF NOT EXIST output (mkdir output)
IF NOT EXIST bin (mkdir bin)
del *.obj bin\*.ilk bin\*.exe bin\*.pdb
g++ src\main.cpp -static-libgcc -static-libstdc++ -lgdi32 -O3 -DPROFILE_TIMINGS -o bin\raytrace.exe
cl /DDEBUG_BUILD /F 4194304 /Zi /MD src\main.cpp /Fe:bin\debug_raytrace.exe user32.lib gdi32.lib /link /debug:full
cl /F 4194304 /Zi /MD /O2it src\main.cpp /Fe:bin\profile_raytrace.exe user32.lib gdi32.lib /link /debug:full
