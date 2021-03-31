@echo off
IF NOT EXIST src (mkdir src)
IF NOT EXIST output (mkdir output)
IF NOT EXIST bin (mkdir bin)
del *.obj bin\*.ilk bin\*.exe bin\*.pdb
g++ src\main.cpp -lgdi32 -o bin\raytrace.exe
cl /DDEBUG_BUILD /F 4194304 /Zi src\main.cpp /Fe:bin\debug_raytrace.exe user32.lib gdi32.lib
