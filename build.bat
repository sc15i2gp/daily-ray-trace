@echo off
IF NOT EXIST src (mkdir src)
IF NOT EXIST output (mkdir output)
IF NOT EXIST bin (mkdir bin)
del *.obj bin\*.ilk bin\*.exe bin\*.pdb
gcc src\main.c -Isrc -o bin\raytrace.exe
cl /DDEBUG_BUILD /Zi /TC src\main.c /Fe:bin\debug_raytrace.exe /link /debug:full
cl /Zi /TC src\main.c /Fe:bin\profile_raytrace.exe /link /debug:full
