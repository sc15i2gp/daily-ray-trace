
@echo off
IF "%~1%" == "D" (
cl /Zi main.cpp user32.lib gdi32.lib
) ELSE (
g++ -g -march=native -mavx -mavx2 -msse -msse2 -msse3 Maths.cpp main.cpp -lgdi32 -o raytrace.exe)

