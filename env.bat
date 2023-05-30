setlocal
SET PATH=%PATH%;C:\cygwin64\bin;%~dp0\bin;
start cmd /k "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
cd src
start cmd /k vim daily_ray_trace.c
