setlocal
SET PATH=%PATH%;C:\Program Files\mingw-w64\x86_64-8.1.0-win32-seh-rt_v6-rev0\mingw64\bin;C:\Program Files (x86)\GnuWin32\bin;%~dp0\bin;
start cmd /k "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
cd src
start cmd /k vim main.cpp
