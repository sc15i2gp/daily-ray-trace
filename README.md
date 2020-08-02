# daily-ray-trace
Raytracer written in "C+".
The specs of the test environment are as follows:
- Windows 10
- make 3.81
- MinGW 8.1.0
- RemedyBG 0.3.0.7

To build, run the build.bat script. If using the command line, the command takes the following form:
- build <target>
Where "target" is either "debug", "all" or nothing (to choose the default option). The default uses MinGW to compile and the debug option uses Cl (MSVC). The debug option is used so that debug information can be generated for use with a debugger. The all option builds both of these.

Currently, the image's scene and parameters to the ray tracer are hard coded. The scene is the Cornell box with a spherical light and a spherical mirror. The image, on the test machine, takes ~20s per pixel sample.
