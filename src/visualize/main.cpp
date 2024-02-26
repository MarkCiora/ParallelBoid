#include <iostream>
#include <fstream>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#ifdef _WIN32
  #include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

const int WIDTH = 800;
const int HEIGHT = 600;
const int FPS = 30;
const char* OUTPUT_FILE = "output.mp4";

int main(int argv, char **argc){

}