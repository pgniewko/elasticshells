#include <fstream>
#include <cstdlib>
#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>
#include <list>
#include <cmath>

using namespace std;

// g++ main.cpp -lm -lGL -lGLU -lglut -lSDL
// can represent vector too
struct Point {
    Point operator +(const Point& p2) {
        return Point(x+p2.x, y+p2.y, z+p2.z);
    }
    Point operator *(float r) {
        return Point(x*r, y*r, z*r);
    }
    Point(float a, float b, float c):x(a),y(b),z(c) {
    }
    void setLength(float r) {
        float rl = r/length();
        x *= rl;
        y *= rl;
        z *= rl;
    }
    float length() {
        return sqrt(x*x+y*y+z*z);
    }
    float x, y, z;
};

struct Triangle {
    Triangle(Point m, Point n, Point o):a(m),b(n),c(o) {};
    Point a, b, c;
};

std::list<Triangle> tris;

float rx = 0;
float ry = 0;

void draw()
{
    glLoadIdentity();

    glTranslatef(0, 0, -5);
    glRotatef(rx, 1, 0, 0);
    glRotatef(ry, 0, 0, 1);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glBegin(GL_TRIANGLES);
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) {
        glVertex3fv(&i->a.x);
        glVertex3fv(&i->b.x);
        glVertex3fv(&i->c.x);
    }
    glEnd();
    SDL_GL_SwapBuffers();
}

void resizeGL(int width, int height)
{
	glViewport(0,0,(GLsizei)(width),(GLsizei)(height));
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0f,(GLfloat)(width)/(GLfloat)(height),1.0f,100.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	return;
}

void subdivide()
{
    std::list<Triangle> newTris;
    float l = tris.begin()->a.length();
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) { // go through all triangles
        Point mid = (i->a + i->b) * 0.5f; // point betwenn points A and B
        mid.setLength(l); // put in on the sphere
        newTris.push_back(Triangle(i->b, i->c, mid)); // remember new triangles
        newTris.push_back(Triangle(i->a, i->c, mid));
    }
    tris.swap(newTris); // use new set of triangles;
}

void writeSurface(const char* filename, bool wrap)
{
    ofstream os(filename);
    
    int index = 1;
    os << 3*(int)tris.size() << "\n" ;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) // go through all triangles
    {
        Point p1 = i->a;
        Point p2 = i->b;
        Point p3 = i->c;

        if (wrap) /* Place all particles in the central box*/
        {
         // do write
        }
        else /* Use particles absolute positions*/
        {
           os << "H" << index++ << " "<< i->a.x << " " << i->a.y << " " << i->a.z << "\n";
//           os << index++ << " " << p1.x << " " << p1.y << " " << p1.z << "\n";
           os << "H" << index++ << " "<<i->b.x << " " << i->b.y << " " << i->b.z << "\n";
//           os << index++ << " " << p2.x << " " << p2.y << " " << p2.z << "\n";
           os << "H" << index++ << " "<< i->c.x << " " << i->c.y << " " << i->c.z << "\n";
//           os << index++ << " " << p3.x << " " << p3.y << " " << p3.z << "\n";
        }

    }

    os.close();
}

void renderingScript(const char* filename)
{
    ofstream os(filename);
    
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    os << "cmd.do(\"load cells.xyz, cells\")\n";
    os << "cmd.do(\"hide all\")\n";
    os << "cmd.do(\"set sphere_color, tv_red\")\n";
    os << "cmd.do(\"set line_color, marine\")\n";
    os << "cmd.do(\"show spheres\")\n";
    os << "cmd.do(\"alter elem h, vdw=0.1\")\n";
    os << "cmd.do(\"rebuild\")\n";


    int index = 1;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        os << "cmd.do(\"bond /cells///UNK`/H"<< index+0 << ", /cells///UNK`/H" << index+1 << "\")\n";
        os << "cmd.do(\"bond /cells///UNK`/H"<< index+1 << ", /cells///UNK`/H" << index+2 << "\")\n";
        os << "cmd.do(\"bond /cells///UNK`/H"<< index+2 << ", /cells///UNK`/H" << index+0 << "\")\n";
        index = index + 3;
    }
    
    os << "cmd.do(\"show lines\")\n";
    os.close();
}

void createCube()
{
    // create cube
    //top
    tris.push_back(Triangle(Point(1, 1, 1), Point(-1, 1, -1), Point(1, 1, -1)));
    tris.push_back(Triangle(Point(-1, 1, -1), Point(1, 1, 1), Point(-1, 1, 1)));
    //bottom
    tris.push_back(Triangle(Point(1, -1, 1), Point(-1, -1, -1), Point(1, -1, -1)));
    tris.push_back(Triangle(Point(-1, -1, -1), Point(1, -1, 1), Point(-1, -1, 1)));

    //right
    tris.push_back(Triangle(Point(1, 1, 1), Point(1, -1, -1), Point(1, 1, -1)));
    tris.push_back(Triangle(Point(1, -1, -1), Point(1, 1, 1), Point(1, -1, 1)));
    //left
    tris.push_back(Triangle(Point(-1, 1, 1), Point(-1, -1, -1), Point(-1, 1, -1)));
    tris.push_back(Triangle(Point(-1, -1, -1), Point(-1, 1, 1), Point(-1, -1, 1)));

    //front
    tris.push_back(Triangle(Point(1, 1, 1), Point(-1, -1, 1), Point(1, -1, 1)));
    tris.push_back(Triangle(Point(-1, -1, 1), Point(1, 1, 1), Point(-1, 1, 1)));
    //rear
    tris.push_back(Triangle(Point(1, 1, -1), Point(-1, -1, -1), Point(1, -1, -1)));
    tris.push_back(Triangle(Point(-1, -1, -1), Point(1, 1, -1), Point(-1, 1, -1)));
    /**/

}



int main ( int argc, char** argv )
{
    if ( SDL_Init( SDL_INIT_VIDEO ) < 0 ) {
        printf( "Unable to init SDL: %s\n", SDL_GetError() );
        return 1;
    }

    atexit(SDL_Quit);

    SDL_Surface* screen = SDL_SetVideoMode(800, 600, 16,
                                           SDL_OPENGL);
    if ( !screen ) {
        printf("Unable to set 640x480 video: %s\n", SDL_GetError());
        return 1;
    }

    resizeGL(800,600);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    createCube();

    bool done = false;
    while (!done) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_QUIT:
                    done = true;
                    break;
                case SDL_MOUSEMOTION:
                    rx = event.motion.x*360.0/800.0;
                    ry = event.motion.y*180.0/600.0;
                    break;
                case SDL_KEYDOWN:
                    if(event.key.keysym.sym == SDLK_w)
                        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    else if(event.key.keysym.sym == SDLK_f)
                        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    else if(event.key.keysym.sym == SDLK_s)
                        subdivide();
                    else if (event.key.keysym.sym == SDLK_ESCAPE)
                        done = true;
                    else if (event.key.keysym.sym == SDLK_n)
                        printf("number of nodes %d \n", (int)tris.size());
                    break;
            }
        }

        draw();
        SDL_GL_SwapBuffers();
    }
    writeSurface("cells.xyz", false);
    renderingScript("render_cells.py");
    return 0;
}
