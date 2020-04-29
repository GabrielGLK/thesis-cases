/**
# A Scattering bview function

The formulation of the function is inspired by the functions under `draw.h` and assumes a `global int` by the name of `n_part` exists. It works in 2D and 3D.
 */
#if dimension == 3 
void glPointParameterfv(GLenum pname, const GLfloat * params);
#endif

struct _scatter {
  coord * loc;    // Coordinates
  float s, pc[3], coefs[3]; // point size, colour and distance attenuation coefs.
};

/**
   The function is inspired by the functions under `draw.h` and assumes a `global int` by the name of `n_part` exists.
 */

trace
void scatter (struct _scatter p){
  bview * view = draw();
#if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.); //from draw_lines()
#else // Dimension == 3
  if (!p.coefs[0]){ // A guess:
    p.coefs[0] = 0.01;
    p.coefs[1] = 0.2;
    p.coefs[2] = 0.5;
  }
  glPointParameterfv(GL_POINT_DISTANCE_ATTENUATION, p.coefs);
#endif
  if (p.pc)
    glColor3f(p.pc[0], p.pc[1], p.pc[2]);
  if (!p.s)
    p.s = 20;
  glPointSize(p.s);
  glBegin (GL_POINTS);
  for (int j = 0; j < n_part; j++){
#if dimension == 2    
    glvertex2d (view, p.loc[j].x, p.loc[j].y);
#else // dimension == 3
    glvertex3d (view, p.loc[j].x, p.loc[j].y, p.loc[j].z);
#endif
  }
  glEnd();
  view->ni++; 
}
/**
## Usage

* [Visualizing flow tracers](particles.h)

## Test

* [Distance attenuation in 3D](distanceat.c)

*/
