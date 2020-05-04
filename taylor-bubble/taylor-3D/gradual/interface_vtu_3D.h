#include "geometry.h"
#include "fractions.h"

#if dimension == 1
coord mycs (Point point, scalar c) {
    coord n = {1.};
    return n;
}
#elif dimension == 2
# include "myc2d.h"
#else // dimension == 3
# include "myc.h"
#endif


void output_vtu (struct OutputFacets p)
{
    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
    scalar c = p.c;
    face vector s = p.s;
    if (!p.fp) p.fp = stdout;

    // print header text
    fputs ("<?xml version=\"1.0\"?>\n", p.fp);
    fputs ("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n", p.fp);
    fputs ("  <UnstructuredGrid>\n", p.fp);

    int nverts = 0;
    int nfacets = 0;

    foreach()
        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
            coord n;
            if (!s.x.i)
                // compute normal from volume fraction
                n = mycs (point, c);
            else {
                // compute normal from face fractions
                double nn = 0.;
                foreach_dimension() {
                    n.x = s.x[] - s.x[1];
                    nn += fabs(n.x);
                }
                assert (nn > 0.);
                foreach_dimension()
                    n.x /= nn;
            }
            double alpha = plane_alpha (c[], n);

            coord v[12];
            int m = facets (n, alpha, v, 1.);
            for (int i = 0; i < m; i++) {
                nverts ++;
            }
            if (m > 0) {
                nfacets ++;
            }
        }

        fprintf (p.fp, "    <Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n", nverts, nfacets);
        fputs ("      <Points>\n", p.fp);
        fputs ("        <DataArray type=\"Float32\" Name=\"vertices\" NumberOfComponents=\"3\" format=\"ascii\">\n", p.fp);

        int offsets[nfacets];

        int ifacet = 0;
        int offset = 0;

        foreach()
            if (c[] > 1e-6 && c[] < 1. - 1e-6) {
                coord n;
                if (!s.x.i)
                    // compute normal from volume fraction
                    n = mycs (point, c);
                else {
                    // compute normal from face fractions
                    double nn = 0.;
                    foreach_dimension() {
                        n.x = s.x[] - s.x[1];
                        nn += fabs(n.x);
                    }
                    assert (nn > 0.);
                    foreach_dimension()
                        n.x /= nn;
                }
                double alpha = plane_alpha (c[], n);

                coord v[12];
                int m = facets (n, alpha, v, 1.);
                for (int i = 0; i < m; i++) {
                    fprintf (p.fp, "%g %g %g ",
                             x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
                }
                if (m > 0) {
                    offset += m;
                    offsets[ifacet] = offset;
                    ifacet ++;
                }
            }


            fputs ("        </DataArray>\n", p.fp);
            fputs ("      </Points>\n", p.fp);
            fputs ("      <Cells>\n", p.fp);

            fputs ("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n", p.fp);

            // print vert numbers
            for (int ivert = 0; ivert < nverts; ivert++)
                fprintf (p.fp, "%i ", ivert);

            fputs ("        </DataArray>\n", p.fp);
            fputs ("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n", p.fp);

            // print offsets
            for (ifacet = 0; ifacet < nfacets; ifacet++)
                fprintf (p.fp, "%i ", offsets[ifacet]);

            fputs ("        </DataArray>\n", p.fp);
            fputs ("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", p.fp);

            // print cell type list
            for (ifacet = 0; ifacet < nfacets; ifacet++)
                fprintf (p.fp, "7 ");

            fputs ("        </DataArray>\n", p.fp);
            fputs ("      </Cells>\n", p.fp);
            fputs ("      <PointData>\n", p.fp);
            fputs ("      </PointData>\n", p.fp);
            fputs ("      <CellData>\n", p.fp);
            fputs ("      </CellData>\n", p.fp);
            fputs ("    </Piece>\n", p.fp);
            fputs ("  </UnstructuredGrid>\n", p.fp);
            fputs ("</VTKFile>\n", p.fp);

            fflush (p.fp);
            #if defined(_OPENMP)
            omp_set_num_threads(num_omp);
            #endif
}

void output_vtu_ascii_foreach (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  vertex scalar marker[];
  int no_points = 0, no_cells=0 ;
  foreach_vertex(){
    marker[] = _k;
    no_points += 1;
  }
  foreach(){
    no_cells += 1;
  }

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
    foreach(){
      fprintf (fp, "%g\n", val(s));
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
    foreach(){
#if dimension == 2
      fprintf (fp, "%g %g 0.\n", val(v.x), val(v.y));
#endif
#if dimension == 3
      fprintf (fp, "%g %g %g\n", val(v.x), val(v.y), val(v.z));
#endif
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
  foreach_vertex(){
#if dimension == 2
    fprintf (fp, "%g %g 0\n", x, y);
#endif
#if dimension == 3
    fprintf (fp, "%g %g %g\n", x, y, z);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  foreach(){
#if dimension == 2
    // Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g \n", marker[], marker[1,0], marker[1,1], marker[0,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0];
    int ape3 = marker[1,1];
    int ape4 = marker[0,1];
    fprintf (fp, "%u %u %u %u \n", ape1, ape2, ape3, ape4);
#endif
#if dimension == 3
    // Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g %g %g %g %g\n", marker[], marker[1,0,0], marker[1,1,0], marker[0,1,0],marker[0,0,1], marker[1,0,1], marker[1,1,1], marker[0,1,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0,0];
    int ape3 = marker[1,1,0];
    int ape4 = marker[0,1,0];
    int ape5 = marker[0,0,1];
    int ape6 = marker[1,0,1];
    int ape7 = marker[1,1,1];
    int ape8 = marker[0,1,1];
    fprintf (fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, ape6, ape7, ape8);

#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

  for (int i = 1; i < no_cells+1; i++){
#if dimension == 2
    fprintf (fp, "%d \n", i*4);
#endif
#if dimension == 3
    fprintf (fp, "%d \n", i*8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  foreach(){
#if dimension == 2
    fputs ("9 \n", fp);
#endif
#if dimension == 3
    fputs ("12 \n", fp);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

