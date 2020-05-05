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
