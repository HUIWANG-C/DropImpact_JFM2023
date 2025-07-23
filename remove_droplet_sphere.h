//Useful functions that uses tag.h. to tag the droplets. 
//If the position of the droplet is outside of a sphere by given Radius, then it will be removed

#include "tag.h"

void remove_droplet_sphere(scalar c, double Radius, double Threshold)
{
    //tag all the droplets
    scalar m[];
    foreach()
      m[] = c[] > Threshold;
    int n = tag(m);

    double v[n];
    coord b[n];

    for (int j=0; j<n; j++)
      v[j] = b[j].x = b[j].y = b[j].z = 0.0;

    //compute volume and position of droplets
    foreach (serial)
      if (m[] > 0) {
        int j = m[] - 1;
        v[j] += dv()*c[];
        coord p = {x,y,z};
        foreach_dimension()
          b[j].x += dv()*c[]*p.x;
      }

  // Reduce for MPI
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  foreach()
    if (m[] > 0) {
      int j = m[] - 1;
      if(sq(b[j].x/v[j])+sq(b[j].y/v[j])+sq(b[j].z/v[j]) > sq(Radius))
	c[] = 0.;
    }
  boundary ({c});
}
