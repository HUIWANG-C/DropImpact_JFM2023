/*In this script, we vary the Re and We, to observe its effect on droplet formation*/
/*Liquid viscosity and surface tension are varied. Air is the surrounding gas*/
/******************************Head Files**************************************************************************************/
#include "grid/octree.h"                          //3D
#include "navier-stokes/centered.h"               //N-S
#include "maxruntime.h"                           //restart  
//For large viscosity and density ratios, the harmonic mean for the viscosity tends to work better than the default arithmetic mean.
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)) 
#include "two-phase.h"                            //2-phase
#include "navier-stokes/conserving.h"             //momentum-conserving VOF advection of the velocity components  
#include "tension.h"                              //surface tension
#include "reduced.h"                              //gravity
#include "tracer.h"                               //passive tracer
#include "view.h"                                 //generate pics&movies online
#include "tag.h"                                  //tag droplets

#include "adapt_wavelet_limited.h"                //apply different max level on regions
#include "remove_droplet_sphere.h"                //remove droplet outside of the sphere

/******************************Macros******************************************************************************************/
//Geometry
#define DB          (4.0e-3)                      //m  drop diameter 4mm
#define y_D         (0.6*DB)                      //m  initial height of the drop 0.1d
#define x_D         (0.)                          //m  center                 
#define z_D         (0.)                          //m  center                                              
#define Drop(x,z,y) (sq((x-x_D)*4/4.3)+sq((z-z_D)*4/4.3)+sq((y-y_D)*4/3.8)) // an oblate shape   
#define DOMAIN      (DB*20.)                      //m  domain size L/DB=20=80mm

//Parameter space
#define Re          (30000.)                      // Reynolds number
#define We          (3000.)                       // Weber number
#define Impact_U    (7.2)                         // m/s     impact speed

//Liquid properties: density, viscosity, surface tension
#define RHO_SW      (1018.3)                      // kg/m^3  liquid density
#define RHO_A       (1.3)                         // kg/m^3  air density 20 degree
#define MU_SW       ((RHO_SW*DB*Impact_U)/Re)     // Pa*s    liquid viscosity
#define MU_A        (18.0e-6)                     // Pa*s    air viscosity
#define SIGMA_SW    ((RHO_SW*DB*Impact_U*Impact_U)/We)   // N/m     surface tension
#define Gravity     (9.8)                         // m/s^2   gravity  

//Running Time
#define MAXTIME     (260.01e-6)                      // s       stop simulation 0.260ms with extra layer
#define Dtdump      (10.0e-6)                     // s       saving dump 10us

//Mesh Refinement
int MINLEVEL = 6;
int MAXLEVEL = 13;

double fErr = 1e-4;              //Error tolerance of VOF
double uErr = 1e-1;              //Error tolerance of velocity

/******************************Boundary conditions*******************************************************************************/
//Top: outflow 
u.n[top]    = neumann (0.);
p[top]      = dirichlet (0.);
pf[top]     = dirichlet (0.);

/******************************Passive tracer*************************************************************************************/
//Passive Tracer for drop
scalar Tdrop[];                                
scalar * tracers = {Tdrop};  

/******************************Main Function*************************************************************************************/
int main(int argc, char * argv[])
{
  maxruntime (&argc, argv);
  if (argc > 1)
    MAXLEVEL = atoi (argv[1]);

  size(DOMAIN);
  origin(-DOMAIN/2., -DOMAIN/2., -DOMAIN/2.);
  init_grid(pow(2.0, MINLEVEL));

  rho1 = RHO_SW;                                     //density of seawater      
  rho2 = RHO_A;                                      //density of air
  mu1 = MU_SW;                                       //viscosicy of sea water
  mu2 = MU_A;                                        //viscosity of air
  f.sigma = SIGMA_SW;                                //surface tension of air-seawater
  G.y = -Gravity;
  
  //Poisson solver tolerance might be modified, to minimise mass conservation errors for very long simulations
  TOLERANCE = 1e-4 [*]; 
  
  fprintf(ferr,"\n\n# CASE SETTING SUMMARY:\n\n"); 
  fprintf(ferr,"# Numerical settings \n");
  fprintf(ferr,"# MAXLEVEL uErr Tolerance minDelta [um]\n");
  fprintf(ferr,"# %d, %g, %g, %g\n\n", MAXLEVEL, uErr, TOLERANCE, DOMAIN/(1 << MAXLEVEL)*1e6);

  fprintf(ferr,"# Physical properties values:\n");
  fprintf(ferr,"# U0 RHO_SW RHO_A MU_SW MU_A SIGMA_SW \n");
  fprintf(ferr,"# %g %g %g %g %g %g \n\n",
  Impact_U, rho1, rho2, mu1, mu2, f.sigma);

  fprintf(ferr,"# Problem parameters:\n");
  fprintf(ferr,"# Re We \n");
  fprintf(ferr,"# %g %g\n\n", Re, We);

  run();
}                   

/******************************Initial Conditions***********************************************************************************/
event init(t = 0)
{
  //restore the simulation from a previous “restart”
  if (!restore (file = "restart", list = (scalar *){all})){
    do{
        fraction (f, sq(DB/2)-Drop(x,z,y)); // The drop
	foreach()    
          u.y[] = (sq(DB/2)-Drop(x,z,y) > 0) ? -Impact_U : 0.0;
    }while (adapt_wavelet ({f, u}, (double[]){1e-6,1e-4,1e-4,1e-4}, MAXLEVEL, MINLEVEL).nf);  
    refine(sq(x)+sq(z)<sq(DB*1.5) && fabs(y) < (0.02*DB) && level < MAXLEVEL);  
      
    scalar m[];
    fraction (m, -y); // The pool
    foreach()
      f[] += m[];  
 
    //Passive tracers
    fraction (Tdrop, sq(DB/2)-Drop(x,z,y));    //Tdrop
    boundary(all);
  }
}

/******************************Adaptive Mesh*******************************************************************************************/
int refRegion(double x, double y, double z)
{
  int lev;  
  if ((sq(x)+sq(z)<sq(DB))&&(y>-0.025*DB)&&(y<0.2*DB))        //inner_cylinder_high level  
    lev = MAXLEVEL+1; 
  else if (sq(x)+sq(z)+sq(y) < sq(39.0e-3))      
    lev = MAXLEVEL;    
  else
    lev = MINLEVEL;

  return lev;
}

event adapt(i++)
{
  adapt_wavelet_limited ({f, u}, (double[]){fErr, uErr, uErr, uErr}, refRegion, MINLEVEL);
}

/******************************Remove top droplets*********************************************************************************************/
//Here we remove the water drop close to top boundary to avoid backflow and spurious current
event removedroplet(t = 0; t += Dtdump; t <= MAXTIME)
{
  remove_droplet_sphere(f, 39.0e-3, 1e-6);  
}

//////////////////////////////////////////////////////////////////////OUTPUT//////////////OUTPUT///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////OUTPUT//////////////OUTPUT///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/******************************Output pics*********************************************************************************************/
/*we output pics every 10us*/
event earlypics(t = 0; t += Dtdump; t <= MAXTIME)  
{ 
  char name[80]; 

  /*f*/
  {
    view (camera = "front", width = 1024, height = 1024, fov = 3); 
    clear();                                                        
    squares ("f", min = 0, max = 1, n = {0,0,1}, alpha = 10e-10 ); 
    sprintf (name, "3f_fov3_%.3f.png", t*1e3);    // the unit is ms
    save (name);
  } 

  /*front view*/
  {
    view (camera = "front", width = 1024, height = 1024, fov = 3);
    clear();
    draw_vof("f");
    sprintf (name, "3vof_fov3_%.3f.png", t*1e3);    // the unit is ms
    save (name);
  }

  /*level*/
  {
    scalar l[];
    foreach()
      l[] = level;
    view (camera = "front", width = 1024, height = 1024, fov = 3); 
    clear();                                                       
    squares ("l", min=MINLEVEL, max=MAXLEVEL+1, n = {0,0,1}, alpha = 10e-10 ); 
    sprintf (name, "3level_fov3_%.3f.png", t*1e3);    // the unit is ms
    save (name);                       
  }
}

/*we output pics every 10us*/
event pics(t = 0; t += Dtdump; t <= MAXTIME)  
{ 
  char name[80]; 

  /*front view*/
  {
    view (camera = "front", width = 1024, height = 1024, fov = 24);
    clear();
    draw_vof("f");
    sprintf (name, "1vof_fov24_%.3f.png", t*1e3);    // the unit is ms
    save (name);
  }  

  /*front view*/
  {
    view (camera = "front", width = 1024, height = 1024, fov = 12);
    clear();
    draw_vof("f");
    sprintf (name, "2vof_fov12_%.3f.png", t*1e3);    // the unit is ms
    save (name);
  } 

  /*level*/
  {
    scalar l[];
    foreach()
      l[] = level;
    view (camera = "front", width = 1024, height = 1024, fov = 12); 
    clear();                                                       
    squares ("l", min=MINLEVEL, max=MAXLEVEL+1, n = {0,0,1}, alpha = 10e-10 ); 
    sprintf (name, "2level_fov12_%.3f.png", t*1e3);    // the unit is ms
    save (name);                       
  }
}

/******************************Perfs*******************************************************************************************/
event perfs (i++) 
{
  static FILE * fp = fopen ("perfs", "w");
  if (i == 0)
    fprintf (fp,
	     "t dt mgp.i mgp.nrelax mgpf.i mgpf.nrelax mgu.i mgu.nrelax "
	     "grid->tn perf.t perf.speed npe\n");
  fprintf (fp, "%g %g %d %d %d %d %d %d %ld %g %g %d\n", 
	   t, dt, mgp.i, mgp.nrelax, mgpf.i, mgpf.nrelax, mgu.i, mgu.nrelax,
	   grid->tn, perf.t, perf.speed, npe());
  fflush (fp);
}

/******************************Output Snapshot************************************************************************************************/
/*we output dump file every 10us*/
event dumpfile(t = 0.; t += Dtdump; t <= MAXTIME) 
{
  char name[80];
  sprintf (name, "dump%.3f", t*1e3);   // the unit is ms
  dump (list = (scalar *){all}, name);
}

/******************************Energy dissipation************************************************************************************************/
int dissipation_rate (double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			     sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    rateWater += mu1/rho1*f[]*sqterm; //water
    rateAir   += mu2/rho2*(1. - f[])*sqterm; //air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}

event energy (i++) {
  static FILE * fpdissi = fopen("dissipation.dat", "w");
  double rates[2];
  dissipation_rate(rates);
  double dissWater_rate = rates[0];      //dissipation rate at this step
  double dissAir_rate   = rates[1];

  if (i == 0) {
    fprintf (fpdissi, "t dt dissAir dissWater\n");
  }

  fprintf (fpdissi, "%g %g %g %g\n", t, dt, dissAir_rate, dissWater_rate);
}

/******************************End Simulation*******************************************************************************************/
event end(t = MAXTIME)                   
{                                         
  printf("i=%d t=%g\n",i,t);
}
