/* file scaled_std_deb.c */
#include <R.h>
#include <math.h>

static double scaled_std_parms[12];
#define L_m   scaled_std_parms[0]
#define p_Am  scaled_std_parms[1]
#define v     scaled_std_parms[2]
#define k_J   scaled_std_parms[3]
#define kap   scaled_std_parms[4]
#define T_A   scaled_std_parms[5]
#define T_ref scaled_std_parms[6]
#define T_b   scaled_std_parms[7]
#define E_G   scaled_std_parms[8]
#define f     scaled_std_parms[9]
#define E_Hb  scaled_std_parms[10]
#define E_Hp  scaled_std_parms[11]
/* #define f_intercept parms[10] */


/* initializer  */
void init_scaled_std(void (* odeparms)(int *, double *))
{
  int N=12;
  odeparms(&N, scaled_std_parms);
}

/* Derivatives and 2 output variable */
void d_scaled_std_deb (int *neq, double *t, double *y, double *ydot,
                    double *yout, int *ip)
{
  if (ip[0] <1) error("nout should be at least 1");
  double E_m;
  double g;
#define e   y[0]
#define l   y[1]
#define uH   y[2]
#define uR y[3]
  
  E_m = p_Am/ v;
  g = E_G/ kap/ E_m;
  double k_M = v / (g * L_m);
  double l_T = 0;
  double k = k_J / k_M;
  double V_m =pow(v/(k_M*g), 3);
  
  double u_Hp = E_Hp/(g*E_m*V_m);
    
  static const double w_E = 23.9; //# molecular weight of reserve g mol^-1
  static const double d_v = 0.16; //# specific density of structure
  static const double mu_E = 550000; //# chemical potential of reserve J / mol 
  double w = p_Am * w_E / (v *d_v * mu_E); //#omega
  
  /* derivatives */
  ydot[0] = k_M * g * (f - y[0]) / y[0];
  ydot[1] = k_M / 3 * (y[0] - y[1] - l_T) / (1 + y[0] / g);
  ydot[2] = k_M * (y[2] < u_Hp) * ((1 - kap) * y[0] * pow(y[1], 2) * (g + y[1] + l_T) / (g + y[0]) - k * y[2]);
  ydot[3] = k_M * (y[2] > u_Hp) * ((1 - kap) * y[0] * pow(y[1], 2) * (g + y[1] + l_T) / (g + y[0]) - k * u_Hp);
  
  /* calculate derived quantities: */
  //wet weight
  yout[0] = pow(y[1]*L_m, 3) * (1 + f * w);
  
  
}