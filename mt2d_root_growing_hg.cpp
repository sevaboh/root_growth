// solver - Crank-Nicholson implicit
// no fractional derivative with respect to the time variable
// changing time steps
// root growing model (rp_form=3)
// specific storage term in left part for the cases with large saturated zone
// surface relief
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <cfloat>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>
#include <float.h>
#include <unistd.h>
#include <sys/time.h>
#include <ctime>
#include <omp.h>
#include "root_growing.h"
////////// auxiliary: Gamma function////////////////////////////

// Note that the functions Gamma and LogGamma are mutually dependent.
double Gamma
(
    double x    // We require x > 0
);

double LogGamma
(
    double x    // x must be positive
)
{
	if (x <= 0.0)
	{
		std::stringstream os;
        os << "Invalid input argument " << x <<  ". Argument must be positive.";
        throw std::invalid_argument( os.str() ); 
	}

    if (x < 12.0)
    {
        return log(fabs(Gamma(x)));
    }

	// Abramowitz and Stegun 6.1.41
    // Asymptotic series should be good to at least 11 or 12 figures
    // For error analysis, see Whittiker and Watson
    // A Course in Modern Analysis (1927), page 252

    static const double c[8] =
    {
		 1.0/12.0,
		-1.0/360.0,
		1.0/1260.0,
		-1.0/1680.0,
		1.0/1188.0,
		-691.0/360360.0,
		1.0/156.0,
		-3617.0/122400.0
    };
    double z = 1.0/(x*x);
    double sum = c[7];
    for (int i=6; i >= 0; i--)
    {
        sum *= z;
        sum += c[i];
    }
    double series = sum/x;

    static const double halfLogTwoPi = 0.91893853320467274178032973640562;
    double logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;    
	return logGamma;
}
double Gamma
(
    double x    // We require x > 0
)
{
	if (x <= 0.0)
	{
		std::stringstream os;
        os << "Invalid input argument " << x <<  ". Argument must be positive.";
        throw std::invalid_argument( os.str() ); 
	}

    // Split the function domain into three intervals:
    // (0, 0.001), [0.001, 12), and (12, infinity)

    ///////////////////////////////////////////////////////////////////////////
    // First interval: (0, 0.001)
	//
	// For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
	// So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
	// The relative error over this interval is less than 6e-7.

	const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

    if (x < 0.001)
        return 1.0/(x*(1.0 + gamma*x));

    ///////////////////////////////////////////////////////////////////////////
    // Second interval: [0.001, 12)
    
	if (x < 12.0)
    {
        // The algorithm directly approximates gamma over (1,2) and uses
        // reduction identities to reduce other arguments to this interval.
		
		double y = x;
        int n = 0;
        bool arg_was_less_than_one = (y < 1.0);

        // Add or subtract integers as necessary to bring y into (1,2)
        // Will correct for this below
        if (arg_was_less_than_one)
        {
            y += 1.0;
        }
        else
        {
            n = static_cast<int> (floor(y)) - 1;  // will use n later
            y -= n;
        }

        // numerator coefficients for approximation over the interval (1,2)
        static const double p[] =
        {
            -1.71618513886549492533811E+0,
             2.47656508055759199108314E+1,
            -3.79804256470945635097577E+2,
             6.29331155312818442661052E+2,
             8.66966202790413211295064E+2,
            -3.14512729688483675254357E+4,
            -3.61444134186911729807069E+4,
             6.64561438202405440627855E+4
        };

        // denominator coefficients for approximation over the interval (1,2)
        static const double q[] =
        {
            -3.08402300119738975254353E+1,
             3.15350626979604161529144E+2,
            -1.01515636749021914166146E+3,
            -3.10777167157231109440444E+3,
             2.25381184209801510330112E+4,
             4.75584627752788110767815E+3,
            -1.34659959864969306392456E+5,
            -1.15132259675553483497211E+5
        };

        double num = 0.0;
        double den = 1.0;
        int i;

        double z = y - 1;
        for (i = 0; i < 8; i++)
        {
            num = (num + p[i])*z;
            den = den*z + q[i];
        }
        double result = num/den + 1.0;

        // Apply correction if argument was not initially in (1,2)
        if (arg_was_less_than_one)
        {
            // Use identity gamma(z) = gamma(z+1)/z
            // The variable "result" now holds gamma of the original y + 1
            // Thus we use y-1 to get back the orginal y.
            result /= (y-1.0);
        }
        else
        {
            // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
            for (i = 0; i < n; i++)
                result *= y++;
        }

		return result;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Third interval: [12, infinity)

    if (x > 171.624)
    {
		// Correct answer too large to display. Force +infinity.
		double temp = DBL_MAX;
		return temp*2.0;
    }

    return exp(LogGamma(x));
}

//////////////////////////////////////
//////////////////////////////////////
//////////////////////////////////////
const double HtoPmult=9.8;

int N;  // number of nodes
int M; // x nodes
int flags=0; 	// bit 0 - no flow bottom condition is boundary_cond=0, otherwise - free flow condition
		// bit 1 - distinguish open and shaded soil
		// bit 3 - inf value particles
		// bit 5 - steady state for initial distribution
		// bit 7 - no js generation
		// bit 8 - always print inputs while PSO
		// bit 9 - full output not only in checking moments of time
		// bit 11 - van Genuchten-Muallem for water permeability coefficient instead of Averianov's model
double Hupper=-10;
int stat_maxsmooth=10000; // maximal number of smooth iteration for calculating initial distribution
double stat_eps=1e-5; // epsilon while smoothing initial conditions
double smooth_f=0.0; // 0 - average, !=0 - diffusion with smooth_f as diff.coef
double smooth_Hb=0; // H value on lower boundary for initial conditions
double LLx=4.0; // default width of the domain
int irr=0; // 1 - use given irrigation data, 0 - simulate irrigation
int fit=0; // 1 - do fitting, 2 - solve direct problem
double klin=0; // linear increase of filtration coefficients
double drip=0.2; // depth of drip irrigation pipeline
double *mask_2d=NULL; // relief map - depth from z=0 at which soil starts
int mask_2d_n;
// linear solver
int ls_max_iter=30;
double ls_eps=1e-12;
double ls_min_tau=0.01;
double ls_mult=1.25;
double ls_percent=0.66;
// plants and sprinklers parameters
int _nplants=0,_nsprinklers=0;
double *_plants_x=NULL;
double *_sprinklers_x=NULL;
double _rootR=0.5, _sprinklerR=0.5;
double max_rld=0.5; // maximal depth (width rootR linearly depends on depth)
double averianov_k=3.5; // base fixed exponent for hydraulic conductivity dependencies (if negative - absolute value is set to all layer of soil model)
// impulse irrigation
double imp_d0=0.15;
double imp_d1=0.3; // depth range
double imp_l0=0.2; // relative uptake level and imp_d0
double max_irr_time=86400;
double min_irr_time=0;
double max_tau_irr=25;
double max_tau=1800;
double zone_threshold=0.1; // threshold for moistened zone determination, kPa
// irrigation on ET
double irr_ET_sum=10; // ET sum in mm to apply irr
double irr_ET_flow=23.8; // in m^3/hour
// root formation
double grad_weight=0.001;
// specific storage
double _vgm_specific_storage=0.1;//0.027;
// H restriction
double glminH=-200.0,glmaxH=5.0;

char *root_systems_file=NULL;
int no_growth=0;
char *root_systems_out_file="root_systems.txt";

int fit_vgm=0; // 0 - only filtration k in one dimension, 1 - filtration k in 2D, 2 - plus theta_r, 3 - plus a and n
int fit_func=0; // 0 - sum of squares of differences, 1 - average relative error
int bounds_size=18;
double bounds[18][2]={{0.8,1.0}, //0 a - don't used
		     {0.8,1.0}, //1 b
		     {1.0,5.0}, //2 L
		     {0.0,1.0}, //3
		     {0.01,0.4}, //4 - a
		     {1.001,3.0}, //5 - n
		     {0.01,0.4}, //6
		     {1e-9,5e-6}, //7 k
		     {0.6,1.0}, //8 min wetness
		     {0.6,1.0}, //9 max wetness
		     {1,100}, //10 irr amount
		     {0.05,0.05}, //11 EVT coefs max, perc_k/irr_k max
		     {0.25,-0.25}, //12
		     {0.05,0}, //13 tikhonov regularization lambda
		     {0,2}, //14 ev_mu
		     {0,10}, //15 averianov_k
		     {0.01,100},//16 Hbottom
		     {0.01,0.01} //17 EVT coefs min, perc_k/irr_k min
		    }; 
#define idx(i,j) ((i)+(j)*(N+2))
/////////////////////////////////////
void grad_func_U(double x,double y,double &d0,double &d1,void *p);
////////////////////////////////////////////
// time-space-fractional solver 1d /////////////
// Da_tH=div(CvDb_zH))-S/C(H)
////////////////////////////////////////////
class H_solver {
public:
	double *b_U; // pressure head
	double *MM,*RP; // multiplication results and right part for implicit scheme

	double *rp_mult;
	// steps
	double tau,tau_m; // current time step length 
	double T; // current time
	double rlw_step; // step on which root layer wetness is computed
	double rlw_value; // value of last computer root layer wetness
	double I_step,I_value;
	double L,dL; // domain length and space variable step length
	// H equation and common
	double gamma,gamma_h; // gamma - space derivative power (vertical), gamma_h - horizontal
	double H0; // initial condition for H
	double Dg,Dg_h,g2b,g2b_h;
	double v; // root layer depth
	double *pr_Tr,*pr_Ev; // transpiration and evaporation rates
	double Hbottom;
	// VGM parameters
	double *vgm_ns,*vgm_s0s,*vgm_s1s,*vgm_as,*vgm_h0,*vgm_h1,*vgm_k,*vgm_kx,*vgm_power,*vgm_specific_storage; // VGM coefs for layered soil
	int vgm_nlayers;
	// EV per day
	double *EV_T;
	double **EV_F;
	double *EV_C;
	double *LAI;
	double ev_mu;
	int nEV_T;
	int nEV_F;
	int isEV;
	int implicit_mode=0;
	// variable root length 
	int nrl;
	double *rlT,*rlV;
	// fixed precipitation
	double *perc_T, *perc_A;
	double perc_k; // multiplier for precipitation (irrigation flow for subsurface DI)
	int nperc;
	// fixed irr
	double *irr_T, *irr_A;
	int nirr;
	double irr_k; // multiplier for irrigation (precipitation for subsurface DI)
	// irrigation scheduling 
	double min_wetness; // min wetness to apply irrigation
	double max_wetness; // wetness to be after irrigation
	double irr_volume; // volume of irrigation water
	double irr_time; // time range to apply irrigation water
	double min_irr_time; // minimal time of watering
	double irr_start_step;
	int irr_start_step0=0; // 1 if last irrigation started on first step
	double last_irr_step;
	double irr_sum_ET=0;
	double *sb_U; //save state before irrigation

	int boundary_cond; // bit 0 - bottom cond type (first or second order), bit 1 - upper cond type
	int rp_form,no_irr;
	double Lx,dLx; // length
	int nplants,nsprinklers;
	double *plants_x;
	double *sprinklers_x;
	double rootR, maxrootR,oldrootR,sprinklerR;
	double *sprinklers_area=NULL;
	double irr_am;
	double _averianov_k=averianov_k;

	double *pr_vgm_ns=NULL,*pr_vgm_s0s,*pr_vgm_s1s,*pr_vgm_as,*pr_vgm_k,*pr_vgm_kx,*pr_vgm_power,*pr_vgm_specific_storage; // VGM coefs per cell
	double *pr_dwdh=NULL;
	double *pr_w=NULL;
	double *pr_K=NULL,*pr_KX=NULL;
	double *pr_A1=NULL,*pr_A2=NULL,*pr_B1=NULL,*pr_B2=NULL,*pr_R=NULL;
	double *pr_D2x=NULL,*pr_D2y=NULL;
	double pr_I,pr_rlw,pr_fw,pr_twc;
	// water balance
	double balance;
	// surface
	double get_depth_at_i(int i)
	{
		if (mask_2d==NULL) return 0.0;
		double i2d=(i/(double)M)*(mask_2d_n-1);
		double ki=i2d-(int)i2d;
		int i1=1;
		if (i2d>=(mask_2d_n-1)) i1=0;
		return mask_2d[(int)i2d]+ ki*(mask_2d[(int)i2d+i1]-mask_2d[(int)i2d]);
	}
	int get_start_z_at_i(int i)
	{
		return get_depth_at_i(i)/dL;
	}
	// root growing model
	root **root_systems=NULL;
	double *root_system_map;
	// calc vgm model coefficient for each layer
	void vgm_calc_coefs()
	{
	    // only once
	    if (pr_vgm_ns==NULL)
	    {
#pragma omp single
	    {
		pr_vgm_ns=new double[N+2];
		pr_vgm_s0s=new double[N+2];
		pr_vgm_s1s=new double[N+2];
		pr_vgm_as=new double[N+2];
		pr_vgm_k=new double[N+2];
		pr_vgm_kx=new double[N+2];
		pr_vgm_power=new double[N+2];
		pr_vgm_specific_storage=new double[N+2];
		memset(pr_vgm_s0s,0,(N+2)*sizeof(double));
	    }
		// calculate
#pragma omp barrier
#pragma omp for
		for (int i=0;i<=N;i++)
		{
		double vgm_s0,vgm_s1,vgm_a,vgm_n,k,kx,avr_power,ss;
	        avr_power=fabs(_averianov_k);
		if (vgm_nlayers!=0)
		{
		    // i inside a layer
		    for (int j=0;j<vgm_nlayers;j++)
		    if (((i*dL)>=vgm_h0[j])&&((i*dL)<vgm_h1[j]))
		    {
			vgm_s0=vgm_s0s[j];
			vgm_s1=vgm_s1s[j];
			vgm_a=vgm_as[j];
			vgm_n=vgm_ns[j];
			k=vgm_k[j];
			kx=vgm_kx[j];
			ss=vgm_specific_storage[j];
			if (_averianov_k>0) avr_power=vgm_power[j];
			goto save;
		    }
		    double h1,h0,minh0=1e300,maxh1=0;
		    int i1,i0;
		    int mi1,mi0;
		    int first=1;
		    for (int j=0;j<vgm_nlayers;j++)
		    {
			if (vgm_h0[j]<minh0)
			{
			    minh0=vgm_h0[j];
			    mi0=j;
			}
			if (vgm_h1[j]>maxh1)
			{
			    maxh1=vgm_h1[j];
			    mi1=j;
			}
			for (int k=0;k<vgm_nlayers;k++)
			if (j!=k)
			if (((i*dL)>=vgm_h1[j])&&((i*dL)<vgm_h0[k]))
			{
			    if (first)
			    {
				h1=vgm_h1[j];
				i1=j;
				h0=vgm_h0[k];
				i0=k;
				first=0;
			    }
			    if (vgm_h1[j]>h1)
			    {
				h1=vgm_h1[j];
				i1=j;
			    }
			    if (vgm_h0[k]<h0)
			    {
				h0=vgm_h0[k];
				i0=k;
			    }
			}
		    }
		    // i between two layers
		    if (first==0)
		    {
			double _k=((i*dL)-h1)/(h0-h1);
			vgm_s0=_k*vgm_s0s[i0]+(1.0-_k)*vgm_s0s[i1];
			vgm_s1=_k*vgm_s1s[i0]+(1.0-_k)*vgm_s1s[i1];
			vgm_a=_k*vgm_as[i0]+(1.0-_k)*vgm_as[i1];
			vgm_n=_k*vgm_ns[i0]+(1.0-_k)*vgm_ns[i1];
			k=_k*vgm_k[i0]+(1.0-_k)*vgm_k[i1];
			kx=_k*vgm_kx[i0]+(1.0-_k)*vgm_kx[i1];
		    }
		    else
		    {
			// i below minimal
			if ((i*dL)<=minh0)
			{
			    vgm_s0=vgm_s0s[mi0];
			    vgm_s1=vgm_s1s[mi0];
			    vgm_a=vgm_as[mi0];
			    vgm_n=vgm_ns[mi0];
			    k=vgm_k[mi0];
			    kx=vgm_kx[mi0];
			    ss=vgm_specific_storage[mi0];
			    if (_averianov_k>0) avr_power=vgm_power[mi0];
			}
			else
			{
			    // i above maximal
			    if ((i*dL)>=maxh1)
			    {
				vgm_s0=vgm_s0s[mi1];
				vgm_s1=vgm_s1s[mi1];
				vgm_a=vgm_as[mi1];
				vgm_n=vgm_ns[mi1];
				k=vgm_k[mi1];
				kx=vgm_kx[mi1];
				ss=vgm_specific_storage[mi1];
				if (_averianov_k>0) avr_power=vgm_power[mi1];
			    }
			}
		    }
	    	}
		// save
save:
		pr_vgm_s0s[i]=vgm_s0;
		pr_vgm_s1s[i]=vgm_s1;
		pr_vgm_as[i]=vgm_a;
		pr_vgm_ns[i]=vgm_n;
		pr_vgm_k[i]=k;
		pr_vgm_kx[i]=kx;
		pr_vgm_specific_storage[i]=ss;
		if (_averianov_k>0) pr_vgm_power[i]=avr_power;
		}
	    }
	    // linear change of k
	    if (klin!=0.0)
#pragma omp for
	    for (int i=0;i<=N;i++)
	    {
		pr_vgm_k[i]+=klin*tau;
		pr_vgm_kx[i]+=klin*tau;
	    }
	}
	/// calc wetness on the base of van Genuchten model ////////
	double wetness(int i,double P)
	{
		P *= HtoPmult;
		if (P <= 0.0)
			return pr_vgm_s1s[i];
		return pr_vgm_s0s[i]+((pr_vgm_s1s[i] - pr_vgm_s0s[i]) / pow(1 + pow(pr_vgm_as[i]*P*10.19716, pr_vgm_ns[i]), (1 - 1 / pr_vgm_ns[i])));
	}
	void pr_wetness()
	{
	    // alloc
#pragma omp single
	    if (pr_w==NULL)
		pr_w=new double[(N+2)*(M+2)];
#pragma omp barrier
#pragma omp for
	    for (int x=0;x<=M;x++)
	    for (int i=0;i<=N;i++)
	    {
		double P;
		int idx1=idx(i, x);
		pr_w[idx1]=wetness(i,-b_U[idx1]);
	    }
	}
	// calc dw/dh 
	void pr_dw_dh()
	{
	    // alloc
#pragma omp single
	    if (pr_dwdh==NULL)
		pr_dwdh=new double[(N+2)*(M+2)];
#pragma omp barrier
#pragma omp for
	    for (int x=0;x<=M;x++)
	    for (int i=0;i<=N;i++)
	    {
		double P;
		double Ch = 0.0;
		int idx=0;
		P = -b_U[idx=idx(i, x)];
		// calculate
		if (P <= 0.0)
			Ch = 0.0;
		else
		{
			P*=10.19716;
			P*=HtoPmult;
			double aPn=pow(pr_vgm_as[i]*P, pr_vgm_ns[i]);
			Ch=-(1.0/(P/(HtoPmult*10.197196)))*(((1.0/pr_vgm_ns[i])-1)*pr_vgm_ns[i]*aPn*pow(aPn+1.0,(1.0/pr_vgm_ns[i])-2.0)*(pr_vgm_s1s[i]-pr_vgm_s0s[i]));
		}
		pr_dwdh[idx]=Ch;
	    }
	}
	// non-linear filtration coefficient 
	void pr_KK()
	{
	    // alloc
#pragma omp single
	    if (pr_K==NULL)
	    {
		pr_K=new double[(N+2)*(M+2)];
		pr_KX=new double[(N+2)*(M+2)];
	    }
#pragma omp barrier
#pragma omp for
	    for (int x=0;x<=M;x++)
	    for (int i=0;i<=N;i++)
	    {
		double Kr=1.0;
		// by Averianov
		if ((flags&2048)==0)
		        Kr = pow(((pr_w[idx(i,x)] - pr_vgm_s0s[i]) / (pr_vgm_s1s[i] - pr_vgm_s0s[i])), pr_vgm_power[i]);
		else // by VGM
			Kr = pow(((pr_w[idx(i,x)] - pr_vgm_s0s[i]) / (pr_vgm_s1s[i] - pr_vgm_s0s[i])), pr_vgm_power[i])
			    *pow(1.0 - pow(1.0 - pow(((pr_w[idx(i,x)] - pr_vgm_s0s[i]) / (pr_vgm_s1s[i] - pr_vgm_s0s[i])), 
				    (1.0 / (1.0 - 1.0 / pr_vgm_ns[i]))), (1.0 - 1.0 / pr_vgm_ns[i])), 2.0);
		pr_K[idx(i,x)]=Kr*pr_vgm_k[i];
		pr_KX[idx(i,x)]=Kr*pr_vgm_kx[i];
	    }
	}
	double avg_root_layer_wetness(int full=0)
	{
		double sum = 0.0;
		int n=0;
		if ((full==0)&&(T==rlw_step)) return rlw_value;
		for (int j = 1;j<M;j++)
			for (int i = 0;i<N;i++)
				if (full==0)
				{
				    if (rp_form!=3)
				    {
					if (i*dL <= v)
					    sum += 0.5*(pr_w[idx(i, j)]*rdf(i)+pr_w[idx(i+1, j)]*rdf(i+1))*rp_mult[j]*dL/v;
				    }
				    else
					if (root_system_map)
						sum += pr_w[idx(i, j)]*root_system_map[idx(i,j)]/nplants;
				}
				else
				{
					sum+=pr_w[idx(i,j)];
					n++;
				}
		if (n!=0) sum/=n;
		if (full==0)
		{
			rlw_step=T;
			rlw_value=sum;
		}
		return sum;
	}
	double total_water_content()
	{
		double sum = 0.0;
		for (int j = 1;j<M-1;j++)
			for (int i = 0;i<N ;i++)
				sum += 0.25*(pr_w[idx(i,j)]+pr_w[idx(i+1,j)]+pr_w[idx(i,j+1)]+pr_w[idx(i+1,j+1)])*dL*dLx;
		return sum;
	}
	// root density function in depth - int( rdf,i=0,v) = v
	double rdf(int si)
	{
		if (rp_form==1)
		// wheat
		//	ret=(2.21-3.72*(si*dL/v)+3.46*(si*dL*si*dL/(v*v))-1.87*(si*dL*si*dL*si*dL/(v*v*v)))*Tr;
		// beans
			return (1.44-0.14*(si*dL/v)-0.61*(si*dL*si*dL/(v*v))-0.69*(si*dL*si*dL*si*dL/(v*v*v)));
		if (rp_form==2) // for impulse irrigation - water uptake from a zone around pipeline
		{
			double a1=-(1-imp_l0)/((drip-imp_d0)*(drip-imp_d0));
			double b1=2.0*(1-imp_l0)/(drip-imp_d0);
			double c1=imp_l0;
			double a2=-1/((imp_d1-drip)*(imp_d1-drip));
			double c2=1.0;
			double imp_int=((1.0/3.0)*a1*((drip-imp_d0)*(drip-imp_d0)*(drip-imp_d0))+0.5*b1*(drip-imp_d0)*(drip-imp_d0)+
					c1*(drip-imp_d0))+((1.0/3.0)*a2*((imp_d1-drip)*(imp_d1-drip)*(imp_d1-drip))+c2*(imp_d1-drip));
			double ret=0;
			if (drip<=imp_d0)
			{
			    imp_int=(1.0/3.0)*a2*((imp_d1-drip)*(imp_d1-drip)*(imp_d1-drip)-(imp_d0-drip)*(imp_d0-drip)*(imp_d0-drip))+c2*(imp_d1-imp_d0);
			    if (((si*dL)>=imp_d0)&&((si*dL)<imp_d1))
				ret=v*(a2*(si*dL-drip)*(si*dL-drip)+c2)/imp_int;
			}
			else
			{
			    if ((si*dL>=drip)&&(si*dL<imp_d1))
				ret=v*(a2*(si*dL-drip)*(si*dL-drip)+c2)/imp_int;
			    if ((si*dL>=imp_d0)&&(si*dL<drip))
				ret=v*(a1*(si*dL-imp_d0)*(si*dL-imp_d0)+b1*(si*dL-imp_d0)+c1)/imp_int;
			}
			//printf("%g (%g %g %g) (%g %g) %g %g\n",si*dL,a1,b1,c1,a2,c2,imp_int,ret);
			return ret;
		}
		return 1.0;
	}
	// linear root density in width
	int Rp2d_mult_first=1;
	double shaded_L=0;
	double Rp2d_mult(int ii)
	{
		if (rootR!=oldrootR) Rp2d_mult_first=1;
		if (Rp2d_mult_first == 1)
#pragma omp critical
		{
			double cr = dLx;
			cr *= 0.5;
			for (int ix = 0;ix < M + 2;ix++)
				rp_mult[ix] =0;
			for (int i = 0;i < nplants;i++)
				{
					double sum=0.0;
					double x, mult;
					x = plants_x[i];
					for (int ix = 1;ix < M ;ix++)
					{
						double x0, r;
						x0 = ix*dLx;
						r = fabs(x - x0);
						if (r <= (rootR+cr))
						{
							if ((r - cr) < 0)
								mult = 1.0;
							else
								mult = (rootR-(r-cr)) / rootR;
							sum+=mult;
						}
					}
					for (int ix = 0;ix < M + 2;ix++)
					{
						double x0, r;
						x0 = ix*dLx;
						r = fabs(x - x0);
						if (r <= (rootR+cr))
						{
							if ((r - cr) < 0)
								mult = 1.0;
							else
								mult = (rootR-(r-cr)) / rootR;
							rp_mult[ix] +=(mult/sum);
						}
					}
				}
			if (nplants) // sum of rp_mult = 1
			for (int ix = 0;ix < M + 2;ix++)
				rp_mult[ix]/=nplants;
			// Ev+Tr in shaded zones, only EV in bare
			shaded_L=0;
			for (int ix = 0;ix < M + 2;ix++)
				if (rp_mult[ix]!=0.0)
					shaded_L+=dLx;
			Rp2d_mult_first = 0;
			oldrootR=rootR;
		}
		return rp_mult[ii];
	}
	// current irrigation flow (m^2/s)
	double _I()
	{
	    double I=0;
	    if (T==I_step) return I_value;
		if ((no_irr==0)||(no_irr==2))
		{
			if ((nirr != 0)&&(no_irr!=2)) // fixed irrigation
			{
				// calc linear combination
				int i;
				for (i = 0;i<nirr;i++)
					if (irr_T[i]>T)
						break;
				if ((i == 0) || (i == nirr))
				{
					if (i == nirr) i--;
					I = irr_A[i];
				}
				else
				{
					double k = (T - irr_T[i - 1]) / (irr_T[i] - irr_T[i - 1]);
					I = k*irr_A[i] + (1 - k)*irr_A[i - 1];
				}
			}
			else
			{
				double aw = pr_rlw;
				if (aw==0) return 0.0;
				if (no_irr==2) // irr per ET
				{
				    if (irr_start_step == -1)
				    {
					if (irr_sum_ET>=irr_ET_sum)
					{
						irr_start_step = T;
						if (T==tau)
							irr_start_step0=1;
						else
							irr_start_step0=0;
						// save state
						memcpy(sb_U,b_U,(N+2)*(M+2)*sizeof(double));
						irr_time=3600*irr_ET_sum*10.0/irr_ET_flow;
						irr_volume=(irr_ET_flow/3600.0)*(Lx/100.0)*0.01;
						irr_sum_ET=0;
					}
					else
					{
						double et=0;
						for (int i=0;i<M;i++)
							et+=pr_Tr[i]+pr_Ev[i];
						et/=M;
						irr_sum_ET+=irr_k*et*1000*tau;
					}
				    }
				    if (irr_start_step != -1) // irrigation is being applied
				    {
					I = irr_volume;
					if ((T - irr_start_step) > irr_time)
					if (irr_am!=0)
					{
						if (fit!=1)
						{
						    if (aw>=max_wetness)
							printf("irrigation amount - %g m^2 (%g-%g)\n", irr_am, irr_start_step, T);
						    else
							printf("irrigation amount - 1000 %g m^2 (%g-%g)\n", irr_am, irr_start_step, T);
						}
						irr_am = 0.0;
						I = 0.0;
						irr_start_step = -1; // stop irrigation	
						last_irr_step=T; // save last irrigation T
					}
				    }
				}
				else
				{
				    if (irr_start_step == -1) // irrigation is not being applied
					if (aw < min_wetness) // start irrigation
					{
						irr_start_step = T;
						if (T==tau)
							irr_start_step0=1;
						else
							irr_start_step0=0;
						// save state
						memcpy(sb_U,b_U,(N+2)*(M+2)*sizeof(double));
					}
				    if (irr_start_step != -1) // irrigation is being applied
				    {
					I = irr_volume;
					if ((((T - irr_start_step) > irr_time) ||
						(aw >= max_wetness)) && ((T-irr_start_step)>min_irr_time))
					if (irr_am!=0)
					{
						if (fit!=1)
						{
						    if (aw>=max_wetness)
							printf("irrigation amount - %g m^2 (%g-%g)\n", irr_am, irr_start_step, T);
						    else
							printf("irrigation amount - 1000 %g m^2 (%g-%g)\n", irr_am, irr_start_step, T);
						}
						irr_am = 0.0;
						I = 0.0;
						irr_start_step = -1; // stop irrigation	
						last_irr_step=T; // save last irrigation T
					}
				    }
				}
			}
		}
		I_step=T;
		I_value=I;
		return I;
	}
	// transpiration in m^2/s 
	double transp(int i,int x)
	{
		if (rp_form==3) // modelling of root system growth
		{
		    if (root_system_map)
			return (root_system_map[idx(i,x)]/nplants)*pr_Tr[x]*Lx;
		    return 0.0;
		}
		if ((i>=0)&&((i*dL)<=v))
			return rdf(i)*Rp2d_mult(x)*pr_Tr[x]*Lx*dL/v; 
		return 0.0;
	}
	// sink term (right part) - roots water uptake and drip irrigation
	double RpF(int i,int x)
	{
		// find nearest plant
		int pl=0;
		double minr=1e300;
		for (int ii = 0;ii < nplants;ii++)
		{
			double r=(dLx*x-plants_x[ii])*(dLx*x-plants_x[ii]);
			if (r<=minr)
			{
			    minr=r;
			    pl=1+ii;
			}
		}
		int upper_point=0;
		if (pl!=0)
			upper_point=get_start_z_at_i(plants_x[pl-1]/dLx);
		double ret=transp(i-upper_point,x);
		int si=i;
		// subsurface DI
		if (sprinklers_area==NULL)
		{
			sprinklers_area=new double[nsprinklers];
			memset(sprinklers_area,0,nsprinklers*sizeof(double));
			for (int ii = 0;ii < nsprinklers;ii++)
				for (int i=0;i<N;i++)
					for (int j=0;j<M;j++)
					{
						double r=(dLx*j-sprinklers_x[ii])*(dLx*j-sprinklers_x[ii])+(i*dL-drip)*(i*dL-drip);
						if (r<=sprinklerR*sprinklerR)
							sprinklers_area[ii]+=dLx*dL;
					}
		}
		int spr=0;
		if (fabs(i*dL-drip)<=sprinklerR)
		for (int ii = 0;ii < nsprinklers;ii++)
		{
			double r=(dLx*x-sprinklers_x[ii])*(dLx*x-sprinklers_x[ii])+(i*dL-drip)*(i*dL-drip);
			if (r<=sprinklerR*sprinklerR)
			{
				spr=1+ii;
				break;
			}
		}
		if (spr)
		{
			    double ii=pr_I;
			    if (ii!=0.0)
			        ret-=ii*irr_k*(dLx*dL/sprinklers_area[spr-1]); // irr amount in m^2/s
		}
		return ret/(dL*dLx);
	}
	// current given precipitation flow (m/s)
	double Perc()
	{
		if (nperc != 0) // percipitation
		{
			// calc linear conbination
			int i;
			double p = 0.0;
			for (i = 0;i<nperc;i++)
				if (perc_T[i]>T)
					break;
			if ((i == 0) || (i == nperc))
			{
				if (i == nperc) i--;
				p= perc_A[i];
			}
			else
			{
				double k = (T - perc_T[i - 1]) / (perc_T[i] - perc_T[i - 1]);
				p= k*perc_A[i] + (1 - k)*perc_A[i - 1];
			}
			return p;
		}
		return 0.0;
	}
	// upper boundary condition (DbU=Uc)
	double Uc(int x)
	{
		double E=pr_Ev[x]; // evaporation
		double I;
		int p0=get_start_z_at_i(x);
		double P=0.0; // precipitation
		if (E<0) E=0;
		I=pr_I;
		P=Perc()*perc_k;
		if (flags&2)
		{
		    if (rp_form!=3)
		    {
			 if (rp_mult[x]!=0) P=0; // no rain under plants
			 else P*=Lx/(Lx-shaded_L); // more rain between plants
		    }
		    else
		    {
			double m=0;
			if (root_system_map)
			{
			    for (int ii=0;ii<N;ii++)
			    {
				m=root_system_map[idx(ii,x)];
				if (m!=0)
				    break;
			    }
			    if (m!=0) P=0; // no rain under plants
			    else P*=Lx/(Lx-shaded_L); // more rain between plants
			}
		    }
		}
		if (b_U[idx(p0,x)]>=glmaxH) P=0;
		double k0=fabs(b_U[idx(p0, x)])/fabs(b_U[idx(p0, x)]+b_U[idx(p0+1, x)]);
		double k1=1.0-k0;
		double kk=(k1*pr_K[idx(p0+1,x)]+k0*pr_K[idx(p0,x)]);
		if (fabs(pr_K[idx(p0+1,x)]/pr_K[idx(p0,x)])>1e10) kk=(pr_K[idx(p0+1,x)]+pr_K[idx(p0,x)])/2;
		double max_flux=-b_U[idx(p0+1,x)]*kk; // restriction on maximal flux
		if ((max_flux>0)&&(E>max_flux)) {E=max_flux;pr_Ev[x]=E;}
		double v=dL*(((E-P)/kk)+1);
		return v;
	}
	// sets current Ev, Tr, root system length(depth) and width
	void set_EV_TR(int xx)
	{
		double v=0.0,k;
		double l=0.0;
		int i=0;
		// ev_t have to be sorted by T
		// find t0<=t<=t1
		double t0=EV_T[0],t1,nx,nr2;
		int it0=0,it1;
		int if0,if1;
		double x=xx*dLx;
		for (i=0;i<nEV_T;i++)
		{
			if (EV_T[3*i+0]>T)
				break;
			if (EV_T[3*i+0]>t0)
			{
			    t0=EV_T[3*i+0];
			    it0=i;
			}
		}
		if (i!=(nEV_T-1))
		{
			t1=EV_T[3*i+0];
			it1=i;
		}
		else
		{
			t1=EV_T[3*(i-1)+0];
			it1=i-1;
		}
		// find nearest EV(x,y) for t0 and t1
		if0=i=it0+1;
		nx=EV_T[3*i+1];
		nr2=(x-nx)*(x-nx);
		while (EV_T[3*i+0]==t0)
		{
			double r2=(x-EV_T[3*i+1])*(x-EV_T[3*i+1]);
			if (r2<nr2)
			{
			    nx=EV_T[3*i+1];
			    if0=i;
			    nr2=r2;
			}
			i++;
		}
		if (t0==t1)
			if1=if0;
		else
		{
			if1=i=it1+1;
			nx=EV_T[3*i+1];
			nr2=(x-nx)*(x-nx);
			while (EV_T[3*i+0]==t0)
			{
			    double r2=(x-EV_T[3*i+1])*(x-EV_T[3*i+1]);
			    if (r2<nr2)
			    {
				nx=EV_T[3*i+1];
				if1=i;
				nr2=r2;
			    }
			    i++;
			}
		}
		// linear combination
		if (t0!=t1)
			k=(T-t0)/(t1-t0);
		else
			k=1.0;
		for (int j=0;j<nEV_F;j++)
			v+=k*EV_F[j][if1]*EV_C[j]+(1-k)*EV_F[j][if0]*EV_C[j];
		l=k*LAI[if1]+(1-k)*LAI[if0];
		// get root length
		if (nrl)
		{
		    for (i=0;i<nrl;i++)
			if (rlT[i]>T)
				break;
		    if ((i==0)||(i==nrl))
		    {
			if (i == nrl) i--;
			this->v=rlV[i];
		    }
		    else
		    {
			k=(T-rlT[i-1])/(rlT[i]-rlT[i-1]);
			this->v=k*rlV[i]+(1-k)*rlV[i-1];
		    }
		    rootR=maxrootR*this->v/max_rld;
		}
		// divide EV into Tr and E
		double M=1.0-exp(-ev_mu*l);
		pr_Tr[xx]=M*v;
		pr_Ev[xx]=(1-M)*v;
	}
	void precalc_Ev_Tr()
	{
		for (int j=0;j<M+1;j++)
		    set_EV_TR(j);
	}	
	// coefficients of five-diagonal linear equations system
	// xy=0 - for (i-1,x), xy=1 - for (i,x-1)
	double A_(int i,int x,int xy)
	{
		if (xy==0) return 0.5*Dg*0.5*(pr_K[idx(i-1,x)]+pr_K[idx(i,x)]);
		return 0.5*Dg_h*0.5*(pr_KX[idx(i,x-1)]+pr_KX[idx(i,x)]);
	}
	// xy=0 - for (i+1,x), xy=1 - for (i,x+1)
	double B(int i,int x,int xy)
	{
		if (xy==0) return 0.5*Dg*0.5*(pr_K[idx(i+1,x)]+pr_K[idx(i,x)]);
		return 0.5*Dg_h*0.5*(pr_KX[idx(i,x+1)]+pr_KX[idx(i,x)]);
	}
	// for (i,x)
	double R(int i,int x)
	{
		return 0.5*0.5*((pr_K[idx(i+1,x)]+pr_K[idx(i-1,x)]+2.0*pr_K[idx(i,x)])*Dg+(pr_KX[idx(i,x+1)]+pr_KX[idx(i,x-1)]+2.0*pr_KX[idx(i,x)])*Dg_h)+
		((pr_dwdh[idx(i,x)]+(pr_w[idx(i,x)]/pr_vgm_s1s[i])*pr_vgm_specific_storage[i])/tau);
	}
	// right part for (ii,x)
	double Rp(int ii,int x,double *b_Uold)
	{
		double ret=0.0;
		int idx1=idx(ii,x);
		int idxm1x=idx(ii-1,x);
		int idxm1y=idx(ii,x-1);
		int idxp1x=idx(ii+1,x);
		int idxp1y=idx(ii,x+1);
		ret = - ((pr_dwdh[idx1]+(pr_w[idx1]/pr_vgm_s1s[ii])*pr_vgm_specific_storage[ii])/tau)*b_Uold[idx1] + RpF(ii,x);
		ret-=pr_A1[idx1]*b_Uold[idxm1x]+pr_A2[idx1]*b_Uold[idxm1y]
		     -(pr_A1[idx1]+pr_A2[idx1]+pr_B1[idx1]+pr_B2[idx1])*b_Uold[idx1]
		     +pr_B1[idx1]*b_Uold[idxp1x]+pr_B2[idx1]*b_Uold[idxp1y];
		ret+=Dg*dL*(pr_K[idxp1x]-pr_K[idx1]);
		return ret;
	}
	// space-fractional derivative non-local part for (ii,x)
	double *sp_precalc=NULL;
	double *sp_precalc_2=NULL;
	double SFD(int ii,int x)
	{
		double nlp=0.0;
		// precalc coefficients
		if (gamma!=1)
		if (sp_precalc==NULL)
#pragma omp critical
		{
			if (sp_precalc) delete [] sp_precalc;
			sp_precalc = new double[N + 2];
			for (int i = 0;i < N + 2;i++)
				sp_precalc[i] = pow((double)i + 1, 1.0-gamma) - pow((double)i , 1-gamma);
		}
		if (gamma_h!=1)
		if (sp_precalc_2==NULL)
#pragma omp critical
		{
			if (sp_precalc_2) delete [] sp_precalc_2;
			sp_precalc_2 = new double[M + 2];
			for (int i = 0;i < M + 2;i++)
				sp_precalc_2[i] = pow((double)i + 1, 1.0-gamma_h) - pow((double)i , 1-gamma_h);
		}
		if (gamma!=1)
		{
		    double m1=g2b/(pow(ii*dL,1.0-gamma)*pow(2.0,gamma));
		    double m2=g2b/(pow(L-ii*dL,1.0-gamma)*pow(2.0,gamma));
		    for (int i = 1;i <N;i++)
			if (i!=ii)
			    nlp+=sp_precalc[abs(ii-i)]*pr_D2x[idx(i,x)]*((i<ii)?m1:m2);
		}
		if (gamma_h!=1)
		{
		    double m1=g2b_h/(pow(x*dLx,1.0-gamma_h)*pow(2.0,gamma_h));
		    double m2=g2b_h/(pow(Lx-x*dLx,1.0-gamma_h)*pow(2.0,gamma_h));
		    for (int i = 1;i <M;i++)
			if (i!=x)
			    nlp+=sp_precalc_2[abs(x-i)]*pr_D2y[idx(ii,i)]*((i<x)?m1:m2);
		}
		return nlp;
	}
	// implicit scheme - matrix multiplication of vector UU
	void mmult(double *UU)
	{
#pragma omp for nowait collapse(2)
		for (int j = 1;j < M ;j++)
			for (int i = 1;i < N ;i++)
			{
				double rr=pr_R[idx(i,j)];
				MM[idx(i,j)]=UU[idx(i-1,j)]*pr_A1[idx(i,j)]+UU[idx(i,j-1)]*pr_A2[idx(i,j)]+
					     UU[idx(i+1,j)]*pr_B1[idx(i,j)]+UU[idx(i,j+1)]*pr_B2[idx(i,j)]-
					     UU[idx(i,j)]*rr;
				// normalization
				if (rr!=0.0)
				    MM[idx(i,j)]/=rr;
				// internal bc on x
				if (get_start_z_at_i(j)<i)
				{
					if (get_start_z_at_i(j-1)>i)
						MM[idx(i,j)]=UU[idx(i,j+1)]-UU[idx(i,j)];
					if (i<get_start_z_at_i(j+1))
						MM[idx(i,j)]=UU[idx(i,j)]-UU[idx(i,j-1)];
				}
			}
#pragma omp for nowait
		for (int j = 1;j < M ;j++)
		{
			// boundary conditions
			if ((boundary_cond&1) == 1)
				MM[idx(N,j)]=UU[idx(N,j)];
			else
				MM[idx(N,j)]=UU[idx(N,j)]-UU[idx(N-1,j)];
			// upper boundary condition
			if (((boundary_cond>>1)&1) == 1)
				MM[idx(get_start_z_at_i(j),j)]=UU[idx(get_start_z_at_i(j),j)];
			else
				MM[idx(get_start_z_at_i(j),j)]=UU[idx(get_start_z_at_i(j),j)]-UU[idx(get_start_z_at_i(j)+1,j)];
			for (int k=0;k<get_start_z_at_i(j);k++)
				MM[idx(k,j)]=0.0;
		}
		// boundary conditions
#pragma omp for
		for (int i = 1;i < N ;i++)
		{
			MM[idx(i,M)]=UU[idx(i,M)]-UU[idx(i,M-1)];
			MM[idx(i,0)]=-UU[idx(i,1)]+UU[idx(i,0)];
		}
	}
	// precalculations for current time step
	void __precalc()
	{
		// H restriction
#pragma omp single
		{
		    for (int j = 0;j <= M ;j++)
			for (int i = 0;i <= N ;i++)
			{
			    if (b_U[idx(i,j)]>=glmaxH) b_U[idx(i,j)]=glmaxH;
			    if (b_U[idx(i,j)]<=glminH) b_U[idx(i,j)]=glminH;
			}
		}
		vgm_calc_coefs();
		pr_wetness();
		pr_dw_dh();
		pr_KK();
		// precalc matrix coefficients
#pragma omp single
	    {
		pr_rlw=avg_root_layer_wetness();
		pr_fw=avg_root_layer_wetness(1);
		pr_twc=total_water_content();
		pr_I=_I();
		// set tau to initial when irrigation starts
		if (irr_start_step==T)
			tau=tau_m;
		if (pr_A1==NULL) pr_A1=new double[(N+2)*(M+2)];
		if (pr_A2==NULL) pr_A2=new double[(N+2)*(M+2)];
		if (pr_B1==NULL) pr_B1=new double[(N+2)*(M+2)];
		if (pr_B2==NULL) pr_B2=new double[(N+2)*(M+2)];
		if (pr_R==NULL) pr_R=new double[(N+2)*(M+2)];
	    }
#pragma omp barrier
		// for space-fractional derivative
		if (gamma!=1)
#pragma omp for nowait collapse(2)
		for (int x = 1;x < M ;x++)
			for (int i = 1;i < N ;i++)
				pr_D2x[idx(i,x)]= 0.5*Dg*(((pr_K[idx(i-1,x)]+pr_K[idx(i,x)])*b_U[idx(i-1,x)] - (pr_K[idx(i-1,x)]+pr_K[idx(i+1,x)]+2*pr_K[idx(i,x)])*b_U[idx(i,x)]+(pr_K[idx(i+1,x)]+pr_K[idx(i,x)])*b_U[idx(i+1,x)])+2.0*dL*(pr_K[idx(i+1,x)]-pr_K[idx(i,x)]));
		if (gamma_h!=1)
#pragma omp for nowait collapse(2)
		for (int i = 1;i < M ;i++)
			for (int ii = 1;ii < N ;ii++)
				pr_D2y[idx(ii,i)] = 0.5*Dg_h*(((pr_KX[idx(ii,i-1)]+pr_KX[idx(ii,i)])*b_U[idx(ii,i - 1)] - (pr_KX[idx(ii,i-1)]+pr_KX[idx(ii,i+1)]+2*pr_KX[idx(ii,i)])*b_U[idx(ii,i)]+(pr_KX[idx(ii,i+1)]+pr_KX[idx(ii,i)])*b_U[idx(ii,i+1)]));
		if (isEV)
			precalc_Ev_Tr();
#pragma omp for nowait collapse(2)
		for (int j = 1;j < M ;j++)
			for (int i = 1;i < N ;i++)
			{
			    pr_A1[idx(i,j)]=A_(i,j,0);
			    pr_A2[idx(i,j)]=A_(i,j,1);
			    pr_B1[idx(i,j)]=B(i,j,0);
			    pr_B2[idx(i,j)]=B(i,j,1);
			    pr_R[idx(i,j)]=R(i,j);
			    // right part
			    RP[idx(i,j)]=Rp(i,j,b_U);
			    if ((gamma!=1)||(gamma_h!=1)) RP[idx(i,j)]-=SFD(i,j);
			    // normalization
			    if (pr_R[idx(i,j)]!=0.0)
				RP[idx(i,j)]/=pr_R[idx(i,j)];
				// internal bc on x
				if (get_start_z_at_i(j)<i)
				{
					if (get_start_z_at_i(j-1)>i)
						RP[idx(i,j)]=0;
					if (i<get_start_z_at_i(j+1))
						RP[idx(i,j)]=0;
				}
			}
#pragma omp for
		for (int j = 1;j < M ;j++)
		{
			// bottom boundary condition
			if ((boundary_cond&1) == 1)
			{
				if (Hbottom!=1e300)
					RP[idx(N,j)]=Hbottom;
				else
					RP[idx(N,j)] = H0 - N*dL;
			}
			else
				RP[idx(N,j)]=(((flags&1)==1)?dL:0.0);
			// upper boundary condition
			if (((boundary_cond>>1)&1) == 1)
				RP[idx(get_start_z_at_i(j),j)]=Hupper;
			else
				RP[idx(get_start_z_at_i(j),j)]=-Uc(j);
			for (int k=0;k<get_start_z_at_i(j);k++)
				RP[idx(k,j)]=0.0;
		}
	}
	void calc_root_systems_map()
	{
		memset(root_system_map,0,(N+2)*(M+2)*sizeof(double));
		for (int i=0;i<nplants;i++)
		{
		    double b0=root_systems[i]->bounds(0);
		    double b1=root_systems[i]->bounds(1);
		    double *d=root_systems[i]->density(M,N);
		    //printf("%g %d %g %g %d\n",T,i,b0,b1,root_systems[i]->size());fflush(stdout);
		    for (int j=0;j<N;j++)
			for (int k=0;k<M;k++)
			{
			    double x=plants_x[i]-b0+k*(2*b0)/M;
			    double y=j*b1/N;
			    int xi=x/dLx, yi=y/dL;
			    if ((xi>=0)&&(xi<=M)&&(yi>=0)&&(yi<=N))
			    {
				double dx=(x/dLx)-xi, dy=(y/dL)-yi;
				root_system_map[idx(yi,xi)]+=(1-dx)*(1-dy)*d[k*N+j];
				root_system_map[idx(yi,xi+1)]+=dx*(1-dy)*d[k*N+j];
				root_system_map[idx(yi+1,xi)]+=(1-dx)*dy*d[k*N+j];
				root_system_map[idx(yi+1,xi+1)]+=dx*dy*d[k*N+j];
				//if (d[k*N+j]!=0.0) printf("%g - %d %d - %g %g - %g - %g %g %g %g\n",T,xi,yi,dx,dy,d[k*N+j],root_system_map[idx(yi,xi)],root_system_map[idx(yi,xi+1)],root_system_map[idx(yi+1,xi)],root_system_map[idx(yi+1,xi+1)]);
			    }
			}
		    delete [] d;
		}
		// normalize 
		double sum=0;
		for (int j=0;j<=M;j++)
		    for (int i=0;i<=N;i++)
			sum+=root_system_map[idx(i,j)];
		if (sum!=0.0)
		for (int j=0;j<=M;j++)
		    for (int i=0;i<=N;i++)
			root_system_map[idx(i,j)]/=(sum/nplants);
		// calc shaded_L
		shaded_L=0;
		for (int j=0;j<M;j++)
		{
		    double m=0;
		    for (int i=0;i<N;i++)
		    {
			m=root_system_map[idx(i,j)];
			if (m!=0)
			    break;
		    }
		    if (m!=0)
			shaded_L+=dLx;
		}
	}
	// tfqmr
	void calc_step(FILE *fi=NULL)
	{
		if (T==tau)
		{
		    balance=total_water_content();
		    for (int i =0;i <= N+1 ;i++)
		    for (int j = 0;j <= M+1 ;j++)
		    {
			MM[idx(i,j)]=0.0;
			RP[idx(i,j)]=0.0;
		    }
		    Rp2d_mult_first = 1;
		    irr_am = 0.0;
		    if (rp_form==3) // initialize root system
		    {
			printf("root generation params: ");
			param_sets[0].write(stdout);
			root_systems=new root*[nplants];
			root_system_map=new double[(N+2)*(M+2)];
			if (root_systems_file==NULL)
			    for (int i=0;i<nplants;i++)
				root_systems[i]=new root(param_sets[0],segment(plants_x[i],0,plants_x[i],tau*param_sets[0].r),grad_func_U,0,this);
			else
			{
			    FILE *fi=fopen(root_systems_file,"rt");
			    for (int i=0;i<nplants;i++)
			    {
				root_systems[i]=new root();
				root_systems[i]->read(fi,grad_func_U,this);
			    }
			    fclose(fi);
			}
			calc_root_systems_map();
		    }
		    if (nirr!=0)
			memcpy(sb_U,b_U,(N+2)*(M+2)*sizeof(double));

		}
		//////////
		double *w=new double[(N+2)*(M+2)];
		double *y[2];
		y[0]=new double[(N+2)*(M+2)];
		y[1]=new double[(N+2)*(M+2)];
		double *rr=new double[(N+2)*(M+2)];
		double *v=new double[(N+2)*(M+2)];
		double *d=new double[(N+2)*(M+2)];
		double *x=new double[(N+2)*(M+2)];
		double theta=0,nu=0,tt=0,ro=0,ro1=0,c=0,b=0;
		double sgm,aa,rv;
		int n,m;
start:
		theta=nu=tt=ro=ro1=c=b=0.0;
		memcpy(x,b_U,(N+2)*(M+2)*sizeof(double));
#pragma omp parallel
	    {
		////////// precalculations
		__precalc();
		mmult(x);
#pragma omp for reduction(+:tt)
		for (int i=0;i<(N+2)*(M+2);i++)
		{
			// w_1=y_1=r_0=b-Ax_0
			// rr0: ro=rr_0*r_0!=0 - rr0=r0
			w[i]=y[0][i]=RP[i]-MM[i];
			// d=0
			d[i]=0;
			// tau=||r0||
			tt+=w[i]*w[i];
		}
#pragma omp single
		tt=sqrt(tt);
#pragma omp barrier
		// random rr0, ro0=rr_0*r_0
#pragma omp for reduction(+:ro)
		for (int i=0;i<(N+2)*(M+2);i++)
		{
			rr[i]=tt*((rand() % 10000) / (10000.0 - 1.0));
			ro+=rr[i]*w[i];
		}
		// v=Ay_1
		mmult(y[0]);
#pragma omp for
		for (int i=0;i<(N+2)*(M+2);i++)
			v[i]=MM[i];
#pragma omp single
		n=1;
#pragma omp barrier
loop1:
		{
			int br=0;
			// sigma_n-1 - rr_0*v_n-1
#pragma omp single
			sgm=0;
#pragma omp barrier
#pragma omp for reduction(+:sgm)
			for (int i=0;i<(N+2)*(M+2);i++)
				sgm+=rr[i]*v[i];
			// a_n-1=po_n-1/sigma_n-1
#pragma omp single
			aa=ro/sgm;
#pragma omp barrier
			// y_2n=y_2n-1 - a_n-1 * v_n-1
#pragma omp for
			for (int i=0;i<(N+2)*(M+2);i++)
				y[1][i]=y[0][i]-aa*v[i];
#pragma omp single
			m=2*n-1;
#pragma omp barrier
loop2:
			{
				double ot=theta,onu=nu;
				// w_m+1=w_m-a_n-1Ay_m
				mmult(y[m-(2*n-1)]);
#pragma omp for
				for (int i=0;i<(N+2)*(M+2);i++)
					w[i]=w[i]-aa*MM[i];
				// theta_m=||w_m+1||/tau_m-1; c_m=1/sqrt(1+theta_m^2)
#pragma omp single
				theta=0;
#pragma omp barrier
#pragma omp for reduction(+:theta)
				for (int i=0;i<(N+2)*(M+2);i++)
					theta+=w[i]*w[i];
#pragma omp single
			    {
				theta=sqrt(theta)/tt;
				c=1.0/sqrt(1.0+theta*theta);
				// tau_m=tau_m-1 * theta_m *c_m; nu_m=c_m^2 *a_n-1
				tt=tt*theta*c;
				nu=c*c*aa;
				rv=0.0;
			    }
#pragma omp barrier
				// d_m = y_m+(theta_m-1^2 nu_m-1 / a_n-1)*d_m-1
				// x_m=x_m-1 + nu_m *d_m
#pragma omp for
				for (int i=0;i<(N+2)*(M+2);i++)
				{
					d[i]=y[m-(2*n-1)][i]+d[i]*(ot*ot*onu/aa);
					x[i]=x[i]+nu*d[i];
				}
				mmult(x);
#pragma omp for reduction(+:rv)
				for (int i=0;i<(N+2)*(M+2);i++)
					rv+=(RP[i]-MM[i])*(RP[i]-MM[i]);
#pragma omp single
				rv=sqrt(rv)/((N+2)*(M+2));
#pragma omp barrier
				if (rv<ls_eps)
				{
				    br=1;
				    goto eloop2;
				}
			}
#pragma omp single
			m++;
#pragma omp barrier
			if (m<=2*n)
			    goto loop2;
eloop2:
			if (br==1)
				goto eloop1;
			// ro_n=rr0*w_2n+1, beta_n=ro_n/ro_n-1
#pragma omp single
			ro1=0;
#pragma omp barrier
#pragma omp for reduction(+:ro1)
			for (int i=0;i<(N+2)*(M+2);i++)
				ro1+=rr[i]*w[i];
#pragma omp single
		    {
			b=ro1/ro;
			ro=ro1;
		    }
#pragma omp barrier
			// y_2n+1 = w_2n+1+beta_n*y_2n
#pragma omp for
			for (int i=0;i<(N+2)*(M+2);i++)
				y[0][i]=w[i]+b*y[1][i];
			// v_n=Ay_2n+1+b*(Ay_2n+b*v_n-1)
			mmult(y[1]);
#pragma omp for
			for (int i=0;i<(N+2)*(M+2);i++)
				v[i]=b*(MM[i]+b*v[i]);
			mmult(y[0]);
#pragma omp for
			for (int i=0;i<(N+2)*(M+2);i++)
				v[i]=MM[i]+v[i];
		}
#pragma omp single
		n++;
#pragma omp barrier
		if (n<ls_max_iter)
		    goto loop1;
eloop1:;
	    }
//		printf("t %g N %d M %d n %d r %g tau %g irr %g tw %g balance %g rlw %g Ev %g Tr %g Perc %g Irrrate %g irram %g outflow %g\n",T,N,M,n,rv,tau,irr_start_step,pr_twc,balance,pr_rlw,pr_Ev[M/2],pr_Tr[M/2],Perc(),pr_I*nsprinklers,irr_am,0);fflush(stdout);
		// change tau and recalc if needed
		if (n==ls_max_iter)
		{
		     if (tau<ls_min_tau)
		     {
			 x[0]=NAN;
			 printf("minimal tau value reached\n");
		     }
		     else
		     {
		        T-=tau;
		        tau/=ls_mult;
		        T+=tau;
		        goto start;
		     }
		}
		// irrigation ammount
		irr_am += pr_I*tau*nsprinklers;
		// debug output
		double outflow=0;
#pragma omp parallel for reduction(+:outflow)
		for (int j=1;j<M;j++)
			outflow+=-pr_K[idx(N,j)]*(((b_U[idx(N,j)]-b_U[idx(N-1,j)])/dL)-1)*dLx;
		double mapsum=1;
		int rss=0;
		if (rp_form==3)
		{
		    rss=0;
		    for (int i=0;i<nplants;i++)
			if (root_systems[i])
			    rss+=root_systems[i]->size();
		    mapsum=0;
		    if (root_system_map)
		    for (int i=0;i<N+1;i++)
			for (int j=0;j<M+1;j++)
			    mapsum+=root_system_map[idx(i,j)];
		}
		if (fi)
		{
			// calculate moistened zone
			double mw=0,md0=drip,md1=drip;
			for (int i=0;i<N;i++)
			    for (int j=0;j<M;j++)
			    {
				// water head increased more than on 0.1 kPa due to irrigation
				if ((b_U[idx(i,j)]-sb_U[idx(i,j)])>zone_threshold/9.81) 
				{
					double dz=i*dL;
					// calc dx to the nearest emitter
					double mindx=1e300;
					for (int k=0;k<nsprinklers;k++)
					    if (fabs(j*dLx-sprinklers_x[k])<mindx)
						mindx=fabs(j*dLx-sprinklers_x[k]);
					// correct zone
					if (dz<md0)
						md0=dz;
					if (dz>md1)
						md1=dz;
					if (2*mindx>mw)
						mw=2*mindx;
				}
			    }
			// output
			fprintf(fi,"t %g N %d M %d n %d r %g tau %g irr %g tw %g balance %g rlw %g Ev %g Tr %g Perc %g outflow %g moistened_zone w %g d0 %g d1 %g ",T,N,M,n,rv,tau,irr_start_step,pr_twc,balance,pr_rlw,pr_Ev[M/2],pr_Tr[M/2],Perc(),outflow,mw,md0,md1);
			if (rp_form==3)
			    fprintf(fi," rs_size %d rs_sum %g\n",rss,mapsum);
			else
			    fprintf(fi,"\n");
			fflush(fi);
		}
		//printf("t %g N %d M %d n %d r %g tau %g irr %g tw %g balance %g rlw %g Ev %g Tr %g Perc %g Irrrate %g irram %g outflow %g\n",T,N,M,n,rv,tau,irr_start_step,pr_twc,balance,pr_rlw,pr_Ev[M/2],pr_Tr[M/2],Perc(),pr_I*nsprinklers,irr_am,outflow);
		// balance
		double sumet=0.0;
		for (int ii=0;ii<M;ii++)
			sumet+=pr_Ev[ii]+pr_Tr[ii]*mapsum;
		balance+=(-sumet*dLx+Perc()*Lx*perc_k+pr_I*nsprinklers*irr_k-outflow)*tau;
		if ((gamma!=1.0)||(gamma_h!=1.0)) balance=0.0; // no balance check for fractional-order model
		// save solution and free memory
		memcpy(b_U,x,(N+2)*(M+2)*sizeof(double));
		delete [] w;
		delete [] y[0];
		delete [] y[1];
		delete [] rr;
		delete [] v;
		delete [] d;
		delete [] x;
		// increase tau if needed and increase T
		if (n<ls_max_iter*ls_percent)
			if ((irr_start_step==-1)||(tau<max_tau_irr)) // maximal tau restriction when irrigating
				 tau*=ls_mult;
		if (tau>max_tau)
			tau=max_tau;
		T+=tau;
		// update root system
		if ((rp_form==3)&&(no_growth==0))
		{
		    for (int i=0;i<nplants;i++)
		    {
			root_systems[i]->growth(tau);
			root_systems[i]->branching(tau);
		    }
		    calc_root_systems_map();
		}
	}
	// constructor
	H_solver(int rpf,int irr,int bc,double _tau_m,double b,double bh,FILE *log,double _L=3.0,double _v=0.5,double Hb=1e300)
	{
		perc_k=1.0;
		irr_k=0.01;
		vgm_ns=vgm_s0s=vgm_s1s=vgm_as=vgm_h0=vgm_h1=vgm_k=vgm_kx=NULL;
		vgm_nlayers=0;
		irr_am=0.0;
		
		rlw_step=-1;
		I_step=-1;

		nrl=0;
		rlT=rlV=NULL;
		Hbottom=Hb;
		b_U=new double[(N+2)*(M+2)];
		sb_U=new double[(N+2)*(M+2)];
		MM=new double[(N+2)*(M+2)];
		RP=new double[(N+2)*(M+2)];
		rp_mult = new double [M + 2];
		pr_D2x = new double[(N+2)*(M+2)];
		pr_D2y = new double[(N+2)*(M+2)];
		pr_Ev = new double[(M+2)];
		pr_Tr = new double[(M+2)];
		
		isEV=0;
		nirr = 0;
		nperc = 0;

		L = _L; // total depth
		H0 = -800.0 / 1000.0; // initial H	
		v = _v; // root layer depth
		
		gamma=b;
		gamma_h=bh;

		boundary_cond=bc; // 0 -dU/dn=0, 1 - U=H0
		rp_form=rpf;
		no_irr=irr;

		min_wetness=0.385 ; // minimal wetness
		max_wetness=0.54; // needed wetness

		irr_volume = (max_wetness - min_wetness)*v ; // in m
		irr_time = 86400; // time range to apply irrigation
		min_irr_time=0;
		irr_start_step = -1;
		last_irr_step=0;

		dL = L/N;
		Lx=LLx; // lengths of x-domain 
		dLx = Lx/M;
		tau = tau_m=_tau_m;
		T=tau;
		Dg=pow(dL,-1-gamma)/Gamma(2.0-gamma);
		Dg_h=pow(dLx,-1-gamma_h)/Gamma(2.0-gamma_h);
		g2b=Gamma(2.0-gamma);
		g2b_h=Gamma(2.0-gamma_h);

		// intialize plants and sprinklers from global variables
		nplants=_nplants;
		nsprinklers=_nsprinklers;
		plants_x=_plants_x;
		sprinklers_x=_sprinklers_x;
		maxrootR=rootR=_rootR;
		oldrootR=-1;
		sprinklerR=_sprinklerR;
		irr_volume *= Lx/((double)nsprinklers);
		// initial conditions
		for (int i = 0;i < N + 1;i++)
		for (int j = 0;j < M + 1;j++)
			b_U[idx(i, j)] = H0-i*dL;
		if (log)
		{
			fprintf(log,"H_solver(right part %d, no_irr %d,boundary cond %d,tau %g,b %g,bh %g)\n",rpf,irr,bc,tau_m,b,bh);
			printf("H_solver(right part %d, no_irr %d,boundary cond %d,tau %g,b %g,bh %g)\n",rpf,irr,bc,tau_m,b,bh);
			fflush(stdout);
		}
	}
	// desctructor
	~H_solver()
	{
		delete [] b_U;
		delete [] sb_U;
		delete [] MM;
		delete [] RP;
		delete[] rp_mult;
		delete [] pr_D2x;
		delete [] pr_D2y;
		delete [] pr_Ev;
		delete [] pr_Tr;
		if (pr_K) delete [] pr_K;
		if (pr_KX) delete [] pr_KX;
		if (pr_w) delete [] pr_w;
		if (pr_dwdh) delete [] pr_dwdh;
		if (pr_A1) delete [] pr_A1;
		if (pr_A2) delete [] pr_A2;
		if (pr_B1) delete [] pr_B1;
		if (pr_B2) delete [] pr_B2;
		if (pr_R) delete [] pr_R;
		if (sp_precalc) delete [] sp_precalc;
		if (sp_precalc_2) delete [] sp_precalc_2;
	}
};
void grad_func_U(double x,double y,double &d0,double &d1,void *p)
{
    H_solver *s=(H_solver *)p;
    int xi=x/s->dLx, yi=y/s->dL;
    if (xi>=(M+1)) xi=M;
    if (yi>=(N+1)) yi=N;
    d0=-grad_weight*((s->b_U[idx(yi,xi)]-s->b_U[idx(yi,xi+1)])/s->dLx);
    d1=-grad_weight*((s->b_U[idx(yi,xi)]-s->b_U[idx(yi+1,xi)])/s->dL);
//    printf("%g - %g %g %g %g - %d %d - %g %g %g %g\n",s->T,x,y,d0,d1,xi,yi,s->b_U[idx(yi,xi-1)],s->b_U[idx(yi,xi)],s->b_U[idx(yi,xi+1)],s->b_U[idx(yi,yi+2)]);fflush(stdout);
}
// fit parameters (by particle swarm) and solve the problem
// values_{T,Z,F} - values to match
// check_{T,Z,F} - values for checking
// init_{Z,F} - initial values
// EV_T,EV_F[] - evapotranspiration values (EV_T must be sorted)
// nEV_F - number of different EV estimations
// to_fit:
//		bit 1 - gamma
//		bit 2 - L
// 		bit 4 - different EV combination weights
//		bit 5 - filtration k
//		bit 6 - percipitation k 
//		bit 7 - ev_mu
//		bit 8 - irrigation flux k
//		bit 9 - averianov power
//		bit 10 - multiplier for bare soil ET
// pso_{n,o,fi_p,fi_g) - classic particle swarm algorithms' parameters
int adaptive_pso=0; // adaptive pso: o - omega change parameter p , fi_p - fi_p,fi_g changing parameter C, fi_g - vmax calculating parameter N
double restart_prox=0.0; // probability of particles "restart" - random initialization
int init_print=1;
double *vgm_ns,*vgm_s0s,*vgm_s1s,*vgm_as,*vgm_h0,*vgm_h1,*vgm_k,*vgm_kx,*vgm_power,*vgm_specific_storage; // VGM coefs for layered soil
int vgm_nlayers=0;
int nrl=0; // variable root length 
double *rlT,*rlV;
int rpf=1;
double rld=0.5;
double Hbottom=1e300;
char *filename=NULL;
int check_as_fit=0;
H_solver *init_solver(double *p,int to_fit,
	double *init_Z,double *init_F,int ninit,
	double *perc_T, double *perc_A, int nperc,
	double *irr_T, double *irr_A, int nirr,
	double *EV_T,double **EV_F,double *LAI,double ev_mu,int nEV_T,int nEV_F,
	int irr,int bc,double tau_m,double *param_values,int nparam_values)
{
	static double *mask=NULL,*sU=NULL;
	double g=1.0,gh=1.0;
	double L=3.0;
	int offs=1;
	if (to_fit&2) { g=p[offs++]; gh=p[offs++]; } else {  g=param_values[1];  gh=param_values[2]; }
	if (to_fit&4) L=p[offs++]; else L=param_values[3];
	if (g>1.0) g=1.0;
	if (gh>1.0) gh=1.0;
	H_solver *ret=new H_solver(rpf,irr,bc,tau_m,g,gh,NULL,L,rld,Hbottom);
	if (vgm_nlayers!=0)
	{
	    double *cp;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_s0s,vgm_nlayers*sizeof(double));
	    ret->vgm_s0s=cp;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_s1s,vgm_nlayers*sizeof(double));
	    ret->vgm_s1s=cp;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_as,vgm_nlayers*sizeof(double));
	    ret->vgm_as=cp;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_ns,vgm_nlayers*sizeof(double));
	    ret->vgm_ns=cp;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_h0,vgm_nlayers*sizeof(double));
	    ret->vgm_h0=cp;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_h1,vgm_nlayers*sizeof(double));
	    ret->vgm_h1=cp;
	    ret->vgm_nlayers=vgm_nlayers;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_k,vgm_nlayers*sizeof(double));
	    ret->vgm_k=cp;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_kx,vgm_nlayers*sizeof(double));
	    ret->vgm_kx=cp;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_power,vgm_nlayers*sizeof(double));
	    ret->vgm_power=cp;
	    cp=new double[vgm_nlayers];
	    memcpy(cp,vgm_specific_storage,vgm_nlayers*sizeof(double));
	    ret->vgm_specific_storage=cp;
	    if (init_print==1)
	    {
		printf("%d layered soil\n",vgm_nlayers);
		for (int i=0;i<vgm_nlayers;i++)
		    printf("%g %g %g %g %g %g %g %g %g\n",vgm_s0s[i],vgm_s1s[i],vgm_as[i],vgm_ns[i],vgm_h0[i],vgm_h1[i],vgm_k[i],vgm_power[i],vgm_specific_storage[i]);
	    }
	}
	ret->isEV=1;
	ret->EV_T=EV_T;
	ret->EV_F=EV_F;
	ret->nEV_T=nEV_T;
	ret->nEV_F=nEV_F;
	ret->LAI=LAI;
	ret->ev_mu=ev_mu;
	ret->EV_C=new double[nEV_F];
	ret->nrl=nrl;
	ret->rlT=rlT;
	ret->rlV=rlV;
	int poffs=nEV_F;
	if (to_fit&16)
    	    for (int i=0;i<nEV_F;i++)
		ret->EV_C[i]=p[offs++];
	else
	    for (int i=0;i<nEV_F;i++)
		 ret->EV_C[i]=param_values[8+i];
	int nfit_vgm=1;
	if (fit_vgm>=1) nfit_vgm+=1;
	if (fit_vgm>=2) nfit_vgm+=1;
	if (fit_vgm>=3) nfit_vgm+=2;
	if (to_fit & 32)
	{
	    if (vgm_nlayers!=0)
	    {
		for (int i=0;i<vgm_nlayers;i++)
		{
		    ret->vgm_k[i]=ret->vgm_kx[i]=p[offs++];
		    if (fit_vgm>=1) ret->vgm_kx[i]=p[offs++];
		    if (fit_vgm>=2) ret->vgm_s0s[i]=ret->vgm_s1s[i]-p[offs++];
		    if (fit_vgm>=3) 
		    {
			ret->vgm_ns[i]=p[offs++];
			ret->vgm_as[i]=p[offs++];
		    }
		    if (ret->vgm_s0s[i]<0) ret->vgm_s0s[i]=0;
		}
		poffs=nEV_F+nfit_vgm*vgm_nlayers-1;
	    }
	}
	else
	{
	    if (vgm_nlayers!=0)
	    {
		for (int i=0;i<vgm_nlayers;i++)
		{
		    ret->vgm_k[i]=ret->vgm_kx[i]=param_values[8+nEV_F+nfit_vgm*i];
		    if (fit_vgm>=1) ret->vgm_kx[i]=param_values[8+nEV_F+nfit_vgm*i+1];
		    if (fit_vgm>=2) ret->vgm_s0s[i]=ret->vgm_s1s[i]-param_values[8+nEV_F+nfit_vgm*i+2];
		    if (fit_vgm>=3)
		    {
			ret->vgm_ns[i]=param_values[8+nEV_F+nfit_vgm*i+3];
			ret->vgm_as[i]=param_values[8+nEV_F+nfit_vgm*i+4];
		    }
		    if (ret->vgm_s0s[i]<0) ret->vgm_s0s[i]=0;
		}
		poffs=nEV_F+nfit_vgm*vgm_nlayers-1;
	    }
	}
	if (to_fit & 64)
	    ret->perc_k = p[offs++];
	else
	     ret->perc_k=param_values[8+nEV_F+nfit_vgm*vgm_nlayers-1+1];
	if ((irr==0)&&(nirr==0))
	{
	    ret->irr_time=max_irr_time;
	    ret->min_irr_time=min_irr_time;
	    if (to_fit & 64)
	    {
	        ret->min_wetness=p[offs++];
	        ret->max_wetness=p[offs++];
	        ret->irr_volume=p[offs++];
	    }
	    else
	    {
	         ret->min_wetness=param_values[8+nEV_F+nfit_vgm*vgm_nlayers-1+2];
	         ret->max_wetness=param_values[8+nEV_F+nfit_vgm*vgm_nlayers-1+3];
	         ret->irr_volume=param_values[8+nEV_F+nfit_vgm*vgm_nlayers-1+4];
	    }
	    if (ret->min_wetness>ret->max_wetness)
	    {
	        double s=ret->min_wetness;
	        ret->min_wetness=ret->max_wetness;
	        ret->max_wetness=s;
	    }
	}
	if (to_fit&128) ret->ev_mu=p[offs++]; else  ret->ev_mu=param_values[8+nEV_F+nfit_vgm*vgm_nlayers-1+5];
	if (to_fit & 256) ret->irr_k = p[offs++]; else  ret->irr_k=param_values[8+nEV_F+nfit_vgm*vgm_nlayers-1+6];
	if (to_fit & 512) averianov_k =ret->_averianov_k= -p[offs++]; else  averianov_k =ret->_averianov_k=-param_values[8+nEV_F+nfit_vgm*vgm_nlayers-1+7];
	if (to_fit & 1024) Hbottom = ret->Hbottom=p[offs++]; else Hbottom=ret->Hbottom=param_values[8+nEV_F+nfit_vgm*vgm_nlayers-1+8];
	// print initial values
	if ((init_print==1)||(flags&256))
	{
	    printf("G[%g %g] L[%g(%g)] EVC[",g,gh,L,ret->Hbottom);
	    for (int i=0;i<nEV_F;i++)
		printf("%g ",ret->EV_C[i]);
	    printf("] KLIN[%g] PK[%g] IK[%g] KpL[",klin,ret->perc_k,ret->irr_k);
	    if (vgm_nlayers!=0)
		for (int i=0;i<vgm_nlayers;i++)
		    printf(" [%g %g %g %g %g]",ret->vgm_k[i],ret->vgm_kx[i],ret->vgm_s0s[i],ret->vgm_ns[i],ret->vgm_as[i]);
	    printf("]");
	    if (irr==0)
		printf(" I[%g %g %g]",ret->min_wetness,ret->max_wetness,ret->irr_volume);
	    printf(" O[%d %g] MU[%g] AVk[%g]",rpf,ret->v,ret->ev_mu,averianov_k);
	    printf("\n");
	}
	 // smooth initial distribution
	if ((ninit>0)&&(flags&32))
	{
	  int NN=(N+2)*(M+2);
	  double diff=0.0;
	  // create mask with fixed points
	  if (to_fit&1024) // clear mask if Hbottom is changing
	  {
		delete [] mask;
		delete [] sU;
		mask=NULL;
		sU=NULL;
	  }
	  if (mask==NULL) 
	  {
		  mask=new double[(N+2)*(M+2)];
		  sU=new double[(N+2)*(M+2)];
		  double avg=0.0;
		  for (int jj=0;jj<ninit;jj++)
			avg+=init_F[jj];
		  avg/=ninit;
		  for (int i = 0;i < N+2 ;i++)
		  for (int j = 0;j < M+2 ;j++)
		  {
			  mask[idx(i,j)]=0;
			  for (int jj=0;jj<ninit;jj++)
				if ((init_Z[3*jj+0]>=ret->dL*i)&&(init_Z[3*jj+0]<ret->dL*(i+1)))
				if ((init_Z[3*jj+1]>=ret->dLx*j)&&(init_Z[3*jj+1]<ret->dLx*(j+1)))
				{
					if (init_F[jj]!=0.0)
						mask[idx(i,j)]=init_F[jj];
					else
						mask[idx(i,j)]=1e-30;
				}
		  }
		  if (smooth_Hb!=0.0)
		  for (int j = 0;j < M+2 ;j++)
			mask[idx(N+1,j)]=smooth_Hb;
		  // smooth initial distribution
		  for (int i = 0;i < NN ;i++)
			ret->b_U[i]=avg;
		  int ii=0;
		  for (ii=0;ii<stat_maxsmooth;ii++)
		  {
			  diff=0.0;
			  for (int i = 0;i <= N+1 ;i++)
			  for (int j = 0;j <= M+1 ;j++)
			  if (mask[idx(i,j)]==0.0)
			  {
				  double v=ret->b_U[idx(i,j)];
				  double vv[4];
				  if (i==0) vv[0]=v; else vv[0]=ret->b_U[idx(i-1,j)];
				  if (i==(N+1)) vv[1]=v; else vv[1]=ret->b_U[idx(i+1,j)];
				  if (j==0) vv[2]=v; else vv[2]=ret->b_U[idx(i,j-1)];
				  if (j==(M+1)) vv[3]=v; else vv[3]=ret->b_U[idx(i,j+1)];

				  int m=0;
				  double av=(vv[0]+vv[1]+vv[2]+vv[3])/4.0;
				  for (int o=1;o<=3;o++) if (fabs(vv[m]-av)<fabs(vv[o]-av)) m=o;
				  double k[4]={1,1,1,1};
				  k[m]=1.5;
				  double sk=k[0]+k[1]+k[2]+k[3];

				  if (((ret->boundary_cond&1)==1)&&(i==(N-1)))
					vv[1]=ret->Hbottom;
				  if (smooth_f==0.0)
					ret->b_U[idx(i,j)]=(k[0]*vv[0]+k[1]*vv[1]+k[2]*vv[2]+k[3]*vv[3])/sk;
				  else
					ret->b_U[idx(i,j)]=v+smooth_f*(vv[0]+vv[1]-2*v+vv[2]+vv[3]-2*v);
				  diff+=(ret->b_U[idx(i,j)]-v)*(ret->b_U[idx(i,j)]-v);
			  }
			  else
				  ret->b_U[idx(i,j)]=mask[idx(i,j)];
			  if (sqrt(diff)<stat_eps)
				break;
		  }
		  printf("smooth %d - %g\n",ii,diff);
		  memcpy(sU,ret->b_U,NN*sizeof(double));
	  }
	  memcpy(ret->b_U,sU,NN*sizeof(double));
	}
	// set percipitation and irrigation
	ret->nirr = nirr;
	ret->nperc = nperc;
	ret->irr_A = irr_A;
	ret->irr_T = irr_T;
	ret->perc_A = perc_A;
	ret->perc_T = perc_T;
	ret->__precalc();
	init_print=0;
	return ret;
}
void clear_solver(H_solver *ss,int to_fit)
{
	if (to_fit&16)
		delete [] ss->EV_C;
	if (vgm_nlayers!=0)
	{
	    delete [] ss->vgm_s0s;
	    delete [] ss->vgm_s1s;
	    delete [] ss->vgm_as;
	    delete [] ss->vgm_ns;
	    delete [] ss->vgm_k;
	    delete [] ss->vgm_kx;
	    delete [] ss->vgm_h0;
	    delete [] ss->vgm_h1;
	    delete [] ss->vgm_power;
	    delete [] ss->vgm_specific_storage;
	}
	delete ss;
}
double solve_and_test(double *p,int to_fit,
	double *init_Z,double *init_F,int ninit,
	double *perc_T, double *perc_A, int nperc,
	double *irr_T, double *irr_A, int nirr,
	double *EV_T,double **EV_F,double *LAI,double ev_mu,int nEV_T,int nEV_F,
	int irr,int bc,double tau_m,FILE *fi,
	double *values_T,double *values_Z,double *values_F,int nvalues,double *param_values,int nparam_values,double t)
{
	double err=0.0;
	double *vs=new double [nvalues];
	H_solver *ss=init_solver(p,to_fit,init_Z,init_F,ninit, perc_T, perc_A, nperc, irr_T, irr_A, nirr, EV_T,EV_F,LAI,ev_mu,nEV_T,nEV_F,irr,bc,tau_m,param_values,nparam_values);

	for (double tt = 0;tt <= t;tt += ss->tau)
	{
		// save old values
		for (int i=0;i<nvalues;i++)
			if ((values_T[i]>=tt)&&(values_T[i]<(tt+ss->tau)))
			{
				for (int ii=0;ii<=N-1;ii++)
				for (int j=0;j<=M-1;j++)
					if (((ss->dL*ii)<=values_Z[3*i+0])&&((ss->dL*(ii+1))>values_Z[3*i+0]))
					if (((ss->dLx*j)<=values_Z[3*i+1])&&((ss->dLx*(j+1))>values_Z[3*i+1]))
					{
						double k2x=(values_Z[3*i+0]-(ss->dL*ii))/ss->dL;
						double k2y=(values_Z[3*i+1]-(ss->dLx*j))/(ss->dLx);
						double v1,v2;
						v1=(1-k2x)*ss->b_U[idx(ii,j)]+k2x*ss->b_U[idx(ii+1,j)];
						v2=(1-k2y)*ss->b_U[idx(ii,j)]+k2y*ss->b_U[idx(ii,j+1)];
						vs[i]=(v1+v2)/2.0;
						break;
					}
			}
		// solve
		ss->calc_step();
		// add to err
		int nn=0;
		for (int i=0;i<nvalues;i++)
			if ((values_T[i]>=tt)&&(values_T[i]<(tt+ss->tau)))
			{
				double k1=(values_T[i]-tt)/ss->tau;
				for (int ii=0;ii<=N-1;ii++)
				for (int j=0;j<=M-1;j++)
					if (((ss->dL*ii)<=values_Z[3*i+0])&&((ss->dL*(ii+1))>values_Z[3*i+0]))
					if (((ss->dLx*j)<=values_Z[3*i+1])&&((ss->dLx*(j+1))>values_Z[3*i+1]))
					{
						double k2x=(values_Z[3*i+0]-(ss->dL*ii))/ss->dL;
						double k2y=(values_Z[3*i+1]-(ss->dLx*j))/(ss->dLx);
						double v1,v2,v;
						v1=(1-k2x)*ss->b_U[idx(ii,j)]+k2x*ss->b_U[idx(ii+1,j)];
						v2=(1-k2y)*ss->b_U[idx(ii,j)]+k2y*ss->b_U[idx(ii,j+1)];
						v=(v1+v2)/2.0;
						v=(1-k1)*vs[i]+k1*v;
					        if (v>0.0)
						    v=0.0;
						if (fit_func==0) err+=(values_F[i]-v)*(values_F[i]-v)+bounds[13][0]*v*v;
						if (fit_func==1)
						{
						    if(values_F[i]!=0.0) 
							err+=fabs((values_F[i]-v)/values_F[i]);
						    if ((values_F[i]==0.0)&&(v!=0.0))
							err+=1.0;
						}
						nn++;
						break;
					}	
			}
		if (nn && (fit_func==1)) err/=nn;
		if (!finite(ss->b_U[0]))
		{
		     err=1e300;
		     break;
		}
	}
	clear_solver(ss,to_fit);
	delete [] vs;
	if (!finite(err)) err = 1e300;
	return err;
}
double rnd()
{
	return ((rand() % 10000) / (10000.0 - 1.0));
}
void init_particle(double *particle, int to_fit,int nEV_F,int idx,int size,int nirr)
{
	int s = 0;	
	// gamma,gamma_h
	if (to_fit & 2)
	{
	    particle[1+s++] = bounds[1][0] + (bounds[1][1]-bounds[1][0])*rnd();
	    particle[1+s++] = bounds[1][0] + (bounds[1][1]-bounds[1][0])*rnd();
	}
	if (to_fit & 4)	particle[1+s++] = bounds[2][0] + (bounds[2][1]-bounds[2][0])*rnd();
	// EVT coeffs from 0 to 1
	if (to_fit & 16)
	{
	    double sum=0;
	    int ss=s;
	    do
	    {
		s=ss;
		sum=0;
		for (int i = 0;i < nEV_F;i++)
		{
			particle[1 + s++] = bounds[17][0]+(bounds[11][0]-bounds[17][0])*rnd();
			sum+=particle[s];			
		}
		if (sum<=bounds[11][0])
		    break;
	    }
	    while (1);
	}
	// filtration k
	if (to_fit & 32)
	{
	    if (vgm_nlayers!=0)
		for (int i=0;i<vgm_nlayers;i++)
		{
		    particle[1 + s++] = bounds[7][0] + (bounds[7][1]-bounds[7][0])*rnd();
		    if (fit_vgm>=1) particle[1 + s++] = bounds[7][0] + (bounds[7][1]-bounds[7][0])*rnd();
		    if (fit_vgm>=2) particle[1 + s++] = vgm_s1s[i]*rnd();
		    if (fit_vgm>=3)
		    {
			particle[1 + s++] = bounds[5][0] + (bounds[5][1]-bounds[5][0])*rnd();
			particle[1 + s++] = bounds[4][0] + (bounds[4][1]-bounds[4][0])*rnd();
		    }
		}
	    else
		particle[1 + s++] = bounds[7][0] + (bounds[7][1]-bounds[7][0])*rnd();
	}
	// percipitation k
	if (to_fit & 64)
	{
	    particle[1 + s++] = bounds[17][1]+(bounds[11][1]-bounds[17][1])*rnd();
	    if ((irr==0)&&(nirr==0))
	    {
		particle[1 + s++] = bounds[8][0] + (bounds[8][1]-bounds[8][0])*rnd();
		particle[1 + s++] = bounds[9][0] + (bounds[9][1]-bounds[9][0])*rnd();
		particle[1 + s++] = bounds[10][0] + (bounds[10][1]-bounds[10][0])*rnd();
	    }
	}
	// ev_mu
	if (to_fit & 128) particle[1+s++] = bounds[14][0] + (bounds[14][1]-bounds[14][0])*rnd();
	// irr_k
	if (to_fit & 256) particle[1+s++] = bounds[11][1]*rnd();
	// averianov k
	if (to_fit & 512) particle[1+s++] = bounds[15][1]*rnd();
	// Hbottom
	if (to_fit & 1024) particle[1+s++] = bounds[16][0]+(bounds[16][1]-bounds[16][0])*rnd();
	for (int i = 0;i < s;i++)
		particle[i + 1 + 2 * s] = particle[i + 1];
	for (int i = 0;i<s;i++)
		particle[i + 1 + s] = 0.0;
}
// saves full solution in some js file
void save_full_solution(H_solver *ss)
{
    char fname2[1024];
    FILE *fi2;
    if (flags&128) return; // no js generation flag
    sprintf(fname2,"res_%9.9d.js",(int)ss->T);
    if((fi2=fopen(fname2,"wt"))!=NULL)
    {
	fprintf(fi2,"data_res.push([%g,[]]);\n",ss->T);
	// process
	for (int i=0;i<N-1;i++)
	for (int j=0;j<M-1;j++)
		fprintf(fi2,"data_res[data_res.length-1][1].push([%g,%g,%g,%g,%g,%g]);\n",i*ss->dL,(j+1)*ss->dLx,ss->b_U[idx(i,j+1)],ss->RpF(i,j+1),((i==0)?ss->Uc(j+1):0.0),ss->transp(i,j+1));
	fclose(fi2);
    }
}
void fit_and_solve(double t,double save_tau,double out_tau,int irr,int bc,double tau_m,
	double *init_Z, double *init_F, int ninit,
	double *perc_T, double *perc_A, int nperc,
	double *irr_T, double *irr_A, int nirr,
	double *values_T,double *values_Z,double *values_F,int nvalues,
	double *check_T,double *check_Z,double *check_F,int ncheck,
	double *EV_T,double **EV_F,double *LAI,double ev_mu,int nEV_T,int nEV_F,
	int to_fit,
	int pso_n,double pso_o,double pso_fi_p,double pso_fi_g,double pso_eps,int pso_max_iter,
	int fit,double *param_values,int nparam_values,char *et_file,char *pv_file)
{
	double ad_pso_p,ad_pso_C; // adaptive PSO parameters
	double ad_v0,*ad_vmax;
	double *ad_pso_fi_p,*ad_pso_fi_g; // per-variable fi_p, fi_g for adaptive PSO
	int size=0;
	double **particles;
	int best;
	double *best_p;
	int iter=0;
	char fn[2048],rfn[2048],bfn[2048];
	if (et_file)
	{
	    char *e=et_file;
	    while (e[0])
	    {
		if (e[0]=='.')
		    e[0]=0;
		e++;
	    }
	}
	if (pv_file)
	{
	    char *e=pv_file;
	    while (e[0])
	    {
		if (e[0]=='.')
		    e[0]=0;
		e++;
	    }
	}	
	sprintf(fn,"log_%s_%s.txt",(et_file?et_file:""),(pv_file?pv_file:""));
	sprintf(rfn,"res_%s_%s.txt",(et_file?et_file:""),(pv_file?pv_file:""));
	sprintf(bfn,"pv_best_%s_%s.txt",(et_file?et_file:""),(pv_file?pv_file:""));
	if (filename)
	{
	    sprintf(fn,"log_%s.txt",filename);
	    sprintf(rfn,"res_%s.txt",filename);
	    sprintf(bfn,"pv_best_%s.txt",filename);
	}
	FILE *fi1 = fopen(fn, "wt");
	// number of variables to optimize
	if (to_fit&2) size+=2;
	if (to_fit&4) size++;
	if (to_fit&16) size+=nEV_F;
	if (to_fit & 32)
	{
	    if (vgm_nlayers!=0)
	    {
		int nfit_vgm=1;
		if (fit_vgm>=1) nfit_vgm+=1;
		if (fit_vgm>=2) nfit_vgm+=1;
		if (fit_vgm>=3) nfit_vgm+=2;
		size+=nfit_vgm*vgm_nlayers;
	    }
	    else	
		size++;
	}
	if (to_fit & 64) { if ((irr==1)||(nirr)) size++; else size+=4; }
	if (to_fit & 128) size++;
	if (to_fit & 256) size++;
	if (to_fit & 512) size++;
	if (to_fit & 1024) size++;
	if (size==0) return;
	best_p=new double[2*size+1];
	if (et_file)
	{
    	    printf("et_file - %s\n",et_file);
    	    fprintf(fi1,"#et_file - %s\n",et_file);
    	}
	if (pv_file)
	{
    	    printf("pv_file - %s\n",pv_file);
    	    fprintf(fi1,"#pv_file - %s\n",pv_file);
    	}
	if (fit==1) // fitting
	{
	    for (int i=0;i<bounds_size;i++)
	    {
		fprintf(fi1,"#bounds[%d]={%g,%g}\n",i,bounds[i][0],bounds[i][1]);
		printf("bounds[%d]={%g,%g}\n",i,bounds[i][0],bounds[i][1]);
	    }
	    particles=new double *[pso_n];
    	    for (int i=0;i<pso_n;i++)
		particles[i]=new double[3*size+2]; // particles[0] contains f value, particles[1+size].. contain velocities, particles[1+2*size]... contain particle's best known position, particles[1+3*size] contains best known f value
    	    // initialize
    	    int nit=0;
    	    do
    	    {
		init_particle(particles[0], to_fit,nEV_F,0,size,nirr);
		particles[0][0]=particles[0][1+3*size]=solve_and_test(particles[0],to_fit,init_Z,init_F,ninit,perc_T,perc_A,nperc,irr_T,irr_A,nirr,EV_T,EV_F,LAI,ev_mu,nEV_T,nEV_F,irr,bc,tau_m,fi1,values_T,values_Z,values_F,nvalues,param_values,nparam_values,t);
		printf(".");fflush(stdout);
		if ((nit++)==100)
		    break;
	    }
	    while ((particles[0][0]>=1e10)&&(!(flags&8)));
	    best=0;
	    fprintf(fi1, "#initial 0 - %g\n", particles[0][0]);
	    printf("initial 0 - %g\n", particles[0][0]);
	    for (int i=1;i<pso_n;i++)
	    {
	        int nit=0;
    	        do
    	        {
		    init_particle(particles[i], to_fit,nEV_F,i,size,nirr);
		    particles[i][0]=particles[i][1+3*size]=solve_and_test(particles[i],to_fit,init_Z,init_F,ninit, perc_T, perc_A, nperc, irr_T, irr_A, nirr, EV_T,EV_F,LAI,ev_mu,nEV_T,nEV_F,irr,bc,tau_m,fi1,values_T,values_Z,values_F,nvalues,param_values,nparam_values,t);
		    printf(".");fflush(stdout);
		    if ((nit++)==100)
			break;
	        }
	        while ((particles[i][0]>=1e10)&&(!(flags&8)));
		if (particles[i][0]<particles[best][0])
			best=i;
		fprintf(fi1, "#initial %d - %g\n",i, particles[i][0]);
		printf("initial %d - %g\n",i, particles[i][0]);
	    }
	    // save best known position
	    for (int j=0;j<=size;j++)
		best_p[j]=particles[best][j];
	    fprintf(fi1, "#initial best: ");
	    printf("initial best: ");
	    for (int j = 0;j <= size;j++)
	    {
		fprintf(fi1, "%2.2g ", best_p[j]);
		printf("%2.2g ", best_p[j]);
	    }
	    fprintf(fi1, "\n");
	    printf("\n");
	    fflush(stdout);
		// adaptive PSO - calc initial and max velocity
		if (adaptive_pso)
		{
			double *xm=new double [2*size]; // min-max x in initial population
			double *ad_v0v=new double[size];
			ad_v0=0;
			for (int i=0;i<size;i++) ad_v0v[i]=0;
			ad_vmax=new double[size];
			ad_pso_fi_p=new double[size];
			ad_pso_fi_g=new double[size];
			for (int i=0;i<pso_n;i++)
				for (int j=0;j<size;j++)
				{
					ad_v0+=particles[i][j+1]*particles[i][j+1];
					ad_v0v[j]+=particles[i][j+1]*particles[i][j+1];
					if (i==0)
						xm[2*j+0]=xm[2*j+1]=particles[i][j+1];
					else
					{
						if (particles[i][j+1]<xm[2*j+0]) xm[2*j+0]=particles[i][j+1];
						if (particles[i][j+1]>xm[2*j+1]) xm[2*j+1]=particles[i][j+1];
					}
				}
			ad_v0=sqrt(ad_v0);
			ad_v0/=pso_n;
			for (int i=0;i<size;i++)
			{
			    ad_vmax[i]=(xm[2*i+1]-xm[2*i+0])/pso_fi_g;
				if (ad_vmax[i]==0.0) ad_vmax[i]=1.0/pso_fi_g;
			    ad_v0v[i]=sqrt(ad_v0v[i]);
			    ad_v0v[i]/=pso_n;
			}
			delete [] xm;
			ad_pso_p=pso_o; ad_pso_C=pso_fi_p;
			pso_o=1.0;
			for (int i=0;i<size;i++)
			{
			    ad_pso_fi_p[i]=ad_pso_C*ad_v0v[i]/ad_vmax[i];
			    ad_pso_fi_g[i]=ad_pso_C*(1-ad_v0v[i]/ad_vmax[i]);
			    if (ad_pso_fi_p[i]>2.0) ad_pso_fi_p[i]=2.0;
			    if (ad_pso_fi_g[i]>2.0) ad_pso_fi_g[i]=2.0;
    			    if (ad_pso_fi_p[i]<0.0) ad_pso_fi_p[i]=0.0;
			    if (ad_pso_fi_g[i]<0.0) ad_pso_fi_g[i]=0.0;
			}
			delete [] ad_v0v;
			printf("adaptive PSO: ||v0||=%g, initial PSO params values %g",ad_v0,pso_o);
			for (int i=0;i<size;i++)
			    printf(",(%g->%g,%g)",ad_vmax[i],ad_pso_fi_p[i],ad_pso_fi_g[i]);
			printf("\n");
		}
	    // process
	    if (pso_max_iter>=1)
	    do
	    {
		// adaptive PSO - change PSO parameters
		if (adaptive_pso)
		{
			double ve=ad_v0*exp(-(2*(iter+1)/(float)(pso_max_iter-iter))*(2*(iter+1)/(float)(pso_max_iter-iter)));
			double vavg=0.0;
			double *vavgc=new double[size];
			for (int j=0;j<size;j++)
			    vavgc[j]=0;
			for (int i=0;i<pso_n;i++)
				for (int j=0;j<size;j++)
				{
					vavg+=particles[i][j+1+size]*particles[i][j+1+size];
					vavgc[j]+=particles[i][j+1+size]*particles[i][j+1+size];
				}
			vavg=sqrt(vavg);
			vavg/=pso_n;
			for (int j=0;j<size;j++)
			{
			    vavgc[j]=sqrt(vavgc[j]);
			    vavgc[j]/=pso_n;
			}
			// change omega
			if (vavg>ve)
				pso_o/=ad_pso_p;
			if (vavg<ve)
				pso_o*=ad_pso_p;
			if (pso_o>2) pso_o=2;
			if (pso_o<0) pso_o=0;
			// change fi_p,fi_g
			for (int i=0;i<size;i++)
			{
			    ad_pso_fi_p[i]=ad_pso_C*vavgc[i]/ad_vmax[i];
			    ad_pso_fi_g[i]=ad_pso_C*(1-vavgc[i]/ad_vmax[i]);
			    if (ad_pso_fi_p[i]>2.0) ad_pso_fi_p[i]=2.0;
			    if (ad_pso_fi_g[i]>2.0) ad_pso_fi_g[i]=2.0;
			    if (ad_pso_fi_p[i]<0.0) ad_pso_fi_p[i]=0.0;
			    if (ad_pso_fi_g[i]<0.0) ad_pso_fi_g[i]=0.0;
			}
			printf("adaptive PSO - vavg %g ve %g params values %g",vavg,ve,pso_o);
			for (int i=0;i<size;i++)
			    printf(",(%g->%g,%g)",vavgc[i],ad_pso_fi_p[i],ad_pso_fi_g[i]);
			printf("\n");
			delete [] vavgc;
		}
		for (int i=0;i<pso_n;i++)
		{
			// update velocity
			for (int j=0;j<size;j++)
			{
				double rp=(rand()%10000)/10000.0;
				double rg=(rand()%10000)/10000.0;
				if (adaptive_pso) // percoordinate fi_p,fi_g
				{
				    pso_fi_p=ad_pso_fi_p[j];
				    pso_fi_g=ad_pso_fi_g[j];
				}
				particles[i][j+1+size]=pso_o*particles[i][j+1+size]+pso_fi_p*rp*(particles[i][j+1+2*size]-particles[i][j+1])+pso_fi_g*rg*(best_p[j+1]-particles[i][j+1]);				    
			}
			// update position
			for (int j=0;j<size;j++)
				particles[i][1+j]+=particles[i][j+1+size];
			// restart particle
			double rp=(rand()%10000)/10000.0;			
			if (rp<((pso_o<1)?(1-pso_o):0.0)*restart_prox)
			{
			    init_particle(particles[i], to_fit,nEV_F,i,size,nirr);
			    printf("r");
			}
			// assure all params are >0
			for (int j=0;j<size;j++)
				if (particles[i][1+j]<0) particles[i][1+j]=0;
			// calc f value
			particles[i][0]=solve_and_test(particles[i],to_fit,init_Z,init_F,ninit, perc_T, perc_A, nperc, irr_T, irr_A, nirr, EV_T,EV_F,LAI,ev_mu,nEV_T,nEV_F,irr,bc,tau_m,fi1,values_T,values_Z,values_F,nvalues,param_values,nparam_values,t);
			// update bests
			if (particles[i][0]<particles[i][1+3*size])
			{
				for (int j=0;j<size;j++)
					particles[i][j+1+2*size]=particles[i][j+1];
				particles[i][1+3*size]=particles[i][0];
			}
			if (particles[i][0]<best_p[0])
			{
				for (int j=0;j<size;j++)
					best_p[j+1]=particles[i][j+1];
				best_p[0]=particles[i][0];
			}
			fflush(stdout);
		}
		// check best-worst difference
		double max = 0.0;
		double avg=0.0;
		for (int i = 0;i < pso_n;i++)
		{
		    if (max < particles[i][0])
			max = particles[i][0];
		    avg+= particles[i][0];
		}
		avg/=pso_n;
		fprintf(fi1, "#%d avg %g best: ", iter,avg);
		printf("%d avg %g best: ", iter,avg);
		for (int j = 0;j <= size;j++)
		{
		    fprintf(fi1, "%g ", best_p[j]);
		    printf("%g ", best_p[j]);
		}
		fprintf(fi1, "\n");
		printf("\n");
		if ((max - best_p[0]) < pso_eps)
		    break;
		iter++;
	    }
	    while ((iter<pso_max_iter)&&(best_p[0]>pso_eps));
	}
	if (fit==2) // only solution
	{
	    for (int i=0;i<size;i++)
		best_p[i+1]=param_values[i];
	    to_fit=0;
	}
	// solve with best parameters values
	H_solver *ss;
	init_print=1;
	if (fit==1)
	{
	    ss=init_solver(best_p,to_fit,init_Z,init_F,ninit, perc_T, perc_A, nperc, irr_T, irr_A, nirr, EV_T,EV_F,LAI,ev_mu,nEV_T,nEV_F,irr,bc,tau_m,param_values,nparam_values);
	    // save results
	    FILE *fi3 = fopen(bfn, "wt");
	    fprintf(fi3,"1 %g %g %g 0 0 0 0 ",ss->gamma,ss->gamma_h,ss->L);
	    for (int i=0;i<nEV_F;i++)
		fprintf(fi3,"%g ",ss->EV_C[i]);
	    if (vgm_nlayers!=0)
	    {
		for (int i=0;i<vgm_nlayers;i++)
		    fprintf(fi3,"%g %g %g ",ss->vgm_k[i],ss->vgm_kx[i],ss->vgm_s1s[i]-ss->vgm_s0s[i]);
		fprintf(fi3,"%g",ss->perc_k);
	    }
	    fprintf(fi3," %g %g %g",ss->min_wetness,ss->max_wetness,ss->irr_volume);
	    fprintf(fi3," %g %g %g %g",ss->ev_mu,ss->irr_k,-averianov_k,Hbottom);
	    fprintf(fi3,"\n");
	    fclose(fi3);
	}
	if (fit==2)
	    ss=init_solver(best_p,to_fit,init_Z,init_F,ninit, perc_T, perc_A, nperc, irr_T, irr_A, nirr, EV_T,EV_F,LAI,ev_mu,nEV_T,nEV_F,irr,bc,tau_m,param_values,nparam_values);
	int nsave = 1;
	int nout=1;
	int d,oldd=-1;
	int logn=0,maxlogn=0;
	double err=0.0,srel_err=0.0;
	std::map<std::pair<double,double>,std::pair<double,int>> errs_per_xy;
	std::map<std::pair<double,double>,std::pair<double,int>>::iterator errs_i;
	int nerr=0;	
	double *vs=new double [ncheck];
	FILE *fi2;
	fi2 = fopen(rfn, "wt");
	if ((fit==1)&&(check_as_fit==0))
	    if (t<check_T[ncheck-1])
		t=check_T[ncheck-1];
	if (check_as_fit>1) t=check_as_fit;
	fit=2;
	for (double tt = 0;tt <= t;tt += ss->tau)
	{
		// save result
		if (tt ==0)
		{
			for (int i = 0;i < N + 1;i++)
				fprintf(fi2, "%g %g %g %g %g %g\n", tt,(double)i*ss->dL,ss->b_U[idx(i,(M+1)/2)],ss->pr_w[idx(i,(M+1)/2)],ss->pr_dwdh[idx(i,(M+1)/2)],ss->pr_K[idx(i,(M+1)/2)]);
			fflush(fi2);
			save_full_solution(ss);
			nsave++;
		}
		// save old 
		int highest_check_T=0;
		for (int i=0;i<ncheck;i++)
			if (check_T[i]<(tt+ss->tau))
			{
				for (int ii=0;ii<=N-1;ii++)
				for (int j=0;j<=M-1;j++)
					if (((ss->dL*ii)<=check_Z[3*i+0])&&((ss->dL*(ii+1))>check_Z[3*i+0]))
					if (((ss->dLx*j)<=check_Z[3*i+1])&&((ss->dLx*(j+1))>check_Z[3*i+1]))
					{
						double k2x=(check_Z[3*i+0]-(ss->dL*ii))/ss->dL;
						double k2y=(check_Z[3*i+1]-(ss->dLx*j))/(ss->dLx);
						double v1,v2;
						v1=(1-k2x)*ss->b_U[idx(ii,j)]+k2x*ss->b_U[idx(ii+1,j)];
						v2=(1-k2y)*ss->b_U[idx(ii,j)]+k2y*ss->b_U[idx(ii,j+1)];
						vs[i]=(v1+v2)/2.0;
						if (check_T[i]>highest_check_T) highest_check_T=check_T[i];
						break;
					}
			}
		// solve	
		ss->calc_step((fit==1)?NULL:((tt>=nsave*save_tau)?fi1:NULL));
		// add to err
		int is_out = 0;
		logn=0;
		    for (int ii=0;ii<=N-1;ii++)
		    {
			for (int i=0;i<ncheck;i++)
			    if ((((flags&512)==0)&&((check_T[i]>=tt)&&(check_T[i]<(tt+ss->tau))))||((flags&512)&&(check_T[i]==highest_check_T)))
			    {
				double k1=(check_T[i]-tt)/ss->tau;
			        for (int j=0;j<=M-1;j++)
			    	   if (((ss->dL*ii)<=check_Z[3*i+0])&&((ss->dL*(ii+1))>check_Z[3*i+0]))
				    if (((ss->dLx*j)<=check_Z[3*i+1])&&((ss->dLx*(j+1))>check_Z[3*i+1]))
				    {
				    		double k2x=(check_Z[3*i+0]-(ss->dL*ii))/ss->dL;
						double k2y=(check_Z[3*i+1]-(ss->dLx*j))/(ss->dLx);
						double v1,v2,v;
						v1=(1-k2x)*ss->b_U[idx(ii,j)]+k2x*ss->b_U[idx(ii+1,j)];
						v2=(1-k2y)*ss->b_U[idx(ii,j)]+k2y*ss->b_U[idx(ii,j+1)];
						v=(v1+v2)/2.0;
						if ((flags&512)==0)
							v=(1-k1)*vs[i]+k1*v;
						if (v>0.0)
							v=0.0;
						err+=(check_F[i]-v)*(check_F[i]-v);
						// save error for xy in H
						if ((errs_i=errs_per_xy.find(std::pair<double,double>(check_Z[3*i+0],check_Z[3*i+1])))!=errs_per_xy.end())
						{
						    errs_i->second.first+=(check_F[i]-v)*(check_F[i]-v);
						    errs_i->second.second++;
						}
						else
						    errs_per_xy.insert(std::pair<std::pair<double,double>,std::pair<double,int> >(
							std::pair<double,double>(check_Z[3*i+0],check_Z[3*i+1]),
							std::pair<double,int>((check_F[i]-v)*(check_F[i]-v),1)));

						if ((check_F[i]==0)&&(v!=0.0))
							srel_err+=1.0;
						if (check_F[i]!=0)
						{
							if (fabs((check_F[i]-v)/check_F[i])<1)
								srel_err+=fabs((check_F[i]-v)/check_F[i]);
							else 
								srel_err+=1.0;
						}
						nerr++;
//						if (fit==1)
						{
						    if (is_out == 0)
							if ((flags&512)==0)
							    fprintf(fi1, "%g ", check_T[i]);
							else
							    fprintf(fi1, "%g ", tt);
						    fprintf(fi1, "%g %g %g %g ", check_Z[3*i+0],  check_Z[3*i+1],check_F[i], v);
						}
						logn++;
						is_out = 1;
				    }
			    }
			}
		if (logn>maxlogn) maxlogn=logn;
//		if (fit==1)
		if (is_out)
			fprintf(fi1,"rlw %g tw %g twc %g sum_err %g avg_err %g avg_rel_err %g\n",
				     ss->pr_rlw,ss->pr_fw,ss->pr_twc,
				     err,(nerr?err/nerr:0),(nerr?srel_err/nerr:0));
		// save result
		if (tt >= nsave*save_tau)
		{
			for (int i = 0;i < N + 1;i++)
				fprintf(fi2, "%g %g %g %g %g %g\n", tt,(double)i*ss->dL,ss->b_U[idx(i,(M+1)/2)],ss->pr_w[idx(i,(M+1)/2)],ss->pr_dwdh[idx(i,(M+1)/2)],ss->pr_K[idx(i,(M+1)/2)]);
			fflush(fi2);
			save_full_solution(ss);
			if ((ncheck==0))//&&(fit==1))
			{
				fprintf(fi1,"%g rlw %g tw %g twc %g sum_err %g avg_err %g avg_rel_err %g\n",
				     tt, ss->pr_rlw,ss->pr_fw,ss->pr_twc,
				     err,(nerr?err/nerr:0),(nerr?srel_err/nerr:0));
			}
			nsave++;
		}
		if (!finite(ss->b_U[0]))
		{
		    err=1e300;
		    break;
		}
	}
	if (ss->irr_am!=0)
		if (fit!=1)
		{
		    if ((ss->irr_start_step0)||(!finite(ss->b_U[0])))
			printf("irrigation amount - 1000 %g m^3 (%g-%g)\n", ss->irr_am, ss->irr_start_step, ss->T);
		    else
			printf("irrigation amount - %g m^3 (%g-%g)\n", ss->irr_am, ss->irr_start_step, ss->T);
		}
//	if (fit==1)
	{
	    fprintf(fi1,"#checking err: %g\n",err);
	    for (errs_i = errs_per_xy.begin(); errs_i != errs_per_xy.end(); errs_i++)
		fprintf(fi1,"#avg abs err in H for (%g,%g) = %g\n",errs_i->first.first,errs_i->first.second,sqrt(errs_i->second.first/errs_i->second.second));
	}
	if (rpf==3)
	{
		FILE *fi=fopen(root_systems_out_file,"wt");
		for (int i=0;i<ss->nplants;i++)
			ss->root_systems[i]->write(fi);
		fclose(fi);
	}
	clear_solver(ss,to_fit);
	fclose(fi1);
}
/// main/////////////////
int main(int argc,char **argv)
{
	double Tm = 64000.0;
	double Sm = 3600.0;
	double Om = 3600.0;
	int bc=0;
	double tau_m=1000.0/20.0;
	char *et_file=NULL;
	char *pv_file=NULL;
	double _ev_mu=1e100;
	int _pso_n=-1;
	double _pso_o=1e100;
	double _pso_fi_p=1e100;
	double _pso_fi_g=1e100;
	double _pso_eps=1e100;
	int _pso_max_iter=-1;
	int _to_fit=-1;
	char *fi_file=NULL;
	double *exclude_depths=NULL;
	int n_ex_depths=0;
	double *fit_list=NULL;
	int n_fit_list=0;
	int start_time=0;
	init_paramsets();
	if (argc==1)
	{
		printf("Tm : end time\n");
		printf("Sm : save time\n");
		printf("Om : log output time\n");
		printf("NB : number of grid cells on depth\n");
		printf("M : number of grid cells on width\n");
		printf("I : 1 - irrigation disabled\n");
		printf("B : boundary condition - 0 -dU/dn=0, 1 - U=H0 - bit 0 - on lower boundary, bit 1 - on upper\n");
		printf("Tau : basic time step\n");
		printf("Fit: 1 - fit and solve, 2 - solve with given values\n");
		printf("fl: additional flags ");
		printf("stMS: maximal number of iterations for initial value computations\n");
		printf("stEPS: eps for steady state initial value computations\n");
		printf("bounds_i_j: bounds for parameter changes in PSO search:\n\t\t i_0 - min, i_1 - max\n\t\t 1 - beta\n\t\t\
			7 - filtration coefficient\n\t\t\
			8,9,10 - min/max wetness and irrigation amount for irrigation application\n\t\t\
			11 - EV coefs max and percipitation coef max\n\t\t\
			14 - Ev_mu bounds\n\t\t");
		printf("vgm: name of file with layers soil VGM coefficients in form [s0,s1,n,a,h0,h1]\n");
		printf("et_file: name of file with ET values\n");
		printf("pv_file: name of file with fixed parameters values\n");
		printf("rpf: root density function form\n");
		printf("rld: root layer depth\n");
		printf("rl_file: file which contains the pairs <time,root length>\n");
		printf("ev_mu,pso_n,pso_o,pso_fi_p,pso_fi_g,pso_eps,pso_max_iter: used to override fitting params\n");
		printf("filename: output filenames prefix\n");
		printf("Hbottom: value of H at the bottom of the domain for boundary condition=1\n");
		printf("fi_file: fit input filename (default: fit_input.txt)\n");
		printf("klin: linear coefficient for Kfilt change in time ((m/s)/s)\n");
		printf("drip: subsurface drip irrigation depth (m)\n");
		printf("StT: starting time, s\n");
		printf("APSO: use adaptive PSO\n");
		printf("restart_prox: possibility of restarting a particle\n");
		printf("exclude_depths: depths of sensors that should not be used\n");
		printf("fit_list: list of time points to use while fitting\n");
		exit(0);
	}
	for (int i=1;i<argc;i+=2)
	{
		if (strcmp(argv[i],"StT")==0)
			start_time=atoi(argv[i+1]);
		if (strcmp(argv[i],"Tm")==0)
			Tm=atof(argv[i+1]);
		if (strcmp(argv[i],"Sm")==0)
			Sm=atof(argv[i+1]);
		if (strcmp(argv[i],"Om")==0)
			Om=atof(argv[i+1]);
		if (strcmp(argv[i],"NB")==0)
			N=atoi(argv[i+1]);
		if (strcmp(argv[i],"I")==0)
			irr=atoi(argv[i+1]);
		if (strcmp(argv[i],"B")==0)
			bc=atoi(argv[i+1]);
		if (strcmp(argv[i],"Tau")==0)
			tau_m=atof(argv[i+1]);
		if (strcmp(argv[i],"M")==0)
			M=atoi(argv[i+1]);
		if (strcmp(argv[i],"Fit")==0)
			fit=atoi(argv[i+1]);
		if (strcmp(argv[i], "fl") == 0)
			flags = atoi(argv[i + 1]);
		if (strcmp(argv[i], "stHb") == 0)
			smooth_Hb = atof(argv[i + 1]);
		if (strcmp(argv[i], "stMS") == 0)
			stat_maxsmooth = atoi(argv[i + 1]);
		if (strcmp(argv[i], "stF") == 0)
			smooth_f = atof(argv[i + 1]);
		if (strcmp(argv[i], "stEPS")== 0)
			stat_eps = atof(argv[i + 1]);
		if (strcmp(argv[i], "rpf")== 0)
			rpf = atoi(argv[i + 1]);
		if (strcmp(argv[i], "rld")== 0)
			max_rld=rld = atof(argv[i + 1]);
		if (strstr(argv[i],"bounds")!=NULL)
		{
		    char *str=argv[i];
		    while ((str[0]!='_')&&(str[0])) str++;
		    if (str[0]=='_')
		    {
			str++;
			char *str2=str;
			while ((str2[0]!='_')&&(str2[0])) str2++;
			if (str2[0]=='_')
			{
			    str2[0]=0;
			    str2++;
			    int ii=atoi(str);
			    int j=atoi(str2);
			    double v=atof(argv[i+1]);
			    if ((j>=0)&&(j<=1))
				if ((ii>=0)&&(ii<=bounds_size-1))
				    bounds[ii][j]=v;
			}
		    }
		}
		if (strstr(argv[i],"vgm")!=NULL)
		{
		    FILE *fi=fopen(argv[i+1],"rt");
		    if (fi)
		    {
			char str[1024];
			int n;
			while (fgets(str,1024,fi)) vgm_nlayers++;
			vgm_ns=new double[vgm_nlayers];
			vgm_as=new double[vgm_nlayers];
			vgm_s0s=new double[vgm_nlayers];
			vgm_s1s=new double[vgm_nlayers];
			vgm_h0=new double[vgm_nlayers];
			vgm_h1=new double[vgm_nlayers];
			vgm_k=new double[vgm_nlayers];
			vgm_kx=new double[vgm_nlayers];
			vgm_power=new double[vgm_nlayers];
			vgm_specific_storage=new double[vgm_nlayers];
			vgm_nlayers=0;
			fseek(fi,0,SEEK_SET);			
			do
			{
			    char str[2048];
			    if (fgets(str,2048,fi)==NULL) break;
			    n=sscanf(str,"%lg %lg %lg %lg %lg %lg %lg %lg %lg",vgm_s0s+vgm_nlayers,vgm_s1s+vgm_nlayers,vgm_ns+vgm_nlayers,vgm_as+vgm_nlayers,vgm_h0+vgm_nlayers,vgm_h1+vgm_nlayers,vgm_k+vgm_nlayers,vgm_power+vgm_nlayers,vgm_specific_storage+vgm_nlayers);
			    vgm_kx[vgm_nlayers]=vgm_k[vgm_nlayers];
			    if (n==8) vgm_specific_storage[vgm_nlayers]=_vgm_specific_storage;
			    if (n>=8) vgm_nlayers++;
			}
			while (1);
		    }		    
		}
		if (strstr(argv[i],"rl_file")!=NULL)
		{
		    FILE *fi=fopen(argv[i+1],"rt");
		    if (fi)
		    {
			char str[1024];
			nrl=0;
			while (fgets(str,1024,fi)) nrl++;
			rlT=new double[nrl];
			rlV=new double[nrl];
			nrl=0;
			fseek(fi,0,SEEK_SET);			
			while (fscanf(fi,"%lg %lg\n",rlT+nrl,rlV+nrl)==2)
			    if ((rlT[nrl]-=start_time)>=0)
				nrl++;
			if (nrl) 
			{
				max_rld=rlV[0];
				for (int i=1;i<nrl;i++) if (max_rld<rlV[i]) max_rld=rlV[i];
				printf("root layer depth %d point read (last %g %g) max %g\n",nrl,rlT[nrl-1],rlV[nrl-1],max_rld);
			}
		    }		    
		}
		if (strstr(argv[i],"et_file")!=NULL)
		    et_file=argv[i+1];
		if (strstr(argv[i],"pv_file")!=NULL)
		    pv_file=argv[i+1];
		if (strstr(argv[i],"ev_mu")!=NULL)
		    _ev_mu=atof(argv[i+1]);
		if (strstr(argv[i],"pso_n")!=NULL)
		    _pso_n=atoi(argv[i+1]);
		if (strstr(argv[i],"pso_o")!=NULL)
		    _pso_o=atof(argv[i+1]);
		if (strstr(argv[i],"pso_fi_p")!=NULL)
		    _pso_fi_p=atof(argv[i+1]);
		if (strstr(argv[i],"pso_fi_g")!=NULL)
		    _pso_fi_g=atof(argv[i+1]);
		if (strstr(argv[i],"pso_eps")!=NULL)
		    _pso_eps=atof(argv[i+1]);
		if (strstr(argv[i],"to_fit")!=NULL)
		    _to_fit=atoi(argv[i+1]);
		if (strstr(argv[i],"fit_vgm")!=NULL)
		    fit_vgm=atoi(argv[i+1]);
		if (strstr(argv[i],"fit_func")!=NULL)
		    fit_func=atoi(argv[i+1]);
		if (strstr(argv[i],"pso_max_iter")!=NULL)
		    _pso_max_iter=atoi(argv[i+1]);
		if (strstr(argv[i],"filename")!=NULL)
		    filename=argv[i+1];
		if (strcmp(argv[i], "Hbottom")== 0)
		    Hbottom = atof(argv[i + 1]);
		if (strcmp(argv[i], "Hupper")== 0)
		    Hupper = atof(argv[i + 1]);
		if (strstr(argv[i],"fi_file")!=NULL)
		    fi_file=argv[i+1];
		if (strstr(argv[i],"klin")!=NULL)
		    klin=atof(argv[i+1]);
		if (strstr(argv[i],"drip")!=NULL)
		    drip=atof(argv[i+1]);
		if (strstr(argv[i],"newrand")!=NULL)
		    srand(time(NULL));
		if (strcmp(argv[i], "APSO") == 0)
			adaptive_pso = atoi(argv[i + 1]);
		if (strcmp(argv[i], "restart_prox") == 0)
			restart_prox = atof(argv[i + 1]);
		if (strstr(argv[i],"averianov_k")!=NULL)
		    averianov_k=atof(argv[i+1]);
		// linear solver max n and eps
		if (strstr(argv[i],"ls_eps")!=NULL)
		    ls_eps=atof(argv[i+1]);
		if (strstr(argv[i],"ls_max_iter")!=NULL)
		    ls_max_iter=atoi(argv[i+1]);
		// linear solver - tau variation
		if (strstr(argv[i],"ls_percent")!=NULL)
		    ls_percent=atof(argv[i+1]);
		if (strstr(argv[i],"ls_min_tau")!=NULL)
		    ls_min_tau=atof(argv[i+1]);
		if (strstr(argv[i],"ls_mult")!=NULL)
		    ls_mult=atof(argv[i+1]);
		// impulse irrigation
		if (strstr(argv[i],"imp_d0")!=NULL)
		    imp_d0=atof(argv[i+1]);
		if (strstr(argv[i],"imp_d1")!=NULL)
		    imp_d1=atof(argv[i+1]);
		if (strstr(argv[i],"imp_l0")!=NULL)
		    imp_l0=atof(argv[i+1]);
		if (strstr(argv[i],"max_irr_time")!=NULL)
		    max_irr_time=atof(argv[i+1]);
		if (strstr(argv[i],"min_irr_time")!=NULL)
		    min_irr_time=atof(argv[i+1]);
		if (strstr(argv[i],"max_tau_irr")!=NULL)
		    max_tau_irr=atof(argv[i+1]);
		if (strstr(argv[i],"max_tau_0")!=NULL)
		    max_tau=atof(argv[i+1]);
		if (strstr(argv[i],"zone_threshold")!=NULL)
		    zone_threshold=atof(argv[i+1]);
		if (strstr(argv[i],"grad_weight")!=NULL)
		    grad_weight=atof(argv[i+1]);
		if (strstr(argv[i],"specific_storage")!=NULL)
		    _vgm_specific_storage=atof(argv[i+1]);
		if (strstr(argv[i],"irr_ET_sum")!=NULL)
		    irr_ET_sum=atof(argv[i+1]);
		if (strstr(argv[i],"irr_ET_flow")!=NULL)
		    irr_ET_flow=atof(argv[i+1]);
		// plants and sprinklers parameters
		if (strstr(argv[i],"rootR")!=NULL)
		    _rootR=atof(argv[i+1]);
		if (strstr(argv[i],"sprinklerR")!=NULL)
		    _sprinklerR=atof(argv[i+1]);
		// threads
		if (strstr(argv[i],"num_threads")!=NULL)
		    omp_set_num_threads(atoi(argv[i+1]));
		// H restriction
		if (strstr(argv[i],"glminH")!=NULL)
		    glminH=atof(argv[i+1]);
		if (strstr(argv[i],"glmaxH")!=NULL)
		    glmaxH=atof(argv[i+1]);
		// root growth
		if (strstr(argv[i],"no_growth")!=NULL)
		    no_growth=atoi(argv[i+1]);
		if (strstr(argv[i],"rg_r")!=NULL)
		    param_sets[2].r=param_sets[1].r=param_sets[0].r=atof(argv[i+1]);
		if (strstr(argv[i],"rg_la")!=NULL)
		    param_sets[2].la=param_sets[1].la=param_sets[0].la=atof(argv[i+1]);
		if (strstr(argv[i],"rg_lb")!=NULL)
		    param_sets[2].lb=param_sets[1].lb=param_sets[0].lb=atof(argv[i+1]);
		if (strstr(argv[i],"rg_ln")!=NULL)
		    param_sets[2].ln=param_sets[1].ln=param_sets[0].ln=atof(argv[i+1]);
		if (strstr(argv[i],"rg_N")!=NULL)
		    param_sets[2].N=param_sets[1].N=param_sets[0].N=atof(argv[i+1]);
		if (strstr(argv[i],"rg_angle")!=NULL)
		    param_sets[2].angle=param_sets[1].angle=param_sets[0].angle=atof(argv[i+1]);
		if (strstr(argv[i],"rg_max_len")!=NULL)
		{
		    param_sets[0].max_len=atof(argv[i+1]);
		    param_sets[1].max_len=param_sets[0].max_len/2.0;
		    param_sets[2].max_len=param_sets[1].max_len/2.0;
		}
		if (strstr(argv[i],"rg_w1")!=NULL)
		    param_sets[2].w1=param_sets[1].w1=param_sets[0].w1=atof(argv[i+1]);
		// plants
		if (strstr(argv[i],"root_systems_file")!=NULL)
		    root_systems_file=argv[i+1];
		if (strstr(argv[i],"root_systems_out_file")!=NULL)
		    root_systems_out_file=argv[i+1];
		if ((strstr(argv[i],"plants")!=NULL)||(strstr(argv[i],"sprinklers")!=NULL)) // x_y_...
		{
			int n=0;
			int n_=1,xy;
			double *x,*y;
			char *s=argv[i+1],*p=s,*p0;
			while (p[0])
			{
				if (p[0]=='_')
					n_++;
				p++;
			}
			n=n_/2;
			x=new double[n];
			y=new double[n];
			p=p0=s;
			n_=0;
			xy=0;
			while (p[0])
			{
				if (p[0]=='_')
				{
					p[0]=0;
					if (xy==0) x[n_]=atof(p0); else y[n_]=atof(p0); 
					xy=1-xy;
					p0=p+1;
					if (xy==0)
						n_++;
					if (n_==n)
						break;
				}
				p++;
			}
			if (xy==1)
				y[n_]=atof(p0);
			if (strstr(argv[i],"plants")!=NULL)
			{
				printf("plants at ");
				for (int i=0;i<n;i++)
					printf("(%g,%g) ",x[i],y[i]);
				printf("\n");
				_nplants=n;
				_plants_x=x;
			}
			else
			{
				printf("sprinklers at ");
				for (int i=0;i<n;i++)
					printf("(%g,%g) ",x[i],y[i]);
				printf("\n");
				_nsprinklers=n;
				_sprinklers_x=x;
			}
		}
		if (strstr(argv[i],"exclude_depths")!=NULL) // z0_z1_...
		{
			int n_=1;
			char *s=argv[i+1],*p=s,*p0;
			while (p[0])
			{
				if (p[0]=='_')
					n_++;
				p++;
			}
			n_ex_depths=n_;
			exclude_depths=new double[n_];
			p=p0=s;
			n_=0;
			while (p[0])
			{
				if (p[0]=='_')
				{
					p[0]=0;
					exclude_depths[n_]=atof(p0);
					p0=p+1;
					n_++;
					if (n_==n_ex_depths)
						break;
				}
				p++;
			}
			exclude_depths[n_]=atof(p0);
			printf("exclude sensors at depths of ");
			for (int i=0;i<n_ex_depths;i++)
				printf("%g ",exclude_depths[i]);
			printf("\n");
		}
		if (strstr(argv[i],"fit_list")!=NULL) // t0_t1_...
		{
			int n_=1;
			char *s=argv[i+1],*p=s,*p0;
			while (p[0])
			{
				if (p[0]=='_')
					n_++;
				p++;
			}
			n_fit_list=n_;
			fit_list=new double[n_];
			p=p0=s;
			n_=0;
			while (p[0])
			{
				if (p[0]=='_')
				{
					p[0]=0;
					fit_list[n_]=atof(p0);
					p0=p+1;
					n_++;
					if (n_==n_fit_list)
						break;
				}
				p++;
			}
			fit_list[n_]=atof(p0);
			printf("points to fit in ");
			for (int i=0;i<n_fit_list;i++)
				printf("%g ",fit_list[i]);
			printf("\n");
		}
		if (strstr(argv[i],"check_as_fit")!=NULL)
		    check_as_fit=atoi(argv[i+1]);
		if (strstr(argv[i],"mask2d")!=NULL)
		{
		    FILE *fi=fopen(argv[i+1],"rt");
			if (fi)
			{
				char str[1024];
				int _mask_3d_m=0;
				while (fgets(str,1024,fi)) _mask_3d_m++;
				fseek(fi,0,SEEK_SET);
				mask_2d_n=1;
				char *s2=str,*s3;
				while (s2[0])
				{
					if (s2[0]==' ')
						mask_2d_n++;
					s2++;
				}
				mask_2d=new double[mask_2d_n];
				double *mp=mask_2d;
				if (fgets(str,1024,fi))
				{
					s2=s3=str;
					while (s2[0])
					{
						if (s2[0]==' ')
						{
							s2[0]=0;
							mp[0]=atof(s3);
							mp++;
							s3=s2+1;	
						}
						s2++;
					}
					mp[0]=atof(s3);
					mp++;
				}
				printf("2d relief map (%d) read\n",mask_2d_n);
			}
		}
	}
	printf("2D - (%d,%d) tend %g tsave %g \n",N,M,Tm,Sm);
	fflush(stdout);		
	{
		FILE *fi1;
		if (fi_file==NULL)
		    fi1=fopen("fit_input.txt", "rt");
		else
		    fi1 = fopen(fi_file, "rt");
		double *values_T,*values_Z,*values_F;
		double *check_T,*check_Z,*check_F;
		double *init_Z,*init_F;		
		double *perc_T, *perc_A;
		double *irr_T, *irr_A;
		double param_values[200];
		int nparam_values=0;
		int nvalues;
		int ninit;
		int ncheck;
		int nperc, nirr;
		double *EV_T,**EV_F,*LAI;
		double ev_mu;
		int nEV_T,nEV_F;
		int to_fit;
		int pso_n;
		double pso_o, pso_fi_p, pso_fi_g, pso_eps;
		int pso_max_iter;
		char str[2048];
		// read data
		if (fscanf(fi1,"%d %d %d %d %d %d %d %d %lg\n",&to_fit,&nvalues,&ninit,&ncheck,&nperc,&nirr,&nEV_T,&nEV_F,&ev_mu)!=9) return 0;
		if (fscanf(fi1,"%d %lg %lg %lg %lg %d\n",&pso_n,&pso_o,&pso_fi_p,&pso_fi_g,&pso_eps,&pso_max_iter)!=6) return 0;
		if (_ev_mu!=1e100) ev_mu=_ev_mu;
		if (_pso_n!=-1) pso_n=_pso_n;
		if (_pso_o!=1e100) pso_o=_pso_o;
		if (_pso_fi_p!=1e100) pso_fi_p=_pso_fi_p;
		if (_pso_fi_g!=1e100) pso_fi_g=_pso_fi_g;
		if (_pso_eps!=1e100) pso_eps=_pso_eps;
		if (_pso_max_iter!=-1) pso_max_iter=_pso_max_iter;
		if (_to_fit!=-1) to_fit=_to_fit;
		printf("%d %d %d %d %d %d %d %d %lg\n",to_fit,nvalues,ninit,ncheck,nperc,nirr,nEV_T,nEV_F,ev_mu);
		printf("%d %lg %lg %lg %lg %d\n",pso_n,pso_o,pso_fi_p,pso_fi_g,pso_eps,pso_max_iter);
		int nf;
		double LLy;
		if (fscanf(fi1,"%d %lg %lg\n",&nf,&LLx,&LLy)!=3) return 0;
		// values to fit -  in 10m
		values_T=new double[nvalues];
		values_Z = new double[3*nvalues]; //  Z for 1d, (x,y,z) for 3d
		values_F=new double[nvalues];
		int j = 0;
		printf("to fit: %d\n",nvalues);
		for (int i = 0;i < nvalues;i++)
			if (fgets(str, 2048, fi1))
			if (sscanf(str, "%lg %lg %lg %lg %lg\n", values_T + j, values_Z + 3 * j + 0, values_Z + 3 * j + 1, values_Z + 3 * j + 2, values_F + j) == 5)
			    if ((values_T[j]-=start_time)>=0)
			    {
				int incl=1;
				if (n_fit_list) incl=0;
				for (int d=0;d<n_fit_list;d++)
					if (values_T[j]==fit_list[d])
						incl=1;
				for (int d=0;d<n_ex_depths;d++)
					if (values_Z[3*j+0]==exclude_depths[d])
						incl=0;
				if (incl) j++;
			    }
		nvalues = j;
		printf("read to fit: %d last:(%g %g %g %g %g)\n",nvalues,values_T[nvalues-1],values_Z[3*(nvalues-1)+0],values_Z[3*(nvalues-1)+1],values_Z[3*(nvalues-1)+2],values_F[nvalues-1]);
		// initial H
		init_Z = new double[3*ninit]; //  Z for 1d, (x,y,z) for 3d
		init_F = new double[ninit];
		j = 0;
		printf("ninit: %d\n",ninit);
		for (int i = 0;i < ninit;i++)
			if (fgets(str, 2048, fi1))
				if (sscanf(str, "%lg %lg %lg %lg\n", init_Z + 3 * j + 0, init_Z + 3 * j + 1, init_Z + 3 * j + 2, init_F + j) == 4)
				    j++;
		ninit = j;
		printf("read init: %d last:(%g %g %g %g)\n",ninit,init_Z[3*(ninit-1)+0],init_Z[3*(ninit-1)+1],init_Z[3*(ninit-1)+2],init_F[ninit-1]);
		// percipitation - fixed on upper plane for 3d in m/s
		perc_T = new double[nperc];
		perc_A = new double[nperc];
		j=0;
		for (int i = 0;i<nperc;i++)
			if (fscanf(fi1, "%lg %lg\n", perc_T + j, perc_A + j) == 2)
			    if ((perc_T[j]-=start_time)>=0)
				j++;
		nperc=j;
		printf("read prcp: %d last:(%g %g)\n",nperc,perc_T[nperc-1],perc_A[nperc-1]);
		// irrigation - fixed on upper plane for 3d in m/s
		irr_T = new double[nirr];
		irr_A = new double[nirr];
		j=0;
		for (int i = 0;i < nirr;i++)
			if (fscanf(fi1, "%lg %lg\n", irr_T + j, irr_A + j) == 2)
			    if ((irr_T[j]-=start_time)>=0)
				j++;
		nirr=j;
		printf("read irr: %d last:(%g %g)\n",nirr,irr_T[nirr-1],irr_A[nirr-1]);
		// values to check in 10m
		check_T=new double[ncheck];
		check_Z = new double[3*ncheck]; //  Z for 1d, (x,y,z) for 3d
		check_F=new double[ncheck];
		j = 0;
		int init_j=0;
		printf("to check: %d\n",ncheck);
		for (int i = 0;i < ncheck;i++)
			if (fgets(str, 2048, fi1))
				if (sscanf(str, "%lg %lg %lg %lg %lg\n", check_T + j, check_Z + 3 * j + 0, check_Z + 3 * j + 1, check_Z + 3 * j + 2, check_F + j) == 5)
					if ((check_T[j]-=start_time)>=0)
					{
					    if ((fabs(check_T[j])<=1)&&(start_time!=0))
					    {	
						if (init_j>=ninit)
						{
						    double *n_init_Z = new double[3*(init_j+1)];
						    double *n_init_F = new double[init_j+1];
						    memcpy(n_init_Z,init_Z,3*ninit*sizeof(double));
						    memcpy(n_init_F,init_F,ninit*sizeof(double));
						    ninit=init_j+1;
						    delete [] init_Z;
						    delete [] init_F;
						    init_Z=n_init_Z;
						    init_F=n_init_F;
						}
						for (int kk=0;kk<3;kk++) init_Z[3*init_j+kk]=check_Z[3*j+kk];
						init_F[init_j]=check_F[j];
						init_j++;
					    }
					    int incl=1;
					    for (int d=0;d<n_ex_depths;d++)
						if (check_Z[3*j+0]==exclude_depths[d])
							incl=0;
					    if (incl) j++;
					}
		ncheck = j;
		if (start_time!=0)
		{
		     ninit=init_j;
		     printf("reread init: %d last:(%g %g %g %g)\n",ninit,init_Z[3*(ninit-1)+0],init_Z[3*(ninit-1)+1],init_Z[3*(ninit-1)+2],init_F[ninit-1]);
		}
		printf("read to check: %d last:(%g %g %g %g %g)\n",ncheck,check_T[ncheck-1],check_Z[3*(ncheck-1)+0],check_Z[3*(ncheck-1)+1],check_Z[3*(ncheck-1)+2],check_F[ncheck-1]);
		// evapotranspiration in m/s
		EV_T = new double[3*nEV_T]; // T for 1d, (T,x,y) for 3d
		EV_F=new double*[nEV_F];
		LAI=new double[nEV_T];
		for (int j = 0;j < nEV_F;j++)
			EV_F[j] = new double[nEV_T];
		int kk = 0;
		FILE *sfi1=fi1,*ff;
		if (et_file)
		    if (ff=fopen(et_file,"rt"))
			fi1=ff;
		for (int i=0;i<nEV_T;i++)
		{
			if (fgets(str, 2048, fi1))
			{
				    int idx=0;
				    char *str1=str;
				    double vs[4];
				    int jj;
				    for (jj=0;jj<4;jj++)
				    {
					idx=0;
					while ((str1[idx]!=' ')&&(str1[idx]!='\n')&&(str1[idx])) idx++;
					if (str1[idx]==' ')
					{
					    str1[idx]=0;
					    vs[jj]=atof(str1);
					    str1=str1+idx+1;
					}
					else 
					    break;
				    }
				    if (jj==4)
				    {
					EV_T[3*kk+0]=vs[0];
					EV_T[3*kk+1]=vs[1];
					EV_T[3*kk+2]=vs[2];
					LAI[kk]=vs[3];
					for (jj=0;jj<nEV_F;jj++)
					{
					    idx=0;
					    while ((str1[idx]!=' ')&&(str1[idx]!='\n')&&(str1[idx])) idx++;
					    if ((str1[idx]==' ')||(jj==(nEV_F-1)))
					    {
						str1[idx]=0;
						EV_F[jj][kk]=atof(str1);
						str1=str1+idx+1;
					    }
					    else 
						break;
					}
					if (jj==nEV_F)
					if ((EV_T[3*kk+0]-=start_time)>=0)
					    kk++;
				    }
			}
			fscanf(fi1, "\n");
		}
		printf("EVn: %d",nEV_T);
		nEV_T = kk;
		printf("<-%d\n", nEV_T);
		printf("ev00 - %g\n",EV_F[0][0]);
		if (et_file)
		{
		    fclose(ff);
		    fi1=sfi1;
		}
		// read fixed parameter values
		if (pv_file)
		    if (ff=fopen(pv_file,"rt"))
			fi1=ff;
		while (fscanf(fi1,"%lg ",&param_values[nparam_values++])==1)
		    if (nparam_values==200)
		        break;
		if (pv_file)
		    fclose(ff);
		if (nparam_values!=200)
		    nparam_values--;
		// run
		fit_and_solve(Tm,Sm,Om,irr,bc,tau_m,init_Z,init_F,ninit,perc_T,perc_A,nperc,irr_T,irr_A,nirr,values_T,values_Z,values_F,nvalues,check_T,check_Z,check_F,ncheck,EV_T,EV_F,LAI,ev_mu,nEV_T,nEV_F,to_fit,pso_n,pso_o,pso_fi_p,pso_fi_g,pso_eps,pso_max_iter,fit,&param_values[0],nparam_values,et_file,pv_file);
	}
	return 0;
}
