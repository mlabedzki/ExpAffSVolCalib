#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <functional>
#include <complex>
#include <Rcpp.h>
#define M_PI 3.14159265358979323846
#define LEN(a) (sizeof(a)/sizeof(a[0]))
#define INIT_FROM_ARRAY(ar) (ar, ar + sizeof(ar) / sizeof(ar[0]))
#include <random>

using namespace Rcpp;

typedef std::complex<double> dcomp;
const dcomp i1(0.0, 1.0);

double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

namespace BT{
  template<class D, class OP, class ARG>
  std::vector<D> Simplex(OP f,                   //target function
			 std::vector<D> init,    //initial guess of the parameters
			 ARG args, //non-optimized arguments for target function
			 D tol=1E8*std::numeric_limits<D>::epsilon(), //termination criteria
			 int iterations=500,
			 std::string logfile = "",
			 std::vector<std::vector<D> > x =  std::vector<std::vector<D> >()
			 //x: The Simplex
			 ){    //iteration step number, e.g. 1E5

    int N=init.size();                         //space dimension
    const double a=1.0, b=1.0, g=0.5, h=0.5;   //coefficients
                                               //a: reflection  -> xr  
                                               //b: expansion   -> xe 
                                               //g: contraction -> xc
                                               //h: full contraction to x1
    std::vector<D> xcentroid_old(N,0);   //simplex center * (N+1)
    std::vector<D> xcentroid_new(N,0);   //simplex center * (N+1)
    std::vector<D> vf(N+1,0);            //f evaluated at simplex vertexes       
    int x1=0, xn=0, xnp1=0;         //x1:   f(x1) = min { f(x1), f(x2)...f(x_{n+1} }
                                    //xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} }
                                    //xn:   f(xn)<f(xnp1) && f(xn)> all other f(x_i)
    int cnt=0; //iteration step number

    if(x.size()==0) //if no initial simplex is specified
      { //construct the trial simplex
	//based upon the initial guess parameters
      std::vector<D> del( init );
      std::transform(del.begin(), del.end(), del.begin(), 
		     std::bind2nd( std::divides<D>() , 20) );//'20' is picked 
                                                             //assuming initial trail close to true
      
      for(int i=0; i<N; ++i){
	std::vector<D> tmp( init );
	tmp[i] +=  del[i];
	x.push_back( tmp );
      }
      x.push_back(init);//x.size()=N+1, x[i].size()=N
      
      //xcentriod
      std::transform(init.begin(), init.end(), 
		xcentroid_old.begin(), std::bind2nd(std::multiplies<D>(), N+1) );
      }//constructing the simplex finished
    
    //optimization begins
    for(cnt=0; cnt<iterations; ++cnt){

	for(int i=0;i<N+1;++i){
	vf[i]= f(x[i],args);
      }
      
      x1=0; xn=0; xnp1=0;//find index of max, second max, min of vf.
      
      for(int i=0;i<vf.size();++i){
	if(vf[i]<vf[x1]){
	  x1=i;
	}
	if(vf[i]>vf[xnp1]){
	  xnp1=i;
	}
      }
      
      xn=x1;
      
      for(int i=0; i<vf.size();++i){ 
	if( vf[i]<vf[xnp1] && vf[i]>vf[xn] )
	  xn=i;
      }
      //x1, xn, xnp1 are found

      std::vector<D> xg(N, 0);//xg: centroid of the N best vertexes
      for(int i=0; i<x.size(); ++i){
	if(i!=xnp1)
	  std::transform(xg.begin(), xg.end(), x[i].begin(), xg.begin(), std::plus<D>() );
      }
      std::transform(xg.begin(), xg.end(), 
		x[xnp1].begin(), xcentroid_new.begin(), std::plus<D>());
      std::transform(xg.begin(), xg.end(), xg.begin(), 
		std::bind2nd(std::divides<D>(), N) );
      //xg found, xcentroid_new updated

      //termination condition
      D diff=0;          //calculate the difference of the simplex centers
                         //see if the difference is less than the termination criteria
      for(int i=0; i<N; ++i)     
	diff += fabs(xcentroid_old[i]-xcentroid_new[i]);

      if (diff/N < tol) break;              //terminate the optimizer
      else xcentroid_old.swap(xcentroid_new); //update simplex center
      
      //reflection:
      std::vector<D> xr(N,0); 
      for( int i=0; i<N; ++i)
	xr[i]=xg[i]+a*(xg[i]-x[xnp1][i]);
      //reflection, xr found
      
      D fxr=f(xr,args);//record function at xr
      
      if(vf[x1]<=fxr && fxr<=vf[xn])
	std::copy(xr.begin(), xr.end(), x[xnp1].begin() );
      
      //expansion:
      else if(fxr<vf[x1]){
	std::vector<D> xe(N,0);
	for( int i=0; i<N; ++i)
	  xe[i]=xr[i]+b*(xr[i]-xg[i]);
	if( f(xe,args) < fxr )
	  std::copy(xe.begin(), xe.end(), x[xnp1].begin() );
	else
	  std::copy(xr.begin(), xr.end(), x[xnp1].begin() );
      }//expansion finished,  xe is not used outside the scope
      
      //contraction:
      else if( fxr > vf[xn] ){
	std::vector<D> xc(N,0);
	for( int i=0; i<N; ++i)
	  xc[i]=xg[i]+g*(x[xnp1][i]-xg[i]);
	if( f(xc,args) < vf[xnp1] )
	  std::copy(xc.begin(), xc.end(), x[xnp1].begin() );
	else{
	  
	  for( int i=0; i<x.size(); ++i ){
	    if( i!=x1 ){ 
	      for(int j=0; j<N; ++j) 
		x[i][j] = x[x1][j] + h * ( x[i][j]-x[x1][j] );
	    }
	  }
	  
	}
      }//contraction finished, xc is not used outside the scope

		//write error to file
		if((cnt+1) % 10 == 0 && logfile!="") {
			//ofstream myfile;
			std::fstream myfile;
			myfile.open (logfile.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
			myfile << sqrt(f(x[x1],args)/30) << ";";
			myfile.close();
		}
		//write error to file
	  
    }//optimization is finished

    if(cnt==iterations){//max number of iteration achieves before tol is satisfied
      std::cout<<"Iteration limit achieves, result may not be optimal"<<std::endl;
    }
	
	if(logfile!="") {
		std::fstream myfile;
		myfile.open (logfile.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		myfile << std::endl;
		myfile.close();
	}
	
    return x[x1];
  }
    
}//BT::

template <typename T>
std::vector<T> rep(std::vector<T> v, int n)
{
	int m = v.size();
    std::vector<T> temp(n);
    for (int i = 0; i < n; i++)
	    for (int j = 0; j < m; j++)
			temp[i*m+j] = v[j];
    return temp;
}

template <typename T>
std::vector<T> subv(std::vector<T> v, int start, int end) 
{
	return std::vector<T>(v.begin() + start,v.begin() + end+1);
}

template <typename T>
std::vector<T> vectorize(T x) 
{
	std::vector<T> v(1);
	v[0] = x;
	return v;
}

template <typename T>
std::vector<T> linspace(T start, T end, int n) 
{
    std::vector<T> temp(n);
    for (int i = 0; i < n; i++)
    {
    	temp.at(i) = start + (end-start)*((double) i/(n-1));
    }
    return temp;	
}

template <typename T>
std::vector<T> linspacing(T start, T end, double spacing) 
{
	int n = (int)((end-start)/spacing) + 2;
    std::vector<T> temp(n);
	temp[0] = start;
    for (int i = 1; i < n; i++)
    {
    	temp[i] = temp[i-1] + spacing;
    }
    return temp;	
}

double trapz(std::vector<double> y, std::vector<double> x) 
{
	int n = x.size();
    double temp=0;
    for (int i = 1; i < n; i++)
    {
    	temp += ( (x[i] - x[i-1]) * (y[i] + y[i-1])) / 2;
    }
    return temp;	
}

std::vector<double> realv(const std::vector<dcomp>& v)
{
            int n = v.size();
            std::vector<double> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = real(v.at(i));
            }
            return temp;
}

std::vector<double> imagv(const std::vector<dcomp>& v)
{
            int n = v.size();
            std::vector<double> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = imag(v.at(i));
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> powv(const std::vector<T>& v, T2 scalar)
{
    int n = v.size();
    std::vector<T> temp(n);
    for (int i = 0; i < n; i++)
    {
		temp.at(i) = pow(v.at(i),scalar);
    }
    return temp;
}

template <typename T>
std::vector<T> sinv(const std::vector<T>& v)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = sin(v.at(i));
            }
            return temp;
}

template <typename T>
std::vector<T> cosv(const std::vector<T>& v)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = cos(v.at(i));
            }
            return temp;
}

template <typename T>
std::vector<T> expv(const std::vector<T>& v)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = exp(v.at(i));
            }
            return temp;
}

template <typename T>
std::vector<T> logv(const std::vector<T>& v)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = log(v.at(i));
            }
            return temp;
}

template <typename T>
std::vector<T> operator-(const std::vector<T>& v)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = -v.at(i);
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator*(const std::vector<T>& v1, const std::vector<T2>& v2)
{
            int n = v1.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v1.at(i)*v2.at(i);
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator*(const std::vector<T>& v, T2 scalar)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v.at(i)*scalar;
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator*(T2 scalar, const std::vector<T>& v)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v.at(i)*scalar;
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator/(const std::vector<T>& v, T2 scalar)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v.at(i)/scalar;
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator/(T2 scalar,const std::vector<T>& v)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = scalar/v.at(i);
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator/(const std::vector<T>& v1, const std::vector<T2>& v2)
{
            int n = v1.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v1.at(i)/v2.at(i);
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T2>& v2)
{
            int n = v1.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v1.at(i) + v2.at(i);
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator+(const std::vector<T>& v, T2 scalar)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v.at(i)+scalar;
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator+(T2 scalar, const std::vector<T>& v)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v.at(i)+scalar;
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T2>& v2)
{
            int n = v1.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v1.at(i) - v2.at(i);
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator-(const std::vector<T>& v, T2 scalar)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = v.at(i)-scalar;
            }
            return temp;
}

template <typename T, typename T2>
std::vector<T> operator-(T2 scalar, const std::vector<T>& v)
{
            int n = v.size();
            std::vector<T> temp(n);
            for (int i = 0; i < n; i++)
            {
                temp.at(i) = scalar-v.at(i);
            }
            return temp;
}

template< typename T>
std::vector<std::vector<T> > makeMatrix(std::vector<T> data, int n, int m)
{
    //vector<vector<T> > mat;
    std::vector<std::vector<T> > mat(n, std::vector<T>(m));
    for (int i=0; i<n; i++) {
    	for (int j=0; j<m; j++) {
    		mat[i][j] = data[i*m+j];
    		//vector<T> tmp;
    		//tmp.push_back(data[i*m+j]);
    	}
    	//mat.push_back(tmp);
	}
	return mat;
}

template< typename T>
std::vector<T> mat2vec(std::vector<std::vector<T> > mat)
{
    //vector<vector<T> > mat;
	int n = mat.size();
	int m = mat[0].size();
    std::vector<T> vec(n*m);
    for (int i=0; i<n; i++) {
    	for (int j=0; j<m; j++) {
    		vec[i*m+j] = mat[i][j];
    		//vector<T> tmp;
    		//tmp.push_back(data[i*m+j]);
    	}
    	//mat.push_back(tmp);
	}
	return vec;
}

double hestonTSerr(
	std::vector<double> params,
	std::vector<std::vector<double> > args)
{
	std::vector<double> tau = args[0];
	std::vector<double> marketTS = args[1];
	double nu0 = exp(params[0]);
	//double nu0 = args[2][0];
	double theta = exp(params[1]);
	double kappa = exp(params[2]);
	//double kappa = args[3][0];
	double err = 0;
	int n = tau.size();
	for (int i = 0; i < n; i++)
	{
		err += fabs(theta + (nu0-theta) * (1-exp(-kappa*tau[i]))/(kappa*tau[i]) - marketTS[i]);
		//err += pow(theta + (nu0-theta) * (1-exp(-kappa*tau[i]))/(kappa*tau[i]) - marketTS[i],2);
	}
	return err;
}

double szTSerr(
	std::vector<double> params,
	std::vector<std::vector<double> > args)
{
	std::vector<double> tau = args[0];
	std::vector<double> marketTS = args[1];
	double omega = args[2][0];
	double nu0 = exp(params[0]);
	double theta = exp(params[1]);
	double kappa = exp(params[2]);
	double x = pow(omega,2)/2/kappa;
	double err = 0;
	int n = tau.size();
	for (int i = 0; i < n; i++)
	{
		err += fabs((exp(-2*kappa*tau[i])*(exp(2*kappa*tau[i])*(pow(theta,2)*(2*kappa*tau[i]-3)+2*theta*nu0+x*(2*kappa*tau[i]-1)+pow(nu0,2))+4*theta*(theta-nu0)*exp(kappa*tau[i])-pow(theta-nu0,2)+x))/(2*kappa*tau[i]) - marketTS[i]);
		//err += pow(theta + (nu0-theta) * (1-exp(-kappa*tau[i]))/(kappa*tau[i]) - marketTS[i],2);
	}
	return err;
}

// [[Rcpp::export]]
std::vector<double> nlVarTS(
	std::vector<double> params,
	std::vector<double> marketTS,
	std::vector<double> tau,
	int modType=1,
	bool err = false,
	int max_it = 8000,
	double omega = 0)
{
	std::vector<std::vector<double> > args;
	args.push_back(tau);
	args.push_back(marketTS);
	//args.push_back(vectorize(params[0]));
	//args.push_back(vectorize(params[2]));
	args.push_back(vectorize(omega));
	std::vector<double> init = logv(params);
	//std::vector<double> init;
	//init.push_back(log(params[1]));
	std::vector<double> out;
	//cout << "The calib Heston params are: " << endl;
	//printcon(BT::Simplex(SSE5, init, args, 1e-11));
	//return SSE5(init,args);
	switch (modType) {
		case 1: {
			out = expv(BT::Simplex(hestonTSerr, init, args, 1e-11, max_it));
			if (err) out.push_back(hestonTSerr(logv(out),args));
			break; }
		case 2: {
			out = expv(BT::Simplex(szTSerr, init, args, 1e-11, max_it));
			if (err) out.push_back(szTSerr(logv(out),args));
			break; }
	}
	//std::vector<double> out3;
	//out3.push_back(params[0]);
	//out3.push_back(out[0]);
	//out3.push_back(params[2]);
	//if (err) out3.push_back(hestonTSerr(logv(out),args));
	return out;
}

// [[Rcpp::export]]
std::vector<double> nlHestonTSv(
	std::vector<double> params1,
	std::vector<double> marketTS1,
	std::vector<double> tau,
	std::vector<double> omega,
	int modType=1,	
	bool err = false,
	int max_it = 8000)
{
	int m = tau.size();
	int nm = marketTS1.size();
	std::vector<std::vector<double> > params = makeMatrix(params1,nm/m,m);
	std::vector<std::vector<double> > marketTS = makeMatrix(marketTS1,nm/m,m);	
	std::vector<std::vector<double> > out;
	for (int i = 0; i < nm/m; i++)
	{
		out.push_back(nlVarTS(params[i],marketTS[i],tau,modType,err,max_it,omega[i]));
	}
	return mat2vec(out);
}


dcomp CFa(dcomp u, std::vector<double> params, double lambda, double tau, bool type = 1)
{
	//#Albrecher et al. (2007)
	double theta = params[0], sigma = params[1], rho = params[2], kappa = params[3], nu0 = params[4];
	double b;
	double zeta;
	dcomp c, d, C, D, h;
	if (type == 1)
		zeta = 1.0;
	if (type == 0)
		zeta = -1.0;
	b = kappa + lambda - (1.0 + zeta) / 2.0 * rho*sigma;
	h = b - u*rho*sigma*i1;
	d = pow(pow(-h,2) - (u*zeta*i1 - pow(u,2))*pow(sigma,2),0.5);
	c = (h - d) / (h + d); //#    == 1 / g;
	if (std::isinf(real(c)))
		C = kappa*theta / pow(sigma,2) * (h - d)*tau;
	else
		C = kappa*theta / pow(sigma,2) * ((h - d)*tau - 2.0 * log((1.0 - c*exp(-d*tau)) / (1.0 - c)));
	D = (h - d) / pow(sigma,2) * ((1.0 - exp(-d*tau)) / (1.0 - c*exp(-d*tau)));
	return exp(C + D*nu0);
}

std::vector<dcomp> CFaV(
	std::vector<dcomp> u, std::vector<double> params, double lambda, double tau, bool type = 1)
{
	//#Albrecher et al. (2007)
	double theta = params[0], sigma = params[1], rho = params[2], kappa = params[3], nu0 = params[4];
	double b;
	double zeta;
	std::vector<dcomp> c, d, C, D, h;
	if (type == 1)
		zeta = 1.0;
	if (type == 0)
		zeta = -1.0;
	b = kappa + lambda - (1.0 + zeta) / 2.0 * rho*sigma;
	h = b - u*rho*sigma*i1;
	d = powv(powv(-h,2) - (u*zeta*i1 - powv(u,2))*pow(sigma,2),0.5);
	c = (h - d) / (h + d); //#    == 1 / g;
	//if (std::isinf(real(c)))
	//	C = (rd - rf)*u*i1*tau + kappa*theta / pow(sigma,2) * (b - rho*sigma*u*i1 - d)*tau;
	//else
	C = kappa*theta / pow(sigma,2) * ((h - d)*tau - 2.0 * logv((1.0 - c*expv(-d*tau)) / (1.0 - c)));
	D = (h - d) / pow(sigma,2) * ((1.0 - expv(-d*tau)) / (1.0 - c*expv(-d*tau)));
	return expv(C + D*nu0);
}

double hestonPre(double F, double K, double rd, double tau,
	std::vector<dcomp> u,
	std::vector<dcomp> I1,
	std::vector<dcomp> I2,
	bool type = 1)
{
	//# type = 1 = > Call
	//# type = -1 = > Put
	double l = log(K/F);
	double P1 = 0.5 + 1 / M_PI * trapz(realv(expv(-i1*u*l)*I1), realv(u));
	double P2 = 0.5 + 1 / M_PI * trapz(realv(expv(-i1*u*l)*I2), realv(u));
	double p1 = (1 - type) + (2 * type - 1)*P1;
	double p2 = (1 - type) + (2 * type - 1)*P2;
	return (2*type-1) * exp(-rd*tau)*(F*p1 - K*p2);
}

double carrPre(double F, double K, double rd, double tau,
	std::vector<dcomp> u,
	std::vector<dcomp> I1,
	double alpha = 5)
{
	double k = log(K);
	return exp(-alpha*k)/M_PI * trapz(realv(expv(-i1*u*(k-log(F)))*exp(-rd*tau)*I1), realv(u));
	//lower limit of integral is problematic, originally was 1e-16 but does not work in some cases, 1e-2 works in all cases
	//upper limit was not checked
}

dcomp CFaSZ(dcomp u, std::vector<double> params, double lambda, double tau, bool type = 1)
{
	double theta = params[0], sigma = params[1], rho = params[2], kappa = params[3], nu0 = params[4];
	dcomp alpha = -0.5*u*(u+i1);
	dcomp beta = 2.0*(kappa-i1*sigma*rho*u);
	double gamma = 2.0*pow(sigma,2);
	dcomp d = sqrt(pow(beta,2)-4.0*alpha*gamma);
	dcomp g = (beta-d)/(beta+d);
	dcomp Bs = 2.0*kappa*theta*(beta-d)*pow(1.0-exp(-d/2.0*tau),2)/d/gamma/(1.0-g*exp(-d*tau));
	dcomp Bv = (beta-d)*(1.0-exp(-d*tau))/2.0/gamma/(1.0-g*exp(-d*tau));
	dcomp A;
	if (std::isinf(real(g)))
		A = 0;
	else
		A = (beta-d)*pow(kappa,2)*pow(theta,2)/pow(d,3)/gamma*(beta*(d*tau-4.0)+d*(d*tau-2.0)+(4.0*exp(-d/2.0*tau)*((pow(d,2)-2.0*pow(beta,2))/(beta+d)*exp(-d/2.0*tau)+2.0*beta))/(1.0-g*exp(-d*tau))) + 0.25*(beta-d)*tau-0.5*log((g*exp(-d*tau)-1.0)/(g-1.0));
		//#A = (beta-d)*kappa^2*theta^2/2/d^3/sigma^2*((d^3 *tau *(e^(d *tau)+1)+2 *d^2 *(beta*tau*e^(d*tau)+2)+8*beta^2*(e^((d*tau)/2)-1)+d*beta*(beta*tau*e^(d*tau)+8*e^((d*tau)/2)-beta*tau))/((beta*(e^(d*tau)-1)+d*(e^(d*tau)+1)))) + 1/4*(beta-d)*tau-1/2*log((g*exp(-d*tau)-1)/(g-1))	#alternative formula
	return exp(A + Bs*nu0 + Bv*pow(nu0,2));
}

dcomp CFaSZ2d(dcomp u, std::vector<double> params, double lambda, double tau, bool type = 1)
{
	return CFaSZ(u,subv(params,0,4),lambda,tau,type)*CFaSZ(u,subv(params,5,9),lambda,tau,type);
}

dcomp CFa2d(dcomp u, std::vector<double> params, double lambda, double tau, bool type = 1)
{
	return CFa(u,subv(params,0,4),lambda,tau,type)*CFa(u,subv(params,5,9),lambda,tau,type);
}

// [[Rcpp::export]]
double sz2dCarr(std::vector<double> paramsPos, double F, double K, double rd, double tau, double alpha=5)
{
	double lambda = 0;
	const int m = 100;
	std::vector<dcomp> I1(m), u = linspace(dcomp(1e-16,0.0), dcomp(100.0,0.0), m);
	for (int k = 0; k < m; k++) {
		I1[k] = CFaSZ2d(u[k]-(alpha+1)*i1, paramsPos, lambda, tau, 0)/(pow(alpha,2) + alpha - pow(u[k],2) + i1*(2*alpha + 1)*u[k]);
	}
	return carrPre(F, K, rd, tau, u, I1, alpha);	
}

// [[Rcpp::export]]
double sz2dHeston(std::vector<double> paramsPos, double F, double K, double rd, double tau)
{
	double lambda = 0;
	const int m = 400;
	std::vector<dcomp> I1(m), I2(m), u = linspace(dcomp(1e-16,0.0), dcomp(100.0,0.0), m);
	dcomp stale = CFaSZ2d(-i1, paramsPos, lambda, tau, 0);
	for (int k = 0; k < m; k++) {
		I1[k] = CFaSZ2d(u[k]-i1, paramsPos, lambda, tau, 0)/(i1 * u[k])/stale;
		I2[k] = CFaSZ2d(u[k], paramsPos, lambda, tau, 0) / (i1 * u[k]);
	}
	return hestonPre(F, K, rd, tau, u, I1, I2);	
}

double attariPre(double F, double K, double rd, double tau,
	std::vector<double> u,
	std::vector<dcomp> I1,
	bool type = 1)
{
	//# type = 1 = > Call
	//# type = -1 = > Put
	int m = u.size();
	std::vector<double> fx(m);
	for (int k = 0; k < m; k++)
		fx[k] = real(exp(-i1*exp(u[k])*log(K/F))*I1[k]);
	double integrand = trapz(fx, u);
	return exp(-rd*tau)*(F - K*(0.5 + 1/M_PI * integrand));
}

double attariPre2(double F, double K, double rd, double tau,
	std::vector<double> u,
	std::vector<dcomp> I1,
	bool type = 1)
{
	//# type = 1 = > Call
	//# type = -1 = > Put
	int m = u.size();
	std::vector<double> fx(m);
	double l = log(K/F);
	for (int k = 0; k < m; k++)
		fx[k] = exp(u[k])*((real(I1[k])+imag(I1[k])/exp(u[k]))*cos(exp(u[k])*l)+(imag(I1[k])-real(I1[k])/exp(u[k]))*sin(exp(u[k])*l))/(1.0+pow(exp(u[k]),2));
	double integrand = trapz(fx, u);
	return exp(-rd*tau)*(F - K*(0.5 + 1/M_PI * integrand));
}

// [[Rcpp::export]]
double hestonAttari(std::vector<double> paramsPos, double F, double K, double rd, double tau)
{
	double lambda = 0;
	std::vector<double> u = linspacing(-17.0, 4.9, .4);
	int m = u.size();
	std::vector<dcomp> I1(m);
	for (int k = 0; k < m; k++) {
		I1[k] = (1.0-i1/exp(u[k]))/(1.0+pow(exp(u[k]),2))*exp(u[k])*CFa(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau, 0);
	}
	return attariPre(F, K, rd, tau, u, I1);	
}

// [[Rcpp::export]]
double hestonAttari2(std::vector<double> paramsPos, double F, double K, double rd, double tau)
{
	double lambda = 0;
	std::vector<double> u = linspacing(-17.0, 4.9, .4);
	int m = u.size();
	std::vector<dcomp> I1(m);
	for (int k = 0; k < m; k++) {
		I1[k] = CFa(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau, 0);
	}
	return attariPre2(F, K, rd, tau, u, I1);	
}

// [[Rcpp::export]]
double batesAttari(std::vector<double> paramsPos, double F, double K, double rd, double tau)
{
	double lambda = 0;
	std::vector<double> u = linspacing(-17.0, 4.9, .4);
	int m = u.size();
	std::vector<dcomp> I1(m);
	for (int k = 0; k < m; k++) {
		I1[k] = (1.0-i1/exp(u[k]))/(1.0+pow(exp(u[k]),2))*exp(u[k])*CFa2d(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau, 0);
	}
	return attariPre(F, K, rd, tau, u, I1);	
}

// [[Rcpp::export]]
double ououAttari(std::vector<double> paramsPos, double F, double K, double rd, double tau)
{
	double lambda = 0;
	std::vector<double> u = linspacing(-17.0, 4.9, .4);
	int m = u.size();
	std::vector<dcomp> I1(m);
	for (int k = 0; k < m; k++) {
		I1[k] = (1.0-i1/exp(u[k]))/(1.0+pow(exp(u[k]),2))*exp(u[k])*CFaSZ2d(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau, 0);
	}
	return attariPre(F, K, rd, tau, u, I1);	
}

// [[Rcpp::export]]
double szAttari(std::vector<double> paramsPos, double F, double K, double rd, double tau)
{
	double lambda = 0;
	std::vector<double> u = linspacing(-17.0, 4.9, .4);
	int m = u.size();
	std::vector<dcomp> I1(m);
	for (int k = 0; k < m; k++) {
		I1[k] = (1.0-i1/exp(u[k]))/(1.0+pow(exp(u[k]),2))*exp(u[k])*CFaSZ(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau, 0);
	}
	return attariPre(F, K, rd, tau, u, I1);	
}

bool feller(std::vector<double> params)
{
	bool x;
	if (params.size()<6) {
		if (2*params[3]*params[0]>pow(params[1],2)) x=true;
		else x=false;
	}
	else {
		if (2*params[3]*params[0]>pow(params[1],2) && 2*params[8]*params[5]>pow(params[6],2)) x=true;
		else x=false;
	}
	return x;
}

bool feller2(std::vector<double> params)
{
	bool x;
	if (params.size()<6) {
		if (params[3]*params[0]>pow(params[1],2)) x=true;
		else x=false;
	}
	else {
		if (params[3]*params[0]>pow(params[1],2) && params[8]*params[5]>pow(params[6],2)) x=true;
		else x=false;
	}
	return x;
}

double ModelErr(
	std::vector<double> params,
	std::vector<std::vector<std::vector<double> > > args)
{
	int quantity = params.size();
	std::vector<std::vector<double> > K = args[0];
	std::vector<std::vector<double> > data = args[1];
	std::vector<std::vector<double> > vega = args[2];
	std::vector<double> tau = args[3][0];
	std::vector<double> rd = args[3][1];
	std::vector<double> F = args[3][2];
	int errType = (int)args[3][3][0];
	int modType = (int)args[3][3][1];
	int fellerMust = (int)args[3][3][2];
	int modsize;
	switch (modType) {
		case 1: modsize = 5; break;
		case 2: modsize = 5; break;
		case 3: modsize = 10; break;
		case 4: modsize = 10; break;
		case 5: modsize = 10; break;
	}	
	std::vector<double> paramsPos(modsize);
	//const int m=200;
	switch (quantity) {
	case 8: {
		paramsPos[0] = exp(params[0]);
		paramsPos[1] = exp(params[1]);
		paramsPos[2] = tanh(params[2]);
		paramsPos[3] = exp(params[3]);
		paramsPos[4] = args[3][4][0];
		paramsPos[5] = exp(params[4]);
		paramsPos[6] = exp(params[5]);
		paramsPos[7] = tanh(params[6]);
		paramsPos[8] = exp(params[7]);
		paramsPos[9] = args[3][4][1];
		break; }
	case 10: {
		paramsPos = expv(params);
		paramsPos[2] = tanh(params[2]);
		paramsPos[7] = tanh(params[7]);
		break; }
	case 6: {
		paramsPos[0] = exp(params[0]);
		paramsPos[1] = args[3][4][0];
		paramsPos[2] = args[3][4][1];
		paramsPos[3] = exp(params[1]);
		paramsPos[4] = exp(params[2]);
		paramsPos[5] = exp(params[5]);
		paramsPos[6] = args[3][4][2];
		paramsPos[7] = args[3][4][3];
		paramsPos[8] = exp(params[1]);
		paramsPos[9] = exp(params[2]);
		break; }
	case 4: {
		paramsPos[0] = args[3][4][0];
		paramsPos[1] = exp(params[0]);
		paramsPos[2] = tanh(params[1]);
		paramsPos[3] = args[3][4][1];
		paramsPos[4] = args[3][4][2];
		paramsPos[5] = args[3][4][3];
		paramsPos[6] = exp(params[2]);
		paramsPos[7] = tanh(params[3]);
		paramsPos[8] = args[3][4][4];
		paramsPos[9] = args[3][4][5];
		break; }		
	case 5: {
		paramsPos = expv(params);
		paramsPos[2] = tanh(params[2]);
		break; }
		//theta,om,rho,ka,nu
	case 3: {
		paramsPos[0] = exp(params[0]);
		paramsPos[1] = args[3][4][0];
		paramsPos[2] = args[3][4][1];
		paramsPos[3] = exp(params[1]);
		paramsPos[4] = exp(params[2]);
		break; }
	/*case 2: {
		paramsPos[0] = args[3][4][0];
		paramsPos[1] = args[3][4][1];
		paramsPos[2] = args[3][4][2];
		paramsPos[3] = args[3][4][3];
		paramsPos[4] = exp(params[0]);
		paramsPos[5] = args[3][4][4];
		paramsPos[6] = args[3][4][5];
		paramsPos[7] = args[3][4][6];
		paramsPos[8] = args[3][4][7];
		paramsPos[9] = exp(params[1]);
		break; }*/
	case 2: {
		paramsPos[0] = args[3][4][0];
		paramsPos[1] = exp(params[0]);
		paramsPos[2] = tanh(params[1]);
		paramsPos[3] = args[3][4][1];
		paramsPos[4] = args[3][4][2];
		break; }
	}
	double lambda = 0, err = 0;
	std::vector<double> paramsPos2(10);
			paramsPos2[0] = paramsPos[0];
			paramsPos2[1] = paramsPos[1];
			paramsPos2[2] = paramsPos[2];
			paramsPos2[3] = paramsPos[3];
			paramsPos2[4] = paramsPos[4];
			paramsPos2[5] = paramsPos[0];
			paramsPos2[6] = paramsPos[1];
			paramsPos2[7] = paramsPos[2];
			paramsPos2[8] = paramsPos[3];
			paramsPos2[9] = paramsPos[4];
	std::vector<double> u = linspacing(-17.0, 4.9, 0.4);
	int m = u.size();
	std::vector<dcomp> I1(m);
	//#if (2*kappa*theta<sigma**2) //#feller cond.
	int n = tau.size();
	if (!feller(paramsPos) && fellerMust==1)
		err=9999;
	else if (!feller2(paramsPos) && fellerMust==2)
		err=9999;	
	else {
	for (int i = 0; i < n; i++)
	{
		switch (modType) {
		case 1://Heston
			for (int k = 0; k < m; k++)
				I1[k] = (1.0-i1/exp(u[k]))/(1.0+pow(exp(u[k]),2))*exp(u[k])*CFa(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau[i], 0);
				//I1[k] = CFa(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau[i], 0);
			break;
		case 2://Schoebel-Zhu
			for (int k = 0; k < m; k++) 
				I1[k] = (1.0-i1/exp(u[k]))/(1.0+pow(exp(u[k]),2))*exp(u[k])*CFaSZ(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau[i], 0);
			break;
		case 3://Bates
			for (int k = 0; k < m; k++)
				I1[k] = (1.0-i1/exp(u[k]))/(1.0+pow(exp(u[k]),2))*exp(u[k])*CFa2d(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau[i], 0);
			break;
		case 4://OUOU
			for (int k = 0; k < m; k++)
				I1[k] = (1.0-i1/exp(u[k]))/(1.0+pow(exp(u[k]),2))*exp(u[k])*CFaSZ2d(dcomp(exp(u[k]),0.0), paramsPos, lambda, tau[i], 0);
			break;			
		case 5: {//OUOU symetrical
			//std::vector<double> paramsPos2 = rep(paramsPos,2);
			for (int k = 0; k < m; k++)
				I1[k] = (1.0-i1/exp(u[k]))/(1.0+pow(exp(u[k]),2))*exp(u[k])*CFaSZ2d(dcomp(exp(u[k]),0.0), paramsPos2, lambda, tau[i], 0);
			break; }		
		}
		switch (errType) {
		case 1://SSE
			for (int j = 0; j < 5; j++)
			{
				err += pow(attariPre(F[i], K[i][j], rd[i], tau[i], u, I1)/vega[i][j] - data[i][j],2);
			}
			break;
		case 2://SAE
			for (int j = 0; j < 5; j++)
			{
				err += fabs(attariPre(F[i], K[i][j], rd[i], tau[i], u, I1) / vega[i][j] - data[i][j]);
			}
			break;
		case 3://SAPE
			for (int j = 0; j < 5; j++)
			{
				err += fabs(attariPre(F[i], K[i][j], rd[i], tau[i], u, I1) / vega[i][j] / data[i][j] - 1);
			}
			break;			
		}
	}}
	return err;
}


// [[Rcpp::export]]
std::vector<double> nlOpter(
	std::vector<double> params,
	std::vector<double> data1,
	std::vector<double> vega1,
	std::vector<double> F, 
	std::vector<double> K1,
	std::vector<double> rd,
	std::vector<double> tau,
	int quantity=5,
	int errType=1,
	int modType=1,
	bool err = false,
	int max_it = 1600,
	int fellerMust = 0,
	std::string logfile = "")
{
	int len = tau.size();
	std::vector<std::vector<double> > K = makeMatrix(K1,len,5);
	std::vector<std::vector<double> > data = makeMatrix(data1,len,5);
	std::vector<std::vector<double> > vega = makeMatrix(vega1,len,5);
	std::vector<std::vector<double> > rates;
	rates.push_back(tau);
	rates.push_back(rd);
	rates.push_back(F);
	std::vector<std::vector<std::vector<double> > > args;
	args.push_back(K);
	args.push_back(data);
	args.push_back(vega);
	std::vector<double> types;
	types.push_back((double)errType);
	types.push_back((double)modType);
	types.push_back((double)fellerMust);
	std::vector<double> other;
	std::vector<double> init(quantity);
	switch (quantity) {
	case 8: {
		init[0] = log(params[0]);
		init[1] = log(params[1]);
		init[2] = atanh(params[2]);
		init[3] = log(params[3]);
		init[4] = log(params[5]);
		init[5] = log(params[6]);
		init[6] = atanh(params[7]);
		init[7] = log(params[8]);
		other.push_back(params[4]);
		other.push_back(params[9]);		
		break; }
	case 10: {
		init = logv(params);		
		init[2] = atanh(params[2]);	
		init[7] = atanh(params[7]);	
		break; }
	case 6: {
		init[0] = log(params[0]);
		init[1] = log(params[3]);
		init[2] = log(params[4]);
		other.push_back(params[1]);
		other.push_back(params[2]);
		init[3] = log(params[5]);
		init[4] = log(params[8]);
		init[5] = log(params[9]);
		other.push_back(params[6]);
		other.push_back(params[7]);		
		break; }
	case 4: {
		init[0] = log(params[1]);
		init[1] = atanh(params[2]);
		other.push_back(params[0]);
		other.push_back(params[3]);
		other.push_back(params[4]);
		init[2] = log(params[6]);
		init[3] = atanh(params[7]);
		other.push_back(params[5]);
		other.push_back(params[8]);
		other.push_back(params[9]);
		break;}		
	case 5: {
		init = logv(params);		
		init[2] = atanh(params[2]);
		break; }
	case 3: {
		init[0] = log(params[0]);
		init[1] = log(params[3]);
		init[2] = log(params[4]);
		other.push_back(params[1]);
		other.push_back(params[2]);
		break; }
	/*case 2: {
		other.push_back(params[0]);
		other.push_back(params[1]);
		other.push_back(params[2]);
		other.push_back(params[3]);
		other.push_back(params[5]);
		other.push_back(params[6]);
		other.push_back(params[7]);
		other.push_back(params[8]);
		init[0] = log(params[4]);
		init[1] = log(params[9]);		
		break; }	*/
	case 2: {
		init[0] = log(params[1]);
		init[1] = atanh(params[2]);
		other.push_back(params[0]);
		other.push_back(params[3]);
		other.push_back(params[4]);
		break;}
	}
	rates.push_back(types);
	rates.push_back(other);
	args.push_back(rates);
	//cout << "The calib Heston params are: " << endl;
	//printcon(BT::Simplex(HestonErr, init, args, 1e-11));
	//return HestonErr(init,args);
	std::vector<double> out_opt = expv(BT::Simplex(ModelErr, init, args, 1e-11,max_it, logfile));
	std::vector<double> out(params.size());
	switch (quantity) {
	case 8: //all but nu0 1&2
		out[0] = out_opt[0];
		out[1] = out_opt[1];
		out[2] = tanh(log(out_opt[2]));
		out[3] = out_opt[3];
		out[4] = params[4];
		out[5] = out_opt[4];
		out[6] = out_opt[5];
		out[7] = tanh(log(out_opt[6]));
		out[8] = out_opt[7];
		out[9] = params[9];
		break;
	case 10:
		out = out_opt;
		out[2] = tanh(log(out_opt[2]));
		out[7] = tanh(log(out_opt[7]));		
		break;
	case 6://all but rho and omega
		out[0] = out_opt[0];
		out[1] = params[1];
		out[2] = params[2];
		out[3] = out_opt[1];
		out[4] = out_opt[2];
		out[5] = out_opt[3];
		out[6] = params[6];
		out[7] = params[7];
		out[8] = out_opt[4];
		out[9] = out_opt[5];
		break;
	case 4:
		out[0] = params[0];
		out[1] = out_opt[0];
		out[2] = tanh(log(out_opt[1]));
		out[3] = params[3];
		out[4] = params[4];
		out[5] = params[5];
		out[6] = out_opt[2];
		out[7] = tanh(log(out_opt[3]));
		out[8] = params[8];
		out[9] = params[9];
		break;
	case 5:
		out = out_opt;
		out[2] = tanh(log(out_opt[2]));
		break;
	case 3://all but rho and omega
		out[0] = out_opt[0];
		out[1] = params[1];
		out[2] = params[2];
		out[3] = out_opt[1];
		out[4] = out_opt[2];
		break;
	/*case 2: // only nu0 1&2
		out[0] = params[0];
		out[1] = params[1];
		out[2] = params[2];
		out[3] = params[3];
		out[4] = out_opt[0];
		out[5] = params[5];
		out[6] = params[6];
		out[7] = params[7];
		out[8] = params[8];
		out[9] = out_opt[1];
		break;	*/
	case 2: //only rho and omega
		out[0] = params[0];
		out[1] = out_opt[0];
		out[2] = tanh(log(out_opt[1]));
		out[3] = params[3];
		out[4] = params[4];
		break;
	}
	if (err) out.push_back(sqrt(ModelErr(logv(out_opt),args)/30));
	return out;
}

double max0(double a) {
	return (a<0)?0:a;
}

// [[Rcpp::export]]
std::vector<double> heston2d_pather(std::vector<double> params, long m=100000)
{
	//const double T[] = {0.08333333, 0.1666667, 0.25, 0.5, 1., 2.};
	const double delta=0.003846154;
	const int times[] = {21,  42,  63, 126, 252, 504};
	std::vector<std::vector<double> > mat(6, std::vector<double>(2*m));
	double theta, sigma, rho, kappa, nu0, theta2, sigma2, rho2, kappa2, nu02;
	double s, sb, v, vb, v2, v2b, W1, W2, W3, W4;
	theta = params[0]; sigma = params[1]; rho = params[2]; kappa = params[3]; nu0 = params[4];
	theta2 = params[5]; sigma2 = params[6]; rho2 = params[7]; kappa2 = params[8]; nu02 = params[9];
	std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<double> dist(0, 1);
	for(int j=0; j<m; j++) {
		s=0, sb=0, v=nu0, vb=nu0, v2=nu02, v2b=nu02;
		for(int i=1; i<times[5]+1; i++)
		{
			W1 = dist(rd);
			W2 = rho*W1 + sqrt(1-pow(rho,2))*dist(rd);
			W3 = dist(rd);
			W4 = rho2*W3 + sqrt(1-pow(rho2,2))*dist(rd);
			s = s - (max0(v)+max0(v2))/2*delta + sqrt(max0(v))*sqrt(delta)*W1 + sqrt(max0(v2))*sqrt(delta)*W3;
			v = v + kappa*(theta-max0(v))*delta +sigma*sqrt(max0(v))*sqrt(delta)*W2 + delta*(W2*W2-1)*pow(sigma,2)/4;
			v2 = v2 + kappa2*(theta2-max0(v2))*delta + sigma2*sqrt(max0(v2))*sqrt(delta)*W4 + delta*(W4*W4-1)*pow(sigma,2)/4;
			sb = sb - (max0(vb)+max0(v2b))/2*delta - sqrt(max0(vb))*sqrt(delta)*W1 - sqrt(max0(v2b))*sqrt(delta)*W3;
			vb = vb + kappa*(theta-max0(vb))*delta - sigma*sqrt(max0(vb))*sqrt(delta)*W2 + delta*(W2*W2-1)*pow(sigma,2)/4;
			v2b = v2b + kappa2*(theta2-max0(v2b))*delta - sigma2*sqrt(max0(v2b))*sqrt(delta)*W4 + delta*(W4*W4-1)*pow(sigma,2)/4;
			if (i==times[0]+1) {
					mat[0][2*j] = s;
					mat[0][2*j+1] = sb;
			}
			else if(i==times[1]+1){
					mat[1][2*j] = s;
					mat[1][2*j+1] = sb;
			}
			else if(i==times[2]+1){
					mat[2][2*j] = s;
					mat[2][2*j+1] = sb;
			}
			else if(i==times[3]+1){
					mat[3][2*j] = s;
					mat[3][2*j+1] = sb;
			}
			else if(i==times[4]+1){
					mat[4][2*j] = s;
					mat[4][2*j+1] = sb;
			}
			else if(i==times[5]){
					mat[5][2*j] = s;
					mat[5][2*j+1] = sb;
			}
		}
	}
	return mat2vec(mat);
}
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
// [[Rcpp::export]]
std::vector<double> ouou2_pather(std::vector<double> params, long m=100000)
{
	//const double T[] = {0.08333333, 0.1666667, 0.25, 0.5, 1., 2.};
	const double delta=0.003846154;
	const int times[] = {21,  42,  63, 126, 252, 504};
	std::vector<std::vector<double> > mat(6, std::vector<double>(2*m));
	double theta, sigma, rho, kappa, nu0, theta2, sigma2, rho2, kappa2, nu02;
	double s, sb, v, vb, v2, v2b, W1, W2, W3, W4, W2b, W4b, tmp;
	theta = params[0]; sigma = params[1]; rho = params[2]; kappa = params[3]; nu0 = params[4];
	theta2 = params[5]; sigma2 = params[6]; rho2 = params[7]; kappa2 = params[8]; nu02 = params[9];
	std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<double> dist(0, 1);
	for(int j=0; j<m; j++) {
		s=0, sb=0, v=nu0, vb=nu0, v2=nu02, v2b=nu02;
		for(int i=1; i<times[5]+1; i++)
		{
			W1 = dist(rd);
			tmp = dist(rd);
			W2 = sgn(v)*rho*W1 + sqrt(1-pow(rho,2))*tmp;
			W2b = sgn(vb)*rho*W1 + sqrt(1-pow(rho,2))*tmp;
			W3 = dist(rd);
			tmp = dist(rd);
			W4 = sgn(v2)*rho2*W3 + sqrt(1-pow(rho2,2))*tmp;
			W4b = sgn(v2b)*rho2*W3 + sqrt(1-pow(rho2,2))*tmp;
			s = s - (v*v+v2*v2)/2*delta + v*sqrt(delta)*W1 + v2*sqrt(delta)*W3;
			v = v + kappa*(theta-v)*delta +sigma*sqrt(delta)*W2 ;
			v2 = v2 + kappa2*(theta2-v2)*delta + sigma2*sqrt(delta)*W4 ;
			sb = sb - (vb*vb+v2b*v2b)/2*delta - vb*sqrt(delta)*W1 - v2b*sqrt(delta)*W3;
			vb = vb + kappa*(theta-vb)*delta - sigma*sqrt(delta)*W2b ;
			v2b = v2b + kappa2*(theta2-v2b)*delta - sigma2*sqrt(delta)*W4b ;
			if (i==times[0]+1) {
					mat[0][2*j] = s;
					mat[0][2*j+1] = sb;
			}
			else if(i==times[1]+1){
					mat[1][2*j] = s;
					mat[1][2*j+1] = sb;
			}
			else if(i==times[2]+1){
					mat[2][2*j] = s;
					mat[2][2*j+1] = sb;
			}
			else if(i==times[3]+1){
					mat[3][2*j] = s;
					mat[3][2*j+1] = sb;
			}
			else if(i==times[4]+1){
					mat[4][2*j] = s;
					mat[4][2*j+1] = sb;
			}
			else if(i==times[5]){
					mat[5][2*j] = s;
					mat[5][2*j+1] = sb;
			}
		}
	}
	return mat2vec(mat);
}
// [[Rcpp::export]]
std::vector<double> ouou_pather(std::vector<double> params, long m=100000)
{
	//const double T[] = {0.08333333, 0.1666667, 0.25, 0.5, 1., 2.};
	const double delta=0.003846154;
	const int times[] = {21,  42,  63, 126, 252, 504};
	std::vector<std::vector<double> > mat(6, std::vector<double>(2*m));
	double theta, sigma, rho, kappa, nu0, theta2, sigma2, rho2, kappa2, nu02;
	double s, sb, v, vb, v2, v2b, W1, W2, W3, W4, W2b, W4b, tmp;
	theta = params[0]; sigma = params[1]; rho = params[2]; kappa = params[3]; nu0 = params[4];
	theta2 = params[5]; sigma2 = params[6]; rho2 = params[7]; kappa2 = params[8]; nu02 = params[9];
	std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<double> dist(0, 1);
	for(int j=0; j<m; j++) {
		s=0, sb=0, v=nu0, vb=nu0, v2=nu02, v2b=nu02;
		for(int i=1; i<times[5]+1; i++)
		{
			W1 = dist(rd);
			tmp = dist(rd);
			W2 = rho*W1 + sqrt(1-pow(rho,2))*tmp;
			W2b = rho*W1 + sqrt(1-pow(rho,2))*tmp;
			W3 = dist(rd);
			tmp = dist(rd);
			W4 = rho2*W3 + sqrt(1-pow(rho2,2))*tmp;
			W4b = rho2*W3 + sqrt(1-pow(rho2,2))*tmp;
			s = s - (v*v+v2*v2)/2*delta + v*sqrt(delta)*W1 + v2*sqrt(delta)*W3;
			v = v + kappa*(theta-v)*delta +sigma*sqrt(delta)*W2 ;
			v2 = v2 + kappa2*(theta2-v2)*delta + sigma2*sqrt(delta)*W4 ;
			sb = sb - (vb*vb+v2b*v2b)/2*delta - vb*sqrt(delta)*W1 - v2b*sqrt(delta)*W3;
			vb = vb + kappa*(theta-vb)*delta - sigma*sqrt(delta)*W2b ;
			v2b = v2b + kappa2*(theta2-v2b)*delta - sigma2*sqrt(delta)*W4b ;
			if (i==times[0]+1) {
					mat[0][2*j] = s;
					mat[0][2*j+1] = sb;
			}
			else if(i==times[1]+1){
					mat[1][2*j] = s;
					mat[1][2*j+1] = sb;
			}
			else if(i==times[2]+1){
					mat[2][2*j] = s;
					mat[2][2*j+1] = sb;
			}
			else if(i==times[3]+1){
					mat[3][2*j] = s;
					mat[3][2*j+1] = sb;
			}
			else if(i==times[4]+1){
					mat[4][2*j] = s;
					mat[4][2*j+1] = sb;
			}
			else if(i==times[5]){
					mat[5][2*j] = s;
					mat[5][2*j+1] = sb;
			}
		}
	}
	return mat2vec(mat);
}

// [[Rcpp::export]]
std::vector<double> test_rnd(long n, long m=100000)
{
	std::vector<double> out;
	double W1,W2,W3,W4;
	std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<double> dist(0, 1);
	for(int j=0; j<m; j++) {
		for(int i=1; i<n; i++)
		{
			W1 = dist(rd);
			W2 = dist(rd);
			W3 = dist(rd);
			W4 = dist(rd);
			out.push_back(W1);
			out.push_back(W2);
			out.push_back(W3);
			out.push_back(W4);
		}
	}
	return out;
}
