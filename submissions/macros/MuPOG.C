//File collecting functions from the Muon POG, used to apply corrections on muons

#include <boost/math/special_functions/erf.hpp>
#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

class SeedSequence {
public:
    explicit SeedSequence(std::initializer_list<uint32_t> seeds)
        : m_seeds(seeds) {}

    template <typename Iter>
    void generate(Iter begin, Iter end) const {
        const size_t n = std::distance(begin, end);
	if (n == 0) return;

	const uint32_t mult = 0x9e3779b9;
	const uint32_t mix_const = 0x85ebca6b;

	std::vector<uint32_t> buffer(n, 0x8b8b8b8b);

	size_t s = m_seeds.size();
	size_t t = (n >= s) ? n - s : 0;

        size_t i = 0;

	for(; i < std::min(n, s); ++i) {
            buffer[i] = buffer[i] ^ (m_seeds[i] + mult * i);
	}
	for(; i < n; ++i) {
            buffer[i] = buffer[i] ^ (mult * i);
        }

	for (size_t k = 0; k < n; ++k) {
            uint32_t z = buffer[(k + n - 1) % n] ^ (buffer[k] >> 27);
	    buffer[k] = (z * mix_const) ^ (buffer[k] << 13);
        }

	std::copy(buffer.begin(), buffer.end(), begin);
    }

private:
    std::vector<uint32_t> m_seeds;
};

struct CrystalBall{
    double pi=3.14159;
    double sqrtPiOver2=sqrt(pi/2.0);
    double sqrt2=sqrt(2.0);
    double m;
    double s;
    double a;
    double n;
    double B;
    double C;
    double D;
    double N;
    double NA;
    double Ns;
    double NC;
    double F;
    double G;
    double k;
    double cdfMa;
    double cdfPa;
CrystalBall():m(0),s(1),a(10),n(10){
    init();
}
CrystalBall(double mean, double sigma, double alpha, double n)
    :m(mean),s(sigma),a(alpha),n(n){
    init();
}
void init(){
    double fa = fabs(a);
    double ex = exp(-fa*fa/2);
    double A  = pow(n/fa, n) * ex;
    double C1 = n/fa/(n-1) * ex; 
    double D1 = 2 * sqrtPiOver2 * erf(fa/sqrt2);
    B = n/fa-fa;
    C = (D1+2*C1)/C1;   
    D = (D1+2*C1)/2;   
    N = 1.0/s/(D1+2*C1); 
    k = 1.0/(n-1);  
    NA = N*A;       
    Ns = N*s;       
    NC = Ns*C1;     
    F = 1-fa*fa/n; 
    G = s*n/fa;    
    cdfMa = cdf(m-a*s);
    cdfPa = cdf(m+a*s);
}
double pdf(double x) const{ 
    double d=(x-m)/s;
    if(d<-a) return NA*pow(B-d, -n);
    if(d>a) return NA*pow(B+d, -n);
    return N*exp(-d*d/2);
}
double pdf(double x, double ks, double dm) const{ 
    double d=(x-m-dm)/(s*ks);
    if(d<-a) return NA/ks*pow(B-d, -n);
    if(d>a) return NA/ks*pow(B+d, -n);
    return N/ks*exp(-d*d/2);

}
double cdf(double x) const{
    double d = (x-m)/s;
    if(d<-a) return NC / pow(F-s*d/G, n-1);
    if(d>a) return NC * (C - pow(F+s*d/G, 1-n) );
    return Ns * (D - sqrtPiOver2 * erf(-d/sqrt2));
}
double invcdf(double u) const{
    if(u<cdfMa) return m + G*(F - pow(NC/u, k));
    if(u>cdfPa) return m - G*(F - pow(C-u/NC, -k) );
    return m - sqrt2 * s * boost::math::erf_inv((D - u/Ns )/sqrtPiOver2);
}
};
