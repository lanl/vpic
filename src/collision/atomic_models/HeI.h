#ifndef _HeI_h_
#define _HeI_h_
#include <float.h>
#include "atom.h"

/* Atomic data routines for collisions with Helium atoms.

  The coefficents used here are based on the methodology developed by
  Khrabrov and Kaganovich and Wang et al. 2017. The basic idea is to
  construct a plausible physical model for the differential scattering
  cross-section and then use experimental measurements of the total,
  momentum, and/or viscous cross-sections to determine free-parameters
  in the model. Specifically, we use the following models:

  en_el    : Scattering in a screened Coulomb potential with an energy
             dependent screening-length.
  et_ex     : Cross-sections are a degenercy averaged cross section for
              excitation from 11S -> nls with n <= 4. Raw cross-sections
              are taken from the analytic fits of Ralchenko 2008. Gamma
              is chosen in a consistent manner. Scattering angle uses
              the elastic angles.
  et_in     : Cross-section is direct ionization from 11S taken from
              Ralchenko 2008.
  ep_rx    : Recombination is tricky. At high energy, resonances due to
             dielectronic recombination dominate the total cross-section,
             but are difficult to include directly since they require very
             high resolution in energy space to capture properly. Here,
             raw data from Nahar 2010 is used and pre-processed to give
             a correct resonace averaged structure.
  in_el/cx : Forward scattering due to a polarization potential and
             backwards scattering due to charge exchange. Forward and
             backward scattering amplitudes are independent with energy
             dependent amplitudes.
  nn_el    : Normally distributed forward scattering with variance
             chosen to match total and viscous cross-sections.

  For details on how the scattering angles are chosen, see Khraborv and
  Kaganovich and Wang et al. However, here we have changed some details
  in the actual implementation. Specifically, the references cited
  specifically compute cos(theta) of the scattering angle, however this
  is not numerically stable for small theta and will not conserve
  momentum. Here we compute delta = tan(theta/2) instead and use that to
  then compute sin(theta) and cos(theta) in the actual collision algorithm.
  This way minimizes the effect of roundoff, but the tradeoff is that we
  have to worry about delta going to infinity for perfect backscattering.
  Also, instead of using look-up tables, we are using polynomial
  interpolants in log2 space.

  References:
    Nahar, S.N., New Astronomy (15), 2010
    Ralchenko, Y. et al., Atomic Data and Nuclear Data Tables (94), 2008

*/


float
HeI_sigma_pt_el(float E){
  E = log2(E + 7.349865e-04f);
  return exp2(((((( 6.012956e-06 *E +2.815147e-05)*E
                   -3.985465e-04)*E -1.499984e-03)*E
                   +9.442074e-03)*E -2.464035e-01)*E +7.177672e+00);
}

float
HeI_sigma_pt_cx(float E){
  E = log2(E + 7.349865e-04f);
  return exp2(((((( 4.581837e-06 *E +1.009797e-05)*E
                   -3.056459e-04)*E -2.674257e-04)*E
                   +4.621865e-04)*E -1.343024e-01)*E +6.043904e+00);
}

float
HeI_sigma_et_el(float E){
  E = log2(E + 7.349865e-04f);
  return exp2(((((( 4.307671e-06 *E +1.273896e-04)*E
                   +8.749181e-04)*E -6.549400e-03)*E
                   -1.136261e-01)*E -4.852011e-01)*E +3.771030e+00);
}

float
HeI_sigma_et_ex(float E){
  float x = E/1.4568;
  if (E < 1) return 0;

  // This is a sum of all excitation cross-sections out of 11S.
  // This may look very long, but is actually shorter than directly
  // computing the first 4 directly (which is needed for 10% accuracy)
  x = log2(x);
  return exp2((((((((((( 6.555543e-05 *x -2.368052e-03)*x
                        +3.769421e-02)*x -3.476701e-01)*x
                        +2.054920e+00)*x -8.124464e+00)*x
                        +2.175721e+01)*x -3.904986e+01)*x
                        +4.534678e+01)*x -3.196479e+01)*x
                        +1.249841e+01)*x -3.129806e+00);

}

float
HeI_sigma_et_in(float E){
  if( E <= 1.807140 ) return 0;

  // Equation 9 from Ralchenko 2008.
  float x = 1.807140/E; // I/E
  float y = 1-x;        // (1-I/E)
  return 5.9071e+00*x*(-5.857e-01*log(x)+(((3.317e+00 *y -2.521e+00)*y
                                           +7.680e-01)*y -4.457e-01)*y);
}

float
HeI_sigma_ep_rx(float E){

  // Radiative recombinaton. Simple but 5% accurate for E >~ 2e-3
  float sigma = 0.9e-5*pow(E+7.349865e-04f, -1.125);
  if( E < 2 || E > 4 ) return sigma;

  // Resonance averaged dielectronic recombination. These resonances
  // are very sharp, but I am averaging them over an energy interval
  // of 0.05 Ry to smear them out. Note when doing averaging, we are
  // really doing <sigma> = <sigma v> / v so that the total rate is
  // conserved (assuming that <.> denotes convolution with a window
  // of unit area). This is a little overkill here, we could get
  // away with only 3 peaks probably, but why not.

  float x;
  x = E - 2.889545;
  sigma += 2.01915e-04 * exp( -96.345*x*x );
  x = E - 2.994697;
  sigma += 4.39237e-03 * exp( -195.95*x*x );
  x = E - 3.703070;
  sigma += 1.55912e-05 * exp( -21.751*x*x );
  x = E - 2.584947;
  sigma += 4.86425e-05 * exp( -70.603*x*x );
  x = E - 3.551628;
  sigma += 1.89497e-04 * exp( -214.397*x*x );

  return sigma;

}

float
HeI_sigma_tt_el(float E){
  E = log2(E + 7.349865e-04f);
  return exp2(((( 1.354527e-05 *E +5.235340e-05)*E
                 -5.428023e-03)*E -1.122875e-01)*E +6.294592e+00);
}

float
HeI_delta_pt_el(rng_t * rng, float E){
  const float sqrt_sqrt_half = 0.8408964152537146f;
  const float one_8          = 0.125f;
  const float five_128       = 0.0390625f;
  const float fifteen_1024   = 0.0146484375f;

  float inv_sqrt_sqrt_a, a, x;
  E = log2(E + 7.349865e-04f);
  a = (((((( -5.692770e-06 *E -1.911994e-04)*E
             +3.192144e-05)*E +2.383839e-02)*E
             -1.522081e-01)*E -3.641591e+00)*E -2.128630e+01);
  inv_sqrt_sqrt_a = exp2(-0.25f*a);
  a = exp2(a);
  x = inv_sqrt_sqrt_a - sqrt_sqrt_half*(1 - a*(one_8 - a*(five_128 - a*fifteen_1024)));
  x = inv_sqrt_sqrt_a - x*frand(rng);
  x *= x;
  x *= x;
  return sqrtf( (1-a*x) / ((2+a)*x - 1) );
}

float
HeI_delta_pt_cx(rng_t * rng, float E){
  const float sqrt_sqrt_half = 0.8408964152537146f;
  const float one_8          = 0.125f;
  const float five_128       = 0.0390625f;
  const float fifteen_1024   = 0.0146484375f;

  float inv_sqrt_sqrt_b, b, x;
  E = log2(E + 7.349865e-04f);
  b = (((((( -1.070780e-05 *E -2.346134e-04)*E
             +6.388344e-04)*E +2.652242e-02)*E
             -1.743078e-01)*E -3.978732e+00)*E -1.685851e+01);
  inv_sqrt_sqrt_b = exp2(-0.25f*b);
  b = exp2(b);
  x = inv_sqrt_sqrt_b - sqrt_sqrt_half*(1 - b*(one_8 - b*(five_128 - b*fifteen_1024)));
  x = inv_sqrt_sqrt_b + x*(frand(rng)-1);
  x *= x;
  x *= x;
  return sqrtf(2*x/(1-b*x) - 1);
}

float
HeI_delta_et_el(rng_t * rng, float E){
  float x,z;
  E *= 13.605693f;
  z = 2.45f * sqrtf(E);
  z = 1 + (z-19.9324f)/(E-2.302041f*z+19.9324f) - z/(E-4.171429f*z+90.1221f);
  x = 1-2*frand(rng);
  x = (z + x)/(1 + z*x);
  return sqrtf((1-x)/(1+x));
}

float
HeI_delta_tt_el(rng_t * rng, float E){
  E = log2(E + 7.349865e-04f);
  return frandn(rng) * exp2((((-8.664093e-05 *E -2.570293e-03)*E
                               -2.534166e-02)*E -2.399119e-01)*E
                               -2.575432e+00);
}

float
HeI_gamma_et_ex(rng_t *rng, float E){
  // At each energy, we should choose the excitation energy E_i with
  // probability sigma_i / sigma_tot, However, this requires evaluating
  // many cross-sections. Instead, we just use <E_i>(E), then compute
  // gamma(E) = sqrt(1-<Ei>/E)
  float x = E/1.4568 - 1;
  if(x <= 0) return x < 0 ;
  x = log2(x);
  x = exp2((2.327251e-02*x +8.498602e-01)*x +1.086187e+00);
  return x/(1+x);
}

float
HeI_gamma_et_in(rng_t *rng, float Eprimary, float *Esecondary){

  // In the non-relativistic binary-encounter-dipole model, the singly
  // differential cross-section is given by
  //
  // d sigma(W,T)/ dW = (S/B)/(t+u+1) [ ((Ni/N)-2)/(t+1) (1/(w+1) + 1/(t-w))
  //                                  + (2-Ni/N) (1/(w+1)^2 + 1/(t-w)^2)
  //                                  + ln(t)/(N (w+1)) df/dw ]
  //
  // where :
  //    W = Energy of secondary electron
  //    B = Binding energy of electron
  //    T = Reduced kinetic energy of primary me v^2 / 2
  //    N = Number of bound electrons in the subshell
  //    U = Average kinetic energy of bound electrons
  //    df/dw = differential oscillator strength
  //    Ni    = int_0^inf df/dw dw = total oscillator strength
  // and lower case letters denote normalization by B.
  //
  // For ground state He:
  //   u        = 1.607
  //   N        = 2
  //   df/dw    = 8.24012e+00 y^3 -1.04769e+01 y^4 +3.96496e+00 y^5 +4.45976e-02 y^6
  //   Ni/N - 2 = -1.18604
  //   y = 1/(1+w)
  //
  // Now integrating the cross-section to obtain the CDF, we find
  //
  //   CDF(w|t) ~   (Ni/N-2) / t+1 [   ln(w+1/t-w) + (t+1) (1/w+1 - 1/t-w)
  //                                 + ln(t)       - (t+1) (1       - 1/t) ]
  //              - ln(t)/N sum_i a_i/i ( (1+w)^-i - 1 )
  //
  // And, the upper bound on w is given by w <= (t-1)/2 so the proportionality
  // coefficient is given by
  //
  //   CDF_norm =  (Ni/N-2) / t+1 [ ln(t) - (t+1) (1 - 1/t) ]
  //             - ln(t)/N sum_i a_i/i ( 2^i (t+1)^-i - 1 )
  //
  // Now this CDF is not easy to invert since it depends on t and w in
  // non-trivial ways, but we can fairly quickly do some Newton iterations.

  int max_iter = 100;
  float CDF, dCDFdw, w, w1, tw, err;
  const float a3   =  2.74671e+00;        // a3 / 3
  const float a4   = -2.61923e+00;        // a4 / 4
  const float a5   =  7.92992e-01;        // a5 / 5
  const float a6   =  7.43293e-03;        // a6 / 6
  const float t    = Eprimary/1.807140;   // t = E/I
  const float t1   = t+1;
  const float c2   = 0.5*log(t);   // ln(t)/N
  const float c1   = -1.18604 / t1;     // (Ni/N - 2)/(t+1)
  const float c0   = c1 * (2*c2 - t1*(1-1/t)) + c2 * (a3 + a4 + a5 + a6);

  w1 = 2/t1;
  const float norm = c0 - c2*((((a6*w1 + a5)*w1 + a4)*w1 + a3)*w1*(w1*w1) );
  const float R    = frand(rng)*norm;
  const float wmax = 0.5*(t-1);

  // SHhould never happen, but test just in case.
  if(t < 1){
    (*Esecondary) = 0;
    return 0;
  }

  // This is a good starting guess.
  w = R*wmax/(t*norm);
  for( ; max_iter ; --max_iter ) {
    w1  = 1/(w+1);
    tw  = 1/(t-w);
    CDF = c0 + c1*( -log(w1/tw) + t1*(w1 - tw) )
             - c2*((((a6*w1 + a5)*w1 + a4)*w1 + a3)*w1*(w1*w1) );
    dCDFdw = c1*( w1*(1-t1*w1) + tw*(1-t1*tw) ) +
             c2*w1*((((6*a6*w1 + 5*a5)*w1 + 4*a4)*w1 + 3*a3)*w1*(w1*w1) );
    err = CDF-R;
    w  -= err/dCDFdw;
    if(w <    0) w = 0;
    if(w > wmax) w = wmax;
    if( err*err <= 1e-6 ) break;
  }

  // Set the secondary energy and compute gamma.
  (*Esecondary) = w*1.807140;
  return sqrtf(1 - (w+1)/t);
}

/* Public interface **********************************************************/

const atomic_properties_t HeI = {
  .alpha_m2    = 0,
  .Nca02       = 0,

  // Cross-sections for the various processes.
  .sigma_pt_el = &HeI_sigma_pt_el ,   // He  + He+ -> He  + He+
  .sigma_pt_cx = &HeI_sigma_pt_cx ,   // He  + He+ -> He+ + He
  .sigma_pt_ex = NULL ,               // He  + He+ -> He* + He+
  .sigma_et_el = &HeI_sigma_et_el ,   // He  +  e  -> He  +  e
  .sigma_et_ex = &HeI_sigma_et_ex ,   // He  +  e  -> He* +  e
  .sigma_et_in = &HeI_sigma_et_in ,   // He  +  e  -> He+ + 2e
  .sigma_ep_rx = &HeI_sigma_ep_rx ,   // He+ +  e  -> He
  .sigma_tt_el = &HeI_sigma_tt_el ,   // He  + He  -> He  + He
  .sigma_tt_ex = NULL ,               // He  + He  -> He* + He

  // Scattering angles for the various elastic processes.
  // Inelastic processes use the corresponding elastic delta.
  .delta_pt_el = &HeI_delta_pt_el ,
  .delta_pt_cx = &HeI_delta_pt_cx ,
  .delta_et_el = &HeI_delta_et_el ,
  .delta_tt_el = &HeI_delta_tt_el ,

  // Coefficient of restitution for the various processes.
  .gamma_pt_ex = NULL ,
  .gamma_et_ex = &HeI_gamma_et_ex ,
  .gamma_tt_ex = NULL ,

  // Ionization coefficient of restitution and secondary energy.
  .gamma_et_in = &HeI_gamma_et_in ,

} ;

#endif /* _HeI_h_ */
