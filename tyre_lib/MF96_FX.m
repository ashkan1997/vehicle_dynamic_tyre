% FX
function [fx] = MF96_FX(fx0,kappa, alpha, phi, Fz, tyre_data)

 % precode

  [Gxa, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  fx = Gxa * fx0;
  
 end
