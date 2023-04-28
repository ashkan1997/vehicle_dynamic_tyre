% FX
function [fx_vec, Gxa] = MF96_FX_vec(fx0_vec, kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

 % precode
 for i = 1:length(fx0_vec)
     [Gxa(i), Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);

 % main code

     fx_vec(i) = Gxa(i) * fx0_vec(i);

 end

end
