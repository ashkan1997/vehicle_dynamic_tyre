function [MZr_vec] = MF96_MZr_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

    MZr_vec = zeros(size(alpha_vec));

    for i= 1:length(alpha_vec)
        
        
        [Br, ~, ~, Dr, ~, ~, alpha__r, ~] = MF96_MZ0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);

        t1 = atan(Br * alpha__r);

        MZr_vec(i) = Dr * cos(t1) * cos(alpha_vec(i));

    end


end