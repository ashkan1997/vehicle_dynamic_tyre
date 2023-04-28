function [t_vec] = MF96_t_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

    t_vec = zeros(size(alpha_vec));

    for i = 1:length(alpha_vec)
        
        [~, Bt, Ct, ~, Dt, Et, ~, alpha__t] = MF96_MZ0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);

        t1 = Bt * alpha__t ;

        %t2 = Ct * atan(t1);
        t2 = Et * (t1 - atan(t1));

        t_vec(i) = Dt * cos(Ct * atan(t1 - t2)) * cos(alpha_vec(i));



    end





end