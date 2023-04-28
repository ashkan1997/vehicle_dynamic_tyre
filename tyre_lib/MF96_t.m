function t = MF96_t(kappa, alpha, phi, Fz, tyre_data)

    [~, Bt, Ct, ~, Dt, Et, ~, alpha__t] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);

    t1 = Bt * alpha__t ;

%     t2 = Ct * atan(t1);
    t2 = Et * (t1 - atan(t1));

    t = Dt * cos(Ct * atan(t1 - t2)) * cos(alpha);





end