function MZr = MF96_MZr(kappa, alpha, phi, Fz, tyre_data)
    
    
    [Br, ~, ~, Dr, ~, ~, alpha__r, ~] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);

    t1 = atan(Br * alpha__r);

    MZr = Dr * cos(t1) * cos(alpha);



end