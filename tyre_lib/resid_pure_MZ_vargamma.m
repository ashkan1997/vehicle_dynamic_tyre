function res = resid_pure_MZ_vargamma(P , MZ , KAPPA , ALPHA , GAMMA_vec , FZ , FY_vec  , tyre_data )

    tmp_tyre_data = tyre_data;

    % [qBz4 , qBz5 , qDz3 , qDz4  ,  qDz8  ,  qDz9 , qEz5 , qHz3 , qHz4]


    tmp_tyre_data.qBz4 = P(1);
    tmp_tyre_data.qBz5 = P(2);
    tmp_tyre_data.qDz3 = P(3);
    tmp_tyre_data.qDz4 = P(4);
    tmp_tyre_data.qDz8 = P(5);
    tmp_tyre_data.qDz9 = P(6);
    tmp_tyre_data.qEz5 = P(7);
    tmp_tyre_data.qHz3 = P(8);
    tmp_tyre_data.qHz4 = P(9);

    res = 0;
    for i = 1:length(ALPHA)

        MZr0 = MF96_MZr(KAPPA , ALPHA(i) , GAMMA_vec(i) , FZ , tmp_tyre_data);
        t    = MF96_t(KAPPA , ALPHA(i) , GAMMA_vec(i) , FZ , tmp_tyre_data);
        MZ0  = MF96_MZ0(t , FY_vec(i) , MZr0);
        res  = res+(MZ0-MZ(i))^2;

    end

    res = res/sum(MZ.^2);

end
