function res = resid_pure_Mz(P , MZ , KAPPA , ALPHA , GAMMA , FZ , FY_vec ,  tyre_data )


    tmp_tyre_data = tyre_data;
     

    tmp_tyre_data.qBz1 =  P(1);
    tmp_tyre_data.qBz9 =  P(2);
    tmp_tyre_data.qBz10 = P(3);
    tmp_tyre_data.qCz1 =  P(4);
    tmp_tyre_data.qDz1 =  P(5);
%     tmp_tyre_data.qDz2 =  P(6);
%     tmp_tyre_data.qDz3 =  P(7);
%     tmp_tyre_data.qDz4 =  P(8);
    tmp_tyre_data.qDz6 =  P(6);
    tmp_tyre_data.qEz1 =  P(7);
    tmp_tyre_data.qEz4 =  P(8);
    tmp_tyre_data.qHz1 =  P(9);


    res = 0;
    for i = 1:length(ALPHA)

        MZr0 = MF96_MZr(KAPPA , ALPHA(i) , GAMMA , FZ , tmp_tyre_data);
        t    = MF96_t(KAPPA , ALPHA(i) , GAMMA , FZ , tmp_tyre_data);
        MZ0  = MF96_MZ0(t , FY_vec(i) , MZr0);
        res  = res+(MZ0-MZ(i))^2;
        
    end

    res = res/sum(MZ.^2);


end