

function res = resid_pure_MZ_varFz(P , MZ , KAPPA , ALPHA , GAMMA , FZ_vec , FY_vec  , tyre_data )
%                                     vec , num   , vec   , num   , vec    , vec , struct
    tmp_tyre_data = tyre_data;

    % [ qBz2 , qBz3 , qDz2 , qDz7 , qEz2 , qEz3 , qHz2 ] 
    tmp_tyre_data.qBz2 = P(1);
    tmp_tyre_data.qBz3 = P(2);
    tmp_tyre_data.qDz2 = P(3);
    tmp_tyre_data.qDz7 = P(4);
    tmp_tyre_data.qEz2 = P(5);
    tmp_tyre_data.qEz3 = P(6);
    tmp_tyre_data.qHz2 = P(7);
    
%     tmp_tyre_data.qHz4 = P(8);




    res = 0;
    for i = 1:length(ALPHA)

        MZr0 = MF96_MZr(KAPPA , ALPHA(i) , GAMMA , FZ_vec(i) , tmp_tyre_data);
        t    = MF96_t(KAPPA , ALPHA(i) , GAMMA , FZ_vec(i) , tmp_tyre_data);
        MZ0  = MF96_MZ0(t , FY_vec(i) , MZr0);
        res  = res+(MZ0-MZ(i))^2;

    end

    res = res/sum(MZ.^2);

end