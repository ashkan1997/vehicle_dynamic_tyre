function [MZ0_vec] = MF96_MZ0_vec(t_vec , FY_vec , MZr_vec)

    MZ0_vec = zeros(size(t_vec));

    for i= 1:length(t_vec)

        MZ0_vec(i) = -t_vec(i) * FY_vec(i) + MZr_vec(i);

    end


end