%% SELF ALLIGNING MOMENT k = 0 , FZ = 220 , variable camber




TData0_gamma = FZ_220;
cut = [103:837 , 990:1718 , 1880:2600 , 2765:3487 , 3653:4376];
TData0_gamma = TData0_gamma(cut , :);
% idx.ALPHA = -0.20 < TData_gamma.SA & TData_gamma.SA < 0.20;
% TData0_gamma = TData_gamma(idx.ALPHA , :);
% plot_selected_data(TData0_gamma);


KAPPA_vec = TData0_gamma.SL;
ALPHA_vec = TData0_gamma.SA;
FY_vec    = TData0_gamma.FY;
FZ_vec    = TData0_gamma.FZ;
MZ_vec    = TData0_gamma.MZ;
GAMMA_vec = TData0_gamma.IA;
FZ0 =  mean(TData0_gamma.FZ);
SA_vec_gamma = linspace(min(ALPHA_vec),max(ALPHA_vec),length(ALPHA_vec));

zeros_vec = zeros(size(ALPHA_vec));
ones_vec  = ones(size(ALPHA_vec));

FY_pred = MF96_FY0_vec(zeros_vec , ALPHA_vec , GAMMA_vec , FZ0*ones_vec , tyre_coeffs);
% Guess values for parameters to be optimised
%       [qBz4 , qBz5 , qDz3 , qDz4  ,  qDz8  ,  qDz9 , qEz5 , qHz3 , qHz4] 
P0_mz_gamma = [ 1 , 1 ,  1  , 1 , 1 ,  1 , 1 , 1 , 1 ];
% P0_mz_fz = [  -3,  2,  10,  -4,  0,  3,  -3 , 5 ];
lb = [];
ub = [];


[P_cam,~,exitflag] = fmincon(@(P)resid_pure_MZ_vargamma(P,MZ_vec, 0, ALPHA_vec, GAMMA_vec , FZ0, FY_pred , tyre_coeffs),...
                               P0_mz_gamma,[],[],[],[],lb,ub);


tyre_coeffs.qBz4 = P_cam(1);
tyre_coeffs.qBz5 = P_cam(2);
tyre_coeffs.qDz3 = P_cam(3);
tyre_coeffs.qDz4 = P_cam(4);
tyre_coeffs.qDz8 = P_cam(5);
tyre_coeffs.qDz9 = P_cam(6);
tyre_coeffs.qEz5 = P_cam(7);
tyre_coeffs.qHz3 = P_cam(8);
tyre_coeffs.qHz4 = P_cam(9);

% GAMMA_0
MZr_vec_1 = MF96_MZr_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_0.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
t_vec_1 = MF96_t_vec(           KAPPA_vec , SA_vec_gamma , mean(GAMMA_0.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
FY_pred_1 = MF96_FY0_vec(zeros_vec , SA_vec_gamma , mean(GAMMA_0.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_1 = MF96_MZ0_vec(t_vec_1 , FY_pred_1 , MZr_vec_1);

% GAMMA_1
MZr_vec_2 = MF96_MZr_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_1.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
t_vec_2 = MF96_t_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_1.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
FY_pred_2 = MF96_FY0_vec(zeros_vec , SA_vec_gamma , mean(GAMMA_1.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_2 = MF96_MZ0_vec(t_vec_2 , FY_pred_2 , MZr_vec_2);

% GAMMA_2
MZr_vec_3 = MF96_MZr_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_2.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
t_vec_3 = MF96_t_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_2.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
FY_pred_3 = MF96_FY0_vec(zeros_vec , SA_vec_gamma , mean(GAMMA_2.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_3 = MF96_MZ0_vec(t_vec_3 , FY_pred_3 , MZr_vec_3);

% GAMMA_3
MZr_vec_4 = MF96_MZr_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_3.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
t_vec_4 = MF96_t_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_3.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
FY_pred_4 = MF96_FY0_vec(zeros_vec , SA_vec_gamma , mean(GAMMA_3.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_4 = MF96_MZ0_vec(t_vec_4 , FY_pred_4 , MZr_vec_4);

% GAMMA_4
MZr_vec_5 = MF96_MZr_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_4.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
t_vec_5 = MF96_t_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_4.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
FY_pred_5 = MF96_FY0_vec(zeros_vec , SA_vec_gamma , mean(GAMMA_4.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_5 = MF96_MZ0_vec(t_vec_5 , FY_pred_5 , MZr_vec_5);

% GAMMA_6
MZr_vec_6 = MF96_MZr_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_5.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
t_vec_6 = MF96_t_vec(KAPPA_vec , SA_vec_gamma , mean(GAMMA_5.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
FY_pred_6 = MF96_FY0_vec(zeros_vec , SA_vec_gamma , mean(GAMMA_5.IA)*ones_vec , FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_6 = MF96_MZ0_vec(t_vec_6 , FY_pred_6 , MZr_vec_6);


% figure

figure('Name','MZ0 vs Gamma')
plot(ALPHA_vec,TData0_gamma.MZ,'.')
hold on
plot(SA_vec_gamma,MZ0_vargamma_vec_1,'-' , 'LineWidth',2 ,'DisplayName','$\gamma_0$')
plot(SA_vec_gamma,MZ0_vargamma_vec_2,'-' , 'LineWidth',2 ,'DisplayName','$\gamma_1$')
plot(SA_vec_gamma,MZ0_vargamma_vec_3,'-' , 'LineWidth',2 ,'DisplayName','$\gamma_2$')
plot(SA_vec_gamma,MZ0_vargamma_vec_4,'-' , 'LineWidth',2 ,'DisplayName','$\gamma_3$')
plot(SA_vec_gamma,MZ0_vargamma_vec_5,'-' , 'LineWidth',2 ,'DisplayName','$\gamma_4$')
plot(SA_vec_gamma,MZ0_vargamma_vec_6,'-' , 'LineWidth',2 ,'DisplayName','$\gamma_5$')
xlabel('$\alpha$ [-]')
ylabel('$M_{z0}$ [Nm/S]')
title('Different camber , $F_{z}$ = 220 [N] , $kappa$ = 0')
legend 

resid_pure_MZ_vargamma = resid_pure_MZ_vargamma(P_cam,MZ_vec, 0, ALPHA_vec, GAMMA_vec , FZ0, FY_pred , tyre_coeffs)
R2 = 1 - resid_pure_MZ_vargamma;
RMSE = sqrt(resid_pure_MZ_vargamma*sum(MZ_vec.^2)/length(ALPHA_vec));
fprintf('R^2 = %6.3f \nRMSE = %6.3f \n', R2, RMSE );

