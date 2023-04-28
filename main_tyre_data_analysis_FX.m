%% Initialisation
clc
clearvars 
close all   

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  16)
set(0,'DefaultLegendFontSize',16)

addpath('tyre_lib/')

to_rad = pi/180;
to_deg = 180/pi;
%% Select tyre dataset
%dataset path
data_set_path = 'TTC_dataset/';
% dataset selection and loading

% Hoostier:----------------------------------------------------------------
%data_set = 'Hoosier_B1464run23'; % pure lateral forces
data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

fprintf('Loading dataset ...')
switch data_set
    case 'Hoosier_B1464run23'
  load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
  cut_start = 27760;
  cut_end   = 54500;
    case 'Hoosier_B1464run30'
  load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  otherwise 
  error('Not found dataset: `%s`\n', data_set) ; 
end


% Goodyear:----------------------------------------------------------------
%data_set = 'Goodyear_B1464run13'; % lateral
%data_set = 'Goodyear_B1464run58'; % long

% switch data_set
%   case 'Goodyear_B1464run13'
%   load ([data_set_path, 'Goodyear_B1464run13.mat']); % pure lateral
%   cut_start = 6295;
%   cut_end   = 18818;
%   case 'Goodyear_B1464run58'
%   load ([data_set_path, 'Goodyear_B1464run58.mat']); % pure longitudinal
%   cut_start = 3008;
%   cut_end   = 18818;
%   otherwise 
%   error('Not found dataset: `%s`\n', data_set) ;  
% end

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')
%% Plot raw data

figure
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')

%figure()
%plot(SA,FY)


%% Select some specific data
% Cut crappy data and select only 80 psi data

vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table
% store raw data in table
tyre_data.SL =  SL(smpl_range);                         % Longitudinal slip
tyre_data.SA =  SA(smpl_range)*to_rad;                  % Side slip
tyre_data.FZ = -FZ(smpl_range);  % 0.453592  lb/kg      % Vertical force
tyre_data.FX =  FX(smpl_range);                         % Longitudinal force
tyre_data.FY =  FY(smpl_range);                         % Lateral force
tyre_data.MZ =  MZ(smpl_range);                         % Self aligning moment 
tyre_data.IA =  IA(smpl_range)*to_rad;                  % Camber angle
            
% Extract points at constant camber inclination angle
GAMMA_tol = 0.05*to_rad;
idx.GAMMA_0 = 0.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 0.0*to_rad+GAMMA_tol;
idx.GAMMA_1 = 1.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 1.0*to_rad+GAMMA_tol;
idx.GAMMA_2 = 2.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 2.0*to_rad+GAMMA_tol;
idx.GAMMA_3 = 3.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 3.0*to_rad+GAMMA_tol;
idx.GAMMA_4 = 4.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 4.0*to_rad+GAMMA_tol;
idx.GAMMA_5 = 5.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 5.0*to_rad+GAMMA_tol;
GAMMA_0  = tyre_data( idx.GAMMA_0, : );
GAMMA_1  = tyre_data( idx.GAMMA_1, : );
GAMMA_2  = tyre_data( idx.GAMMA_2, : );
GAMMA_3  = tyre_data( idx.GAMMA_3, : );
GAMMA_4  = tyre_data( idx.GAMMA_4, : );
GAMMA_5  = tyre_data( idx.GAMMA_5, : );

% Extract points at constant vertical load
% 200, 150, 50, 250, 100
FZ_tol = 100;
idx.FZ_220  = 220-FZ_tol < tyre_data.FZ & tyre_data.FZ < 220+FZ_tol;
idx.FZ_440  = 440-FZ_tol < tyre_data.FZ & tyre_data.FZ < 440+FZ_tol;
idx.FZ_700  = 700-FZ_tol < tyre_data.FZ & tyre_data.FZ < 700+FZ_tol;
idx.FZ_900  = 900-FZ_tol < tyre_data.FZ & tyre_data.FZ < 900+FZ_tol;
idx.FZ_1120 = 1120-FZ_tol < tyre_data.FZ & tyre_data.FZ < 1120+FZ_tol;
FZ_220  = tyre_data( idx.FZ_220, : );
FZ_440  = tyre_data( idx.FZ_440, : );
FZ_700  = tyre_data( idx.FZ_700, : );
FZ_900  = tyre_data( idx.FZ_900, : );
FZ_1120 = tyre_data( idx.FZ_1120, : );

% The slip angle is varied continuously between -4 and +12° and then
% between -12° and +4° for the pure slip case

% The slip angle is varied step wise for longitudinal slip tests
% 0° , - 3° , -6 °
SA_tol = 0.5*to_rad;
idx.SA_0    =  0-SA_tol          < tyre_data.SA & tyre_data.SA < 0+SA_tol;
idx.SA_3neg = -(3*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -3*to_rad+SA_tol;
idx.SA_6neg = -(6*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -6*to_rad+SA_tol;
SA_0     = tyre_data( idx.SA_0, : );
SA_3neg  = tyre_data( idx.SA_3neg, : );
SA_6neg  = tyre_data( idx.SA_6neg, : );

figure()
tiledlayout(3,1)

ax_list(1) = nexttile;
plot(tyre_data.IA*to_deg)
hold on
plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(2) = nexttile;
plot(tyre_data.FZ)
hold on
plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
plot(vec_samples(idx.FZ_440),FZ_440.FZ,'.');
plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')


ax_list(3) = nexttile;
plot(tyre_data.SA*to_deg)
hold on
plot(vec_samples(idx.SA_0),   SA_0.SA*to_deg,'.');
plot(vec_samples(idx.SA_3neg),SA_3neg.SA*to_deg,'.');
plot(vec_samples(idx.SA_6neg),SA_6neg.SA*to_deg,'.');
title('Slide slip')
xlabel('Samples [-]')
ylabel('[rad]')


%% Intersect tables to obtain specific sub-datasets

% Dataset 30: -- Longitudinal
[TData0, ~] = intersect_table_data( SA_0, GAMMA_0, FZ_220 );
[TData3, ~] = intersect_table_data( SA_3neg, GAMMA_0, FZ_220 );
[TData6, ~] = intersect_table_data( SA_6neg, GAMMA_0, FZ_220 );

% Dataset 23: -- Lateral
%[TData0, ~] = intersect_table_data( GAMMA_0, FZ_220 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10 );
%[TDataSub, ~] = intersectTableData( ALPHA_00, GAMMA_00, VX_10, FZ_6500 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10, FZ_6500 );

% get data for tyre deflection (radius) versus speed
%[TDataSubRho, ~] = intersectTableData( KAPPA_00, ALPHA_00, GAMMA_00, FZ_6500 );


%% Plot_selected_data

figure('Name','Selected-data')
plot_selected_data(TData0);

%% INITIALISE FITTING DATA
% initialise tyre data
Fz0 = 220;   % [N] nominal load is given
R0  = 0.327; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***

tyre_coeffs = initialise_tyre_data(R0, Fz0);

%% PURE LONGITUDINAL SLIP FITTING with Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10 -----> DONE

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TData0.SL));
ones_vec  = ones(size(TData0.SL));

FX0_guess = MF96_FX0_vec(TData0.SL,zeros_vec , zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
%figure()
%plot(TData0.SL,TData0.FX,'.')
%hold on
%plot(TData0.SL,FX0_guess,'-')

% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
%    [pCx1 pDx1 pEx1 pEx4 pHx1 pKx1  pVx1] 
P0 = [  1,   2,   1,   0,   0,    1,    0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1] 
lb = [1,   0.1,   0,   0,  -10,    0,   -10];
ub = [2,    4,   1,   1,   10,   100,  10];


KAPPA_vec = TData0.SL;
FX_vec    = TData0.FX;


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fx(P,FX_vec, KAPPA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.pCx1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDx1 = P_fz_nom(2) ;  
tyre_coeffs.pEx1 = P_fz_nom(3) ;
tyre_coeffs.pEx4 = P_fz_nom(4) ;
tyre_coeffs.pHx1 = P_fz_nom(5) ; 
tyre_coeffs.pKx1 = P_fz_nom(6) ;
tyre_coeffs.pVx1 = P_fz_nom(7) ;

SL_vec = -0.3:0.001:0.3;
FX0_fz_nom_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                              FZ0.*ones(size(SL_vec)),tyre_coeffs);

figure()
plot(TData0.SL,TData0.FX,'.')
hold on
plot(SL_vec,FX0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')
legend('Raw data','Fitting',Location='best')
title('Fit pure longitudinal force $F_{x0}$')

%% Fit LONGITUDINAL coeefficient with variable load ---------------------------------------> DONE

% extract data with variable load
[TDataDFz, ~] = intersect_table_data( SA_0, GAMMA_0 );

% long slip

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
%FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TDataDFz.SL));
ones_vec  = ones(size(TDataDFz.SL));

FX0_guess = MF96_FX0_vec(TDataDFz.SL,zeros_vec , zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
figure()
plot(TDataDFz.SL,TDataDFz.FX,'.')
hold on
plot(TDataDFz.SL,FX0_guess,'-')

% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
%    [pDx2 pEx2 pEx3 pHx2  pKx2  pKx3  pVx2] 
P0 = [  0,   0,   0,  0,   0,   0,   0]; 


% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1 
lb = [];
ub = [];


KAPPA_vec = TDataDFz.SL;
FX_vec    = TDataDFz.FX;
FZ_vec    = TDataDFz.FZ;

% check guess
SL_vec = -0.3:0.001:0.3;
FX0_dfz_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                           TDataDFz.FZ,tyre_coeffs);
% 
 figure
 plot(KAPPA_vec,FX_vec,'.')
 hold on
 plot(SL_vec,FX0_dfz_vec,'.')


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_dfz,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varFz(P,FX_vec, KAPPA_vec,0,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

disp(exitflag)
% Change tyre data with new optimal values                             
tyre_coeffs.pDx2 = P_dfz(1) ; % 1
tyre_coeffs.pEx2 = P_dfz(2) ;  
tyre_coeffs.pEx3 = P_dfz(3) ;
tyre_coeffs.pHx2 = P_dfz(4) ;
tyre_coeffs.pKx2 = P_dfz(5) ; 
tyre_coeffs.pKx3 = P_dfz(6) ;
tyre_coeffs.pVx2 = P_dfz(7) ;


res_FX0_dfz_vec = resid_pure_Fx_varFz(P_dfz,FX_vec,SL_vec,0 , FZ_vec,tyre_coeffs);

tmp_zeros = zeros(size(SL_vec));
tmp_ones = ones(size(SL_vec));

% Visulisation--- use mean(Fz to avoid oscilltion)
FX0_fz_var_vec1 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec2 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec3 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec4 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','Fx0(Fz0)')
plot(TDataDFz.SL,TDataDFz.FX,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
%plot(SL_vec,FX0_dfz_vec,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec1,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec2,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec3,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec4,'-','LineWidth',2)

xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);

Calfa_vec1 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','C_alpha')
subplot(2,1,1)
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(mean(FZ_220.FZ),Calfa_vec1_0,'+','LineWidth',2)
plot(mean(FZ_700.FZ),Calfa_vec3_0,'+','LineWidth',2)
plot(mean(FZ_900.FZ),Calfa_vec4_0,'+','LineWidth',2)
plot(mean(FZ_1120.FZ),Calfa_vec2_0,'+','LineWidth',2)
legend({'Fz_{220}','Fz_{700}','Fz_{900}','Fz_{1120}'})

subplot(2,1,2)
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SL_vec,Calfa_vec1,'-','LineWidth',2)
plot(SL_vec,Calfa_vec2,'-','LineWidth',2)
plot(SL_vec,Calfa_vec3,'-','LineWidth',2)
plot(SL_vec,Calfa_vec4,'-','LineWidth',2)
legend({'Fz_{220}','Fz_{700}','Fz_{900}','Fz_{1120}'})


%% Fit LONGITUDINAL coefficient with variable camber --------------------------------------> DONE

% extract data with variable load
[TDataGamma, ~] = intersect_table_data( SA_0, FZ_220 );

% Fit the coeffs { pDx3}

% Guess values for parameters to be optimised
P0 = [0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%lb = [0, 0,  0, 0,  0,  0,  0];
%ub = [2, 1e6,1, 1,1e1,1e2,1e2];
lb = [];
ub = [];


zeros_vec = zeros(size(TDataGamma.SL));
ones_vec  = ones(size(TDataGamma.SL));

KAPPA_vec = TDataGamma.SL;
GAMMA_vec = TDataGamma.IA; 
FX_vec    = TDataGamma.FX;
FZ_vec    = TDataGamma.FZ;

figure()
plot(KAPPA_vec,FX_vec);


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varGamma(P,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pDx3 = P_varGamma(1) ; % 1

FX0_varGamma_vec = MF96_FX0_vec(KAPPA_vec,zeros_vec , GAMMA_vec, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

figure('Name','Fx0 vs Gamma')
plot(KAPPA_vec,TDataGamma.FX,'.')
hold on
plot(KAPPA_vec,FX0_varGamma_vec,'-')
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')
% Calculate the residuals with the optimal solution found above
res_Fx0_varGamma  = resid_pure_Fx_varGamma(P_varGamma,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fx0_varGamma);


[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
% 
fprintf('Bx      = %6.3f\n',Bx);
fprintf('Cx      = %6.3f\n',Cx);
fprintf('mux     = %6.3f\n',Dx/tyre_coeffs.FZ0);
fprintf('Ex      = %6.3f\n',Ex);
fprintf('SVx     = %6.3f\n',SVx);
fprintf('kappa_x = %6.3f\n',kappa__x);
fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_coeffs.FZ0);

% % Longitudinal stiffness
% Kx_vec = zeros(size(load_vec));
% for i = 1:length(load_vec)
%   [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, load_vec(i), tyre_data);
%   Kx_vec(i) = Bx*Cx*Dx/tyre_data.Fz0;
% end
% 
% figure('Name','Kx vs Fz')
% plot(load_vec,Kx_vec,'o-')
%% COMBINED LONGITUDINAL FX - alpha = 0 ---------------------------------------------------> DONE

% Initialise:--------------------------------------------------------------
FZ0_0       = mean(TData0.FZ);
zeros_vec_0 = zeros(size(TData0.SL));
ones_vec_0  = ones(size(TData0.SL));

% Optimisation:------------------------------------------------------------

% Guess for parameters to be optimised:
%         [pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1]
P0_pure_0 = [   1,    2,    1,    0,    0,    1,    0]; 
lb_pure_0 = [   1,  0.1,    0,    0,  -10,    0,  -10];
ub_pure_0 = [   2,    4,    1,    1,   10,  100,   10];

KAPPA_vec_0 = TData0.SL;
ALPHA_vec_0 = TData0.SA;
FX_vec_0    = TData0.FX;

% Pure parameters optimisation:

[P_opt_pure_0,fval,exitflag] = fmincon(@(P)resid_pure_Fx(P,FX_vec_0, KAPPA_vec_0,0,FZ0_0, tyre_coeffs),...
                                     P0_pure_0,[],[],[],[],lb_pure_0,ub_pure_0);

% Update tyre data with new optimal values                             
tyre_coeffs.pCx1 = P_opt_pure_0(1) ;
tyre_coeffs.pDx1 = P_opt_pure_0(2) ;  
tyre_coeffs.pEx1 = P_opt_pure_0(3) ;
tyre_coeffs.pEx4 = P_opt_pure_0(4) ;
tyre_coeffs.pHx1 = P_opt_pure_0(5) ; 
tyre_coeffs.pKx1 = P_opt_pure_0(6) ;
tyre_coeffs.pVx1 = P_opt_pure_0(7) ;

% Calculate FX0_pure:------------------------------------------------------
FX0_pure_0 = MF96_FX0_vec(TData0.SL,zeros_vec_0 , zeros_vec_0, ...
                              FZ0_0.*ones_vec_0,tyre_coeffs);

% Combined parameters optimisation:----------------------------------------

%           [rBx1, rBx2, rCx1, rHx1]
P0_comb_0 = [   1,    1,    1,    1]; 
lb_comb_0 = [   0,    0,    0,    0];
ub_comb_0 = [   1,    1,    1,    1];

[P_opt_comb_0,fval,exitflag] = fmincon(@(P)resid_comb_Fx(P,FX_vec_0, KAPPA_vec_0,ALPHA_vec_0,0,FZ0_0, tyre_coeffs,FX0_pure_0),...
                                     P0_comb_0,[],[],[],[],lb_comb_0,ub_comb_0);

% Update coeffs:
tyre_coeffs.rBx1 = P_opt_comb_0(1) ;
tyre_coeffs.rBx2 = P_opt_comb_0(2) ;  
tyre_coeffs.rCx1 = P_opt_comb_0(3) ;
tyre_coeffs.rHx1 = P_opt_pure_0(4) ;

% Calculate FX Combined:---------------------------------------------------
[FX_0 , Gxa_0] = MF96_FX_vec(KAPPA_vec_0,ALPHA_vec_0 , zeros_vec_0, ...
                             FZ0_0.*ones_vec_0,tyre_coeffs , FX0_pure_0);

%% COMBINED LONGITUDINAL FX - alpha = -3 --------------------------------------------------> DONE

% Initialise:--------------------------------------------------------------
FZ0_3       = mean(TData3.FZ);
zeros_vec_3 = zeros(size(TData3.SL));
ones_vec_3  = ones(size(TData3.SL));

% Optimisation ------------------------------------------------------------

% Guess for parameters to be optimised:
%         [pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1]
P0_pure_3 = [   1,    2,    1,    0,    0,    1,    0]; 
lb_pure_3 = [   1,  0.1,    0,    0,  -10,    0,  -10];
ub_pure_3 = [   2,    4,    1,    1,   10,  100,   10];

KAPPA_vec_3 = TData3.SL;
FX_vec_3    = TData3.FX;
ALPHA_vec_3 = TData3.SA;

% Pure parameters optimisation:

[P_opt_pure_3,fval,exitflag] = fmincon(@(P)resid_pure_Fx(P,FX_vec_3, KAPPA_vec_3,0,FZ0_3, tyre_coeffs),...
                                     P0_pure_3,[],[],[],[],lb_pure_3,ub_pure_3);

% Update tyre data with new optimal values                             
tyre_coeffs.pCx1 = P_opt_pure_3(1) ;
tyre_coeffs.pDx1 = P_opt_pure_3(2) ;  
tyre_coeffs.pEx1 = P_opt_pure_3(3) ;
tyre_coeffs.pEx4 = P_opt_pure_3(4) ;
tyre_coeffs.pHx1 = P_opt_pure_3(5) ; 
tyre_coeffs.pKx1 = P_opt_pure_3(6) ;
tyre_coeffs.pVx1 = P_opt_pure_3(7) ;

% Calculate FX0_pure:------------------------------------------------------
FX0_pure_3 = MF96_FX0_vec(TData3.SL,zeros_vec_3, zeros_vec_3, ...
                              FZ0_3.*ones_vec_3,tyre_coeffs);

% Combined parameters optimisation:

%           [rBx1, rBx2, rCx1, rHx1]
P0_comb_3 = [   1,    1,    1,    1]; 
lb_comb_3 = [   0,    0,    0,    0];
ub_comb_3 = [   1,    1,    1,    1];

[P_opt_comb_3,fval,exitflag] = fmincon(@(P)resid_comb_Fx(P,FX_vec_3, KAPPA_vec_3,ALPHA_vec_3,0,FZ0_3, tyre_coeffs,FX0_pure_3),...
                                     P0_comb_3,[],[],[],[],lb_comb_3,ub_comb_3);

% Update coeffs:
tyre_coeffs.rBx1 = P_opt_comb_3(1) ;
tyre_coeffs.rBx2 = P_opt_comb_3(2) ;  
tyre_coeffs.rCx1 = P_opt_comb_3(3) ;
tyre_coeffs.rHx1 = P_opt_pure_3(4) ;

% Calculate FX Combined:
[FX_3 , Gxa_3] = MF96_FX_vec(TData3.SL,zeros_vec_3 , zeros_vec_3, ...
                              FZ0_3.*ones_vec_3,tyre_coeffs , FX0_pure_3);


%% COMBINED LONGITUDINAL FX - alpha = -6 --------------------------------------------------> DONE

% Initialise:
FZ0_6       = mean(TData6.FZ);
zeros_vec_6 = zeros(size(TData6.SL));
ones_vec_6  = ones(size(TData6.SL));

% Optimisation ----------

% Guess for parameters to be optimised:
%         [pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1]
P0_pure_6 = [   1,    2,    1,    0,    0,    1,    0]; 
lb_pure_6 = [   1,  0.1,    0,    0,  -10,    0,  -10];
ub_pure_6 = [   2,    4,    1,    1,   10,  100,   10];

KAPPA_vec_6 = TData6.SL;
FX_vec_6    = TData6.FX;
ALPHA_vec_6 = TData6.SA;

% Pure parameters optimisation:

[P_opt_pure_6,fval,exitflag] = fmincon(@(P)resid_pure_Fx(P,FX_vec_6, KAPPA_vec_6,0,FZ0_6, tyre_coeffs),...
                                     P0_pure_6,[],[],[],[],lb_pure_6,ub_pure_6);

% Update tyre data with new optimal values                             
tyre_coeffs.pCx1 = P_opt_pure_6(1) ;
tyre_coeffs.pDx1 = P_opt_pure_6(2) ;  
tyre_coeffs.pEx1 = P_opt_pure_6(3) ;
tyre_coeffs.pEx4 = P_opt_pure_6(4) ;
tyre_coeffs.pHx1 = P_opt_pure_6(5) ; 
tyre_coeffs.pKx1 = P_opt_pure_6(6) ;
tyre_coeffs.pVx1 = P_opt_pure_6(7) ;

% Calculate FX0_pure:
FX0_pure_6 = MF96_FX0_vec(TData6.SL,zeros_vec_6 , zeros_vec_6, ...
                              FZ0_6.*ones_vec_6,tyre_coeffs);

% Combined parameters optimisation:

%           [rBx1, rBx2, rCx1, rHx1]
P0_comb_6 = [   1,    1,    1,    1]; 
lb_comb_6 = [   0,    0,    0,    0];
ub_comb_6 = [   1,    1,    1,    1];

[P_opt_comb_6,fval,exitflag] = fmincon(@(P)resid_comb_Fx(P,FX_vec_6, KAPPA_vec_6,ALPHA_vec_6,0,FZ0_6, tyre_coeffs,FX0_pure_6),...
                                     P0_comb_6,[],[],[],[],lb_comb_6,ub_comb_6);

% Update coeffs:
tyre_coeffs.rBx1 = P_opt_comb_6(1) ;
tyre_coeffs.rBx2 = P_opt_comb_6(2) ;  
tyre_coeffs.rCx1 = P_opt_comb_6(3) ;
tyre_coeffs.rHx1 = P_opt_pure_6(4) ;

% Calculate FX Combined:
[FX_6 , Gxa_6] = MF96_FX_vec(TData6.SL,ALPHA_vec_6 , zeros_vec_6, ...
                              FZ0_6.*ones_vec_6,tyre_coeffs , FX0_pure_6);
%% PLOT combined FX for all alphas: -------------------------------------------------------> DONE
% Plot Fx:
figure()
plot(TData0.SL,TData0.FX,'.' , Color = [0 0.4470 0.7410])
hold on
plot(TData0.SL,FX_0,'-','LineWidth',2 , Color = [0.8500 0.3250 0.0980] )
plot(TData3.SL,TData3.FX,'.', Color = [0 0.4470 0.7410])
plot(TData3.SL,FX_3,'-','LineWidth',2 , Color = [0.9290 0.6940 0.1250])
plot(TData6.SL,TData6.FX,'.', Color = [0 0.4470 0.7410])
plot(TData6.SL,FX_6,'-','LineWidth',2 , Color = [0.4660 0.6740 0.1880])
title('Combined Longitudinal Force $F_x$ for $\alpha = \{0,-3,-6\}$')
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')
legend('', '$F_{x}$ combined with $\alpha = 0\deg$'  , ...
       '', '$F_{x}$ combined with $\alpha = -3\deg$' , ...
       '', '$F_{x}$ combined with $\alpha = -6\deg$' , ...
       Location='best')
hold off

%% PLOT Gx --------------------------------------------------------------------------------> DONE

sa = [0,3,6,10,20]; % side slip in radians
sl = linspace(-1,1,1e4);   % longitudinal slip

Gxa_k = zeros(length(sa),length(sl));
for i = 1:length(sa)
    for j = 1:length(sl)
        Gxa_k(i,j) = MF96_FXFYCOMB_coeffs(sl(j), sa(i)*to_rad, 0, FZ0_0, tyre_coeffs); %alpha row, k column
    end
end

figure, grid on, hold on;
plot(sl,Gxa_k , LineWidth=1)
xlabel('longitudinal slip $k$(-)')
ylabel('$G_{xa}(-)$')
ylim('padded')
leg = cell(size(sa));
for i = 1:length(sa)
    leg{i} = ['$\alpha$ = ',num2str(sa(i)),' deg'];
end
legend(leg,Location="best")
title('Weighting function $G_{xa}$ as a function of $k$')
hold off

sa = linspace(-20,20,1e4);
sl = [0,0.1,0.2,0.5,0.8,1];
Gxa_a = zeros(length(sl),length(sa));
for i = 1:length(sl)
    for j = 1:length(sa)
        Gxa_a(i,j) = MF96_FXFYCOMB_coeffs(sl(i), sa(j), 0, FZ0_0, tyre_coeffs); % k row, alpha column
    end 
end

figure, grid on, hold on;
for i = 1:length(sl)
    plot(sa,Gxa_a(i,:), LineWidth=1)
end
xlabel('side slip angle $\alpha$(deg)')
ylabel('$G_{xa}(-)$')
ylim('padded')
leg = cell(size(sl));
for i = 1:length(sl)
    leg{i} = ['$k$ = ',num2str(sl(i))];
end
legend(leg,Location="best")
title('Weighting function $G_{xa}$ as a function of $\alpha$')
hold off

%% SELF ALLIGNING MOMENT FITTING FX: ------------------------------------------------------> TO DO

% GUESS:-------------------------------------------------------------------
% Initialise:
zeros_vec = zeros(size(TData0.SL));
ones_vec  = ones(size(TData0.SL));
FZ0       = mean(TData0.FZ);
KAPPA_vec = TData0.SL;
ALPHA_vec = TData0.SA;
FY_vec    = TData0.FY;
MZ_vec    = TData0.MZ;

MZR_guess = MF96_MZr_vec(KAPPA_vec,ALPHA_vec,zeros_vec,FZ0.*ones_vec,tyre_coeffs);
t_guess   = MF96_t_vec(KAPPA_vec,ALPHA_vec,zeros_vec,FZ0.*ones_vec,tyre_coeffs);
MZ0_guess = MF96_MZ0_vec(t_guess,FY_vec,MZR_guess);

figure()
plot(ALPHA_vec,TData0.MZ,'.')
hold on
plot(ALPHA_vec,MZ0_guess)
hold off

% FITTING:-----------------------------------------------------------------

% Guess for parameters to be optimised:
%         [ qBz1, qBz9, qBz10, qCz1, qDz1, qDz2, qDz3, qDz4, qDz6, qEz1, qEz4, qHz1]
%P0_pure = [    1,    6,     0,    1,    0,    0,    0,    4,    0,   10,    0,    1]; 
P0_pure = [    1,    0,     3,    3,    0,    3,    6,    0,    3,   3,    0,    3];
lb_pure = [];
ub_pure = [];


% Pure parameters optimisation:
[P_opt_pure,fval,exitflag] = fmincon(@(P)resid_pure_Mz(P,MZ_vec,0, ALPHA_vec,0,FZ0, tyre_coeffs,FY_vec),...
                                     P0_pure,[],[],[],[],lb_pure,ub_pure);

% Update tyre data with new optimal values                             
tyre_coeffs.qBz1  = P_opt_pure(1) ;
tyre_coeffs.qBz9  = P_opt_pure(2) ;  
tyre_coeffs.qBz10 = P_opt_pure(3) ;
tyre_coeffs.qCz1  = P_opt_pure(4) ;
tyre_coeffs.qDz1  = P_opt_pure(5) ; 
tyre_coeffs.qDz2  = P_opt_pure(6) ;
tyre_coeffs.qDz3  = P_opt_pure(7) ;
tyre_coeffs.qDz4  = P_opt_pure(8) ;
tyre_coeffs.qDz6  = P_opt_pure(9) ;  
tyre_coeffs.qEz1  = P_opt_pure(10) ;
tyre_coeffs.qEz4  = P_opt_pure(11) ;
tyre_coeffs.qHz1  = P_opt_pure(12) ; 

% Calculate MZ0_pure:------------------------------------------------------

SA_vec = -0.3:0.001:0.3; % use so dont f**k up plot
MZR0 = MF96_MZr_vec(KAPPA_vec,ALPHA_vec,zeros_vec,FZ0.*ones_vec,tyre_coeffs);
t0   = MF96_t_vec(KAPPA_vec,ALPHA_vec,zeros_vec,FZ0.*ones_vec,tyre_coeffs);
MZ0  = MF96_MZ0_vec(t0,ALPHA_vec,MZR0);

figure()
plot(ALPHA_vec,TData0.MZ,'.')
hold on
plot(ALPHA_vec,MZ0_guess)
hold off


%% Save tyre data structure to mat file

save(['tyre_' data_set,'.mat'],'tyre_coeffs');




