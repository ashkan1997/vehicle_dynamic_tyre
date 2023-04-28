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
data_set = 'Hoosier_B1464run23'; % pure lateral forces
%data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

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
%[TData0, ~] = intersect_table_data( SA_0, GAMMA_0, FZ_220 );
%[TData3, ~] = intersect_table_data( SA_3neg, GAMMA_0, FZ_220 );
%[TData6, ~] = intersect_table_data( SA_6neg, GAMMA_0, FZ_220 );

% Dataset 23: -- Lateral
[TData0, ~] = intersect_table_data( GAMMA_0, FZ_220 );

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

%% PURE LATERAL SLIP FITTING with Fz=Fz_nom= 220N and camber=0  kappa = 0 VX= 10 ----------> DONE

% Fit the coeffs {pCy1, pDy1, pEy1, pEy4, pKy1, pHy1, pVy1}
FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TData0.SA));
ones_vec  = ones(size(TData0.SA));

%FY0_guess = MF96_FY0_vec(zeros_vec,TData0.SA, zeros_vec , tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
%figure()
%plot(TData0.SA,TData0.FY,'.')
%hold on
%plot(TData0.SA,FY0_guess,'-')

% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
%    [pCy1, pDy1, pEy1, pHy1,  pKy1,  pKy2, pVy1] 
P0 = [   1,    2,   -1,    1,     1,     1,    0];

% NOTE: many local minima => limits on parameters are fundamental
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [pCy1, pDy1, pEy1, pHy1,  pKy1,  pKy2, pVy1] 
lb = [];
ub = [];

ALPHA_vec = TData0.SA;
FY_vec    = TData0.FY;

% LSM_pure_Fy returns the residual, so minimize the residual varying Y. It
% is an unconstrained minimization problem 

[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fy(P,FY_vec, ALPHA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.pCy1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDy1 = P_fz_nom(2) ;  
tyre_coeffs.pEy1 = P_fz_nom(3) ;
tyre_coeffs.pHy1 = P_fz_nom(4) ;
tyre_coeffs.pKy1 = P_fz_nom(5) ; 
tyre_coeffs.pKy2 = P_fz_nom(6) ;
tyre_coeffs.pVy1 = P_fz_nom(7) ;

SA_vec = -0.3:0.001:0.3;
FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(ALPHA_vec)), SA_vec, zeros(size(ALPHA_vec)), ...
                              FZ0.*ones(size(ALPHA_vec)),tyre_coeffs);

figure('Name','Fy0(Fz0)')
plot(TData0.SA,-TData0.FY,'.')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,-FY0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}$ [N]')
legend('Raw data','Fitting',Location='best')
%title('Fit pure lateral force $F_{y0}$')

% Calculate the residuals:

FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(ALPHA_vec)), ALPHA_vec, zeros(size(ALPHA_vec)), ...
                              FZ0.*ones(size(ALPHA_vec)),tyre_coeffs);

res_FY0_nom = resid_pure_Fy(P_fz_nom,FY_vec,ALPHA_vec,0,FZ0,tyre_coeffs);


% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_FY0_nom);

[alpha__y,By,Cy,Dy,Ey,SVy,~, ~, ~] = MF96_FY0_coeffs(0, 0, 0, tyre_coeffs.FZ0, tyre_coeffs);

fprintf('By      = %6.3f\n',By);
fprintf('Cy      = %6.3f\n',Cy);
fprintf('Ey      = %6.3f\n',Ey);
fprintf('SVy     = %6.3f\n',SVy);
fprintf('alpha__y = %6.3f\n',alpha__y);
fprintf('Ky      = %6.3f\n',By*Cy*Dy/tyre_coeffs.FZ0);

%% PURE LATERAL SLIP FITTING with Fz=Fz_nom= 440N and camber=0  kappa = 0 VX= 10 ----------> DONE

tyre_coeffs = initialise_tyre_data(R0, Fz0);

% Fit the coeffs {pCy1, pDy1, pEy1, pEy4, pKy1, pHy1, pVy1}
[TData440, ~] = intersect_table_data( GAMMA_0, FZ_440 );
FZ0_440 = mean(TData440.FZ);

zeros_vec_440 = zeros(size(TData440.SA));
ones_vec_440  = ones(size(TData440.SA));

% Guess values for parameters to be optimised
%    [pCy1, pDy1, pEy1, pHy1,  pKy1,  pKy2, pVy1] 
P0 = [   1,    2,   -1,    1,     1,     1,    0];
lb = [   1, -500,    0, -500,  -500,  -500, -500];
ub = [   2,  500,    1,  500,   500,   500,  500];

ALPHA_vec_440 = TData440.SA;
FY_vec_440    = TData440.FY;

% LSM_pure_Fy returns the residual, so minimize the residual varying Y. It
% is an unconstrained minimization problem 

[P_fz_nom_440,fval,exitflag] = fmincon(@(P)resid_pure_Fy(P,FY_vec_440, ALPHA_vec_440,0,FZ0_440, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.pCy1 = P_fz_nom_440(1) ; % 1
tyre_coeffs.pDy1 = P_fz_nom_440(2) ;  
tyre_coeffs.pEy1 = P_fz_nom_440(3) ;
tyre_coeffs.pHy1 = P_fz_nom_440(4) ;
tyre_coeffs.pKy1 = P_fz_nom_440(5) ; 
tyre_coeffs.pKy2 = P_fz_nom_440(6) ;
tyre_coeffs.pVy1 = P_fz_nom_440(7) ;

SA_vec = -0.3:0.001:0.3;
FY0_fz_nom_vec_440 = MF96_FY0_vec(zeros(size(ALPHA_vec_440)), SA_vec, zeros(size(ALPHA_vec_440)), ...
                              FZ0.*ones(size(ALPHA_vec_440)),tyre_coeffs);

figure('Name','Fy0(Fz0)')
plot(TData440.SA,-TData440.FY,'.')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,-FY0_fz_nom_vec_440,'-','LineWidth',2)
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}$ [N]')
legend('Raw data','Fitting',Location='best')
title('Fit pure lateral force $F_{y0}$ @ $F_z = 440$N')

%% FIT WITH VARIABLE LOAD -----------------------------------------------------------------> CHECK

% extract data with variable load
TDataDFz_tmp = [FZ_220; FZ_440; FZ_700; FZ_900; FZ_1120];
% Remember fitting var load must use zero camber:
TDataDFz     = intersect_table_data(GAMMA_0 , TDataDFz_tmp);




% Initialise:
zeros_vec = zeros(size(TDataDFz.SA));
ones_vec  = ones(size(TDataDFz.SA));
KAPPA_vec = TDataDFz.SL;
ALPHA_vec = TDataDFz.SA;
FY_vec    = TDataDFz.FY;
FZ_vec    = TDataDFz.FZ;

% Guess for parameters to be optimised:
%               [pDy2, pEy2, pHy2, pVy2]
P0_varFZ = [   -5,    -7,    0,    0]; 
lb_varFZ = [];
ub_varFZ = [];

% Optimisation:
[P_varFZ,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varFz(P,FY_vec, ALPHA_vec,0,FZ_vec, tyre_coeffs),...
                                     P0_varFZ,[],[],[],[],lb_varFZ,ub_varFZ);

% Update tyre data with new optimal values                             
tyre_coeffs.pDy2 = P_varFZ(1) ;
tyre_coeffs.pEy2 = P_varFZ(2) ;
tyre_coeffs.pHy2 = P_varFZ(3) ;
tyre_coeffs.pVy2 = P_varFZ(4) ;

% check guess
SA_vec = -0.3:0.001:0.3;
% FY0_dfz_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
%                            TDataDFz.FZ,tyre_coeffs);
% figure()
% plot(ALPHA_vec,FY_vec,'.')
% hold on
% plot(SA_vec,FY0_dfz_vec,'.')

res_FY0_dfz_vec = resid_pure_Fy_varFz(P_varFZ,FY_vec,SA_vec,0 , FZ_vec,tyre_coeffs);

tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));

% Visulisation--- use mean(Fz to avoid oscilltion)
FY0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec5 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','Fy0(Fz0)')
plot(TDataDFz.SA,-TDataDFz.FY,'.')
hold on
plot(SA_vec,-FY0_fz_var_vec1,'-','LineWidth',2,Color = [0.8500 0.3250 0.0980])
plot(SA_vec,-FY0_fz_var_vec2,'-','LineWidth',2,Color = [0.9290 0.6940 0.1250])
plot(SA_vec,-FY0_fz_var_vec3,'-','LineWidth',2,Color = [0.4940 0.1840 0.5560])
plot(SA_vec,-FY0_fz_var_vec4,'-','LineWidth',2,Color = [0.4660 0.6740 0.1880])
plot(SA_vec,-FY0_fz_var_vec5,'-','LineWidth',2,Color = [0.6350 0.0780 0.1840])
legend('','Fit $F_z = 220$ N', 'Fit $F_z = 440$ N', 'Fit $F_z = 700$ N', 'Fit $F_z = 900$ N', 'Fit $F_z = 1120$ N' , Loc ...
    = 'best');
xlabel('$\alpha$ [-]');
ylabel('$F_{y0}$ [N]');


[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_440.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec5_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);

Calfa_vec1 = MF96_CorneringStiffnessFY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffnessFY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffnessFY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffnessFY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec5 = MF96_CorneringStiffnessFY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','C_alpha')
subplot(2,1,1)
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(mean(FZ_220.FZ) ,Calfa_vec1_0,'+','LineWidth',2)
plot(mean(FZ_440.FZ) ,Calfa_vec2_0,'+','LineWidth',2)
plot(mean(FZ_700.FZ) ,Calfa_vec3_0,'+','LineWidth',2)
plot(mean(FZ_900.FZ) ,Calfa_vec4_0,'+','LineWidth',2)
plot(mean(FZ_1120.FZ),Calfa_vec5_0,'+','LineWidth',2)
legend({'Fz$_{220}$','Fz$_{440}$','Fz$_{700}$','Fz$_{900}$','Fz$_{1120}$'})

subplot(2,1,2)
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,Calfa_vec1,'-','LineWidth',2)
plot(SA_vec,Calfa_vec2,'-','LineWidth',2)
plot(SA_vec,Calfa_vec3,'-','LineWidth',2)
plot(SA_vec,Calfa_vec4,'-','LineWidth',2)
plot(SA_vec,Calfa_vec5,'-','LineWidth',2)
legend({'Fz$_{220}$','Fz$_{440}$','Fz$_{700}$','Fz$_{900}$','Fz$_{1120}$'})

% Calculate the residuals:

res_FYdz_nom = resid_pure_Fy_varFz(P_varFZ,FY_vec,ALPHA_vec,0,FZ_vec,tyre_coeffs);

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_FY0_nom);

[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] = MF96_FY0_coeffs(0, 0, FZ_vec(4), tyre_coeffs.FZ0, tyre_coeffs);
% 
fprintf('By      = %6.3f\n',By);
fprintf('Cy      = %6.3f\n',Cy);
fprintf('mu_y    = %6.3f\n',Dy/tyre_coeffs.FZ0);
fprintf('Ey      = %6.3f\n',Ey);
fprintf('SVy     = %6.3f\n',SVy);
fprintf('alpha_y = %6.3f\n',alpha__y);
fprintf('Ky      = %6.3f\n',By*Cy*Dy/tyre_coeffs.FZ0);

%% FIT WITH VARIABLE CAMBER, Fz = Fz_nom --------------------------------------------------> CHECK

% extract data with variable load
TDataGamma = FZ_220;

% Guess values for parameters to be optimised
%    [pDy3, pEy3, pEy4, pHy3, pKy3, pVy3, pVy4]
P0 = [   1,    1,  0.5,    1,  0.2,  0.1,    1]; 
lb = [];
ub = [];

zeros_vec = zeros(size(TDataGamma.SA));
ones_vec  = ones(size(TDataGamma.SA));
KAPPA_vec = TDataGamma.SL;
ALPHA_vec = TDataGamma.SA;
GAMMA_vec = TDataGamma.IA; 
FY_vec    = TDataGamma.FY;
FZ_vec    = TDataGamma.FZ;

%figure()
%plot(ALPHA_vec,FY_vec);

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pDy3 = P_varGamma(1);
tyre_coeffs.pEy3 = P_varGamma(2);
tyre_coeffs.pEy4 = P_varGamma(3);
tyre_coeffs.pHy3 = P_varGamma(4);
tyre_coeffs.pKy3 = P_varGamma(5);
tyre_coeffs.pVy3 = P_varGamma(6);
tyre_coeffs.pVy4 = P_varGamma(7);

SA_vec = linspace(-0.3,0.3,length(GAMMA_vec));

FY0_varGamma_vec = MF96_FY0_vec(zeros_vec,SA_vec , GAMMA_vec, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

figure('Name','Fy0 vs Gamma')
plot(ALPHA_vec,-TDataGamma.FY,'.')
hold on
plot(SA_vec,-FY0_varGamma_vec,'-')
xlabel('$\alpha$ [-]')
ylabel('$F_{y0}$ [N]')

% Calculate the residuals:

res_Fy0_varGamma  = resid_pure_Fy_varGamma(P_varGamma,FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fy0_varGamma);


[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] = MF96_FY0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
% 
fprintf('By      = %6.3f\n',By);
fprintf('Cy      = %6.3f\n',Cy);
fprintf('mu_y    = %6.3f\n',Dy/tyre_coeffs.FZ0);
fprintf('Ey      = %6.3f\n',Ey);
fprintf('SVy     = %6.3f\n',SVy);
fprintf('alpha_y = %6.3f\n',alpha__y);
fprintf('Ky      = %6.3f\n',By*Cy*Dy/tyre_coeffs.FZ0);


%% COMBINED LATERAL FY - kappa = 0 --------------------------------------------------------> CHECK

% Initialise:
% zeros_vec_0 = zeros(size(TData0.SA));
% ones_vec_0  = ones(size(TData0.SA));
% FZ0_0       = mean(TData0.FZ);
% KAPPA_vec_0 = TData0.SL;
% ALPHA_vec_0 = TData0.SA;
% FY_vec_0    = TData0.FY;

TDataDFz = [FZ_220; FZ_440; FZ_700; FZ_900; FZ_1120];
zeros_vec_0 = zeros(size(TDataDFz.SA));
ones_vec_0  = ones(size(TDataDFz.SA));
FZ0_0       = mean(TDataDFz.FZ);
KAPPA_vec_0 = TDataDFz.SL;
ALPHA_vec_0 = TDataDFz.SA;
FY_vec_0    = TDataDFz.FY;
GAMMA_vec   = TDataDFz.IA;
FZ_vec      = TDataDFz.FZ;

% % PURE PARAMETERS----------------------------------------------------------
% % Guess for parameters to be optimised:
% %           [pCy1, pDy1, pEy1, pHy1, pKy1, pKy2, pVy1]
% P0_pure_0 = [   1,    2,   -1,    1,     1,     1,    0]; 
% %lb_pure_0 = [   1,  0.1,    0,  -10,     0,    -2,  -10];
% %ub_pure_0 = [   2,    4,    1,   10,   100,     2,   10];
% lb_pure_0 = [];
% ub_pure_0 = [];
% 
% 
% % Pure parameters optimisation:
% [P_opt_pure_0,fval,exitflag] = fmincon(@(P)resid_pure_Fy(P,FY_vec_0, ALPHA_vec_0,0,FZ0_0, tyre_coeffs),...
%                                      P0_pure_0,[],[],[],[],lb_pure_0,ub_pure_0);
% 
% % Update tyre data with new optimal values                             
% tyre_coeffs.pCy1 = P_opt_pure_0(1) ;
% tyre_coeffs.pDy1 = P_opt_pure_0(2) ;  
% tyre_coeffs.pEy1 = P_opt_pure_0(3) ;
% tyre_coeffs.pHy1 = P_opt_pure_0(4) ;
% tyre_coeffs.pKy1 = P_opt_pure_0(5) ; 
% tyre_coeffs.pKy2 = P_opt_pure_0(6) ;
% tyre_coeffs.pVy1 = P_opt_pure_0(7) ;

% % VARIABLE FZ: ------------------------------------------------------------
% FZ_vec = TData0.FZ;
% % Guess for parameters to be optimised:
% %                 [pDy2, pEy2, pHy2, pVy2]
% P0_pure_varFZ_0 = [   1,    2,    1,    0]; 
% lb_pure_varFZ_0 = [   0,  0.1,   -1,   -2];
% ub_pure_varFZ_0 = [   2,    4,    1,    1];
% 
% % Optimisation:
% [P_opt_pure_varFZ_0,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varFz(P,FY_vec_0, ALPHA_vec_0,0,FZ_vec, tyre_coeffs),...
%                                      P0_pure_varFZ_0,[],[],[],[],lb_pure_varFZ_0,ub_pure_varFZ_0);
% 
% % Update tyre data with new optimal values                             
% tyre_coeffs.pDy2 = P_opt_pure_varFZ_0(1) ;
% tyre_coeffs.pEy2 = P_opt_pure_varFZ_0(2) ;
% tyre_coeffs.pHy2 = P_opt_pure_varFZ_0(3) ;
% tyre_coeffs.pVy2 = P_opt_pure_varFZ_0(4) ;
% 
% % VARIABLE GAMMA: ---------------------------------------------------------
% GAMMA_vec = TData0.IA;
% % Guess for parameters to be optimised:
% %                    [pDy3, pEy3, pEy4, pHy3, pKy3, pVy3, pVy4]
% P0_pure_varGAMMA_0 =  [   1,    2,    1,    0,    1,    1,    1]; 
% lb_pure_varGAMMA_0 = [   1,  0.1,    0,    0,    0,    0,    0];
% ub_pure_varGAMMA_0 = [   2,    4,    1,    1,    0,    0,    0];
% 
% % Optimisation:
% [P_opt_pure_varGAMMA_0,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec_0, ALPHA_vec_0,GAMMA_vec,FZ0_0, tyre_coeffs),...
%                                         P0_pure_varGAMMA_0,[],[],[],[],lb_pure_varGAMMA_0,ub_pure_varGAMMA_0);
% 
% % Update tyre data with new optimal values                             
% tyre_coeffs.pDy2 = P_opt_pure_varGAMMA_0(1) ;
% tyre_coeffs.pEy2 = P_opt_pure_varGAMMA_0(2) ;  
% tyre_coeffs.pHy2 = P_opt_pure_varGAMMA_0(3) ;
% tyre_coeffs.pVy2 = P_opt_pure_varGAMMA_0(4) ;
% 
 % Calculate FY0_pure:------------------------------------------------------
 
SA_vec = -0.3:0.001:0.3; % use so dont f**k up plot
FY0_pure_0 = MF96_FY0_vec(zeros_vec_0, SA_vec, zeros_vec_0,tyre_coeffs.FZ0*ones_vec_0,tyre_coeffs);

% COMBINED:----------------------------------------------------------------

%           [ rBy1, rBy2, rBy3, rCy1, rHy1, rVy1, rVy4, rVy5, rVy6]
P0_comb = [ 1,   1,   2,   0,   1,   1,   3,   4,   1]; 
lb_comb = [];
ub_comb = [];

[P_opt_comb,fval,exitflag] = fmincon(@(P)resid_comb_Fy(P,FY_vec_0, KAPPA_vec_0,ALPHA_vec_0,GAMMA_vec,tyre_coeffs.FZ0*ones_vec_0, tyre_coeffs,FY0_pure_0),...
                                     P0_comb,[],[],[],[],lb_comb,ub_comb);

% Update coeffs:
tyre_coeffs.rBy1 = P_opt_comb(1) ;
tyre_coeffs.rBy2 = P_opt_comb(2) ;  
tyre_coeffs.rBy3 = P_opt_comb(3) ;
tyre_coeffs.rCy1 = P_opt_comb(4) ;
tyre_coeffs.rHy1 = P_opt_comb(5) ;
tyre_coeffs.rVy1 = P_opt_comb(6) ;  
tyre_coeffs.rVy4 = P_opt_comb(7) ;
tyre_coeffs.rVy5 = P_opt_comb(8) ;
tyre_coeffs.rVy6 = P_opt_comb(9) ;

% Calculate FY Combined:---------------------------------------------------
[FY_0 , Gya_0] = MF96_FY_vec(KAPPA_vec_0,SA_vec , GAMMA_vec, FZ_vec,tyre_coeffs , FY0_pure_0);

% Plot combined Fy:
figure()
plot(TData0.SA,-TData0.FY,'b.')
hold on
plot(SA_vec,-FY_0,'r-','LineWidth',2)
title('Combined Lateral Force Fy for $\kappa = 0$')
xlabel('$\alpha$ [-]')
ylabel('$F_{y}$ [N]')
legend('', '$F_{y}$ with $\kappa = 0$ deg',Location='best')
hold off

%% SELF ALLIGNING MOMENT FIT FY: ----------------------------------------------------------> BROKEN

% GUESS:-------------------------------------------------------------------
% Initialise:

[TDataMz, ~] = intersect_table_data( GAMMA_0, FZ_220 );

zeros_vec = zeros(size(TDataMz.SA));
ones_vec  = ones(size(TDataMz.SA));
FZ0       = mean(TDataMz.FZ);
KAPPA_vec = TDataMz.SL;
ALPHA_vec = TDataMz.SA;
FY_vec    = TDataMz.FY;
MZ_vec    = TDataMz.MZ;

MZR_guess = MF96_MZr_vec(zeros_vec,ALPHA_vec,zeros_vec,FZ0.*ones_vec,tyre_coeffs);
t_guess   = MF96_t_vec(zeros_vec,ALPHA_vec,zeros_vec,FZ0.*ones_vec,tyre_coeffs);
MZ0_guess = MF96_MZ0_vec(t_guess,FY_vec,MZR_guess);

figure()
plot(ALPHA_vec,TData0.MZ,'.')
hold on
plot(ALPHA_vec,MZ0_guess)
hold off

%%

% FITTING:-----------------------------------------------------------------

% Guess for parameters to be optimised:
%         [ qBz1, qBz9, qBz10, qCz1, qDz1, qDz2, qDz3, qDz4, qDz6, qEz1, qEz4, qHz1]
%P0_pure = [    1,    6,     0,    1,    0,    0,    0,    4,    0,   10,    0,    1]; 
%P0_Mz = [    1,    0,     3,    3,    0,    3,    6,    0,    3,   3,    0,    3];
P0_Mz = [    0,    0,     0,    0,    0,    0,    0,    0,    0,   0,    0,    0];
lb_Mz = [   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1];
ub_Mz = [  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100];


% Pure parameters optimisation:
[P_opt_Mz,fval,exitflag] = fmincon(@(P)resid_pure_Mz(P,MZ_vec,0, ALPHA_vec,0,FZ0, tyre_coeffs,FY_vec),...
                                     P0_Mz,[],[],[],[],lb_Mz,ub_Mz);

% Update tyre data with new optimal values                             
tyre_coeffs.qBz1  = P_opt_Mz(1) ;
tyre_coeffs.qBz9  = P_opt_Mz(2) ;  
tyre_coeffs.qBz10 = P_opt_Mz(3) ;
tyre_coeffs.qCz1  = P_opt_Mz(4) ;
tyre_coeffs.qDz1  = P_opt_Mz(5) ; 
tyre_coeffs.qDz2  = P_opt_Mz(6) ;
tyre_coeffs.qDz3  = P_opt_Mz(7) ;
tyre_coeffs.qDz4  = P_opt_Mz(8) ;
tyre_coeffs.qDz6  = P_opt_Mz(9) ;  
tyre_coeffs.qEz1  = P_opt_Mz(10) ;
tyre_coeffs.qEz4  = P_opt_Mz(11) ;
tyre_coeffs.qHz1  = P_opt_Mz(12) ; 

% Calculate MZ0_pure:------------------------------------------------------

SA_vec = -0.3:0.001:0.3; % use so dont f**k up plot
MZR0 = MF96_MZr_vec(KAPPA_vec,ALPHA_vec,zeros_vec,FZ0.*ones_vec,tyre_coeffs);
t0   = MF96_t_vec(KAPPA_vec,ALPHA_vec,zeros_vec,FZ0.*ones_vec,tyre_coeffs);
MZ0  = MF96_MZ0_vec(t0,ALPHA_vec,MZR0);

figure()
plot(ALPHA_vec,TData0.MZ,'.')
hold on
plot(ALPHA_vec,MZ0_guess)
title('Aligning moment $M_z$ with $F_z = 220$N, zero camber, long. slip $\kappa = 0$');
legend('Raw' , 'Fit', Location = 'best');
xlabel('Side slip angle $\alpha$ (deg)');
ylabel('$M_z$(Nm)');
hold off

%% PLOT Gyk -------------------------------------------------------------------------------> DONE

sa = [0,3,6,10,20];        % side slip in radians
sl = linspace(-50,50,1e4);   % longitudinal slip

Gyk_k = zeros(length(sa),length(sl)); %alpha row, k column
for i = 1:length(sa)
    for j = 1:length(sl)
        [~,Gyk_k(i,j),~] = MF96_FXFYCOMB_coeffs(sl(j), sa(i), 0, FZ0_0, tyre_coeffs);
    end
end

figure, grid on, hold on;
plot(sl,Gyk_k , LineWidth=1.5)
xlabel('longitudinal slip $k$(-)')
ylabel('$G_{yk}(-)$')
ylim('padded')
leg = cell(size(sa));
for i = 1:length(sa)
    leg{i} = ['$\alpha$ = ',num2str(sa(i)),' deg'];
end
legend(leg,Location="best")
title('Weighting function $G_{yk}$ as a function of $k$')
hold off


sa = linspace(-10,10,1e4);
sl = [0,0.1,0.2,0.5,0.8,1];
Gya = zeros(length(sl),length(sa));
for i = 1:length(sl)
    for j = 1:length(sa)
        [~,Gya(i,j),~] = MF96_FXFYCOMB_coeffs(sl(i), sa(j), 0, FZ0_0, tyre_coeffs); % k row, alpha column
    end 
end

figure, grid on, hold on;
for i = 1:length(sl)
    plot(sa,Gya(i,:), LineWidth=1.5)
end
xlabel('side slip angle $\alpha$(deg)')
ylabel('$G_{ya}(-)$')
ylim('padded')
leg = cell(size(sl));
for i = 1:length(sl)
    leg{i} = ['$k$ = ',num2str(sl(i))];
end
legend(leg,Location="best")
title('Weighting function $G_{ya}$ as a function of $\alpha$')
hold off



%% Save tyre data structure to mat file
%
save(['tyre_' data_set,'.mat'],'tyre_coeffs');


