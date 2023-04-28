%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               % FY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISE: ------------------------------------------------------------
close all;
clear all;
clc;

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
%% SELECT TYRE DATA -------------------------------------------------------
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


diameter = 18*2.56;   % Converting inches to cm
Fz0 = 220;            % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m)

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')
%% PLOT RAW DATA ----------------------------------------------------------

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


%% SELECT SPECIFIC DATA ---------------------------------------------------
% Cut crappy data and select only 80 psi data

vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table
% store raw data in table
tyre_data.SL =  SL(smpl_range);                         % Longitudinal slip
tyre_data.SA =  SA(smpl_range)*to_rad;                  % Side slip
tyre_data.FZ = -FZ(smpl_range);  % 0.453592  lb/kg      % Vertical force
tyre_data.FX =  FX(smpl_range);                         % Longitudinal force
tyre_data.FY = -FY(smpl_range);                         % Lateral force
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





%% PURE LATERAL SLIP FITTING with Fz=Fz_nom= 220N and camber=0  kappa = 0 VX= 10 ----------> DONE

% Initialise:
tyre_coeffs = initialise_tyre_data(R0, Fz0);
[TData0, ~] = intersect_table_data( GAMMA_0, FZ_220 );

FZ0 = mean(TData0.FZ);
zeros_vec = zeros(size(TData0.SA));
ones_vec  = ones(size(TData0.SA));
ALPHA_vec = TData0.SA;
FY_vec    = TData0.FY;

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


% LSM_pure_Fy returns the residual, so minimize the residual varying Y. It
% is an unconstrained minimization problem 

[P_fz_nom,~,~] = fmincon(@(P)resid_pure_Fy(P,FY_vec, ALPHA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.pCy1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDy1 = P_fz_nom(2) ;  
tyre_coeffs.pEy1 = P_fz_nom(3) ;
tyre_coeffs.pHy1 = P_fz_nom(4) ;
tyre_coeffs.pKy1 = P_fz_nom(5) ; 
tyre_coeffs.pKy2 = P_fz_nom(6) ;
tyre_coeffs.pVy1 = P_fz_nom(7) ;

SA_vec = min(ALPHA_vec):0.001:max(ALPHA_vec);

FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs);

figure('Name','Fy0(Fz0)')
plot(TData0.SA,TData0.FY,'.')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,FY0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}$ [N]')
legend('Raw data','Fitting',Location='best')
title('Fit pure lateral force $F_{y0}$')

%% VARIABLE LOAD LATERAL FY ---------------------------------------------------------------> DONE

% extract data with variable load
TDataDFz = GAMMA_0;

% Initialise:
zeros_vec = zeros(size(TDataDFz.SA));
ones_vec  = ones(size(TDataDFz.SA));
KAPPA_vec = TDataDFz.SL;
ALPHA_vec = TDataDFz.SA;
FY_vec    = TDataDFz.FY;
FZ_vec    = TDataDFz.FZ;

% Guess for parameters to be optimised:
%               [                  pDy2,                  pEy2,                    pHy2,                  pVy2]
P0_pure_varFZ = [    -0.274390567365254,    -0.278086233494965,    -0.00117998967899828,    -0.040168008694034]; 
lb_pure_varFZ = [];
ub_pure_varFZ = [];

% Optimisation:
[P_dfz,~,~] = fmincon(@(P)resid_pure_Fy_varFz(P,FY_vec, ALPHA_vec,0,FZ_vec, tyre_coeffs),...
                                     P0_pure_varFZ,[],[],[],[],lb_pure_varFZ,ub_pure_varFZ);

% Update tyre data with new optimal values                             
tyre_coeffs.pDy2 = P_dfz(1) ;
tyre_coeffs.pEy2 = P_dfz(2) ;
tyre_coeffs.pHy2 = P_dfz(3) ;
tyre_coeffs.pVy2 = P_dfz(4) ;

SA_vec = min(ALPHA_vec):1e-4:max(ALPHA_vec);
tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));

% Visulisation--- use mean(Fz to avoid oscilltion)
FY0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec5 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','Fy0(Fz0)')
plot(TDataDFz.SA,TDataDFz.FY,'.')
hold on
plot(SA_vec,FY0_fz_var_vec1,'-','LineWidth',2,Color = [0.8500 0.3250 0.0980])
plot(SA_vec,FY0_fz_var_vec2,'-','LineWidth',2,Color = [0.9290 0.6940 0.1250])
plot(SA_vec,FY0_fz_var_vec3,'-','LineWidth',2,Color = [0.4940 0.1840 0.5560])
plot(SA_vec,FY0_fz_var_vec4,'-','LineWidth',2,Color = [0.4660 0.6740 0.1880])
plot(SA_vec,FY0_fz_var_vec5,'-','LineWidth',2,Color = [0.6350 0.0780 0.1840])
legend('','Fit $F_z = 220$ N', 'Fit $F_z = 440$ N', 'Fit $F_z = 700$ N', 'Fit $F_z = 900$ N', 'Fit $F_z = 1120$ N' , Location='best');
xlabel('$\alpha$ [-]');
ylabel('$F_{y0}$ [N]');
title('Pure lateral slip at different vertical loads')


res_FY0_dfz_vec = resid_pure_Fy_varFz(P_dfz,FY_vec,SA_vec,0 , FZ_vec,tyre_coeffs);


%% CORNERING STIFFNESS --------------------------------------------------------------------> DONE

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



%% VARIABLE CAMBER LATERAL FY -------------------------------------------------------------> CHECK

% extract data with variable load
TDataGamma = FZ_220;

% Guess values for parameters to be optimised
%    [pDy3, pEy3, pEy4, pHy3, pKy3, pVy3, pVy4]
%P0 = [   5.16349361588200, 0.883205306970767, -4.20062483943437, -0.0313517996563904, 1.34159743861844, -2.92370451981899, -2.81666715709586]; 
P0 = [ 1,   1,   1,   1,   1,   1,   1  ]; 
lb = [-1e4,    0,    0, -1e4, -1e4, -1e4, -1e4];
ub = [1e4,    1,    1, 1e4, 1e4, 1e4, 1e4];

KAPPA_vec = TDataGamma.SL;
ALPHA_vec = TDataGamma.SA;
GAMMA_vec = TDataGamma.IA; 
FY_vec    = TDataGamma.FY;
FZ_vec    = TDataGamma.FZ;
SA_vec    = min(ALPHA_vec):0.001:max(ALPHA_vec);
zeros_vec = zeros(size(SA_vec));
ones_vec  = ones(size(SA_vec));
FZ0       = mean(TDataGamma.FZ);

[P_varGamma,~,~] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec, ALPHA_vec,GAMMA_vec,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pDy3 = P_varGamma(1);
tyre_coeffs.pEy3 = P_varGamma(2);
tyre_coeffs.pEy4 = P_varGamma(3);
tyre_coeffs.pHy3 = P_varGamma(4);
tyre_coeffs.pKy3 = P_varGamma(5);
tyre_coeffs.pVy3 = P_varGamma(6);
tyre_coeffs.pVy4 = P_varGamma(7);

FY0_varGamma_vec1 = MF96_FY0_vec(zeros_vec, SA_vec, mean(GAMMA_0.IA).*ones_vec, FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec2 = MF96_FY0_vec(zeros_vec, SA_vec, mean(GAMMA_1.IA).*ones_vec, FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec3 = MF96_FY0_vec(zeros_vec, SA_vec, mean(GAMMA_2.IA).*ones_vec, FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec4 = MF96_FY0_vec(zeros_vec, SA_vec, mean(GAMMA_3.IA).*ones_vec, FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec5 = MF96_FY0_vec(zeros_vec, SA_vec, mean(GAMMA_4.IA).*ones_vec, FZ0*ones_vec,tyre_coeffs);



figure('Name','Fy0 vs Gamma')
plot(ALPHA_vec,TDataGamma.FY,'.')
hold on
plot(SA_vec,FY0_varGamma_vec1,'-', LineWidth=1.5)
plot(SA_vec,FY0_varGamma_vec2,'-', LineWidth=1.5)
plot(SA_vec,FY0_varGamma_vec3,'-', LineWidth=1.5)
plot(SA_vec,FY0_varGamma_vec4,'-', LineWidth=1.5)
plot(SA_vec,FY0_varGamma_vec5,'-', LineWidth=1.5)
xlabel('$\alpha$ [rad]')
ylabel('$F_{y0}$ [N]')
legend('','$F_y$ with $\gamma=0$ deg', ...
          '$F_y$ with $\gamma=1$ deg', ...
          '$F_y$ with $\gamma=2$ deg', ...
          '$F_y$ with $\gamma=3$ deg', ...
          '$F_y$ with $\gamma=4$ deg', ...
          '$F_y$ with $\gamma=5$ deg', ...
          Location='best');
title('Pure lateral slip at different camber angles')


% Calculate the residuals with the optimal solution found above
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

%% MZ(alpha)

[TData0, ~] = intersect_table_data( GAMMA_3, FZ_220 );
% cut the tails
idx.ALPHA = -0.2 < TData0.SA & TData0.SA < 0.2;
ALPHA_vec_1 = TData0(idx.ALPHA , :);

FY_vec    = ALPHA_vec_1.FY;
MZ_vec    = ALPHA_vec_1.MZ;
ALPHA_VEC = ALPHA_vec_1.SA;
ALPHA_vec = ALPHA_vec_1.SA;
zeros_vec = zeros(size(ALPHA_vec_1.SA));
ones_vec  = ones(size(ALPHA_vec_1.SA));

% MZr_vec = MF96_MZr_vec(zeros_vec , ALPHA_vec , zeros_vec , ...
%                             FZ0*ones_vec, tyre_coeffs);
% t_vec = MF96_t_vec(zeros_vec , ALPHA_vec , zeros_vec , ...
%                             FZ0*ones_vec, tyre_coeffs);
% MZ0_vec = MF96_MZ0_vec(t_vec , FY_vec , MZr_vec);
% 
% figure
% plot( KAPPA_vec,MZ_vec , 'o')
% hold on
% plot( KAPPA_vec,MZ0_vec , '-')

% Guess values for parameters to be optimised
%       [qBz1 , qBz9 , qBz10 , qCz1 , qDz1 , qDz2 , qDz3 , qDz4 , qDz6 , qEz1 , qEz4 , qHz1] 
% P0_mz = [   0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
P0_mz = [   1, -1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0];
lb_mz = [];
ub_mz = []

SA_vec = linspace(-0.3,0.3,length(ALPHA_vec))
[P_fz_nom_mz , fval , exitflag] = fmincon(@(P)resid_pure_Mz(P,MZ_vec , 0 , ALPHA_vec , 0 , FZ0 ,tyre_coeffs , FY_vec) , ...
                                        P0_mz , [],[],[],[],lb_mz,ub_mz);

tyre_coeffs.qBz1 = P_fz_nom_mz(1) ; % 1
tyre_coeffs.qBz9 = P_fz_nom_mz(2) ;  
tyre_coeffs.qBz10 = P_fz_nom_mz(3) ;
tyre_coeffs.qCz1 = P_fz_nom_mz(4) ; 
tyre_coeffs.qDz1 = P_fz_nom_mz(5) ;
tyre_coeffs.qDz2 = P_fz_nom_mz(6) ;
tyre_coeffs.qDz3 = P_fz_nom_mz(7) ;
tyre_coeffs.qDz4 = P_fz_nom_mz(8) ;
tyre_coeffs.qDz6 = P_fz_nom_mz(9) ;
tyre_coeffs.qEz1 = P_fz_nom_mz(10) ;
tyre_coeffs.qEz4 = P_fz_nom_mz(11) ;
tyre_coeffs.qHz1 = P_fz_nom_mz(12) ;

SA_zeros = zeros(size(SA_vec));
SA_ones = ones(size(SA_vec));


MZr_vec = MF96_MZr_vec(SA_zeros , SA_vec , SA_zeros , ...
                            FZ0*SA_ones, tyre_coeffs);
t_vec = MF96_t_vec(SA_zeros , SA_vec , SA_zeros , ...
                            FZ0*SA_ones, tyre_coeffs);
MZ0_vec = MF96_MZ0_vec(t_vec , FY_vec , MZr_vec);


figure('Name','MZ0(Fz0)')
plot(ALPHA_vec_1.SA,ALPHA_vec_1.MZ,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,MZ0_vec,'-','LineWidth',2)
xlabel('$\alpha$ [deg]')
ylabel('$M_{z0}$ [Nm/s]')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             % FX & FY COMB %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                               
%% INITIALISATION ---------------------------------------------------------
clc
clearvars -except tyre_coeffs

%% SELECT TYRE DATA -------------------------------------------------------
%dataset path
data_set_path = 'TTC_dataset/';
% dataset selection and loading

% Hoostier:----------------------------------------------------------------
%data_set = 'Hoosier_B1464run23'; % pure lateral forces
data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

diameter = 18*2.56;   % Converting inches to cm
Fz0 = 220;            % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m)
to_rad = pi/180;
to_deg = 180/pi; 


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

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')
%% PLOT RAW DATA ----------------------------------------------------------

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


%% SELECT SPECIFIC DATA ---------------------------------------------------
% Cut crappy data and select only 80 psi data

vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table
% store raw data in table
tyre_data.SL =  SL(smpl_range);                         % Longitudinal slip
tyre_data.SA =  SA(smpl_range)*to_rad;                  % Side slip
tyre_data.FZ = -FZ(smpl_range);  % 0.453592  lb/kg      % Vertical force
tyre_data.FX =  FX(smpl_range);                         % Longitudinal force
tyre_data.FY =  -FY(smpl_range);                         % Lateral force
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



%% PLOT SELECTED DATA -----------------------------------------------------

[TData0, ~] = intersect_table_data( SA_0, GAMMA_0, FZ_220 );
figure('Name','Selected-data')
plot_selected_data(TData0);


%% PURE LONGITUDINAL SLIP FX, Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10 --------------> DONE

[TData0, ~] = intersect_table_data( SA_0, GAMMA_0, FZ_220 );

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
title('Fit pure longitudinal force $F_{x0}$ with $F_z=F_{znom}= 220$N, $\gamma=0$ , $\alpha=0$')

%% VARIABLE LOAD LONGITUDINAL FX ----------------------------------------------------------> DONE

% extract data with variable load
[TDataDFz, ~] = intersect_table_data( SA_0, GAMMA_0 );


KAPPA_vec = TDataDFz.SL;
FX_vec    = TDataDFz.FX;
FZ_vec    = TDataDFz.FZ;
zeros_vec = zeros(size(TDataDFz.SL));
ones_vec  = ones(size(TDataDFz.SL));

FX0_guess = MF96_FX0_vec(KAPPA_vec,zeros_vec , zeros_vec, FZ_vec, tyre_coeffs);

SL_vec = min(KAPPA_vec):1e-4:max(KAPPA_vec);
tmp_zeros = zeros(size(SL_vec));
tmp_ones = ones(size(SL_vec));

% Guess values for parameters to be optimised
%    [pDx2 pEx2 pEx3 pHx2  pKx2  pKx3  pVx2] 
P0 = [  0,   0,   0,  0,   0,   0,   0]; 
lb = [];
ub = [];

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_dfz,~,~] = fmincon(@(P)resid_pure_Fx_varFz(P,FX_vec, KAPPA_vec,0,FZ_vec, tyre_coeffs),...
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


res_FX0_dfz_vec = resid_pure_Fx_varFz(P_dfz,FX_vec,KAPPA_vec,0 , FZ_vec,tyre_coeffs);

% Visulisation--- use mean(Fz to avoid oscilltion)
FX0_fz_var_vec1 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec2 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec3 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec4 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','Fx0(Fz0)')
plot(TDataDFz.SL,TDataDFz.FX,'.')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
%plot(SL_vec,FX0_dfz_vec,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec1,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec2,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec3,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec4,'-','LineWidth',2)
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

%% CORNERING STIFFNESS

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


%% VARIABLE CAMBER LONGITUDINAL -----------------------------------------------------------> DONE

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
[P_varGamma,~,~] = fmincon(@(P)resid_pure_Fx_varGamma(P,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pDx3 = P_varGamma(1) ; % 1

FX0_varGamma_vec0 = MF96_FX0_vec(KAPPA_vec,zeros_vec , GAMMA_0.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FX0_varGamma_vec2 = MF96_FX0_vec(KAPPA_vec, zeros_vec , GAMMA_2.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FX0_varGamma_vec4 = MF96_FX0_vec(KAPPA_vec, zeros_vec , GAMMA_4.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

figure('Name','Fx0(Gamma)')
plot(KAPPA_vec,TDataGamma.FX,'.', 'DisplayName','Raw')
hold on
plot(KAPPA_vec,FX0_varGamma_vec0,'-', 'LineWidth',1.5, 'DisplayName',  'Fit $\gamma = 0$')
plot(KAPPA_vec,FX0_varGamma_vec2,'-', 'LineWidth',1.5, 'DisplayName',  'Fit $\gamma = 2$')
plot(KAPPA_vec,FX0_varGamma_vec4,'-', 'LineWidth',1.5, 'DisplayName',  'Fit $\gamma = 4$')
title('Pure Longitudinal slip at verticle load $F_{z}$ = 220[N], $\alpha = 0$ ')
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')
legend('Location','best')

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


%% COMBINED LONGITUDINAL FX ----------------------------------------------------------------> DONE

[TDataComb, ~] = intersect_table_data(GAMMA_0, FZ_220);

% Initialise:--------------------------------------------------------------
FZ0            = mean(TDataComb.FZ);
zeros_vec      = zeros(size(TDataComb.SL));
ones_vec       = ones(size(TDataComb.SL));
ALPHA_comb_vec = TDataComb.SA;
FX_comb_vec    = TDataComb.FX;
KAPPA_comb_vec = TDataComb.SL;

% Optimisation:------------------------------------------------------------

% Guess for parameters to be optimised:
%    [rBx1, rBx2, rCx1, rHx1]
P0 = [  10,    5,    1,    0];
lb = 1;
ub = 10;

% Pure parameters optimisation:

FX0_vec =  MF96_FX0_vec(KAPPA_comb_vec, zeros_vec, zeros_vec, FZ0.*ones_vec, tyre_coeffs);



[P_comb,~,~] = fmincon(@(P)resid_comb_Fx(P,FX_comb_vec, KAPPA_comb_vec,ALPHA_comb_vec,0,FZ0,tyre_coeffs,FX0_vec),...
                                     P0,[],[],[],[],lb,ub);


% Update tyre data with new optimal values                             
tyre_coeffs.rBx1 = P_comb(1) ; 
tyre_coeffs.rBx2 = P_comb(2) ;  
tyre_coeffs.rCx1 = P_comb(3) ;
tyre_coeffs.rHx1 = P_comb(4) ;

fx_SA0 = MF96_FX_vec(FX0_vec, KAPPA_comb_vec, mean(SA_0.SA).*ones_vec   , zeros_vec, FZ0.*ones_vec, tyre_coeffs);
fx_SA3 = MF96_FX_vec(FX0_vec, KAPPA_comb_vec, mean(SA_3neg.SA).*ones_vec, zeros_vec, FZ0.*ones_vec, tyre_coeffs);
fx_SA6 = MF96_FX_vec(FX0_vec, KAPPA_comb_vec, mean(SA_6neg.SA).*ones_vec, zeros_vec, FZ0.*ones_vec, tyre_coeffs);

figure, grid on, hold on;
plot(KAPPA_comb_vec*to_deg,TDataComb.FX,'.')
plot(KAPPA_comb_vec*to_deg,fx_SA0,'LineWidth',2,Color = [0.9290 0.6940 0.1250])
plot(KAPPA_comb_vec*to_deg,fx_SA3,'LineWidth',2,Color = [0.4940 0.1840 0.5560])
plot(KAPPA_comb_vec*to_deg,fx_SA6,'LineWidth',2,Color = [0.4660 0.6740 0.1880])
xlabel('$k(-)$')
ylabel('$F_x[N]$')
ylim('padded')
legend('Raw Data','Fit $\alpha$=0 deg',...
    'Fit $\alpha$=-3 deg','Fit $\alpha$=-6 deg',Location='southeast')
title('Combined Slip Longitudinal Force')


%% PLOT Gx --------------------------------------------------------------------------------> CHECK

sa = [0,3,6,10,20]; % side slip in radians
sl = linspace(-1,1,1e4);   % longitudinal slip

Gxa_k = zeros(length(sa),length(sl));
for i = 1:length(sa)
    for j = 1:length(sl)
        [Gxa_k(i,j),~,~] = MF96_FXFYCOMB_coeffs(sl(j), sa(i)*to_rad, 0, FZ0, tyre_coeffs); %alpha row, k column
    end
end

figure, grid on, hold on;
plot(sl,Gxa_k , LineWidth=1.5)
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
        Gxa_a(i,j) = MF96_FXFYCOMB_coeffs(sl(i), sa(j), 0, FZ0, tyre_coeffs); % k row, alpha column
    end 
end

figure, grid on, hold on;
for i = 1:length(sl)
    plot(sa,Gxa_a(i,:), LineWidth=1.5)
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


%% COMBINED LATERAL FY --------------------------------------------------------------------> CHECK

[TDataComb,~] = intersect_table_data( GAMMA_0, FZ_220 );

% Initialise:
FZ0       = mean(FZ_220.FZ);
KAPPA_vec = TDataComb.SL;
ALPHA_vec = TDataComb.SA;
FY_vec    = TDataComb.FY;

% Guess for parameters to be optimised:
%    [ rBy1, rBy2, rBy3, rCy1, rHy1, rVy1, rVy4, rVy5, rVy6]
P0 = [    2,    3,0.002,    2, 0.04, -0.2,    1, -0.2, -0.2];  
lb = [];
ub = [];

[P_comb,~,~] = fmincon(@(P)resid_comb_Fy(P,FY_vec, KAPPA_vec,ALPHA_vec,0,FZ0, tyre_coeffs),...
                                     P0,[],[],[],[],lb,ub);

% Update coeffs:
tyre_coeffs.rBy1 = P_comb(1) ;
tyre_coeffs.rBy2 = P_comb(2) ;  
tyre_coeffs.rBy3 = P_comb(3) ;
tyre_coeffs.rCy1 = P_comb(4) ;
tyre_coeffs.rHy1 = P_comb(5) ;
tyre_coeffs.rVy1 = P_comb(6) ;  
tyre_coeffs.rVy4 = P_comb(7) ;
tyre_coeffs.rVy5 = P_comb(8) ;
tyre_coeffs.rVy6 = P_comb(9) ;


SL_vec    = min(KAPPA_vec):1e-4:max(KAPPA_vec);
ones_vec  = ones(size(SL_vec));
zeros_vec = zeros(size(SL_vec));

FY_comb_sa0 = MF96_FY_vec(SL_vec, mean(SA_0.SA   ).*ones_vec, zeros_vec, FZ0.*ones_vec, tyre_coeffs);
FY_comb_sa3 = MF96_FY_vec(SL_vec, mean(SA_3neg.SA).*ones_vec, zeros_vec, FZ0.*ones_vec, tyre_coeffs);
FY_comb_sa6 = MF96_FY_vec(SL_vec, mean(SA_6neg.SA).*ones_vec, zeros_vec, FZ0.*ones_vec, tyre_coeffs);

figure()
plot(KAPPA_vec,FY_vec,'.')
hold on
plot(SL_vec,FY_comb_sa0, LineWidth=1.5,Color = [0.9290 0.6940 0.1250])
plot(SL_vec,FY_comb_sa3, LineWidth=1.5,Color = [0.4940 0.1840 0.5560])
plot(SL_vec,FY_comb_sa6, LineWidth=1.5,Color = [0.4660 0.6740 0.1880])
xlabel('$k[-]$ ')
ylabel('$F_y[N]$')
legend('','alpha = 0 deg','alpha = -3 deg','alpha = -6 deg',Location='best')
title('Combined lateral slip')
hold off

%%

Gyk = zeros(3,length(KAPPA_vec));
for i=1:length(SL_vec)
    [~,Gyk(1,i),~] = MF96_FXFYXCOMB_coeffs(SL_vec(i), mean(SA_0.SA   ), 0, FZ0, tyre_coeffs);
    [~,Gyk(2,i),~] = MF96_FXFYXCOMB_coeffs(SL_vec(i), mean(SA_3neg.SA), 0, FZ0, tyre_coeffs);
    [~,Gyk(3,i),~] = MF96_FXFYXCOMB_coeffs(SL_vec(i), mean(SA_6neg.SA), 0, FZ0, tyre_coeffs);
end
figure()
plot(SL_vec,Gyk)
xlabel('k')
ylabel('Gyk')
title('Gyk(k)')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             % MZ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% Save tyre data structure to mat file
%
save(['tyre_' data_set,'.mat'],'tyre_coeffs');


