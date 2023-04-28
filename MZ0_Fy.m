%% INITIALISATION ---------------------------------------------------------
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
 %% SELECT TYRE DATA -------------------------------------------------------
% %dataset path
% data_set_path = 'TTC_dataset/';
% % dataset selection and loading
% 
% % Hoostier:----------------------------------------------------------------
% data_set = 'Hoosier_B1464run23'; % pure lateral forces
% % data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined
% 
% fprintf('Loading dataset ...')
% switch data_set
%     case 'Hoosier_B1464run23'
%   load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
%   cut_start = 27760;
%   cut_end   = 54500;
%     case 'Hoosier_B1464run30'
%   load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
%   cut_start = 19028;
%   cut_end   = 37643;
%   otherwise 
%   error('Not found dataset: `%s`\n', data_set) ; 
% end
% 
% 
% % Goodyear:----------------------------------------------------------------
% %data_set = 'Goodyear_B1464run13'; % lateral
% %data_set = 'Goodyear_B1464run58'; % long
% 
% % switch data_set
% %   case 'Goodyear_B1464run13'
% %   load ([data_set_path, 'Goodyear_B1464run13.mat']); % pure lateral
% %   cut_start = 6295;
% %   cut_end   = 18818;
% %   case 'Goodyear_B1464run58'
% %   load ([data_set_path, 'Goodyear_B1464run58.mat']); % pure longitudinal
% %   cut_start = 3008;
% %   cut_end   = 18818;
% %   otherwise 
% %   error('Not found dataset: `%s`\n', data_set) ;  
% % end
% 
% % select dataset portion
% smpl_range = cut_start:cut_end;
% 
% fprintf('completed!\n')

%% GOODYEAR TYRE DATASET
%dataset path
data_set_path = 'TTC_dataset/';
% dataset selection and loading

% GOODYEAR 
data_set = 'Goodyear_B1464run13'; % pure lateral + combined
% data_set = 'Goodyear_B1464run58';  % pure longitudinal + combined


% tyre geometric data:
% Goodyear D2704 20.0x7.0-13
% 20 diameter in inches
% 7.0 section width in inches
% tread width in inches
diameter = 20*2.56; % convert inch to cm
Fz0 = 220;   % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m)




fprintf('Loading dataset ...')
switch data_set
  case 'Goodyear_B1464run13'
  load ([data_set_path, 'Goodyear_B1464run13.mat']); % pure lateral
  cut_start = 6295;
  cut_end   = 18818;
  case 'Goodyear_B1464run58'
  load ([data_set_path, 'Goodyear_B1464run58.mat']); % pure longitudinal
  cut_start = 3008;
  cut_end   = 18818;
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


%% INTERSECT TABLES -------------------------------------------------------

% % Dataset 30: -- Longitudinal
% [TData0, ~] = intersect_table_data( SA_0, GAMMA_0, FZ_220 );
% [TData3, ~] = intersect_table_data( SA_3neg, GAMMA_0, FZ_220 );
% [TData6, ~] = intersect_table_data( SA_6neg, GAMMA_0, FZ_220 );

% Dataset 23: -- Lateral
[TData0, ~] = intersect_table_data( GAMMA_0, FZ_1120 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10 );
%[TDataSub, ~] = intersectTableData( ALPHA_00, GAMMA_00, VX_10, FZ_6500 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10, FZ_6500 );

% get data for tyre deflection (radius) versus speed
%[TDataSubRho, ~] = intersectTableData( KAPPA_00, ALPHA_00, GAMMA_00, FZ_6500 );


%% PLOT SELECTED DATA -----------------------------------------------------

% figure('Name','Selected-data')
% plot_selected_data(TData0);

%% INITIALISE FITTING DATA ------------------------------------------------
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

SA_vec = linspace(-0.2,0.2,128);
FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(ALPHA_vec)), SA_vec, zeros(size(ALPHA_vec)), ...
                              FZ0.*ones(size(ALPHA_vec)),tyre_coeffs);

figure('Name','Fy0(Fz0)')
plot(TData0.SA,TData0.FY,'.')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,FY0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}$ [N]')
legend('Raw data','Fitting',Location='best')
%title('Fit pure lateral force $F_{y0}$')

%% SELF ALLIGNING MOMENT k = 0 , Fz = 220  [ashkan]

TData = intersect_table_data( GAMMA_0, FZ_220 );
idx.ALPHA = -0.2 < TData.SA & TData.SA < 0.2;
TData0 = TData(idx.ALPHA , :);
% TData0 = smoothdata(TData0);
plot_selected_data(TData0);
% Guess values for parameters to be optimised
%       [qBz1 , qBz9 , qBz10 , qCz1 , qDz1 , qDz2 , qDz3 , qDz4 , qDz6 , qEz1 , qEz4 , qHz1] 
% P0_mz = [   0.1, 0.8, 0.5, 0.2, 0.1, 0.1, 0.1, -0.1, -0.1, 0.1, 0.1, 0.3];
P0_mz = [   1, -1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0];
lb_mz = [];
ub_mz = [];

zeros_vec = zeros(size(TData0));
ones_vec = ones(size(TData0));

KAPPA_VEC = TData0.SL;
ALPHA_vec = TData0.SA;
FZ_vec = TData0.FZ;
FZ0 = mean(FZ_vec);
MZ_vec = TData0.MZ;
FY_vec = TData0.FY;
SA_vec = linspace(max(ALPHA_vec),min(ALPHA_vec),length(ALPHA_vec));

figure()
plot(ALPHA_vec , MZ_vec , '.');

[P_fz_nom_mz , fval , exitflag] = fmincon(@(P)resid_pure_Mz(P , MZ_vec , 0 , ALPHA_vec , 0 , FZ0 ,FY_vec ,  tyre_coeffs ), ...
    P0_mz , [],[],[],[],lb_mz , ub_mz);

tyre_coeffs.qBz1 =  P_fz_nom_mz(1);
tyre_coeffs.qBz9 =  P_fz_nom_mz(2);
tyre_coeffs.qBz10 =  P_fz_nom_mz(3);
tyre_coeffs.qCz1 =  P_fz_nom_mz(4);
tyre_coeffs.qDz1 =  P_fz_nom_mz(5);
tyre_coeffs.qDz2 =  P_fz_nom_mz(6);
% tyre_coeffs.qDz3 =  P_fz_nom_mz(7);
% tyre_coeffs.qDz4 =  P_fz_nom_mz(8);
tyre_coeffs.qDz6 =  P_fz_nom_mz(7);
tyre_coeffs.qEz1 =  P_fz_nom_mz(8);
tyre_coeffs.qEz4 =  P_fz_nom_mz(9);
tyre_coeffs.qHz1 =  P_fz_nom_mz(10);

MZr_vec = MF96_MZr_vec(zeros_vec , SA_vec , zeros_vec , FZ0*ones_vec , tyre_coeffs); 
t_vec = MF96_t_vec(zeros_vec , SA_vec , zeros_vec , FZ0*ones_vec , tyre_coeffs );
MZ0_nom = MF96_MZ0_vec(t_vec , FY_vec , MZr_vec);

figure('Name','MZ0(Fz0)')
plot(TData0.SA,TData0.MZ,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,MZ0_nom,'-','LineWidth',2)
xlabel('$\alpha$ [deg]')
ylabel('$M_{z0}$ [Nm/s]')

















