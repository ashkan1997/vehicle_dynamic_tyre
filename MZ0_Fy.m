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
% data_set_path = 'TTC_dataset/';
% % dataset selection and loading
% 
% % GOODYEAR 
% data_set = 'Goodyear_B1464run13'; % pure lateral + combined
% % data_set = 'Goodyear_B1464run58';  % pure longitudinal + combined
% 
% 
% % tyre geometric data:
% % Goodyear D2704 20.0x7.0-13
% % 20 diameter in inches
% % 7.0 section width in inches
% % tread width in inches
% diameter = 20*2.56; % convert inch to cm
% Fz0 = 220;   % [N] nominal load is given
% R0  = diameter/2/100; % [m] get from nominal load R0 (m)
% 
% 
% 
% 
% fprintf('Loading dataset ...')
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
%   
% end
% 
% % select dataset portion
% smpl_range = cut_start:cut_end;
% 
% fprintf('completed!\n')

%% hoosier dataset
% dataset path
data_set_path = 'TTC_dataset/';

% dataset selection and loading
data_set = 'Hoosier_B1464run23'; % pure lateral forces
% data_set = 'Hoosier_B1464run30';  % braking/traction (pure long. force) + combined

% tyre geometric data:
% Hoosier	18.0x6.0-10
% 18 diameter in inches
% 6.0 section width in inches
% tread width in inches
diameter = 18*2.56;   % Converting inches to cm
Fz0 = 220;            % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m)


fprintf('Loading dataset ...')
switch data_set
    case 'Hoosier_B1464run23'
       load ([data_set_path,data_set, '.mat']); % pure lateral
       cut_start = 32380;
       cut_end   = 54500;
    case 'Hoosier_B1464run30'
        load ([data_set_path,data_set,'.mat']); % pure longitudinal
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


% The longitudinal slip varied between  



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
[TData0, ~] = intersect_table_data( GAMMA_0, FZ_220 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10 );
%[TDataSub, ~] = intersectTableData( ALPHA_00, GAMMA_00, VX_10, FZ_6500 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10, FZ_6500 );

% get data for tyre deflection (radius) versus speed
%[TDataSubRho, ~] = intersectTableData( KAPPA_00, ALPHA_00, GAMMA_00, FZ_6500 );


%% PLOT SELECTED DATA -----------------------------------------------------

figure('Name','Selected-data')
plot_selected_data(TData0);

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

[TData , ~] = intersect_table_data( GAMMA_0, FZ_220 );
% idx.ALPHA = -0.2 < TData.SA & TData.SA < 0.2;
% TData0 = TData(idx.ALPHA , :);
% TData0 = smoothdata(TData0);
cut = 104:840;
TData0 = TData(cut , :);
plot_selected_data(TData0);


% Guess values for parameters to be optimised
%       [qBz1 , qBz9 , qBz10 , qCz1 , qDz1  , qDz6 , qEz1 , qEz4 , qHz1] 
% P0_mz = [   0.1, 0.8, 0.5, 0.2, 0.1, 0.1, 0.1, -0.1, -0.1, 0.1, 0.1, 0.3];
P0_mz = [      1, 1, 1, 1, 1, 1, 1, 1, 1];
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

[P_fz_nom_mz , ~ , ~] = fmincon(@(P)resid_pure_Mz(P , MZ_vec , 0 , ALPHA_vec , 0 , FZ0 ,FY_vec ,  tyre_coeffs ), ...
    P0_mz , [],[],[],[],lb_mz , ub_mz);

tyre_coeffs.qBz1 =  P_fz_nom_mz(1);
tyre_coeffs.qBz9 =  P_fz_nom_mz(2);
tyre_coeffs.qBz10 =  P_fz_nom_mz(3);
tyre_coeffs.qCz1 =  P_fz_nom_mz(4);
tyre_coeffs.qDz1 =  P_fz_nom_mz(5);
% tyre_coeffs.qDz2 =  P_fz_nom_mz(6);
% tyre_coeffs.qDz3 =  P_fz_nom_mz(7);
% tyre_coeffs.qDz4 =  P_fz_nom_mz(8);
tyre_coeffs.qDz6 =  P_fz_nom_mz(6);
tyre_coeffs.qEz1 =  P_fz_nom_mz(7);
tyre_coeffs.qEz4 =  P_fz_nom_mz(8);
tyre_coeffs.qHz1 =  P_fz_nom_mz(9);

% tmp_fy = MF96_FY0_vec(zeros_vec , SA_vec , zeros_vec , FZ_vec , tyre_coeffs);

MZr_vec = MF96_MZr_vec(zeros_vec , SA_vec , zeros_vec , FZ0*ones_vec , tyre_coeffs); 
t_vec = MF96_t_vec(zeros_vec , SA_vec , zeros_vec , FZ0*ones_vec , tyre_coeffs );
MZ0_nom = MF96_MZ0_vec(t_vec , FY_vec , MZr_vec);

%%
figure('Name','MZ0(Fz0)')
plot(TData0.SA,TData0.MZ,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,MZ0_nom,'-','LineWidth',2)
xlabel('$\alpha$ [deg]')
ylabel('$M_{z0}$ [Nm/s]')
title(' Aligning moment $M_{z}$ , $F_{z}$ = 220 [N] , $kappa$ = 0 [-] , $gamma$ = 0 [deg]')

%%
res_Mz0 = resid_pure_Mz(P , MZ_vec , 0 , ALPHA_vec , 0 , FZ0 ,FY_vec ,  tyre_coeffs );
R2 = 1 - res_Mz0;
RMSE = sqrt(res_Mz0*sum(MZ_vec.^2)/length(ALPHA_vec));
fprintf('R^2 = %6.3f \nRMSE = %6.3f \n', R2, RMSE );
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

SA_vec = linspace(min(ALPHA_vec) , max(ALPHA_vec) , 3525);
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




%% SELF ALLIGNING MOMENT k = 0 , variable Load

TDatadfz = GAMMA_0;
% idx.ALPHA = -0.20 < TDatadfz.SA & TDatadfz.SA < 0.20;
% TData0dfz = TDatadfz(idx.ALPHA , :);
% TData0 = smoothdata(TData0);

% figure()
% plot(TData0dfz.FZ)

cut = [95:840 , 950:1730 ,1770:2630 , 2758:3493 , 3550:4380];
TData0dfz = TDatadfz(cut , :);
idx.ALPHA = -0.20 < TData0dfz.SA & TData0dfz.SA < 0.20;
TData0dfz = TData0dfz(idx.ALPHA , :);
plot_selected_data(TData0dfz);

% Guess values for parameters to be optimised
%       [qBz2 , qBz3 , qDz7 , qEz2 , qEz3 , qHz2 , qDz2] 
P0_mz_fz = [ 1 ,   1   ,  1  ,    1 ,   1 ,  1 , 1  ];
% P0_mz_fz = [  -3,  2,  10,  -4,  0,  3,  -3 , 5 ];
lb = [];
ub = [];

KAPPA_VEC = TData0dfz.SL;  % zero
ALPHA_vec = TData0dfz.SA;
FZ_vec = TData0dfz.FZ;
% FZ0 = mean(FZ_vec);
MZ_vec = TData0dfz.MZ;
FY_vec = TData0dfz.FY;
SA_vec_mz = linspace(max(ALPHA_vec),min(ALPHA_vec),length(ALPHA_vec));
% SA_vec_mz = linspace(0.3 , -0.3 , length(TData0dfz.SA));
% SA_vec_mz = linspace(0.3 , -0.3 , 4023);
% figure()
% plot(ALPHA_vec , MZ_vec , '.');


[P_dfz,~,exitflag] = fmincon(@(P)resid_pure_MZ_varFz(P,MZ_vec, 0, ALPHA_vec, 0 , FZ_vec, FY_vec , tyre_coeffs),...
                               P0_mz_fz,[],[],[],[],lb,ub);


tyre_coeffs.qBz2 = P_dfz(1);
tyre_coeffs.qBz3 = P_dfz(2);
tyre_coeffs.qDz7 = P_dfz(3);
tyre_coeffs.qEz2 = P_dfz(4);
tyre_coeffs.qEz3 = P_dfz(5);
tyre_coeffs.qHz2 = P_dfz(6);
tyre_coeffs.qDz2 = P_dfz(7);
% tyre_coeffs.qHz4 = P(8);

tmp_zeros = zeros(size(SA_vec_mz));
tmp_ones  = ones(size(SA_vec_mz));

FY_vec_pred = MF96_FY0_vec(tmp_zeros , SA_vec_mz , tmp_zeros , FZ_vec , tyre_coeffs);

res_MZ0_fz_vec = resid_pure_MZ_varFz(P_dfz , MZ_vec , 0 , SA_vec_mz , 0 , FZ_vec , FY_vec , tyre_coeffs);


% FZ_220
cut_1 = [1381:2166];
tmp_data_1 = TData0dfz(cut_1 , :);
tmp_sa_z_1 = linspace(max(ALPHA_vec),min(ALPHA_vec) , length(tmp_data_1.SA));
tmp_zeros_1 = zeros(length(tmp_data_1.SA));
ones_vec_1 = ones(length(tmp_data_1.SA));

MZr_fz_1 = MF96_MZr_vec(tmp_zeros_1 , tmp_sa_z_1 , tmp_zeros_1 , mean(FZ_220.FZ)*ones_vec_1 , tyre_coeffs); 
t_fz_1 = MF96_t_vec(tmp_zeros_1 , tmp_sa_z_1 , tmp_zeros_1 , mean(FZ_220.FZ)*ones_vec_1 , tyre_coeffs );
FY_vec_pred = tmp_data_1.FY;
MZ0_fz_1 = MF96_MZ0_vec(t_fz_1 , FY_vec_pred , MZr_fz_1);


% FZ_440
cut_2 = [2829:3585];
tmp_data_2 = TData0dfz(cut_2 , :);
tmp_sa_z_2 = linspace(max(ALPHA_vec),min(ALPHA_vec) , length(tmp_data_2.SA));
tmp_zeros_2 = zeros(length(tmp_data_2.SA));
ones_vec_2 = ones(length(tmp_data_2.SA));

MZr_fz_2 = MF96_MZr_vec(tmp_zeros_2 , tmp_sa_z_2 , tmp_zeros_2 , mean(FZ_440.FZ)*ones_vec_2 , tyre_coeffs); 
t_fz_2 = MF96_t_vec(tmp_zeros_2 , tmp_sa_z_2 , tmp_zeros_2 , mean(FZ_440.FZ)*ones_vec_2 , tyre_coeffs );
FY_vec_pred = tmp_data_2.FY;
MZ0_fz_2 = MF96_MZ0_vec(t_fz_2 , FY_vec_pred , MZr_fz_2);

% FZ_700
cut_3 = [673:1375];
tmp_data_3 = TData0dfz(cut_3 , :);
tmp_sa_z_3 = linspace(max(ALPHA_vec),min(ALPHA_vec) , length(tmp_data_3.SA));
tmp_zeros_3 = zeros(length(tmp_data_3.SA));
ones_vec_3 = ones(length(tmp_data_3.SA));

MZr_fz_3 = MF96_MZr_vec(tmp_zeros_3 , tmp_sa_z_3 , tmp_zeros_3 , mean(FZ_700.FZ)*ones_vec_3 , tyre_coeffs); 
t_fz_3 = MF96_t_vec(tmp_zeros_3 , tmp_sa_z_3 , tmp_zeros_3 , mean(FZ_700.FZ)*ones_vec_3 , tyre_coeffs );
FY_vec_pred = tmp_data_3.FY;
MZ0_fz_3 = MF96_MZ0_vec(t_fz_3 , FY_vec_pred , MZr_fz_3);

% FZ_900
cut_4 = [1:670];
tmp_data_4 = TData0dfz(cut_4 , :);
tmp_sa_z_4 = linspace(max(ALPHA_vec),min(ALPHA_vec) , length(tmp_data_4.SA));
tmp_zeros_4 = zeros(length(tmp_data_4.SA));
ones_vec_4 = ones(length(tmp_data_4.SA));

MZr_fz_4 = MF96_MZr_vec(tmp_zeros_4 , tmp_sa_z_4 , tmp_zeros_4 , mean(FZ_700.FZ)*ones_vec_4 , tyre_coeffs); 
t_fz_4 = MF96_t_vec(tmp_zeros_4 , tmp_sa_z_4 , tmp_zeros_4 , mean(FZ_700.FZ)*ones_vec_4 , tyre_coeffs );
FY_vec_pred = tmp_data_4.FY;
MZ0_fz_4 = MF96_MZ0_vec(t_fz_4 , FY_vec_pred , MZr_fz_4);

% FZ_1120
cut_5 = [2167:2827];
tmp_data_5 = TData0dfz(cut_5 , :);
tmp_sa_z_5 = linspace(max(ALPHA_vec),min(ALPHA_vec) , length(tmp_data_5.SA));
tmp_zeros_5 = zeros(length(tmp_data_5.SA));
ones_vec_5 = ones(length(tmp_data_5.SA));

MZr_fz_5 = MF96_MZr_vec(tmp_zeros_5 , tmp_sa_z_5 , tmp_zeros_5 , mean(FZ_700.FZ)*ones_vec_5 , tyre_coeffs); 
t_fz_5 = MF96_t_vec(tmp_zeros_5 , tmp_sa_z_5 , tmp_zeros_5 , mean(FZ_700.FZ)*ones_vec_5 , tyre_coeffs );
FY_vec_pred = tmp_data_5.FY;
MZ0_fz_5 = MF96_MZ0_vec(t_fz_5 , FY_vec_pred , MZr_fz_5);


figure('Name','MZ0(Fz)')
plot(TData0dfz.SA , TData0dfz.MZ , '.')
hold on
plot(tmp_sa_z_1 , MZ0_fz_1 , '-' , 'LineWidth',1 , DisplayName='$F_{Z220}$')
plot(tmp_sa_z_2 , MZ0_fz_2 , '-' , 'LineWidth',1 , DisplayName='$F_{Z440}$')
plot(tmp_sa_z_3 , MZ0_fz_3 , '-' , 'LineWidth',1 , DisplayName='$F_{Z700}$')
plot(tmp_sa_z_4 , MZ0_fz_4 , '-' , 'LineWidth',1 , DisplayName='$F_{Z900}$')
plot(tmp_sa_z_5 , MZ0_fz_5 , '-' , 'LineWidth',1 , DisplayName='$F_{Z1120}$')

xlabel('$\alpha$ [-]')
ylabel('$M_{Z0}$ [N]')
legend

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


%% SELF ALLIGNING MOMENT k = 0 , FZ = 220 , variable camber




TData0_gamma = FZ_220;
% idx.ALPHA = -0.20 < TData_gamma.SA & TData_gamma.SA < 0.20;
% TData0_gamma = TData_gamma(idx.ALPHA , :);
plot_selected_data(TData0_gamma);


KAPPA_vec = TData0_gamma.SL;
ALPHA_vec = TData0_gamma.SA;
FY_vec    = TData0_gamma.FY;
FZ_vec    = TData0_gamma.FZ;
MZ_vec    = TData0_gamma.MZ;
GAMMA_vec = TData0_gamma.IA;

zeros_vec = zeros(size(ALPHA_vec));
ones_vec  = ones(size(ALPHA_vec));


% Guess values for parameters to be optimised
%       [qBz4 , qBz5 , qDz3 , qDz4  ,  qDz8  ,  qDz9 , qEz5 , qHz3 , qHz4] 
P0_mz_gamma = [ 1 , 1 ,  1  , 1 , 1 ,  1 , 1 , 1 , 1 ];
% P0_mz_fz = [  -3,  2,  10,  -4,  0,  3,  -3 , 5 ];
lb = [];
ub = [];


[P_cam,~,exitflag] = fmincon(@(P)resid_pure_MZ_vargamma(P,MZ_vec, 0, ALPHA_vec, GAMMA_vec , tyre_coeffs.FZ0, FY_vec , tyre_coeffs),...
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
MZr_vec_1 = MF96_MZr_vec(KAPPA_vec , ALPHA_vec , mean(GAMMA_0.IA)*ones_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
t_vec_1 = MF96_t_vec(KAPPA_vec , ALPHA_vec , mean(GAMMA_0.IA)*ones_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_1 = MF96_MZ0_vec(t_vec_1 , FY_vec , MZr_vec_1);

% GAMMA_1
MZr_vec_2 = MF96_MZr_vec(KAPPA_vec , ALPHA_vec , mean(GAMMA_1.IA)*ones_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
t_vec_2 = MF96_t_vec(KAPPA_vec , ALPHA_vec , mean(GAMMA_1.IA)*ones_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_2 = MF96_MZ0_vec(t_vec_2 , FY_vec , MZr_vec_2);

% GAMMA_2
MZr_vec_3 = MF96_MZr_vec(KAPPA_vec , ALPHA_vec , mean(GAMMA_2.IA)*ones_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
t_vec_3 = MF96_t_vec(KAPPA_vec , ALPHA_vec , mean(GAMMA_2.IA)*ones_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_3 = MF96_MZ0_vec(t_vec_3 , FY_vec , MZr_vec_3);

% GAMMA_3
MZr_vec_4 = MF96_MZr_vec(KAPPA_vec , ALPHA_vec , mean(GAMMA_3.IA)*ones_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
t_vec_4 = MF96_t_vec(KAPPA_vec , ALPHA_vec , mean(GAMMA_3.IA)*ones_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec_4 = MF96_MZ0_vec(t_vec_4 , FY_vec , MZr_vec_4);

% figure

figure('Name','MZ0 vs Gamma')
plot(ALPHA_vec,TData0_gamma.MZ,'.')
hold on
plot(ALPHA_vec,MZ0_vargamma_vec_1,'-' , 'LineWidth',1 ,'DisplayName','$\gamma_0$')
plot(ALPHA_vec,MZ0_vargamma_vec_2,'-' , 'LineWidth',1 ,'DisplayName','$\gamma_1$')
plot(ALPHA_vec,MZ0_vargamma_vec_3,'-' , 'LineWidth',1 ,'DisplayName','$\gamma_2$')
plot(ALPHA_vec,MZ0_vargamma_vec_4,'-' , 'LineWidth',1 ,'DisplayName','$\gamma_3$')
xlabel('$\alpha$ [-]')
ylabel('$M_{z0}$ [Nm/S]')
title('Different camber , $F_{z}$ = 220 [N] , $kappa$ = 0')
legend 



%%
data_gamma_0 = intersect_table_data(GAMMA_0 , FZ_220);

MZr_vec = MF96_MZr_vec(KAPPA_vec , ALPHA_vec , GAMMA_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
t_vec = MF96_t_vec(KAPPA_vec , ALPHA_vec , GAMMA_vec , tyre_coeffs.FZ0*ones_vec , tyre_coeffs);
MZ0_vargamma_vec = MF96_MZ0_vec(t_vec , FY_vec , MZr_vec);

