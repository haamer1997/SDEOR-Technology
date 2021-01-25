clear 
clc
close all
opt = struct('nkr',        2, ...
    'shouldPlot', 1 );
%     opt = merge_options(opt, varargin{:});

%% Load necessary modules, etc
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties


%% Basic Parameters
tol=1e-5;
physdim = [1990 1990 150]*ft;

wellRadius = 0.25*ft;

matrixPoro = 0.07;
matrixPerm = 500*nano*darcy;

% realFracPoro = 1-1e-8;
realw_f = 0.01*ft;
% realk_f = (0.5*2000)*darcy;
% realk_f = 0.5*darcy;
% realk_f = ((0.5*10)/28.571429)*darcy;
% realk_f = (0.5/28.571429)*darcy;
% realk_f = (0.5/(10*28.571429))*darcy;

% fracPoro = realFracPoro/wfMultiplier;
% fracPerm = realk_f/wfMultiplier;
% w_f = realw_f*wfMultiplier;    %2*centi*meter
x_f = 350*ft;
% 
% C_fD_real = (realk_f*realw_f) / (matrixPerm*x_f);
% C_fD_model = (fracPerm*w_f) / (matrixPerm*x_f);
% fprintf("Dimensionless fracture conductivity is %f \n", C_fD_real, C_fD_model)
% fprintf("modified K_f = %e;    modified frac poro = %f;  modified w_f = %f\n", fracPerm, fracPoro, w_f)

griddim = [91 19 7];
cellsize = physdim./griddim;

%% Create Box Reservoir Grid
%make ny, nz even
[G_matrix,frac_ids,well_ids]=explicitStencil(physdim, 'aperture',realw_f,...
      'numStages',1,'fracSpacing',  physdim(1), 'fracHalfLength'  , x_f,...
      'nxRefine',16,'tipNY',10,'ny',15, 'nz',2, 'gridType','cartesian'); 

G_matrix = computeGeometry(G_matrix);

if (opt.shouldPlot)
    figure,
    plotGrid(G_matrix)%, view(5,45)
end

G_matrix.rock = makeRock(G_matrix, matrixPerm, matrixPoro);

%% Create Fracture System

 %-Hydraulic Fracture Creation
fracplanes = struct;
fracplanes(1).points = [940 645 0;
                        940 1345 0;
                        940 1345 150;
                        940 645 150]*ft;

fracplanes(2).points = [995 645 0;
                        995 1345 0;
                        995 1345 150;
                        995 645 150]*ft;
                    
fracplanes(3).points = [1055 645 0;
                        1055 1345 0;
                        1055 1345 150;
                        1055 645 150]*ft;

fracplanes(1).aperture = 0.0001*ft;
fracplanes(1).poro= 1-1e-8;
fracplanes(1).perm= 0.005*nano*darcy;

fracplanes(2).aperture =  0.01*ft;
fracplanes(2).poro= 1-1e-8;
fracplanes(2).perm= 0.5*darcy;

fracplanes(3).aperture = 0.0001*ft;
fracplanes(3).poro= 1-1e-8;
fracplanes(3).perm= 0.005*nano*darcy;
% %     fracplanes(1).epsL=1.0;
% 
% fem = iscoplanar(fracplanes(1).points);
% if ~all(fem) 
%     fprintf('The frac is not coplanar');
% end
%     
% checkIfCoplanar(fracplanes)

if (opt.shouldPlot)
    plotfracongrid(G_matrix,fracplanes); % visualize to check before pre-process
end

%% Create Wells
% wells = struct;
% 
% wells(1).points=[20,70,35; 100,70,35]*meter;
% wells(1).radius=0.01*meter;

wells = struct;

%We use a short horizontal well to allow us find fracture/well intersection
halfToe2HeelDist = 50;
midpoint = physdim/2;
wells(1).points=[midpoint(1)-halfToe2HeelDist,midpoint(2),midpoint(3); midpoint(1)+halfToe2HeelDist,midpoint(2),midpoint(3)]*meter;
wells(1).radius=wellRadius;


%% visualize to check before pre-process
G=G_matrix;


%% EDFM PreProcessing

[G,fracplanes]=EDFMgrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
  %-Frac-Matrix NNCs
G = fracturematrixNNC3D(G,tol);
  %-Frac-Frac NNCs
[G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol);

%% Well fracture nnc
%-Well-Fracs NNCs
[G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol);
% Ge=G;%edfm processing end here
G = pMatFracNNCs3D(G,tol)
%% Set up EDFM operators
TPFAoperators = setupEDFMOperatorsTPFA(G, G.rock, tol);
%% fluid

fluid = initSimpleADIFluid('phases', 'W',...
                           'mu' , 1 .* centi*poise     , ...    
                           'rho', (1000) .* kilogram/meter^3, ...
                           'n'  , (2));

% Add compressibility to fluid
pRef = 5*barsa;   %bakken initial rsv pressure
c_w = 1e-5/psia;
c_o = 2.03e-5/psia;
c_g = 1e-3/psia;

fluid.bW = @(p) exp((p - pRef)*c_w);
fluid.bO = @(p) exp((p - pRef)*c_o);
fluid.bG = @(p) exp((p - pRef)*c_g);

%Pc = 170*barsa;
%alpha = 0.5;
%Pmax = 200*barsa;
%m = 0.5;
%fluid.KGangiFn = @(p) power((1-power(((Pc - alpha.*p)./Pmax),m)),3);

%% Define three-phase compressible flow model
% We define a three-phase black-oil model without dissolved gas or vaporized
% oil. This is done by first instantiating the blackoil model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid.

gravity reset off
model = WaterModel(G,[], fluid);
model.operators = TPFAoperators;

%Surface Conditions
% p_sc = 101325; %atmospheric pressure
%T_sc = 288.706;% 60 Farenheit
%[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

bc = [];
bc  = pside(bc, G, 'LEFT', 10*barsa);
bc  = pside(bc, G, 'RIGHT', 1*barsa);

% bce = [];
% bce  = pside(bce, Ge, 'LEFT', 10*barsa);
% bce  = pside(bce, Ge, 'RIGHT', 1*barsa);
% 
% bcr = [];
% bcr  = pside(bcr, Gr, 'LEFT', 10*barsa);
% bcr  = pside(bcr, Gr, 'RIGHT', 1*barsa);

%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water and oil initially. We also set up a simple-time step strategy that
% ramps up gradually towards 30 day time-steps.

s0 = (1);

state  = initResSol(G, pRef, s0);
% state_e = initResSol(Ge, pRef, s0);
% state_r = initResSol(Gr, pRef, s0);

totTime = 0.5*year;
nSteps = 50;
dt = rampupTimesteps(totTime, 30*day, nSteps);

schedule = simpleSchedule(dt, 'bc', bc);
% schedule_e = simpleSchedule(dt, 'bc', bce);
% schedule_r = simpleSchedule(dt, 'bc', bcr);

%% Simulate problem
[ws, states, reports] = simulateScheduleAD(state, model, schedule, ...
   'afterStepFn', getPlotAfterStep(state, model, schedule));

% [ws_e, states_e, reports_e] = simulateScheduleAD(state_e, model_e, schedule_e, ...
%    'afterStepFn', getPlotAfterStep(state_e, model_e, schedule_e));
% 
% [ws_r, states_r, reports_r] = simulateScheduleAD(state_r, model_r, schedule_r, ...
%    'afterStepFn', getPlotAfterStep(state_r, model_r, schedule_r));

%% plotting
figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;


%Pressure contour comparsion
[X,Y]=meshgrid(cellsize(1)/2:cellsize(1):physdim(1),cellsize(2)/2:cellsize(2):physdim(2));
nLayer=1;
NumGridLayer=griddim(1)*griddim(2);
startId=(nLayer-1)*NumGridLayer+1;
endId=(nLayer)*NumGridLayer;
% 
pres=reshape(states{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);
% 
figure, 
title('Pressure contour pEDFM')
surf(X,Y,pres,'EdgeAlpha',0.3), view(130,20)
colorbar

pres_e=reshape(states_e{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

figure, 
title('Pressure contour EDFM')
surf(X,Y,pres_e,'EdgeAlpha',0.3), view(130,20)
colorbar

% [X,Y]=meshgrid(cellsize_r(1)/2:cellsize_r(1):physdim(1),cellsize_r(2)/2:cellsize_r(2):physdim(2));
% nLayer=1;
% NumGridLayer=griddim_r(1)*griddim_r(2);
% startId=(nLayer-1)*NumGridLayer+1;
% endId=(nLayer)*NumGridLayer;
% 
% pres_r=reshape(states_r{end,1}.pressure(startId:endId),[griddim_r(1),griddim_r(2)]);
% 
% figure, 
% title('Pressure contour Explicit')
% surf(X,Y,pres_r,'EdgeAlpha',0.3), view(130,20)
% colorbar

%Pressure overline comparsion
% IJ_hFrac=ones(11,2)*6; IJ_hFrac(:,1)=1:11;
% CenterCellIDs = sub2ind(G.cartDims(1:2), IJ_hFrac(:,1), IJ_hFrac(:,2));
% 
% IJ_hFrac_r=ones(225,2)*113; IJ_hFrac_r(:,1)=1:225;
% CenterCellIDs_r = sub2ind(Gr.cartDims(1:2), IJ_hFrac_r(:,1), IJ_hFrac_r(:,2));
% 
% x_line=IJ_hFrac(:,1)*cellsize(1)-cellsize(1)/2;
% 
% pres_line=states{end,1}.pressure(CenterCellIDs);
% pres_line_e=states_e{end,1}.pressure(CenterCellIDs);
% 
% 
% x_line_r=IJ_hFrac_r(:,1)*cellsize_r(1)-cellsize_r(1)/2;
% pres_line_r=states_r{end,1}.pressure(CenterCellIDs_r);
% 
% 
% figure('rend','painters','pos',[10 10 800 600]);
% plot(x_line, pres_line,'r-', 'LineWidth', 2,'DisplayName','pEDFM');
% hold on;
% plot(x_line, pres_line_e,'b--', 'LineWidth', 2,'DisplayName','EDFM');
% plot(x_line_r, pres_line_r,'k--', 'LineWidth', 2,'DisplayName','Explicit');
% hold off;
% 
% xlim([0 9]);
% ylim([1*barsa 10*barsa]);
% 
% set(gca,'FontSize',25);
% xlabel('X [m]');
% ylabel('Pressure [Pa]');
% legend;
