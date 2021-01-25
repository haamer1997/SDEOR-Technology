%% Clearing MATLAB File
clear 
clc
close all
opt = struct('nkr',        1, ...
    'shouldPlot', 0 );
%% Load necessary modules, etc
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core;
mrstModule add ad-props 
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;
%% Basic Parameters
physdim = [200 150 120]*meter;
wellRadius=0.1*meter;
%Rock Properties
matrix_perm = 1*nano*darcy; %100*micro
matrix_poro= 0.09;
%HF Parameters
w_f = 0.0030;
x_f = 74.99;
frac_poro = 0.5;
frac_perm = 1*milli*darcy;
%% Create Box Reservoir Grid
load('G.mat');G_matrix = G;
if (opt.shouldPlot)
    figure,
    plotGrid(G_matrix), view(5,45)
end
G_matrix.rock=makeShaleRock(G_matrix, matrix_perm, matrix_poro);
%% Create Fracture System
[fracplanes]=createMultiStageHFs('numStages',1,'fracSpacing', 110,...
    'numFracsPerStage', 1,'fracHalfLength', 75,'fracHeight',120,...
    'clusterSpacing', 0,'heelCoord', [100.0 75, 60.0]);

for i=1:numel(fracplanes)
     fracplanes(i).aperture = w_f;
     fracplanes(i).poro = frac_poro;
     fracplanes(i).perm = frac_perm;
end

if (opt.shouldPlot)
    figure,
    plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
end
%% Create Wells
wells = struct;
wells(1).points=[50 75, 60.0; 150,75,60]*meter;
wells(1).radius=0.1*meter;
%% EDFM PreProcessing
G = G_matrix;
tol=1e-5;
[G,fracplanes]=EDFMshalegrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
  %-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
  %-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
  %-Well-Fracs NNCs
[G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol);
if (opt.shouldPlot)
    figure,
    plotFracSystem(G,fracplanes,wells,'label',false)
end
TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);
%% Diffusion Modeling
%Fick's Diffusion
% G.rock.shaleMechanisms.diffusion = 1;
% G.rock.Di=[2.8,2.5,1.9]*10^-7;
% G.rock.Di_o=[0,0,0];%G.rock.Di/100; %[0,0,0];
% G.rock.tau = 1.7;

G.rock.shaleMechanisms.dispersion = 1;
G.rock.Di=[2.8,2.5,1.9]*10^-7;
G.rock.Di_o=[0,0,0];%G.rock.Di/100; %[0,0,0];
G.rock.tau = 1.7;

%% Define three-phase compositional flow model
useNatural =  true;
casename = 'barnett3comps';
pwf = 500*psia;
[fluid, info] = getShaleCompFluidCase(casename);
eosname = 'prcorr';  

%Unit Cell Approach
[fluidMixture, info_G1cell] = getShaleCompFluidCase('barnett3comps_modified');
G1cell = cartGrid([1, 1, 1]);
G1cell = computeGeometry(G1cell);
cmodel = GenericOverallCompositionModel(G1cell, makeRock(G1cell, 1, 1), initSimpleADIFluid(), fluidMixture);
cmodel = cmodel.validateModel(); % Set up groups
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);
G1cell.rock=makeShaleRock(G1cell, matrix_perm, matrix_poro);

%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]); 

gravity reset on

arg_MS_Diffusion = {G, [], flowfluid, fluid, 'water', false};
arg_G1cell = {G1cell, G1cell.rock, flowfluid, fluid, 'water', false};

diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);

if useNatural
    constructor = @NatVarsShaleModel;
    constructor_G1cell = @GenericNaturalVariablesModel;
else
    model = OverallCompShaleModel(G_matrix, [], flowfluid, fluid, 'water', false);
end
modelDiagonalAD_MS_Diffusion = constructor(arg_MS_Diffusion{:}, 'AutoDiffBackend', diagonal_backend);
modelDiagonalAD_MS_Diffusion.operators = TPFAoperators;
ncomp = fluid.getNumberOfComponents();
s0 = [0.1, 0.9];
%% Set up initial state
statec = initCompositionalState(G1cell, info_G1cell.pressure, info_G1cell.temp, [0.3, 0.4, 0.3], info_G1cell.initial, cmodel.EOSModel);
state_MS_Diffusion = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelDiagonalAD_MS_Diffusion.EOSModel);
%% Pick linear solver
linsolve = selectLinearSolverAD(modelDiagonalAD_MS_Diffusion);
disp(linsolve)
nls = NonLinearSolver('LinearSolver', linsolve);
%% Add Well
W1 = [];
for wi=1:numel(wells)
    W1 = addWellEDFMshale(W1, G, G.Matrix.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir','x');
end
for wi=1:numel(wells)
    W1(wi).components = info.initial;
end
%% Set up Schedule
totTime = 15*year;
nSteps = 15;
dt = rampupTimesteps(totTime, 30*day, nSteps);
schedule = simpleSchedule(dt, 'W', W1);
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
%% Simulation Run
[ws_MS_Diffusion, states_MS_Diffusion, reports_MS_Diffusion] = simulateScheduleAD(state_MS_Diffusion, modelDiagonalAD_MS_Diffusion, schedule, 'nonlinearsolver', nls,'Verbose', true);
%% plotting
if (opt.shouldPlot)
    figure, 
    plotToolbar(G_matrix, states_MS_Diffusion)
    view(40,30);
    axis tight equal;
    plotWellSols(ws_MS_Diffusion, cumsum(schedule.step.val))
end
%% Extracting cumulative gas production (field units)
q_gs_MS_Diffusion = zeros(numel(ws_MS_Diffusion),1);
for i = 1:numel(ws_MS_Diffusion)
    q_gs_MS_Diffusion(i)=-3.051e+6*ws_MS_Diffusion{i}.qGs;
end
Q_gs_MS_Diffusion = cumtrapz(tinDays,q_gs_MS_Diffusion)/1e6; %cum gas produced in MMscf/D
%% Save variables in HPC
% if ~opt.shouldPlot
%     fpath =  '/scratch/ahass16/';
%     fullFinalOut = [fpath, 'Diffusion_1nD.mat'];
%     save(fullFinalOut,'Q_gs_diffusion','Q_gs_dispersion','Q_gs_ECL');
% end