clear 
clc
close all
opt = struct('nkr',        1, ...
    'shouldPlot', 0 );
%% Load necessary modules, etc
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core;
mrstModule add ad-props ad-core
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;
%% Basic Parameters
physdim = [200 150 120]*meter;
wellRadius=0.1*meter;
%HF Parameters
w_f = 0.0030;
x_f = 74.99;
%% Create Box Reservoir Grid
[G,frac_ids,well_ids]=explicitStencil(physdim, 'aperture',w_f,...
     'numStages',1,'fracSpacing', 200, 'fracHalfLength'  , x_f,      ...
     'nxRefine',512,'tipNY',2,'ny',8, 'nz',8, 'gridType','geomspacing'); 
G = computeGeometry(G);
if (opt.shouldPlot)
    figure,
    plotGrid(G), view(5,45)
end
%Rock Properties
matrix_perm = 1*nano*darcy; %100*micro
matrix_poro= 0.09;
%HF Properties
frac_poro = 0.5;
frac_perm = 1*milli*darcy;
G.rock=makeShaleRock(G, matrix_perm, matrix_poro);
%% Create Fracture System
numHFplanes = size(frac_ids,1);
for i=1:numHFplanes
     frac_id=frac_ids(i,:);
     G.rock.poro(frac_id)=frac_poro;
     G.rock.perm(frac_id)=frac_perm;     
end
%% Create Wells
wells = struct;
wells(1).XFracCellIDs=well_ids;
wells(1).radius=wellRadius;
%% Plot well and permeability
% Since the well is inside the reservoir, we remove a section around the
% well so that we can see the well path

% figure
% %Plot matrix
% show = true([G.cells.num, 1]);
% show(frac_ids(:)) = false;% Hide well cell
% plotCellData (G , convertTo(G.rock.perm,milli*darcy),show, ...
%     'EdgeColor', 'k','facealpha',0.0);
% colorbar ('horiz'); view(5,45); axis equal tight;
% hold on;
% 
% %Plot frac plane
% G.rock.perm(well_ids)=2.0*darcy; %Highlight the well cell
% show = false([G.cells.num, 1]);
% show(frac_ids(:)) = true;% Hide well cell
% plotCellData (G , convertTo(G.rock.perm,milli*darcy),show, ...
%     'EdgeColor', 'k');
% G.rock.perm(well_ids)=1.0*darcy;
% hold off;
%%
useNatural =  true;
% Name of problem and pressure range
casename = 'barnett3comps';
pwf = 500*psia;
pinj = 8000*psia;
G2 = G;
G3 = G;
G4 = G;
%Diffusion
G2.rock.shaleMechanisms.diffusion = 1;
G2.rock.Di=[2.8,2.5,1.9]*10^-7;
G2.rock.Di_o=[0,0,0];%G.rock.Di/100; %[0,0,0];
G2.rock.tau = 1.7;
%Dispersion
G3.rock.shaleMechanisms.dispersion = 1;
G3.rock.Di=[2.8,2.5,1.9]*10^-7;
G3.rock.Di_o=[0,0,0];%G.rock.Di/100; %[0,0,0];
G3.rock.tau = 1.7;
%Eclipse Diffusion
fluidname = 'barnett3comps_modified_2'; 
G4.rock.shaleMechanisms.ECLdiffusion = 1;
G4.rock.Di=[2.8,2.5,1.9]*10^-7;
G4.rock.Di_o=[0,0,0];
G4.rock.tau = 1.7;
[G4.rock.Di,G4.rock.Di_o] = HassanD_i_a(fluidname,G4.rock.Di,G4.rock.Di_o);
%% Define three-phase compositional flow model
[fluid, info] = getShaleCompFluidCase(casename);
eosname = 'prcorr';  %'srk','rk','prcorr'
[fluidMixture, info_G1cell] = getShaleCompFluidCase('barnett3comps_modified');

%Unit Cell Approach
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

arg_no_diffusion = {G, G.rock, flowfluid, fluid, 'water', false};
arg_diffusion = {G2, G2.rock, flowfluid, fluid, 'water', false};
arg_dispersion = {G3, G3.rock, flowfluid, fluid, 'water', false};
arg_ECL = {G4, G4.rock, flowfluid, fluid, 'water', false};
arg_G1cell = {G1cell, G1cell.rock, flowfluid, fluid, 'water', false};

diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);

if useNatural
    constructor = @NatVarsShaleModel;
    constructor_G1cell = @GenericNaturalVariablesModel;
else
    model = OverallCompShaleModel(G, [], flowfluid, fluid, 'water', false);
end
modelDiagonalAD_no_diffusion= constructor(arg_no_diffusion{:}, 'AutoDiffBackend', diagonal_backend);
modelDiagonalAD_diffusion = constructor(arg_diffusion{:}, 'AutoDiffBackend', diagonal_backend);
modelDiagonalAD_dispersion = constructor(arg_dispersion{:}, 'AutoDiffBackend', diagonal_backend);
modelDiagonalAD_ECL = constructor(arg_ECL{:}, 'AutoDiffBackend', diagonal_backend);

% modelSparseAD_G1cell = constructor_G1cell(arg_G1cell{:}, 'AutoDiffBackend', sparse_backend);
% modelDiagonalAD_G1cell = constructor_G1cell(arg_G1cell{:}, 'AutoDiffBackend', diagonal_backend);
% modelMexDiagonalAD_G1cell = constructor_G1cell(arg_G1cell{:}, 'AutoDiffBackend', mex_backend);

ncomp = fluid.getNumberOfComponents();
s0 = [0.1, 0.9];
%% Calculate Nabla operator's denominator for ECL diffusion calculations
dd = computeNabla(modelDiagonalAD_ECL.G, modelDiagonalAD_ECL.G.rock);
modelDiagonalAD_ECL.operators.dd = dd(modelDiagonalAD_ECL.operators.internalConn);
modelDiagonalAD_ECL.operators.g_gradz = modelDiagonalAD_ECL.operators.Grad(G.cells.centroids(:,3))*9.81;
%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water, oil and gas initially. We also set up a simple-time step strategy 
% that ramps up gradually towards 30 day time-steps.
statec = initCompositionalState(G1cell, info_G1cell.pressure, info_G1cell.temp, [0.3, 0.4, 0.3], info_G1cell.initial, cmodel.EOSModel);
state_no_diffusion = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelDiagonalAD_no_diffusion.EOSModel);
state_diffusion = initCompositionalState(G2, info.pressure, info.temp, s0, info.initial, modelDiagonalAD_diffusion.EOSModel);
state_dispersion = initCompositionalState(G3, info.pressure, info.temp, s0, info.initial, modelDiagonalAD_dispersion.EOSModel);
state_ECL = initCompositionalState(G4, info.pressure, info.temp, s0, info.initial, modelDiagonalAD_ECL.EOSModel);
%% Pick linear solver
% The AMGCL library is one possible solver option for MRST. It is fairly
% easy to write interfaces to other solvers using MEX files and/or the
% LinearSolverAD class. We call the routine for automatically selecting a
% reasonable linear solver for the specific model.
linsolve = selectLinearSolverAD(modelDiagonalAD_diffusion);
disp(linsolve)
nls = NonLinearSolver('LinearSolver', linsolve);
%% Add Well
W1 = [];
for wi=1:numel(wells)
    W1 = addWell(W1, G, G.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir','x');
end
for wi=1:numel(wells)
    W1(wi).components = info.initial;
end

totTime = 15*year;
nSteps = 15;
dt = rampupTimesteps(totTime, 30*day, nSteps);
schedule = simpleSchedule(dt, 'W', W1);
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
%% Simulation
% [ws_no_diffusion, states_no_diffusion, reports_no_diffusion] = simulateScheduleAD(state_no_diffusion, modelDiagonalAD_diffusion, schedule, 'nonlinearsolver', nls,'Verbose', true);
[ws_diffusion, states_diffusion, reports_diffusion] = simulateScheduleAD(state_diffusion, modelDiagonalAD_diffusion, schedule, 'nonlinearsolver', nls,'Verbose', true);
[ws_dispersion, states_dispersion, reports_dispersion] = simulateScheduleAD(state_dispersion, modelDiagonalAD_dispersion, schedule, 'nonlinearsolver', nls,'Verbose', true);
[ws_ECL, states_ECL, reports_ECL] = simulateScheduleAD(state_dispersion, modelDiagonalAD_ECL, schedule, 'nonlinearsolver', nls,'Verbose', true);
%% plotting
if (opt.shouldPlot)
    ws = {ws_no_diffusion,ws_diffusion, ws_dispersion,ws_ECL};
    names = {'no_diffusion','Difussion', 'Dispersion','ECL_Diffusion'};
    plotWellSols(ws, cumsum(schedule.step.val), 'datasetnames', names,'field','qGs','linestyles',{'r','b','g','y'})
    figure, 
    plotToolbar(G, states1)
    view(40,30);
    axis tight equal;
    figure,
    plotWellSols(ws1,cumsum(schedule1.step.val))
end
%% Extracting cumulative gas production (field units)
q_gs_diffusion = zeros(numel(ws_diffusion),1);
q_gs_dispersion = zeros(numel(ws_dispersion),1);
q_gs_ECL = zeros(numel(ws_ECL),1);
for i = 1:numel(ws_diffusion)
    q_gs_diffusion(i)=-3.051e+6*ws_diffusion{i}.qGs;
    q_gs_dispersion(i)=-3.051e+6*ws_dispersion{i}.qGs;
    q_gs_ECL(i)=-3.051e+6*ws_ECL{i}.qGs;
end
Q_gs_diffusion = cumtrapz(tinDays,q_gs_diffusion)/1e6; %cum gas produced in MMscf/D
Q_gs_dispersion = cumtrapz(tinDays,q_gs_dispersion)/1e6; %cum gas produced in MMscf/D
Q_gs_ECL = cumtrapz(tinDays,q_gs_ECL)/1e6; %cum gas produced in MMscf/D
%% Save variables in HPC
if ~opt.shouldPlot
    fpath =  '/scratch/ahass16/';
    fullFinalOut = [fpath, 'Diffusion_x512.mat'];
    save(fullFinalOut,'Q_gs_diffusion','Q_gs_dispersion','Q_gs_ECL');
end