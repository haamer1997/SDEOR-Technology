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
[G,frac_ids,well_ids]=explicitStencil(physdim, 'aperture',w_f,...
     'numStages',1,'fracSpacing', 200, 'fracHalfLength'  , x_f,      ...
     'nxRefine',64,'tipNY',2,'ny',8, 'nz',8, 'gridType','geomspacing'); 
G = computeGeometry(G);
if (opt.shouldPlot)
    figure,
    plotGrid(G), view(5,45)
end
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
%% Diffusion Modeling
%MS Diffusion
G.rock.fluidname = 'barnett3comps_nat_vars'; 
G.rock.shaleMechanisms.MS_ECL_Diffusion = 1;
G.rock.tau = 1.7;
% [G.rock.Dg_MS,G.rock.Do_MS] = MSdiffNatVars(fluidname); %computes MS Diffusion Coefficients for the reservoir fluid at initial conditions
% [G.rock.Dg_MS,G.rock.Do_MS] = HassanD_i_a(fluidname,Dg_MS,Do_MS);
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

arg_Diffusion = {G, G.rock, flowfluid, fluid, 'water', false};
arg_G1cell = {G1cell, G1cell.rock, flowfluid, fluid, 'water', false};

diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);

if useNatural
    constructor = @NatVarsShaleModel;
    constructor_G1cell = @GenericNaturalVariablesModel;
else
    model = OverallCompShaleModel(G, [], flowfluid, fluid, 'water', false);
end
modelDiagonalAD_Diffusion = constructor(arg_Diffusion{:}, 'AutoDiffBackend', diagonal_backend);

ncomp = fluid.getNumberOfComponents();
s0 = [0.1, 0.9];
%% Calculate Nabla operator's denominator for MS diffusion calculations
dd = computeNabla(modelDiagonalAD_Diffusion.G, modelDiagonalAD_Diffusion.G.rock);
modelDiagonalAD_Diffusion.operators.dd = dd(modelDiagonalAD_Diffusion.operators.internalConn);
modelDiagonalAD_Diffusion.operators.g_gradz = modelDiagonalAD_Diffusion.operators.Grad(G.cells.centroids(:,3))*9.81;
%% Set up initial state
statec = initCompositionalState(G1cell, info_G1cell.pressure, info_G1cell.temp, [0.3, 0.4, 0.3], info_G1cell.initial, cmodel.EOSModel);
state_Diffusion = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelDiagonalAD_Diffusion.EOSModel);
%% Pick linear solver
linsolve = selectLinearSolverAD(modelDiagonalAD_Diffusion);
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
%% Set up Schedule
totTime = 15*year;
nSteps = 15;
dt = rampupTimesteps(totTime, 30*day, nSteps);
schedule = simpleSchedule(dt, 'W', W1);
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
%% Simulation Run
[ws_Diffusion, states_Diffusion, reports_Diffusion] = simulateScheduleAD(state_Diffusion, modelDiagonalAD_Diffusion, schedule, 'nonlinearsolver', nls,'Verbose', true);
%% plotting
if (opt.shouldPlot)
    figure, 
    plotToolbar(G, states_Diffusion)
    view(40,30);
    axis tight equal;
    plotWellSols(ws_Diffusion, cumsum(schedule.step.val))
end
%% Extracting cumulative gas production (field units)
q_gs_Diffusion = zeros(numel(ws_Diffusion),1);
for i = 1:numel(ws_Diffusion)
    q_gs_Diffusion(i)=-3.051e+6*ws_Diffusion{i}.qGs;
end
Q_gs_Diffusion = cumtrapz(tinDays,q_gs_Diffusion)/1e6; %cum gas produced in MMscf/D
%% Save variables in HPC
% if ~opt.shouldPlot
%     fpath =  '/scratch/ahass16/';
%     fullFinalOut = [fpath, 'Diffusion_1nD.mat'];
%     save(fullFinalOut,'Q_gs_diffusion','Q_gs_dispersion','Q_gs_ECL');
% end