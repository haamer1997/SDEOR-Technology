%% File Header
clear 
clc
close all
Globals
tol=1e-5;
opt = struct('nkr',        1, ...
             'shouldPlot', 0);
%% Load necessary modules, etc
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-props ad-core
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;
%% Basic Parameters
physdim = [200 150 120]*meter;
wellRadius=0.1*meter;
%HF Parameters
w_f = 0.1*ft;
x_f = 74.99;
%% Create Box Reservoir Grid
[G_matrix,frac_ids,well_ids]=explicitStencil(physdim, 'aperture',w_f,...
                             'numStages',1,'fracSpacing', 200, 'fracHalfLength'  , x_f,      ...
                             'nxRefine',501,'tipNY',2,'ny',10, 'nz',10, 'gridType','geomspacing'); 
G_matrix = computeGeometry(G_matrix);
if (opt.shouldPlot)
    figure,
    plotGrid(G_matrix), view(5,45)
end
%Rock Parameters
matrix_perm = 10*nano*darcy;
matrix_poro= 0.09;

%HF Parameters
frac_poro = 0.5;
frac_perm = 1*darcy;
G_matrix.rock=makeShaleRock(G_matrix, matrix_perm, matrix_poro);
%% Compute MS diffusion coefficients to populate into matrix and fracture cells
% Molecular Diffusion Modeling
fluidname = 'barnett3comps_modified'; %this is a modified fluid such that a dummy component is added to the end (NatVars formulation)
G_matrix.rock.shaleMechanisms.Diffusion = true; 
G_matrix.rock.shaleMechanisms.MS_Diffusion = 1;
[Dg_MS,Do_MS] = MSdiffNatVars(fluidname); %MS static implementation
Dg_MS = reshape(Dg_MS,[1 numel(Dg_MS)]); Do_MS = reshape(Do_MS,[1 numel(Do_MS)]);
G_matrix.rock.Dg = repmat(Dg_MS,G_matrix.cells.num,1); G_matrix.rock.Do = repmat(Do_MS,G_matrix.cells.num,1);
G_matrix.rock.tau = 1.5;
%% Create Fracture System
numHFplanes = size(frac_ids,1);
for i=1:numHFplanes
     frac_id=frac_ids(i,:);
     G_matrix.rock.poro(frac_id)=frac_poro;
     G_matrix.rock.perm(frac_id)=frac_perm;  
     %Hassan Edit:
     G_matrix.rock.Dg(frac_id,:) = repmat(Dg_MS,size(frac_id,2),1).*1e2; 
     G_matrix.rock.Do(frac_id,:) = repmat(Do_MS,size(frac_id,2),1).*1e2; 
end
%% Create Wells
wells = struct;
wells(1).XFracCellIDs=well_ids;
wells(1).radius=wellRadius;
G=G_matrix;
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
% plotfracongrid(G_matrix,fracplanes,'label',false);
% G.rock.perm(well_ids)=1.0*darcy;
% hold off;
%%
useNatural =  true;
% Name of problem and pressure range
casename = 'barnett3comps'; 
pwf = 500*psia;
pinj = 8000*psia;
G2=G;
%% Define three-phase compositional flow model
% We define a three-phase compositional model with or without water.
% This is done by first instantiating the compositional model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid.

[fluid, info] = getShaleCompFluidCase(casename);

eosname = 'prcorr';  %'srk','rk','prcorr'
G1cell = cartGrid([1, 1, 1]);
G1cell = computeGeometry(G1cell);
[fluidMixture, info_G1cell] = getShaleCompFluidCase('barnett3comps_modified');

EOSModel = EquationOfStateModel(G1cell, fluid, eosname);
G1cell.rock=makeShaleRock(G1cell, matrix_perm, matrix_poro);
%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]); 

gravity reset on

arg = {G, G.rock, flowfluid, fluid, 'water', false};
arg_G1cell = {G1cell, G1cell.rock, flowfluid, fluid, 'water', false};

diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);
sparse_backend = SparseAutoDiffBackend();

if useNatural
    constructor = @NatVarsShaleModel;
    constructor_G1cell = @GenericNaturalVariablesModel;
%     model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', true);
else
%     model = OverallCompositionShaleCompositionalModel(G, [], flowfluid, fluid, 'water', true);
    model = OverallCompShaleModel(G, [], flowfluid, fluid, 'water', false);
end
%%
modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
modelDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', diagonal_backend);
modelMexDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', mex_backend);

ncomp = fluid.getNumberOfComponents();
s0 = [0.1, 0.9];
%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water, oil and gas initially. We also set up a simple-time step strategy 
% that ramps up gradually towards 30 day time-steps.
state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelSparseAD.EOSModel);
%% Calculate Nabla operator's denominator for ECL diffusion calculations
dd = computeNabla(modelDiagonalAD.G, modelDiagonalAD.G.rock);
modelDiagonalAD.operators.dd = dd(modelDiagonalAD.operators.internalConn);
%% Pick linear solver
% The AMGCL library is one possible solver option for MRST. It is fairly
% easy to write interfaces to other solvers using MEX files and/or the
% LinearSolverAD class. We call the routine for automatically selecting a
% reasonable linear solver for the specific model.
linsolve = selectLinearSolverAD(modelMexDiagonalAD);
disp(linsolve)
nls = NonLinearSolver('LinearSolver', linsolve);
%%
W = [];
for wi=1:numel(wells)
    W = addWell(W, G, G.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir','x');
end
for wi=1:numel(wells)
    W(wi).components = info.initial;
end

totTime = 15*year;
nSteps = 15;
% plotWell(G,W1);
dt = rampupTimesteps(totTime, 30*day, nSteps);
% schedule = simpleSchedule(dt, 'W', W);
%% Equilibration (Establish vertical equilibrium between gravity and capillary forces)
schedule = struct();
[W, W_equil] = deal(W);

W_equil(1).status = false; 

schedule.control = [struct('W', W);...  % Normal Schedule
                    struct('W', W_equil)];... % Equilibration
dt_equil = rampupTimesteps(40*day, 5*day, 7);
schedule.step.val = [dt_equil;dt];
schedule.step.control  = [2*ones(numel(dt_equil),1);1*ones(numel(dt),1)];
%% Simulate Problem
[ws, states, reports1] = simulateScheduleAD(state, modelDiagonalAD, schedule,'Verbose', true, 'nonlinearsolver', nls);
%% plotting
if (opt.shouldPlot)
    figure, 
    plotToolbar(G, states)
    view(40,30);
    axis tight equal;
    figure,
    plotWellSols(ws,cumsum(schedule.step.val))
end
%% Extracting cumulative gas production (field units)
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
q_gs = zeros(numel(ws),1);
for i = 1:numel(ws)
    q_gs(i)=-3.051e+6*ws{i}.qGs;
end
Q_gs = cumtrapz(tinDays,q_gs)/1e6; %cum gas produced in MMscf/D
%% Save variables in HPC
if ~opt.shouldPlot
    fpath =  '/scratch/ahass16/';
    fullFinalOut = [fpath, 'Verification1_Diffusion_Reference_x1024.mat'];
    save(fullFinalOut,'Q_gs','ws','G','schedule');
end