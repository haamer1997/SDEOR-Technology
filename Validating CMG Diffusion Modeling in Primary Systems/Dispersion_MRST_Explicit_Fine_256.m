clear 
clc
close all
opt = struct('nkr',        1, ...
    'shouldPlot', 1 );

%% Load necessary modules, etc
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core;
mrstModule add ad-props ad-core
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;

%% Basic Parameters
%physdim = [200 100 120]*meter;
physdim = [200 150 120]*meter;
%cellsize = physdim./griddim;
wellRadius=0.1*meter;
%HF Parameters
w_f = 0.0030;
x_f = 74.99;
%% Create Box Reservoir Grid
% G_matrix = tensorGrid(0:cellsize(1):physdim(1),...
%                       0:cellsize(2):physdim(2),...
%                       0:cellsize(3):physdim(3));

% [G_matrix,frac_ids,well_ids]=explicitStencilSingleFrac(physdim, 'aperture',w_f,...
%       'numStages',1,'fracSpacing', 200, 'fracHalfLength'  , x_f,      ...
%        'nxRefine',11,'tipNX',8,'nx',10, 'nz',8, 'gridType','linspacing');  

[G_matrix,frac_ids,well_ids]=explicitStencil(physdim, 'aperture',w_f,...
     'numStages',1,'fracSpacing', 200, 'fracHalfLength'  , x_f,      ...
     'nxRefine',250,'tipNY',2,'ny',8, 'nz',8, 'gridType','geomspacing'); 

%  [G_matrix,frac_ids,well_ids]=explicitStencil(physdim, 'aperture',w_f,...
%      'numStages',1,'fracSpacing', 200, 'fracHalfLength'  , x_f,      ...
%      'nxRefine',24,'tipNY',2,'ny',10, 'nz',10, 'gridType','linspacing');
 
G_matrix = computeGeometry(G_matrix);
% save(pwd,'G_matrix');
if (opt.shouldPlot)
    figure,
    plotGrid(G_matrix), view(5,45)
end

%Rock Properties
matrix_perm = 1*nano*darcy; %100*micro
matrix_poro= 0.063;

frac_poro = 0.5;
frac_perm = 1*milli*darcy;


G_matrix.rock=makeShaleRock(G_matrix, matrix_perm, matrix_poro);

%% Create Fracture System
numHFplanes = size(frac_ids,1);
for i=1:numHFplanes
     frac_id=frac_ids(i,:);
     G_matrix.rock.poro(frac_id)=frac_poro;
     G_matrix.rock.perm(frac_id)=frac_perm;     
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
% G.rock.perm(well_ids)=1.0*darcy;
% hold off;

%%
useNatural =  true;

% Name of problem and pressure range
casename = 'barnett3comps';
pwf = 500*psia;
pinj = 8000*psia;
G2=G;

%Use this to turn sorption on/off
% G.rock.shaleMechanisms.sorption = 1;


% %dispersion
% G.rock.shaleMechanisms.dispersion = 1;
% G.rock.Di=[2.8,2.5,1.9]*10^-7;
% G.rock.Di_o=[0,0,0];
% G.rock.tau = 2;

%ECLIPSE diffusion
G.rock.shaleMechanisms.ECLdiffusion = 1;
G.rock.Di=[2.8,2.5,1.9]*10^-7;
G.rock.Di_o=[0,0,0];

G1=G;

%% Define three-phase compositional flow model
% We define a three-phase compositional model with or without water.
% This is done by first instantiating the compositional model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid.

[fluid, info] = getShaleCompFluidCase(casename);

eosname = 'prcorr';  %'srk','rk','prcorr'
G1cell = cartGrid([1 1],[1 1]);
G1cell = computeGeometry(G1cell);
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]); 

gravity reset on

arg = {G, G.rock, flowfluid, fluid, 'water', false};
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);
sparse_backend = SparseAutoDiffBackend();

if useNatural
    constructor = @NatVarsShaleModel;
%     model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', true);
else
%     model = OverallCompositionShaleCompositionalModel(G, [], flowfluid, fluid, 'water', true);
    model = OverallCompShaleModel(G, [], flowfluid, fluid, 'water', false);
end

modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
modelDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', diagonal_backend);
modelMexDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', mex_backend);
ncomp = fluid.getNumberOfComponents();
s0 = [0.1, 0.9];

%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water, oil and gas initially. We also set up a simple-time step strategy 
% that ramps up gradually towards 30 day time-steps.

state1 = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelSparseAD.EOSModel);
%% Pick linear solver
% The AMGCL library is one possible solver option for MRST. It is fairly
% easy to write interfaces to other solvers using MEX files and/or the
% LinearSolverAD class. We call the routine for automatically selecting a
% reasonable linear solver for the specific model.
linsolve = selectLinearSolverAD(modelMexDiagonalAD);
disp(linsolve)
nls = NonLinearSolver('LinearSolver', linsolve);
%%
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
% plotWell(G,W1);
dt = rampupTimesteps(totTime, 30*day, nSteps);


schedule1 = simpleSchedule(dt, 'W', W1);
tinSecs = cumsum(schedule1.step.val);
tinDays = tinSecs./86400;

% [ws1, states1, reports1] = simulateScheduleAD(state1, model, schedule1, 'Verbose', true);
[ws1, states1, reports1] = simulateScheduleAD(state1, modelDiagonalAD, schedule1, 'nonlinearsolver', nls,'Verbose', true);
%% plotting
figure, 
plotToolbar(G, states1)
view(40,30);
axis tight equal;
figure,
plotWellSols(ws1,cumsum(schedule1.step.val))
%% Extracting cumulative gas production (field units)
q_gs = zeros(numel(ws1),1);
for i = 1:numel(ws1)
    q_gs(i)=-3.051e+6*ws1{i}.qGs;
end
Q_gs = cumtrapz(tinDays,q_gs)/1e6; %cum gas produced in MMscf/D
%% Modeling system without sorption
% 
% clear model,
% clear model.operators;
% clear G;
% 
% G=G2;
% if  strcmp(fractureModel,'pEDFM')
%     TPFAoperators= setupPEDFMOpsTPFA(G, G.rock, tol);
% else
%     TPFAoperators= setupEDFMOpsTPFA(G, G.rock, tol);
% end
% 
% if useNatural
%     model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', true);
% else
%     model = OverallCompositionCompositionalModel(G, [], flowfluid, fluid, 'water', true);
% end
% model.operators = TPFAoperators;
% 
% state2 = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);
% 
% W2 = [];
% for wi=1:numel(wells)
%     W2 = addWellEDFMshale(W2, G, G.Matrix.rock, wells(wi).XFracCellIDs, ...
%         'comp_i', [0.23, 0.76, 0.01],'Name', ['Prod',num2str(wi)], 'Val',...
%         pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir','x');
% end
% for wi=1:numel(wells)
%     W2(wi).components = info.initial;
% end
% 
% schedule2= simpleSchedule(dt, 'W', W2);
% 
% [ws2, states2, reports2] = simulateScheduleAD(state2, model, schedule2, 'Verbose', true);

%% plotting

% figure, 
% plotToolbar(G1, states1)
% view(40,30);
% axis tight equal;
% 
% figure, 
% plotToolbar(G2, states2)
% view(40,30);
% axis tight equal;
% 
% ws = {ws1, ws2};
% states = {states1, states2};
% 
% names = {'With Difussion', 'Without Difussion'};
% plotWellSols(ws, cumsum(schedule2.step.val), 'datasetnames', names,'field','qGs','linestyles',{'r','b'})
%%
% DX=diff(facesXcoords)';
% DY=diff(facesYcoords)';
% DZ=diff(facesZcoords)';
% fid = fopen('C:\Users\ahass16\Desktop\Research Work\mrst-2020a\modules\shale_chapter\examples\Diffusion_CMG_Validation\DI_DJ_DK_Fine.txt', 'w'); 
% fprintf(fid, 'DI IVAR     ** %d float\n', numel(DX));
% fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DX(:),ft));fprintf(fid, '\n');
% fprintf(fid, 'DJ JVAR     ** %d float\n', numel(DY));
% fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DY(:),ft)); fprintf(fid, '\n');
% fprintf(fid, 'DK KVAR     ** %d float\n', numel(DZ));
% fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DZ(:),ft)); fprintf(fid, '\n');
% fclose(fid);
