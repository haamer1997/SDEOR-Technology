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
%Rock Properties
matrix_perm = 1*nano*darcy;
matrix_poro= 0.063;
%% Create Box Reservoir Grid
G1cell = cartGrid([1 1],[1 1]);
G1cell = computeGeometry(G1cell);
G1cell.rock=makeShaleRock(G1cell, matrix_perm, matrix_poro);
%%
useNatural =  true;

% Name of problem and pressure range
casename = 'simple';
pwf = 500*psia;
%% Define three-phase compositional flow model
[fluid, info] = getCompositionalFluidCase(casename);

eosname = 'prcorr';  %'srk','rk','prcorr'

EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]); 

gravity reset on

arg_G1cell = {G1cell, G1cell.rock, flowfluid, fluid, 'water', false};

diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);
sparse_backend = SparseAutoDiffBackend();

if useNatural
    constructor_G1cell = @GenericNaturalVariablesModel;
else
    constructor_G1cell = OverallCompShaleModel(G, [], flowfluid, fluid, 'water', false); %not used
end
%%
modelSparseAD_G1cell = constructor_G1cell(arg_G1cell{:}, 'AutoDiffBackend', sparse_backend);
modelDiagonalAD_G1cell = constructor_G1cell(arg_G1cell{:}, 'AutoDiffBackend', diagonal_backend);
modelMexDiagonalAD_G1cell = constructor_G1cell(arg_G1cell{:}, 'AutoDiffBackend', mex_backend);

ncomp = fluid.getNumberOfComponents();
s0 = [0.1, 0.9];
%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water, oil and gas initially. We also set up a simple-time step strategy 
% that ramps up gradually towards 30 day time-steps.
state0 = initCompositionalState(G1cell, info.pressure, info.temp, s0, info.initial, modelSparseAD_G1cell.EOSModel);
%% We can also initialize AD-variables in state to get derivatives
% plot state functions
modelDiagonalAD_G1cell = modelDiagonalAD_G1cell.validateModel();
groups = modelDiagonalAD_G1cell.getStateFunctionGroupings();
figure
plotStateFunctionGroupings(groups,'Stop','ComponentTotalMass')
figure
plotStateFunctionGroupings(groups,'Stop','ComponentTotalFlux')
% Initialize pressure as a AD-variable
stateAD = state0;
% stateAD.pressure = initVariablesADI(state0.pressure);
stateAD.x = initVariablesADI(state0.x);
stateAD.y = initVariablesADI(state0.y);
% Outputs are now ADI, with differentiation with respect to pressure
% mob = modelDiagonalAD_G1cell.getProp(state0, 'Mobility');
fug =  modelDiagonalAD_G1cell.getProp(state0, 'Fugacity'); %this works as values not AD variable
fugAD = modelDiagonalAD_G1cell.getProp(stateAD, 'Fugacity'); %I want to get fugacity derivative wrt liquid and gas composition of each component
disp(fugAD)