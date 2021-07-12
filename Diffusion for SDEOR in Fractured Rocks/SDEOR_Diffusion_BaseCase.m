%% File Header
clear; 
clc;
close all;
Globals;
case2run ='ProdBot_InjTop';
opt = struct('nkr',        1, ... %keep at 1 for simplicity
             'shouldPlot', 0 ); %change to 0 if running on HPC
%% Load necessary modules
mrstModule add hfm;             % hybrid fracture module        
mrstModule add ad-props ad-core 
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;
%% Basic SDEOR Parameters
tol=1e-5;    
physdim = [200 100 80]*meter; %sets grid's physical dimensions (x,y,z)
wf = 10*milli*meter; %SlotFrac aperture
fracloc = [10 70]*meter; %location of SlotFrac
layer_z = 20; %No. of layers in z-direction, excluding Slotfrac layers
layer_frac = fix(layer_z/(numel(fracloc)*2)); %calculate no. of layers above or under frac plan   
%% Create a Box Reservoir Grid Using Geometric Grid (OMO):    
% space on top of top frac
endGeom = 10;
numPts = 20; %20
dist2grid = fracloc(1);
ptZ = geomspace(wf/2, dist2grid, numPts, wf);
refZ = fracloc(1);
facesZcoords_aboveTop = refZ - fliplr(ptZ);
facesZcoords_aboveTop(1)=0;
% space below bottom frac
refZ = fracloc(2);
facesZcoords_belowBottom = refZ + ptZ;
% space above bottom frac
numPts = 15; %15
dist2grid = endGeom;%0.5*sum(fracloc)-fracloc(1); 
ptZ = geomspace(wf/2, dist2grid, numPts, wf);
facesZcoords_aboveBottom = refZ - fliplr(ptZ);
% space below top frac
refZ = fracloc(1);
facesZcoords_belowTop = refZ + ptZ;
% space in-between two SD fractures starting from endGeom point 1 and 2
last_log_interval = facesZcoords_belowTop(end)-facesZcoords_belowTop(end-1);
linear_points = ceil(((fracloc(2)-endGeom)-(fracloc(1)+endGeom))/last_log_interval);
facesZcoords_between = linspace(fracloc(1)+endGeom,fracloc(2)-endGeom,linear_points);
z = [facesZcoords_aboveTop, facesZcoords_belowTop(1:end-1),...
                facesZcoords_between,facesZcoords_aboveBottom(2:end), facesZcoords_belowBottom];

frac_z=[];
for ii = 1:length(fracloc) %this nested loop locate fracture z-layer index from fracture location.
    n_layer = 0;
    for jj = 1:length(z)-1
        if ~((fracloc(ii) > z(jj)) && (fracloc(ii) < z(jj+1)))
            n_layer = n_layer+1;
        else
            frac_z = [frac_z,n_layer+1];
            break;
        end
    end
end
%% Assign Parameters to Matrix and SD Fractures
% Matrix
G_matrix = tensorGrid(0:5:physdim(1), 0:5:physdim(2),z);
G_matrix = computeGeometry(G_matrix);    
G_matrix.rock=makeShaleRock(G_matrix,10*micro*darcy,0.07);
    
% SD Fractures
[fracIndx_X,fracIndx_Y,fracIndx_Z] = meshgrid(1:G_matrix.cartDims(1), 1:G_matrix.cartDims(2), frac_z);
fraccells = sub2ind(G_matrix.cartDims, reshape(fracIndx_X,numel(fracIndx_X),1),reshape(fracIndx_Y,...
            numel(fracIndx_X),1), reshape(fracIndx_Z,numel(fracIndx_X),1));
G_matrix.rock.poro(fraccells) = 0.33;
G_matrix.rock.perm(fraccells) = 10*darcy; 
%% Compute MS diffusion coefficients to populate into matrix and fracture cells
% Molecular Diffusion Modeling
fluidname = 'oil_1_modified'; %this is a modified fluid such that a dummy component is added to the end (NatVars formulation)
G_matrix.rock.shaleMechanisms.Diffusion = true; 
G_matrix.rock.shaleMechanisms.MS_ECL_NNC_Static = 1; 
[Dg_MS,Do_MS] = MSdiffNatVars(fluidname); %MS static implementation
G_matrix.rock.Dg = repmat(Dg_MS,G_matrix.cells.num,1); G_matrix.rock.Do = repmat(Do_MS,G_matrix.cells.num,1);
G_matrix.rock.tau = 1.5;
%% Create 2 Intersecting NFs
fracplanes = struct;
fracplanes(1).points=[75,40,60; 75,60,60; 125,60,40; 125,40,40]; %[1,0,0; 1,1,0; 2,1,1; 2,0,1];
fracplanes(2).points=[75,60,40; 75,40,40; 125,40,60; 125,60,60];

for i=1:numel(fracplanes)
    fracplanes(i).aperture = 0.1*ft;
    fracplanes(i).poro=0.5;
    fracplanes(i).perm=10*milli*darcy;
    fracplanes(i).Dg = Dg_MS*1e2; %assume MS diffusion coefficeints are 2 orders of magnitude higher in NF planes
    fracplanes(i).Do = Do_MS*1e2;
end 
if (opt.shouldPlot) % visualize to check before pre-process
    figure,
    plotfracongrid(G_matrix,fracplanes,'label',false); 
end
G=G_matrix;
%% Process fracture(s)
[G,fracplanes]=EDFMshalegrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
%% Fracture-Matrix NNCs
G=fracturematrixShaleNNC3D(G,tol);
%% Fracture-Fracture NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
%% OMO: Projection-based NNCs
G = pMatFracNNCs3D(G,tol);
% Set up EDFM operators
% TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);
% Set up pEDFM operators
 TPFAoperators = setupPEDFMOpsTPFA(G, G.rock, tol);
%% Define three-phase compressible flow model
useNatural = true;
casename = 'oil_1';
pwf = 2500*psia;
rate = 0.003277;


[fluid, info] = getShaleCompFluidCase(casename);
eosname = 'prcorr'; 
G1cell = cartGrid([1 1],[1 1]);
G1cell = computeGeometry(G1cell);
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[L, x, y, Z_L, Z_V, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);
flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]);    % flowfluid.KGangiFn = @(p) power((1-power(((Pc - alpha.*p)./Pmax),m)),3);
gravity reset on

arg = {G, [], flowfluid, fluid, 'water', true};

diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);
sparse_backend = SparseAutoDiffBackend();

if useNatural
%     constructor_generic = @GenericNaturalVariablesModel; %GenericNatVarsShaleModel, GenericNaturalVariablesModel
    constructor = @NatVarsShaleModel;
else
    constructor = @GenericOverallCompositionModel;
end

modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
modelDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', diagonal_backend);
modelMexDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', mex_backend);

modelSparseAD.operators = TPFAoperators;
modelDiagonalAD.operators = TPFAoperators;
modelMexDiagonalAD.operators = TPFAoperators;
modelMexDiagonalAD.extraStateOutput = true;
%% Calculate Nabla operator's denominator for ECL diffusion calculations
% dd = computeNabla(modelDiagonalAD.G, modelDiagonalAD.G.rock);
% % dd = computeNabla(G_matrix, G_matrix.rock);
% modelDiagonalAD.operators.dd = dd(modelDiagonalAD.operators.internalConn);
% modelDiagonalAD.operators.g_gradz = modelDiagonalAD.operators.Grad(G.cells.centroids(:,3))*9.81;
% 
% dd = computeNabla(modelMexDiagonalAD.G, modelDiagonalAD.G.rock);
% % dd = computeNabla(G_matrix, G_matrix.rock);
% modelMexDiagonalAD.operators.dd = dd(modelMexDiagonalAD.operators.internalConn);
% modelMexDiagonalAD.operators.g_gradz = modelMexDiagonalAD.operators.Grad(G.cells.centroids(:,3))*9.81;
%% Set up initial state
totTime = 30*year;
nSteps =15;
ncomp = fluid.getNumberOfComponents();
s0 = [0.23, 0.70, 0.07]; 
state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelSparseAD.EOSModel);
%% Pick linear solver
% The AMGCL library is one possible solver option for MRST. It is fairly
% easy to write interfaces to other solvers using MEX files and/or the
% LinearSolverAD class. We call the routine for automatically selecting a
% reasonable linear solver for the specific model.
linsolve = selectLinearSolverAD(modelDiagonalAD,'useAMGCL',true,'useCPR',true);
disp(linsolve)
nls = NonLinearSolver('LinearSolver', linsolve);
%% Schedule
wellRadius = 0.01;
W = [];
switch case2run
case 'ProdBot' 
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 5, 3, frac_z(2), ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
case 'ProdBot_InjTop'
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 20, 10, frac_z(2), ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius); 
%     W = verticalWell(W, G.Matrix, G.Matrix.rock, 20, 10, frac_z(1), ...
%         'comp_i', [0 0 1],'Name', 'Inj_Top', 'Val', pinj, 'sign', 1, 'Type', 'bhp','Radius', wellRadius); %control by injection pressure
    % Injector
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 20, 10, frac_z(1), ...
        'comp_i', [0 0 1],'Name', 'Inj_Top', 'Val', rate, 'sign', 1, 'Type', 'rate','Radius', wellRadius); %control by injection rate
    W(1).components = info.initial;
    W(2).components = info.injection;
otherwise
    warning('Case Does Not Exist. Running case with Prod Only at Bottom')
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 80, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
end 
dt = rampupTimesteps(totTime, 30*day, nSteps); %20*day
% schedule = simpleSchedule(dt, 'W', W);
%% Equilibration (Establish vertical equilibrium between gravity and capillary forces)
schedule = struct();
[W, W_equil] = deal(W);

W_equil(1).status = false; W_equil(2).status = false;

schedule.control = [struct('W', W);...  % Normal Schedule
                    struct('W', W_equil)];... % Equilibration
dt_equil = rampupTimesteps(40*day, 5*day, 7);
schedule.step.val = [dt_equil;dt];
schedule.step.control  = [2*ones(numel(dt_equil),1);ones(numel(dt),1)]; 
%% Simulate problem
[ws, states, reports] = simulateScheduleAD(state, modelDiagonalAD, schedule, 'nonlinearsolver', nls, 'Verbose', true);
% [ws, states, reports] = simulateScheduleAD(state, modelMexDiagonalAD, schedule, 'nonlinearsolver', nls, 'Verbose', true);
%% plotting
if (opt.shouldPlot)
    figure, 
    plotToolbar(G, states)
    view(40,30);
    axis tight equal;
    plotWellSols(ws,cumsum(schedule.step.val))
end
%% Calculate RF
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
Boi = ws{1,1}(1).qOr/ws{1,1}(1).qOs; %calculate the initial oil formation volume factor from production rates
STOIIP = 6.289811*sum((G_matrix.cells.volumes .* G_matrix.rock.poro))*(1-s0(1))/Boi; %in STB

qO = [];
for ii=1:size(ws)
    qO = [qO,-543439.65*ws{ii,1}(1).qOs]; %covnert from m^3/s to stb/d
end
Np = trapz(tinDays,qO);
RF = Np/STOIIP;
%% Get cumulative reservoir fluid withdrawn
x = zeros(length(ws),1);
for i = 1:length(ws)
   x(i)= -543439.7*ws{i}(1).qTr;
end
QTr = trapz(tinDays,x);
%% Save Output Variables (Used in HPC).
if ~opt.shouldPlot
    fpath =  '/scratch/ahass16/';
    fullFinalOut = [fpath, 'EOR_512NF_SemiCircle_30years.mat'];
    save(fullFinalOut,'-v7.3');
end