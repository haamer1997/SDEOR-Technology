%% Loading Modules
clear; 
clc;
close all;
Globals;
case2run = 'ProdBot_InjTop'; %'ProdBot';
opt = struct('nkr',        2, ...
             'shouldPlot', 0 ); %change to 0 if running on HPC
 %% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core;         
mrstModule add ad-props ad-core 
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;
%% Basic Parameters
tol=1e-5;    
physdim = [200 100 80]*meter; %sets grid's physical dimensions (x,y,z)
wf = 10*milli*meter; %SlotFrac aperture
fracloc = [10 70]*meter; %location of SlotFrac
layer_z = 20; %No. of layers in z-direction, excluding Slotfrac layers
layer_frac = fix(layer_z/(numel(fracloc)*2)); %calculate no. of layers above or under frac plan   
%% Create a Box Reservoir Grid Using Geometric Grid (OMO):    
% space on top of top frac
numPts = 20;
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
dist2grid = 0.5*sum(fracloc)-fracloc(1); 
ptZ = geomspace(wf/2, dist2grid, numPts, wf);
facesZcoords_aboveBottom = refZ - fliplr(ptZ);
% space below top frac
refZ = fracloc(1);
facesZcoords_belowTop = refZ + ptZ;
z = [facesZcoords_aboveTop, facesZcoords_belowTop(1:end-1)...
                facesZcoords_aboveBottom, facesZcoords_belowBottom];

G_matrix = tensorGrid(0:20:physdim(1), 0:20:physdim(2),z);
G_matrix = computeGeometry(G_matrix);    
G_matrix.rock=makeShaleRock(G_matrix,1*micro*darcy,0.07); %20*micro*darcy
    
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
[fracIndx_X,fracIndx_Y,fracIndx_Z] = meshgrid(1:G_matrix.cartDims(1), 1:G_matrix.cartDims(2), frac_z);
fraccells = sub2ind(G_matrix.cartDims, reshape(fracIndx_X,numel(fracIndx_X),1),reshape(fracIndx_Y,...
            numel(fracIndx_X),1), reshape(fracIndx_Z,numel(fracIndx_X),1));
G_matrix.rock.poro(fraccells) = 0.33;% 0.33*3/100;  %should be related to size of w_f wrt frac cell size
G_matrix.rock.perm(fraccells) = 10*darcy; %5*darcy

if (opt.shouldPlot)
    figure,
    plotCellData(G_matrix,convertTo(G_matrix.rock.perm,milli*darcy));
    colorbar('horiz'); axis equal tight; view(3);
end
%% Microfractures parallel to bedding planes
rng(12345); 
tol4domain = tol*1e3;

% set1 = Field(DFN('dim',3,'n',30,'dir',20,'ddir',0,'minl',10,...
%            'mu',15,'maxl',30,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',-45,'ddip',0,...
%            'shape','l','q',4),'Poly'); 
% set2 = Field(DFN('dim',3,'n',30,'dir',45,'ddir',-1e9,'minl',10,...
%             'mu',15,'maxl',30,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',45,'ddip',-1e9,...
%             'shape','l','q',4),'Poly');  
% set3 = Field(DFN('dim',3,'n',40,'dir',-45,'ddir',-1e9,'minl',10,...
%             'mu',15,'maxl',30,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',135,'ddip',-1e9,...
%             'shape','s'),'Poly'); %'l','q',4

set1 = Field(DFN('dim',3,'n',10,'dir',20,'ddir',0,'minl',0.3,...
           'mu',5,'maxl',10,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',-45,'ddip',0,...
           'shape','l','q',4),'Poly'); 
set2 = Field(DFN('dim',3,'n',10,'dir',45,'ddir',-1e9,'minl',0.6,...
            'mu',5,'maxl',10,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',45,'ddip',-1e9,...
            'shape','l','q',4),'Poly');  
set3 = Field(DFN('dim',3,'n',10,'dir',-45,'ddir',-1e9,'minl',0.6,...
            'mu',5,'maxl',10,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',135,'ddip',-1e9,...
            'shape','s'),'Poly'); %'l','q',4

set2 = [set2;set3];
[set1_,nonPlanarSets1,fracArea1] = processStochFracs(set1);
[set2_,nonPlanarSets2,fracArea2] = processStochFracs(set2); 
fracArea = [fracArea1;fracArea2];

fprintf('%d of %d set1 fracs were OK while %d of %d set2 fracs were OK \n',...
    numel(set1_),numel(set1),numel(set2_),numel(set2));

fprintf('Number of nonplanar sets in sets 1 and 2 are : %d and %d respectively\n',...
    numel(nonPlanarSets1),numel(nonPlanarSets2));
if (opt.shouldPlot)
    figure;
    Draw('ply',set1_);
    Draw('ply',set2_);view(45,30)
end 
fracSet = [set1_ ;set2_]; 
fracplanes = struct;   
for i=1:numel(fracSet)
    fracplanes(i).points = fracSet{i}(1:end-1,:);
    fracplanes(i).aperture = 0.1*ft;% 100*nano*meter; %1*micro*meter;
    fracplanes(i).poro=0.5;
    fracplanes(i).perm=1*milli*darcy;%100*micro*darcy;%0.01*darcy;%0.5*milli*darcy
end 
if (opt.shouldPlot)
    figure,
    plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
end
G=G_matrix;
%% Process fracture(s)
[G,fracplanes]=EDFMshalegrid(G,fracplanes,...
    'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracSet));
%% Fracture-Matrix NNCs
G=fracturematrixNNC3D(G,tol);
%% Fracture-Fracture NNCs
[G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol);
%% OMO: Projection-based NNCs
%G = pMatFracNNCs3D(G,tol);
% Set up EDFM operators
TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);
% Set up pEDFM operators
 %TPFAoperators = setupPEDFMOpsTPFA(G, G.rock, tol);
%% Define three-phase compressible flow model
useNatural = true;
% casename = 'bakken_light';
% pwf = 1000*psia;
% pinj = 1600*psia;
casename = 'oil_1' %'bakken_light';
pwf = 2500*psia;
pinj = 5700*psia;

% % diffusion
% G.rock.shaleMechanisms.diffusion = 1;
% %G.rock.Di=[2.8,2.0,1.2,1.2,1,0.95]*10^-7; %diffusivity coefficients for Bakken Light fluid
% %G.rock.Di=[0.28,0.187,0.95]*10^-7; %diffusivity coefficients for oil_1
% %fluid in oil phase
% G.rock.Di=[2.8,1.87,0.95]*10^-7; %diffusivity coefficients for oil_1 fluid in gas phase
% G.rock.tau = 2;

% %Dispersion
G.rock.shaleMechanisms.dispersion = 1;
G.rock.Di=[2.8,1.87,0.95]*10^-7;
G.rock.Di_o=G.rock.Di./100;
G.rock.tau = 2;

% % ECLIPSE diffusion
% G.rock.shaleMechanisms.ECLdiffusion = 1;
% G.rock.Di=[2.8,1.87,0.95]*10^-7;
% G.rock.Di_o=G.rock.Di./100;

[fluid, info] = getShaleCompFluidCase(casename);
eosname = 'prcorr';  %'srk','rk','prcorr'
G1cell = cartGrid([1 1],[1 1]);
G1cell = computeGeometry(G1cell);
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[L, x, y, Z_L, Z_V, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);
%     nkr = 1;
flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]);    % flowfluid.KGangiFn = @(p) power((1-power(((Pc - alpha.*p)./Pmax),m)),3);
gravity reset on

arg = {G, [], flowfluid, fluid, 'water', true};
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);
sparse_backend = SparseAutoDiffBackend();

if useNatural
    constructor = @GenericNaturalVariablesModel;
else
    constructor = @GenericOverallCompositionModel;
end

modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
modelDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', diagonal_backend);
modelMexDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', mex_backend);

modelSparseAD.operators = TPFAoperators;
modelDiagonalAD.operators = TPFAoperators;
modelMexDiagonalAD.operators = TPFAoperators;

%% Set up initial state
totTime = 8*year;
nSteps =12;
ncomp = fluid.getNumberOfComponents();
s0 = [0.23, 0.70, 0.07];   %s0 = [0.23, 0.77, 0.07];
%                                 (G, p, T, s0, z0, eos)
% state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);
state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelSparseAD.EOSModel);
%% Pick linear solver
% The AMGCL library is one possible solver option for MRST. It is fairly
% easy to write interfaces to other solvers using MEX files and/or the
% LinearSolverAD class. We call the routine for automatically selecting a
% reasonable linear solver for the specific model.
linsolve = selectLinearSolverAD(modelDiagonalAD);
disp(linsolve)
nls = NonLinearSolver('LinearSolver', linsolve);
%%
wellRadius = 0.01;
W = [];
%% Schedule
switch case2run
case 'ProdTop' 
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 10, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
case 'ProdTop_InjBot'
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 2, 2, 10, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    % Injector
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 70, ...
        'comp_i', [0 0 1],'Name', 'Inj_Bot', 'Val', pinj, 'sign', 1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
    W(2).components = info.injection;
case 'ProdBot' 
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 5, 3, frac_z(2), ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
case 'ProdBot_InjTop'
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 5, 3, frac_z(2), ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius); 
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 5, 3, frac_z(1), ...
        'comp_i', [0 0 1],'Name', 'Inj_Top', 'Val', pinj, 'sign', 1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
    W(2).components = info.injection;
otherwise
    warning('Case Does Not Exist. Running case with Prod Only at Bottom')
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 80, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
end
plotWell(G,W); 
dt = rampupTimesteps(totTime, 30*day, nSteps); %20*day
schedule = simpleSchedule(dt, 'W', W);  
%% Simulate problem
% [ws, states, reports] = simulateScheduleAD(state, model, schedule, 'Verbose', true);
[ws, states, reports] = simulateScheduleAD(state, modelDiagonalAD, schedule, 'nonlinearsolver', nls,'Verbose', true);
%% plotting
figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;
plotWellSols(ws,cumsum(schedule.step.val))
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
Qo = cumtrapz(tinDays,qO);
RF = Np/STOIIP
%% Save Output Variables (Used in HPC).
if ~opt.shouldPlot
    fpath =  '/scratch/ahass16/';
    fullFinalOut = [fpath, 'Diffusion_30NF_1nD_100psi_Oil.mat'];
    save(fullFinalOut,'ws','RF','states','G','schedule');
end