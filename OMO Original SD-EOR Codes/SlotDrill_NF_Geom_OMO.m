	%% Loading Modules
    clear 
    clc
    close all
    
	Globals;
    case2run = 'ProdBot_InjTop';
    %case2run = 'ProdBot';
    opt = struct('nkr',        2, ...
                 'shouldPlot', 1 ); %change to 0 if running on HPC
    %opt = merge_options(opt, varargin{:});
    
    % Load necessary modules, etc 
    mrstModule add hfm;             % hybrid fracture module
    mrstModule add ad-core;         
    mrstModule add ad-props ad-core 
    mrstModule add mrst-gui;        % plotting routines
    mrstModule add compositional; %mrstModule add ad-blackoil;
    mrstModule add upr;
    tol=1e-5;
%% Basic Parameters
    physdim = [300 100 80]*meter; %sets grid's physical dimensions (x,y,z)
    wf = 10*milli*meter; %SlotFrac aperture
    fracloc = [10 70]*meter; %location of SlotFrac
    layer_z = 20; %No. of layers in z-direction, excluding Slotfrac layers
    layer_frac = fix(layer_z/(numel(fracloc)*2)); %calculate no. of layers above or under frac plan
    
%% Create a Box Reservoir Grid Using Geometric Grid (OMO):    
    % space on top of top frac
    numPts = 15;
    dist2grid = fracloc(1);
    ptZ = geomspace(wf/2, dist2grid, numPts, wf);
    refZ = fracloc(1);
    facesZcoords_aboveTop = refZ - fliplr(ptZ);
    facesZcoords_aboveTop(1)=0;
    % space below bottom frac
    refZ = fracloc(2);
    facesZcoords_belowBottom = refZ + ptZ;
    % space above bottom frac
    numPts = 6;
    dist2grid = 0.5*sum(fracloc)-fracloc(1); 
    ptZ = geomspace(wf/2, dist2grid, numPts, wf);
    facesZcoords_aboveBottom = refZ - fliplr(ptZ);
    
    % space below top frac
    refZ = fracloc(1);
    facesZcoords_belowTop = refZ + ptZ;
    z = [facesZcoords_aboveTop, facesZcoords_belowTop(1:end-1)...
                    facesZcoords_aboveBottom, facesZcoords_belowBottom];
   
%     DZ = diff(facesZcoords); 
%     z_1 = fliplr(fracloc(1) - geomspace(wf/2,fracloc(1),layer_frac+1,10*wf)); z_1(1)=0;
%     z_2 = geomspace(fracloc(1) + (wf/2),fracloc(1)+((fracloc(2)-fracloc(1))/2),layer_frac+10,10*wf);
%     z_11 = fliplr(fracloc(2)-geomspace(wf/2,((fracloc(2)-fracloc(1))/2),layer_frac+10,10*wf));
%     z_22 = geomspace(fracloc(2)+wf/2,physdim(3),layer_frac+1,10*wf);
%     z = [z_1, z_2(1:end-1), z_11, z_22];
    
    %G_matrix = tensorGrid(0:15:physdim(1), 0:20:physdim(2),z);
    G_matrix = tensorGrid(0:60:physdim(1), 0:25:physdim(2),z);
    G_matrix = computeGeometry(G_matrix);    
    G_matrix.rock=makeRock(G_matrix,20*micro*darcy,0.07);
    
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
    G_matrix.rock.perm(fraccells) = 10*darcy;
    
    if (opt.shouldPlot)
        %Plot triangular grid
        figure,
        %plotGrid(G_matrix), view(5,45)
        plotCellData(G_matrix,convertTo(G_matrix.rock.perm,milli*darcy));
        colorbar('horiz'); axis equal tight; view(3);
    end
%% 
%   Microfractures parallel to bedding planes
    rng(1234567); 
   
    tol4domain = tol*1e3;
    set1 = Field(DFN('dim',3,'n',200,'dir',0,'ddir',0,'minl',0.3,...
               'mu',5,'maxl',10,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',90,'ddip',0,...
               'shape','l','q',100),'Poly'); 
    set2 = Field(DFN('dim',3,'n',200,'dir',45,'ddir',-1e9,'minl',0.6,...
                'mu',5,'maxl',10,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',45,'ddip',-1e9,...
                'shape','l','q',4),'Poly');  
    set3 = Field(DFN('dim',3,'n',100,'dir',-45,'ddir',-1e9,'minl',0.6,...
                'mu',5,'maxl',10,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',45,'ddip',-1e9,...
                'shape','l','q',4),'Poly');
 
    set2 = [set2;set3];
    [set1_,nonPlanarSets1,fracArea1] = processStochFracs(set1);
    [set2_,nonPlanarSets2,fracArea2] = processStochFracs(set2); 
    fracArea = [fracArea1;fracArea2];
    
    fprintf('%d of %d set1 fracs were OK while %d of %d set2 fracs were OK \n',...
        numel(set1_),numel(set1),numel(set2_),numel(set2));
    
    fprintf('Number of nonplanar sets in sets 1 and 2 are : %d and %d respectively\n',...
        numel(nonPlanarSets1),numel(nonPlanarSets2));
    
    if (opt.shouldPlot)
%         figure,
%         Draw('ply',nonPlanarSets1);
%         Draw('ply',nonPlanarSets2);view(45,30)
        
        
        Draw('ply',set1_);
        Draw('ply',set2_);view(45,30)
    end         
%     Draw('ply',badSet);view(45,30)
% plotfracongrid(G,fracplanes(279))
    
    fracSet = [set1_ ;set2_]; 
    
    fracplanes = struct;
    
%     for i=1:numel(nonPlanarSets2)
%         fracplanes(i).points = nonPlanarSets2{i}(1:end-1,:);
%         fracplanes(i).aperture = 1*micro*meter; %1*micro*meter;
%         fracplanes(i).poro=0.5;
%         fracplanes(i).perm=100*micro*darcy;%0.01*darcy;
%     end
    
    for i=1:numel(fracSet)
        fracplanes(i).points = fracSet{i}(1:end-1,:);
        fracplanes(i).aperture = 0.5*micro*meter;% 100*nano*meter; %1*micro*meter;
        fracplanes(i).poro=0.5;
        fracplanes(i).perm=0.5*milli*darcy;%100*micro*darcy;%0.01*darcy;
    end
    
    
%     checkIfCoplanar(fracplanes)
    if (opt.shouldPlot)
        figure,
        plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
    end
    G=G_matrix;
    
%     Boi = 1.274424026342976;
%     STOIIP_NF = 6.289811*sum(G_matrix.cells.volumes .* G_matrix.rock.poro)*(1-0.23)/Boi %in stb
%     STOIIP = 6.289811*sum(G.cells.volumes .* G.rock.poro)*(1-0.23)/Boi %in stb
    % GlobTri = globalTriangulation(G_matrix);
    % [G,fracplanes] = preProcessingFractures(G, fracplanes, ...
    %                  'GlobTri', GlobTri);
    %% Process fracture(s)
    [G,fracplanes]=EDFMgrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
    
    
    %% Fracture-Matrix NNCs
    G=fracturematrixNNC3D(G,tol);
    %%
    [G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol);
    
    %% OMO: Projection-based NNCs
%     G = pMatFracNNCs3D(G,tol);
    %%
    % MRST includes both natural variables and overall composition. This toggle
    % can switch between the modes.
    useNatural = true;
    % Name of problem and pressure range
    casename = 'eagleford';%'onlydecane'; %'bakken';%'testFluid';
    pwf = 3500*psia;
    pinj = 8200*psia;
    % Set up EDFM operators
    TPFAoperators = setupEDFMOperatorsTPFA(G, G.rock, tol);
%     TPFAoperators = setupPEDFMOperatorsTPFA(G, G.rock, tol);
    %% Define fluid properties
%     Pc = 170*barsa;
%     alpha = 0.5;
%     Pmax = 200*barsa;
%     m = 0.5;
    %% Define three-phase compressible flow model
    % We define a three-phase black-oil model without dissolved gas or vaporized
    % oil. This is done by first instantiating the blackoil model, and then
    % manually passing in the internal transmissibilities and the topological
    % neighborship from the embedded fracture grid.
    % gravity reset off
    % model = ThreePhaseBlackOilModel(G, [], fluid, 'disgas', false, 'vapoil', false);
    % model.operators = TPFAoperators;
    [fluid, info] = getShaleCompositionalFluidCase(casename);
    
    eosname = 'prcorr';  %'srk','rk','prcorr'
    G1cell = cartGrid([1 1],[1 1]);
    G1cell = computeGeometry(G1cell);
    EOSModel = EquationOfStateModel(G1cell, fluid, eosname);
%     flowfluid = initSimpleADIFluid('n', [nkr, nkr, nkr], 'rho', [1000, 800, 10]);
%     model = OverallCompositionCompositionalModel(G, [], flowfluid, fluid, 'water', false);
    %Surface Conditions
    p_sc = 101325; %atmospheric pressure
    T_sc = 288.706;% 60 Farenheit
    [L, x, y, Z_L, Z_V, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);
%     nkr = 1;
    flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]);    % flowfluid.KGangiFn = @(p) power((1-power(((Pc - alpha.*p)./Pmax),m)),3);
    gravity reset on
    if useNatural
        model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', true);
    else
        model = OverallCompositionCompositionalModel(G, [], flowfluid, fluid, 'water', true);
    end
    model.operators = TPFAoperators;
    %% Set up initial state and schedule
    % We set up a initial state with the reference pressure and a mixture of
    % water and oil initially. We also set up a simple-time step strategy that
    % ramps up gradually towards 30 day time-steps.
    %% Set up initial state and schedule
    % We set up a initial state with the reference pressure and a mixture of
    % water and oil initially. We also set up a simple-time step strategy that
    % ramps up gradually towards 30 day time-steps.
    totTime = 7*year;
    nSteps =15;
    ncomp = fluid.getNumberOfComponents();
    s0 = [0.23, 0.70, 0.07];   %s0 = [0.23, 0.77, 0.07];
    %                                 (G, p, T, s0, z0, eos)
    state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);
    wellRadius = 0.01;
    W = [];
    
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
        W = verticalWell(W, G.Matrix, G.Matrix.rock, 3, 3, frac_z(2), ...
            'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
        W(1).components = info.initial;
    case 'ProdBot_InjTop'
        % Producer
        W = verticalWell(W, G.Matrix, G.Matrix.rock, 3, 3, frac_z(2), ...
            'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius); 
        W = verticalWell(W, G.Matrix, G.Matrix.rock, 3, 3, frac_z(1), ...
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
    % plotWell(G,W);
    
    
    % bc = [];
    % bc = pside(bc, G, 'LEFT', 90.0*barsa, 'sat', [0.3, 0.3, 0.4]);
    % bc.components = info.initial;
    % state  = initResSol(G, pRef, s0);
    dt = rampupTimesteps(totTime, 15*day, nSteps);
    % dt = dt(1:25);
    % dt = repmat(30*day,61,1);
    % schedule = simpleSchedule(dt, 'bc', bc);
    schedule = simpleSchedule(dt, 'W', W);
    % % Manually set the injection composition
    % [schedule.control.W.components] = deal([1     0     0     0     0     0     0     0]);
    % % Injection is pure gas
    % [schedule.control.W.compi] = deal([1, 0]);
    
    %% Simulate problem
%     fn = getPlotAfterStep(state, model, schedule);
%     [ws, states, report] = simulateScheduleAD(state, model, schedule, ...
%        'afterStepFn', getPlotAfterStep(state, model, schedule));
   
   
    [ws, states, reports] = simulateScheduleAD(state, model, schedule, 'Verbose', true);
%     numIters = zeros(25,1);    
%     for ii=1:25
%         numIters(ii) = reports{ii}.Iterations;
%         
%         spy(reports{ii}.StepReports{1}.NonlinearReport{1}.LinearSolver.A);
%         saveas(gcf,['/Users/folorode/Documents/MATLAB/mrst-2018b/OMOresults/','spyEDFM_',int2str(ii),'.png']);
%     end
%     simtimeEDFM = simtime(1:25);
    %% plotting
    % Plot the results in the interactive viewer
    figure, 
    plotToolbar(G, states)
    view(40,30);
    axis tight equal;
    
    plotWellSols(ws,cumsum(schedule.step.val))
    tinSecs = cumsum(schedule.step.val);
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
    
    %% Save Output Variables (Used in HPC).
    
    if ~opt.shouldPlot
        fpath =  '/scratch/ahass16/';
        
        fullFinalOut = [fpath, 'OutputVars2'];
        mkdir(fullFinalOut)
        save(fullFinalOut,'ws');
        save(fullFinalOut,'RF');
        save(fullFinalOut,'states');
        save(fullFinalOut','G');
    end