function[ws,states, schedule,G,reports] = SlotDrill(case2run,varargin)

%     clear
%     close all

    opt = struct('nkr',        1, ...
                 'shouldPlot', 1 );
    opt = merge_options(opt, varargin{:});

    % Load necessary modules, etc 
    % mrstModule add hfm;             % hybrid fracture module
    mrstModule add ad-core;         
    mrstModule add ad-props ad-core 
    mrstModule add mrst-gui;        % plotting routines
    mrstModule add compositional; %mrstModule add ad-blackoil;
    mrstModule add upr;


%     nx = 10;
%     ny = 5;
%     Lx = 400;
%     Ly = 200;

    % Set up grid and rock
    x_max = 120;
    x_min = 0;
    max_nLogGrids  = 100;
    wf = 10*milli*meter;
    deltaR  = wf/2;  % Frac width will be 2*deltaR    0.001meters = 1mm 
    factor  = 1.3;

    fracLoc = 10;

    pt = zeros(1,10);
    pt(1) = deltaR;
    fprintf('Pt[%d]= %24.16f \n',1,pt(1));
    for ilogSpacnNum = 2:max_nLogGrids 
       deltaR = deltaR * factor;
       pt(ilogSpacnNum) = pt(ilogSpacnNum-1) + deltaR;
       if(pt(ilogSpacnNum) >= ((x_max-x_min)/2.0))
          fprintf('num of log grids is: %d \n', ilogSpacnNum);
          break
       end
       fprintf('Pt[%d]= %24.16f \n',ilogSpacnNum,pt(ilogSpacnNum));
    end

    ilogSpacnNum = ilogSpacnNum - 1;
    pt = pt(1:ilogSpacnNum);

    pt_under = pt(1:24);
    pt_above = pt(1:28);
    myX1 = [0, fracLoc + [-fliplr(pt_under) ,pt_above] ];
    myX = [myX1, fliplr(80-myX1)]; 

    % G = cartGrid([nx ny],[Lx Ly]);
    % G = makeLayeredGrid(G,diff(myX));

    G = tensorGrid(0:200:400, 0:100:200,myX);

    G = computeGeometry(G);
    
    if (opt.shouldPlot)
        %Plot triangular grid
        figure(1)
        plotGrid(G), view(5,45)
    end

    %%
    % MRST includes both natural variables and overall composition. This toggle
    % can switch between the modes.
    useNatural = true;


    % Name of problem and pressure range
    casename = 'bakken';
    pwf = 2000*psia;
    pinj = 8000*psia;

     rock = makeRock(G, 10*micro*darcy, 0.063);
    % rock = makeRock(G, 100*nano*darcy, 0.063);
    % rock = makeRock(G, 10*micro*darcy, 0.063);

    % fracIndx_X = 1:G.cartDims(1);
    % fracIndx_Y = 1:G.cartDims(2);
    % numFracCells = length(fracIndx_Y);
    % fracIndx_Y = ones(numFracCells,1)*ceil(G.cartDims(1)/2);

    [fracIndx_X,fracIndx_Y,fracIndx_Z] = meshgrid(1:G.cartDims(1), 1:G.cartDims(2), [26,80]);
    fraccells = sub2ind(G.cartDims, reshape(fracIndx_X,numel(fracIndx_X),1),reshape(fracIndx_Y,...
        numel(fracIndx_X),1), reshape(fracIndx_Z,numel(fracIndx_X),1));
    rock.poro(fraccells) = 0.33*3/100;
    rock.perm(fraccells) = 50*darcy;


    % Set up EDFM operators
    % TPFAoperators = setupEDFMOperatorsTPFA(G, rock, tol);

    %% Define three-phase compressible flow model
    % We define a three-phase black-oil model without dissolved gas or vaporized
    % oil. This is done by first instantiating the blackoil model, and then
    % manually passing in the internal transmissibilities and the topological
    % neighborship from the embedded fracture grid.
    % gravity reset off
    % model = ThreePhaseBlackOilModel(G, [], fluid, 'disgas', false, 'vapoil', false);
    % model.operators = TPFAoperators;

    
    [fluid, info] = getCompositionalFluidCase(casename);
    
    eosname = 'pr';%'prcorr';  %'srk','rk','prcorr'
    G1cell = cartGrid([1 1],[1 1]);
    G1cell = computeGeometry(G1cell);
    EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

    %Surface Conditions
    p_sc = 101325; %atmospheric pressure
    T_sc = 288.706;% 60 Farenheit
    [~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);
 
    flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]);
    
    % flowfluid.KGangiFn = @(p) power((1-power(((Pc - alpha.*p)./Pmax),m)),3);

    gravity reset on

    if useNatural
    % The version without rock as 2nd argument is used when TPFA operators will
    % be set up differently (eg in EDFM)
    %     model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', true);
        model = NaturalVariablesCompositionalModel(G, rock, flowfluid, fluid, 'water', true);
    else
    %     model = OverallCompositionCompositionalModel(G, [], flowfluid, fluid, 'water', true);
        model = OverallCompositionCompositionalModel(G, rock, flowfluid, fluid, 'water', true);
    end
    % model.operators = TPFAoperators;


    %% Set up initial state and schedule
    % We set up a initial state with the reference pressure and a mixture of
    % water and oil initially. We also set up a simple-time step strategy that
    % ramps up gradually towards 30 day time-steps.

    totTime = 20*year; 
    nSteps =25;

    ncomp = fluid.getNumberOfComponents();
    s0 = [0.23, 0.70, 0.07];   %s0 = [0.23, 0.77, 0.07];

    %                                 (G, p, T, s0, z0, eos)
    state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);


    wellRadius = 0.01;
    W = [];
    
    switch case2run
    case 'ProdTop' 
        % Producer
        W = verticalWell(W, G, rock, 1, 1, 26, ...
            'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
        W(1).components = info.initial;
    case 'ProdTop_InjBot'
        % Producer
        W = verticalWell(W, G, rock, 1, 1, 26, ...
            'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
        % Injector
        W = verticalWell(W, G, rock, 1, 1, 80, ...
            'comp_i', [0 0 1],'Name', 'Inj_Bot', 'Val', pinj, 'sign', 1, 'Type', 'bhp','Radius', wellRadius);
        W(1).components = info.initial;
        W(2).components = info.injection;
    case 'ProdBot' 
        % Producer
        W = verticalWell(W, G, rock, 1, 1, 80, ...
            'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
        W(1).components = info.initial;
    case 'ProdBot_InjTop'
        % Producer
        W = verticalWell(W, G, rock, 1, 1, 80, ...
            'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius); 
        % Injector
        W = verticalWell(W, G, rock, 1, 1, 26, ...
            'comp_i', [0 0 1],'Name', 'Inj_Top', 'Val', pinj, 'sign', 1, 'Type', 'bhp','Radius', wellRadius);
        W(1).components = info.initial;
        W(2).components = info.injection;
    otherwise
        warning('Case Does Not Exist. Running case with Prod Only at Bottom')
        % Producer
        W = verticalWell(W, G, rock, 1, 1, 80, ...
            'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
        W(1).components = info.initial;
    end
    % plotWell(G,W);
    


    % bc = [];
    % bc = pside(bc, G, 'LEFT', 90.0*barsa, 'sat', [0.3, 0.3, 0.4]);
    % bc.components = info.initial;

    % state  = initResSol(G, pRef, s0);
    dt = rampupTimesteps(totTime, 30*day, nSteps);
%     dt = dt(1:25);
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

    G.rock = rock;
    %% Plot the results in the interactive viewer
    figure(2); 
    plotToolbar(G, states)
    view(40,30);
    axis tight equal;
    
    %figure(3);
    plotWellSols(ws,cumsum(schedule.step.val))
    %plotWellSols(ws,cumsum(schedule.step.val)/86400)
    tinSecs = cumsum(schedule.step.val);
    tinDays = tinSecs./86400; 
    Boi = ws{1,1}(1).qOr/ws{1,1}(1).qOs; %calculate the initial oil formation volume factor from production rates
    % Boi = 1.274424026342976;  % Bo = qOr(1)/qOs_ProdTop_InjDown(1)
    STOIIP = 6.289811*sum((G.cells.volumes .* rock.poro))*(1-s0(1))/Boi; %in STB
    % STOIIP = 6.289811*sum(G.cells.volumes .* rock.poro)*(1-s0(1))/Boi; %in stb
    qO = [];
    for ii=1:size(ws)
        qO = [qO,-543439.65*ws{ii,1}(1).qOs]; %covnert from m^3/s to stb/d
    end
    Np = trapz(tinDays,qO);
    RF = Np/STOIIP
    % Np = trapz(tinDays,qOs_ProdTop_InjDown)
    % Np = trapz(tinDays,qOs_ProdBot_InjTop)
    % Np = trapz(tinDays,qOs_ProdatTop)
    % RF = Np/STOIIP
    % figure(1)
    % plot(tinDays,qOs_ProdatTop,tinDays,qOs_ProdTop_InjDown,tinDays,qOs_ProdBot_InjTop,tinDays,qOs_ProdBot)
    % figure(2)
    % plot(tinDays,cumtrapz(tinDays,qOs_ProdatTop),tinDays,cumtrapz(tinDays,qOs_ProdTop_InjDown),tinDays,cumtrapz(tinDays,qOs_ProdBot_InjTop),tinDays,cumtrapz(tinDays,qOs_ProdBot))

end
