close all
% clear
clc
addpath('../ADFNE15')
addpath('modules/UnconvResourceSim/')
names = {'Producer_NF', 'Producer','EOR_NF', 'EOR'};

nkr = 1;
shouldPlot = 1;

[ws_NF_PI, states_NF_PI, schedule_NF_PI,G_NF_PI,reports_NF_PI] =  SlotDrill_NF('ProdBot_InjTop','nkr',nkr,'shouldPlot',shouldPlot); 
[ws_PI, states_PI, schedule_PI, G_PI,reports_PI] = SlotDrill('ProdBot_InjTop','nkr',nkr,'shouldPlot',shouldPlot);
[ws_P, states_P, schedule_P, G_P,reports_P] = SlotDrill('ProdBot','nkr',nkr,'shouldPlot',shouldPlot); 

[ws_NF, states_NF, schedule_NF,G_NF,reports_NF] = SlotDrill_NF('ProdBot','nkr',nkr,'shouldPlot',shouldPlot);
[ws_, states_, schedule_, G_,reports_] = SlotDrill('ProdBot','nkr',nkr,'shouldPlot',shouldPlot);

% spy(reports_PI.ControlstepReports{269}.StepReports{1}.NonlinearReport{1}.LinearSolver.A)

numIters = zeros(25,1);    
for ii=1:25
    numIters(ii) = reports_NF_PI.ControlstepReports{ii}.Iterations;
    spy(reports_NF_PI.ControlstepReports{ii}.StepReports{1}.NonlinearReport{1}.LinearSolver.A);
    saveas(gcf,['/Users/folorode/Documents/MATLAB/mrst-2018b/OMOresults/','spyEDFM_',int2str(ii),'.png']);
end
simtimeEDFM = reports.SimulationTime(1:25);


tinSecs = cumsum(schedule_.step.val);
tinSecs_NF = cumsum(schedule_NF.step.val);

ws = {ws_NF, ws_ , ws_NF_PI, ws_PI };
states = {states_NF, states_ , states_NF_PI, states_PI};
G = {G_NF, G_ , G_NF_PI, G_PI};

shortname = {'ProdOnly_NF', 'ProdOnly','ProdnInj_NF', 'ProdnInj'};

if ~shouldPlot
    fpath =  '/scratch/folorode/';

    fullFinalOut = [fpath, 'OutputVars2'];
    mkdir(fullFinalOut)
    save(fullFinalOut,'ws');
    save(fullFinalOut,'states');
    save(fullFinalOut','G');
    
else
    plotWellSols(ws, tinSecs, 'datasetnames', names,'field','qOs','linestyles',{'r','b','r--','b--'})

    figure(6)
    plotToolbar(G_, states_);
    axis equal tight off
    daspect([1 1 0.2])
    view(85, 20);

    figure(7)
    plotToolbar(G_NF, states_NF);
    axis equal tight off
    daspect([1 1 0.2])
    view(85, 20);

    figure(8)
    plotToolbar(G_PI, states_PI);
    axis equal tight off
    daspect([1 1 0.2])
    view(85, 20);

    figure(9)
    plotToolbar(G_NF_PI, states_NF_PI);
    axis equal tight off
    daspect([1 1 0.2])
    view(85, 20);

    % plotWell(G_NF, W);
    % title(names{2});
    % colorbar('horiz')

    qOs_ = zeros(numel(ws_),1);
    for i=1:numel(ws_)
        qOs_(i) = ws_{i}.qOs;
    end 
    qOs_ = -qOs_;

    % plot(tinSecs/86400,qOs_)

    Boi = 1.274424026342976;  % Bo = qOr(1)/qOs_ProdTop_InjDown(1)
    STOIIP_ = 6.289811*sum(G_.cells.volumes .* G_.rock.poro)*(1-0.23)/Boi %in stb   
    STOIIP_NF = 6.289811*sum(G_NF.cells.volumes .* G_NF.rock.poro)*(1-0.23)/Boi %in stb
    % STOIIP_PI = 6.289811*sum(G_PI.cells.volumes .* G_PI.rock.poro)*(1-0.23)/Boi %in stb   
    % STOIIP_NF_PI = 6.289811*sum(G_NF_PI.cells.volumes .* G_NF_PI.rock.poro)*(1-0.23)/Boi %in stb 
    % Np = trapz(tinSecs,qOs_);
    % Np = trapz(tinDays,qOs_ProdBot_InjTop)
    % Np = trapz(tinDays,qOs_ProdatTop)
    % RF = Np/STOIIP
    % figure(1)
    % plot(tinDays,qOs_ProdatTop,tinDays,qOs_ProdTop_InjDown,tinDays,qOs_ProdBot_InjTop,tinDays,qOs_ProdBot)
    % figure(2)
    % plot(tinDays,cumtrapz(tinDays,qOs_ProdatTop),tinDays,cumtrapz(tinDays,qOs_ProdTop_InjDown),tinDays,cumtrapz(tinDays,qOs_ProdBot_InjTop),tinDays,cumtrapz(tinDays,qOs_ProdBot))
end