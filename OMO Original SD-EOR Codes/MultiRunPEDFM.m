close all
% clear
clc
addpath('../ADFNE15')
addpath('modules/UnconvResourceSim/')
names = {'Producer_NF', 'Producer','EOR_NF', 'EOR'};

nkr = 1;
shouldPlot = 1;

[ws_NF_PIE, states_NF_PIE, schedule_NF_PIE, G_NF_PIE, reports_NF_PIE] =  SlotDrill_NF('ProdBot_InjTop','nkr',nkr,'shouldPlot',shouldPlot);
[ws_NF_PI, states_NF_PI, schedule_NF_PI, G_NF_PI, reports_NF_PI] =  SlotDrill_NFpEDFM('ProdBot_InjTop','nkr',nkr,'shouldPlot',shouldPlot); 
% [ws_PI, states_PI, schedule_PI, G_PI, reports_PI] = SlotDrill('ProdBot_InjTop','nkr',nkr,'shouldPlot',shouldPlot);
% [ws_NF, states_NF, schedule_NF, G_NF, reports_NF] = SlotDrill_NF('ProdBot','nkr',nkr,'shouldPlot',shouldPlot);
% [ws_, states_, schedule_, G_, reports_] = SlotDrill('ProdBot','nkr',nkr,'shouldPlot',shouldPlot);

numItersPEDFMlowCf = zeros(25,1);    
for ii=1:25
    numItersPEDFMlowCf(ii) = reports_NF_PI.ControlstepReports{ii}.Iterations;
    spy(reports_NF_PI.ControlstepReports{ii}.StepReports{1}.NonlinearReport{1}.LinearSolver.A);
    saveas(gcf,['/Users/folorode/Documents/MATLAB/mrst-2018b/OMOresults/','spyPEDFMlowCf_',int2str(ii),'.png']);
end
simtimePEDFMlowCf = reports_NF_PI.SimulationTime(1:25);

numItersEDFMlowCf = zeros(25,1);    
for ii=1:25
    numItersEDFMlowCf(ii) = reports_NF_PIE.ControlstepReports{ii}.Iterations;
    spy(reports_NF_PIE.ControlstepReports{ii}.StepReports{1}.NonlinearReport{1}.LinearSolver.A);
    saveas(gcf,['/Users/folorode/Documents/MATLAB/mrst-2018b/OMOresults/','spyEDFMlowCf_',int2str(ii),'.png']);
end
simtimeEDFMlowCf = reports_NF_PIE.SimulationTime(1:25);

% tinSecs = cumsum(schedule_.step.val);
tinSecs_NF_PI_EDFM = cumsum(schedule_NF_PIE.step.val);
tinSecs_NF_PI_PEDFM = cumsum(schedule_NF_PI.step.val);
tinSecs = {tinSecs_NF_PI_EDFM, tinSecs_NF_PI_PEDFM};
% ws = {ws_NF, ws_ , ws_NF_PI, ws_PI };
% states = {states_NF, states_ , states_NF_PI, states_PI};
% G = {G_NF, G_ , G_NF_PI, G_PI};
ws = {ws_NF_PIE, ws_NF_PI};
states = {states_NF_PIE, states_NF_PI};
G = {G_NF_PIE, G_NF_PI};

shortname = {'ProdnInj_NF_EDFM', 'ProdnInj_NF_PEDFM'};

if ~shouldPlot
    fpath =  '/Users/folorode/Documents/MATLAB/mrst-2018b/';

    fullFinalOut = [fpath, 'OMOresults/'];
%     mkdir(fullFinalOut)
    save(fullFinalOut,'ws');
    save(fullFinalOut,'states');
    save(fullFinalOut','G');
    
else
    plotWellSols(ws, tinSecs, 'datasetnames', shortname,'field','qOs','linestyles',{'r','b'})

    figure,
    plotToolbar(G_NF_PIE, states_NF_PIE);
    axis equal tight off
    daspect([1 1 0.2])
    view(85, 20);
    
    figure,
    plotToolbar(G_NF_PI, states_NF_PI);
    axis equal tight off
    daspect([1 1 0.2])
    view(85, 20);


%     qOs_ = zeros(numel(ws_),1);
%     for i=1:numel(ws_)
%         qOs_(i) = ws_{i}.qOs;
%     end
%     qOs_ = -qOs_;
% 
%     % plot(tinSecs/86400,qOs_)
% 
%     Boi = 1.274424026342976;  % Bo = qOr(1)/qOs_ProdTop_InjDown(1)
%     STOIIP_ = 6.289811*sum(G_.cells.volumes .* G_.rock.poro)*(1-0.23)/Boi %in stb   
%     STOIIP_NF = 6.289811*sum(G_NF.cells.volumes .* G_NF.rock.poro)*(1-0.23)/Boi %in stb
%     % STOIIP_PI = 6.289811*sum(G_PI.cells.volumes .* G_PI.rock.poro)*(1-0.23)/Boi %in stb   
%     % STOIIP_NF_PI = 6.289811*sum(G_NF_PI.cells.volumes .* G_NF_PI.rock.poro)*(1-0.23)/Boi %in stb 
%     % Np = trapz(tinSecs,qOs_);
%     % Np = trapz(tinDays,qOs_ProdBot_InjTop)
%     % Np = trapz(tinDays,qOs_ProdatTop)
%     % RF = Np/STOIIP
%     % figure(1)
%     % plot(tinDays,qOs_ProdatTop,tinDays,qOs_ProdTop_InjDown,tinDays,qOs_ProdBot_InjTop,tinDays,qOs_ProdBot)
%     % figure(2)
%     % plot(tinDays,cumtrapz(tinDays,qOs_ProdatTop),tinDays,cumtrapz(tinDays,qOs_ProdTop_InjDown),tinDays,cumtrapz(tinDays,qOs_ProdBot_InjTop),tinDays,cumtrapz(tinDays,qOs_ProdBot))
end