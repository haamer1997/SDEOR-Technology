%% PEDFM Validation case of + fracture of single phase flow
% Sec3.1 https://www.sciencedirect.com/science/article/pii/S0309170817300994
% water injection in a 3-dimensional fractured porous media using the pEDFM
% method
% This case regenerates Figure 3 in this paper
clear 
clc
close all
opt = struct('shouldPlot', 1 );
%     opt = merge_options(opt, varargin{:});
%% Load necessary modules
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties

%% Basic Parameters for the EDFM and pEDFM grids
tol=1e-5;
physdim = [3 1 1]*meter;
griddim = [3 1 1];
cellsize = physdim./griddim;

%% Create Box Reservoir Grid
G = tensorGrid(0:cellsize(1):physdim(1),...
               0:cellsize(2):physdim(2),...
               0:cellsize(3):physdim(3));
G = computeGeometry(G);
G.rock=makeRock(G,1*darcy,0.3);

%% Create Fracture System
fracplanes = struct;
fracplanes(1).points=[1,0,0; 1,1,0; 2,1,1; 2,0,1];
% fracplanes(2).points=[1.2,0,0; 1.2,1,0; 2,1,0.8; 2,0,0.8];
% fracplanes(3).points=[1,0,0.2; 1,1,0.2; 1.8,1,1; 1.8,0,1];


for i=1:numel(fracplanes)
     fracplanes(i).aperture = 1/25;     % unit is meters
     fracplanes(i).poro = 0.3;
     fracplanes(i).perm = 1000*darcy; % high-perm fracture case
%      fracplanes(i).perm = 1e-8*darcy;   % flow barrier case
end

%% EDFM PreProcessing 
[G,fracplanes]=EDFMshalegrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
%-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
%-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
G_edfm = G; %Store EDFM grid
% OMO: Projection-based NNCs % Switch from EDFM to pEDFM
G_pedfm = pMatFracNNCs3D(G,tol); %pEDFM

if (opt.shouldPlot)
    figure,
    plotfracongrid(G_pedfm,fracplanes);
%     plotfracongrid(G_edfm,fracplanes); % both give identical plots
end

%% Define resolved explict grid
% griddim_r = [225 225 3];
% cellsize_r = physdim./griddim_r;
% 
% G_r = tensorGrid(0:cellsize_r(1):physdim(1),...
%                 0:cellsize_r(2):physdim(2),...
%                 0:cellsize_r(3):physdim(3));
% G_r = computeGeometry(G_r);
% 
% G_r.rock=makeRock(G_r,1*darcy,0.3);
% 
% %Set fracture perm
% % Compute the I and J indices for the horizontal fracture
% IJ_hFrac = ones(125,2)*113; % There are 125 horizontal frac cells with J=113
% IJ_hFrac(:,1) = 51:175;     % horizontal frac starts at I=51; ends at I=175
% % Compute the I and J indices for the vertical fracture
% IJ_vFrac = ones(125,2)*113; % There are 125 vertical frac cells with I=113
% IJ_vFrac(:,2) = 51:175;     % vertical frac starts at J=51; ends at J=175
% 
% IJ_Fracs = [IJ_hFrac; IJ_vFrac]; % concatenate the IJ coords of both fractures
% globFcellIdx = sub2ind(G_r.cartDims(1:2), IJ_Fracs(:,1), IJ_Fracs(:,2)); % obtain global Fcell indices
% 
% fractureCells = [];
% numCellsPerLayer = griddim_r(1)*griddim_r(2);
% for layerID = 1:griddim_r(3)
%     fractureCells = [fractureCells; globFcellIdx+(layerID-1)*numCellsPerLayer];
% end
% G_r.rock.perm(fractureCells) = fracplanes(1).perm; %Frac perm
% 
% 
% figure,
% show = false([G_r.cells.num, 1]); % Hide matrix cells
% show(fractureCells) = true;      % Show Fracture cells
% clf; 
% plotCellData(G_r,convertTo(G_r.rock.perm,darcy),show,'EdgeColor', 'None','facealpha',0.5);
% plotCellData(G_r,convertTo(G_r.rock.perm,darcy),~show,'EdgeColor', 'None','facealpha',0.0);
% view(3); colorbar; axis equal tight


%% Set up TPFA operators for PEDFM, EDFM and the Refined Explicit Grid solution
TPFAoperators_pedfm = setupPEDFMOpsTPFA(G_pedfm, G_pedfm.rock, tol);
TPFAoperators_edfm = setupShaleEDFMOpsTPFA(G_edfm, G_edfm.rock, tol);
% TPFAoperators_r = setupOperatorsTPFA(G_r, G_r.rock);
    
%% Define fluid properties
% Define a three-phase fluid model without capillarity.
pRef = 5*barsa;  % specify reference/initial pressure
fluid = initSimpleADIFluid('phases','W',       ... % Fluid phase: water
                           'mu',  1*centi*poise, ... % Viscosity
                           'rho', 1000,          ... % Surface density [kg/m^3]
                           'c',   1e-4/barsa,    ... % Fluid compressibility
                           'cR',  1e-5/barsa     ... % Rock compressibility
                           );
                       
%% Define water flow model
gravity reset off
model_pedfm = WaterModel(G_pedfm,[], fluid);
model_pedfm.operators = TPFAoperators_pedfm;
% 
model_edfm = WaterModel(G_edfm,[], fluid);
model_edfm.operators = TPFAoperators_edfm;
% 
% model_r = WaterModel(G_r,[], fluid);
% model_r.operators = TPFAoperators_r;


s0 = (1);
bc_pedfm = [];
bc_pedfm = pside(bc_pedfm, G_pedfm, 'LEFT', 10*barsa,'sat', s0);
bc_pedfm = pside(bc_pedfm, G_pedfm, 'RIGHT', 1*barsa,'sat', s0);

bc_edfm = [];
bc_edfm = pside(bc_edfm, G_edfm, 'LEFT', 10*barsa,'sat', s0);
bc_edfm = pside(bc_edfm, G_edfm, 'RIGHT', 1*barsa,'sat', s0);
% 
% bc_r = [];
% bc_r  = pside(bc_r, G_r, 'LEFT', 10*barsa,'sat', s0);
% bc_r  = pside(bc_r, G_r, 'RIGHT', 1*barsa,'sat', s0);


clf, figure
plotGrid(G,'FaceColor', 'none'); view(3);
plotFaces(G, bc_pedfm.face(strcmp(bc_pedfm.type,'pressure')), 'r','facealpha',0.7);
plotFaces(G, bc_pedfm.face(strcmp(bc_pedfm.type,'pressure')), 'r','facealpha',0.7);
hold on;
hf={};
for i = 1:length(fracplanes)
    X=fracplanes(i).points(:,1);
    Y=fracplanes(i).points(:,2);
    Z=fracplanes(i).points(:,3);
    hf{end+1}=fill3(X,Y,Z,'b');
    set(hf{end},'facealpha',0.7);
end
set(gca, 'CameraPosition', [-24.1107  -33.7559  -19.0289]);


%% Set up initial state and schedule
state_pedfm  = initResSol(G_pedfm, pRef, s0);
state_edfm = initResSol(G_edfm, pRef, s0);
% state_r = initResSol(G_r, pRef, s0);

totTime = 60*day;
nSteps = 10;
dt = rampupTimesteps(totTime, 30*day, nSteps);

schedule_pedfm = simpleSchedule(dt, 'bc', bc_pedfm);
schedule_edfm = simpleSchedule(dt, 'bc', bc_edfm);
% schedule_r = simpleSchedule(dt, 'bc', bc_r);

%% Simulate problem
[ws_pedfm, states_pedfm, reports_pedfm] = simulateScheduleAD(state_pedfm, model_pedfm, schedule_pedfm, ...
   'afterStepFn', getPlotAfterStep(state_pedfm, model_pedfm, schedule_pedfm));

[ws_edfm, states_edfm, reports_edfm] = simulateScheduleAD(state_edfm, model_edfm, schedule_edfm, ...
   'afterStepFn', getPlotAfterStep(state_edfm, model_edfm, schedule_edfm));
% 
% [ws_r, states_r, reports_r] = simulateScheduleAD(state_r, model_r, schedule_r, ...
%    'afterStepFn', getPlotAfterStep(state_r, model_r, schedule_r));

%% Calculating the error for the simulation run
%L-2 norm
tinSecs = cumsum(schedule_edfm.step.val);
tinDays = tinSecs./86400;
l2_error = zeros(length(ws_edfm),1);
max_error_percentage = zeros(length(ws_edfm),1);
for ii=1:size(states_edfm)
    error_pressure = states_edfm{ii}.pressure - states_pedfm{ii}.pressure; %calculate pressure difference beteween EDFM and pEDFM
    [max_error,max_id] = max(abs(error_pressure)); %calculate and locate the maximum error difference 
    max_error_percentage(ii) = max_error / states_edfm{ii}.pressure(max_id)*100; %calculate the maximum error percentage assuming EDFM is the correct value
    l2_error(ii) = norm(error_pressure); %calculate the L2-norm for the error difference between EDFM and pEDFM pressures
end
%% plotting
% figure, 
% plotToolbar(G_pedfm, states_pedfm)
% view(40,30);
% axis tight equal;


%% Pressrue plot
% figure,
% plotCellData(G_pedfm, convertTo(states_pedfm{end,1}.pressure, barsa()), 'EdgeColor', 'k');
% xlabel('x'), ylabel('y'), zlabel('Depth');
% view(3); axis tight;colormap('jet');
% h=colorbar; caxis([0 10]); 
% h.Label.String = 'Pressure (bar)';
% set(gca, 'CameraPosition', [-24.1107  -33.7559  -19.0289]);


%% Pressure contour comparsion
% Figure and figure settings
figure('Position',[200 460 1200 400]);

% [X,Y]=meshgrid(cellsize_r(1)/2:cellsize_r(1):physdim(1),cellsize_r(2)/2:cellsize_r(2):physdim(2));
% nLayer=1;
% NumGridLayer=griddim_r(1)*griddim_r(2);
% startId=(nLayer-1)*NumGridLayer+1;
% endId=(nLayer)*NumGridLayer;
% 
% pres_r=reshape(states_r{end,1}.pressure(startId:endId),[griddim_r(1),griddim_r(2)]);
% 
% subplot(1,3,1);
% surf(X,Y,convertTo(pres_r, barsa()),'EdgeAlpha',0.1), view(130,20)
% colormap('jet');axis tight;
% title('Explicit solution pressure')

% [X,Y]=meshgrid(cellsize(1)/2:cellsize(1):physdim(1),cellsize(2)/2:cellsize(2):physdim(2));
% nLayer=1;
% NumGridLayer=griddim(1)*griddim(2);
% startId=(nLayer-1)*NumGridLayer+1;
% endId=(nLayer)*NumGridLayer;
% 
% pres_pedfm=reshape(states_pedfm{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);
% 
% subplot(1,3,2);
% surf(X,Y,convertTo(pres_pedfm, barsa()),'EdgeAlpha',0.3), view(130,20)
% colormap('jet');axis tight;
% title('pEDFM pressure')
% 
% pres_edfm=reshape(states_edfm{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);
% 
% subplot(1,3,3);
% surf(X,Y,convertTo(pres_edfm, barsa()),'EdgeAlpha',0.3), view(130,20)
% colormap('jet');axis tight;
% title('EDFM pressure')

%% Pressure overline comparsion
% IJ_hFrac=ones(griddim(1),2)*((griddim(1)-1)/2+1); IJ_hFrac(:,1)=1:griddim(1);
% CenterCellIDs = sub2ind(G.cartDims(1:2), IJ_hFrac(:,1), IJ_hFrac(:,2));
% 
% IJ_hFrac_r=ones(225,2)*113; IJ_hFrac_r(:,1)=1:225;
% CenterCellIDs_r = sub2ind(G_r.cartDims(1:2), IJ_hFrac_r(:,1), IJ_hFrac_r(:,2));
% 
% x_line=IJ_hFrac(:,1)*cellsize(1)-cellsize(1)/2;
% 
% pres_line_pedfm=states_pedfm{end,1}.pressure(CenterCellIDs);
% % pres_line_edfm=states_edfm{end,1}.pressure(CenterCellIDs);
% 
% 
% % x_line_r=IJ_hFrac_r(:,1)*cellsize_r(1)-cellsize_r(1)/2;
% % pres_line_r=states_r{end,1}.pressure(CenterCellIDs_r);
% 
% 
% figure('rend','painters','pos',[10 10 640 500]);
% title('Pressure solution over the horizontal center line');
% plot(x_line_r, convertTo(pres_line_r, barsa()),'k-', 'LineWidth', 2,'DisplayName','Explicit');
% hold on;
% plot(x_line, convertTo(pres_line_pedfm, barsa()),'o','Color','r', ...
%     'LineWidth', 2,'DisplayName','pEDFM');
% plot(x_line, convertTo(pres_line_edfm, barsa()),'^','Color','b', ...
%     'LineWidth', 2,'DisplayName','EDFM');
% hold off;
% 
% xlim([-1e-3 9]);
% ylim([-1e-3 10]);
% 
% set(gca,'FontSize',15);
% xlabel('X [m]');
% ylabel('Pressure [Bar]');
% legend; axis tight;
