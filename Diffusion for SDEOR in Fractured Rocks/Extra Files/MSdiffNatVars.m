function [Dg,Do] = MSdiffNatVars(fluidname) %static mode implementation
%% Set up Unit Cell
G = cartGrid([1, 1, 1]);
G = computeGeometry(G);
G.rock.shaleMechanisms.Diffusion = false; 
[fluidMixture, info] = getShaleCompFluidCase(fluidname);
cmodel = GenericNaturalVariablesModel(G, makeRock(G, 1, 1), initSimpleADIFluid(), fluidMixture);
cmodel = cmodel.validateModel(); % Set up groups
%% Set up a compositional state
useAD = true;
statec = initCompositionalState(G, info.pressure,  info.temp, [0.01, 0.01, 0.98], info.initial, cmodel.EOSModel); %G, p, T, saturationVector, zi,eosModel
stateAD = cmodel.validateState(statec);
stateAD = cmodel.getStateAD(stateAD, useAD); % No AD initialization
[L, ~, ~, ~, ~, ~, ~] = standaloneFlash(info.pressure, info.temp, info.initial, cmodel.EOSModel); %check phases from flash calculations
%% Evaluate fugacity
fug = cmodel.getProp(stateAD, 'Fugacity');
ncomp = numel(info.injection)-1; %gets actual number of components ignoring dummy component
%% Check phases in Fluid
if isempty(fug{1,1}.jac{end-1}) %check if dfL/dx exists
    phase = 'g';
elseif isempty(fug{1,2}.jac{2}) %check if dfV/dy exists
    phase = 'o';
else
    phase ='og';
end
% % a different approach for phase identification using flash
% if L == 1 %MATLAB checks for integer 1 without issues
%     phase = 'o';
% elseif L == 0 
%     phase = 'g';
% else
%     phase ='og';
% end
%% Compute thermodynamic factor
dfv_dy = zeros(ncomp,ncomp);
dfl_dx = zeros(ncomp,ncomp);
z_x_dlnPhil_dx = zeros(ncomp,ncomp);
z_x_dlnPhiv_dy = zeros(ncomp,ncomp);
for i = 1:ncomp
   for j = 1:ncomp
       switch phase
           case 'o'
               dfl_dx(i,j) = fug{i,1}.jac{j+ncomp+3};
               z_x_dlnPhil_dx(i,j) = (statec.x(i)*dfl_dx(i,j))/fug{i,1}.val;
               z_x_dlnPhiv_dy(i,j) = 0;
           case 'g'
               dfv_dy(i,j) = fug{i,2}.jac{j+1};
               z_x_dlnPhil_dx(i,j) = 0;
               z_x_dlnPhiv_dy(i,j) = (statec.y(i)*dfv_dy(i,j))/fug{i,2}.val;               
           case 'og'
               dfl_dx(i,j) = fug{i,1}.jac{j+ncomp+3};
               dfv_dy(i,j) = fug{i,2}.jac{j+1};
               z_x_dlnPhil_dx(i,j) = (statec.x(i)*dfl_dx(i,j))/fug{i,1}.val;
               z_x_dlnPhiv_dy(i,j) = (statec.y(i)*dfv_dy(i,j))/fug{i,2}.val;              
       end
   end
end

%% Compute MS Diffusion Coefficients 
rel2totalMass = 1;

Fij = z_x_dlnPhiv_dy;
Dideal = dIdealNatVars(statec,statec.y,fluidMixture,rel2totalMass);
% Dg = Dideal*Fij;
% Dg = Dideal;

% Dg = reshape(Dideal*Fij,[1 numel(Dideal)]); %MS diffusion matrix is now a row vector
Dg = reshape(Dideal,[1 numel(Dideal)]); %MS diffusion matrix is now a row vector without \Gamma_{ij}
% Dg = reshape(Dg,[sqrt(numel(Dg)),sqrt(numel(Dg))]); %This line reshapes Dg vector into an n-by-n matrix

Fij = z_x_dlnPhil_dx;
Dideal = dIdealNatVars(statec,statec.x,fluidMixture,rel2totalMass);
% Do = Dideal*Fij;
% Do = Dideal;

% Do = reshape(Dideal*Fij,[1 numel(Dideal)]);
Do = reshape(Dideal,[1 numel(Dideal)]); %MS diffusion matrix is now a row vector without \Gamma_{ij}
% Do = reshape(Do,[sqrt(numel(Do)),sqrt(numel(Do))]);
end