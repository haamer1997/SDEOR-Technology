mrstModule add ad-core ad-props compositional
G = cartGrid([1, 1, 1]);
G = computeGeometry(G);
% [fluidMixture, info] = getBenchmarkMixture('spe5');
[fluidMixture, info] = getShaleCompFluidCase('barnett3comps_nat_vars');
cmodel = GenericOverallCompositionModel(G, makeRock(G, 1, 1), initSimpleADIFluid(), fluidMixture);
% cmodel = GenericNaturalVariablesModel(G, makeRock(G, 1, 1), initSimpleADIFluid(), fluidMixture);
cmodel = cmodel.validateModel(); % Set up groups
%% Set up a compositional state
useAD = true;
statec = initCompositionalState(G, info.pressure,  info.temp, [0.01, 0.01, 0.98], info.initial, cmodel.EOSModel); %G, p, T, saturationVector, zi,eosModel
stateAD = cmodel.validateState(statec);
stateAD = cmodel.getStateAD(stateAD, useAD); % No AD initialization
% Show the initialized container
% disp(stateAD.PVTProps)

%% Evaluate fugacity
fug = cmodel.getProp(stateAD, 'Fugacity');
component = cmodel.getProp(stateAD, 'ComponentPhaseMoleFractions');
ncomp = numel(info.injection)-1; %gets actual number of components ignoring dummy component
%% Show the PVT property cache after evaluation

%% Check phases in Fluid

%% Compute thermodynamic factor
dfv_dz = zeros(ncomp,ncomp);
dfl_dz = zeros(ncomp,ncomp);
z_x_dlnPhil_dz = zeros(ncomp,ncomp);
z_x_dlnPhiv_dz = zeros(ncomp,ncomp);
for i = 1:ncomp
   for j = 1:ncomp
       dfl_dz(i,j) = fug{i+1,1}.jac{j+1};
       dfv_dz(i,j) = fug{i+1,2}.jac{j+1};
       z_x_dlnPhil_dz(i,j) = (component{i+1,2}.val*dfl_dz(i,j))/fug{i+1,1}.val;
       z_x_dlnPhiv_dz(i,j) = (component{i+1,3}.val*dfv_dz(i,j))/fug{i+1,2}.val;
   end
end

rel2totalMass = 1;

Fij = z_x_dlnPhiv_dz;
Dideal = dIdealNatVars(statec,state.y,fluidMixture,rel2totalMass);
D = Dideal*Fij
