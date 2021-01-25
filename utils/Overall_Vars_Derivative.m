mrstModule add ad-core ad-props compositional
G = cartGrid([1, 1, 1]);
G = computeGeometry(G);
% [fluidMixture, info] = getBenchmarkMixture('spe5');
[fluidMixture, info] = getShaleCompFluidCase('oil_1_modified_2');
cmodel = GenericOverallCompositionModel(G, makeRock(G, 1, 1), initSimpleADIFluid(), fluidMixture);
% cmodel = GenericNaturalVariablesModel(G, makeRock(G, 1, 1), initSimpleADIFluid(), fluidMixture);
cmodel = cmodel.validateModel(); % Set up groups
%% Set up a compositional state
useAD = true;
statec = initCompositionalState(G, info.pressure,  info.temp, [0.3, 0.4, 0.3], info.initial, cmodel.EOSModel); %G, p, T, saturationVector, zi,eosModel
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

%% Calculate df_i/dz_i
df_V_by_dz = zeros(1,ncomp);
df_L_by_dz = zeros(1,ncomp); 
for i = 1:ncomp
   df_L_by_dz(i) = fug{i+1,1}.jac{i+1};
   df_V_by_dz(i) = fug{i+1,2}.jac{i+1};
end