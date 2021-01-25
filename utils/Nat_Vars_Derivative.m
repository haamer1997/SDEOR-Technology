mrstModule add ad-core ad-props compositional
G = cartGrid([1, 1, 1]);
G = computeGeometry(G);
[fluidMixture, info] = getShaleCompFluidCase('bakken_light_dummy');
cmodel = GenericNaturalVariablesModel(G, makeRock(G, 1, 1), initSimpleADIFluid(), fluidMixture);
cmodel = cmodel.validateModel(); % Set up groups
%% Set up a compositional state
useAD = true;
statec = initCompositionalState(G, info.pressure,  info.temp, [0.3, 0.4, 0.3], info.initial, cmodel.EOSModel); %G, p, T, saturationVector, zi,eosModel
stateAD = cmodel.validateState(statec);
stateAD = cmodel.getStateAD(stateAD, useAD); % No AD initialization
%% Evaluate fugacity
fug = cmodel.getProp(stateAD, 'Fugacity');
component = cmodel.getProp(stateAD, 'ComponentPhaseMoleFractions');
ncomp = numel(info.injection)-1; %gets actual number of components ignoring dummy component
%% Check phases in Fluid
if isempty(fug{1,1}.jac{2})
    phase = 'o';
elseif isempty(fug{1,2}.jac{end})
    phase = 'g';
else
    phase ='og';
end
%% Calculate df_i/dx_i
df_V_by_dy = zeros(1,ncomp);
df_L_by_dx = zeros(1,ncomp); 
Da_V = zeros(1,ncomp);
Da_L = zeros(1,ncomp);
for i = 1:ncomp
   switch phase
       case 'o'
           df_L_by_dx(i) = fug{i,2}.jac{i+ncomp+3};
           Da_V(i) = Di(i); 
           Da_L(i) = fug{i,2}.val*Di_o(i)/(component{i,3}.val*df_L_by_dx(i));             
       case 'g'
           df_V_by_dy(i) = fug{i,1}.jac{i+1};
           %df_L_by_dx(i) = 1;
%            Da_V(i) = fug{i,1}.val*Di(i)/(component{i,2}.val*df_V_by_dy(i)); 
%            Da_L(i) = Di_o(i);  
       case 'og'
           df_V_by_dy(i) = fug{i,1}.jac{i+1};
           df_L_by_dx(i) = fug{i,2}.jac{i+ncomp+3};
%            Da_V(i) = fug{i,1}.val*Di(i)/(component{i,2}.val*df_V_by_dy(i)); 
%            Da_L(i) = fug{i,2}.val*Di_o(i)/(component{i,3}.val*df_L_by_dx(i)); 
   end
end