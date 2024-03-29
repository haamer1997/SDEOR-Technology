function [problem, state] = eqnsShaleNaturalVars(state0, state, model, dt, drivingForces, varargin)
% Equations for natural variables formulation

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
%%
opt = struct('Verbose',     mrstVerbose,...
            'reverseMode', false,...
            'resOnly',     false,...
            'staticWells',  false, ...
            'propsPressure', [], ...
            'reduceToPressure', false, ...
            'iteration',   -1);

opt = merge_options(opt, varargin{:});
% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
f = model.fluid;
W = drivingForces.W;


fluid = model.fluid;
compFluid = model.EOSModel.fluid;
% Properties at current timestep
[p, sW, sO, sG, x, y, z, temp, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'so', 'sg', 'x', 'y', 'z', 'T', 'wellSol');
z = expandMatrixToCell(z);
[p0, sW0, sO0, sG0, x0, y0, temp0, wellSol0] = model.getProps(state0, ...
    'pressure', 'water', 'so', 'sg', 'x', 'y', 'T', 'wellSol');
[pureLiquid, pureVapor, twoPhase] = model.getFlag(state);

if 1
    stol = 1e-6;
    pureWater = sO + sG < stol;
    sO(~pureVapor & pureWater) = stol;
    sG(~pureLiquid & pureWater) = stol;
    
    [pureLiquid0, pureVapor0, twoPhase0] = model.getFlag(state0);
    pureWater0 = sO0 + sG0 < stol;
    sO0(~pureVapor0 & pureWater0) = stol;
    sG0(~pureLiquid0 & pureWater0) = stol;
end

multiPhase = twoPhase;
freeOil = twoPhase;
freeGas = twoPhase;

z_tol = model.EOSModel.minimumComposition;
state0.x = ensureMinimumFraction(state0.x, z_tol);
state0.y = ensureMinimumFraction(state0.y, z_tol);
state0.components = ensureMinimumFraction(state0.components, z_tol);

x = ensureMinimumFraction(x, z_tol);
y = ensureMinimumFraction(y, z_tol);
x = expandMatrixToCell(x);
y = expandMatrixToCell(y);

[wellvars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
if opt.staticWells
    wellvars = {[]};
end
ncomp = compFluid.getNumberOfComponents();
[xnames, ynames, cnames] = deal(model.EOSModel.fluid.names);
for i = 1:ncomp
    xnames{i} = ['v_', cnames{i}];
    ynames{i} = ['w_', cnames{i}];
end

twoPhaseIx = find(twoPhase);

wtmp = ones(nnz(twoPhase), 1);

w = cell(ncomp, 1);
[w{:}] = deal(wtmp);

nc = model.G.cells.num;
for i = 1:(ncomp-1)
    w{i} = y{i}(twoPhase);
end

so = sO(freeOil);

nwellvar = sum(cellfun(@numel, wellvars));
nwelleqs = numel(wellvars);


if opt.resOnly
    initfn = @deal;
else
    initfn = @(varargin) model.AutoDiffBackend.initVariablesAD(varargin{:});
end
sg = sG(freeGas);

if model.water
    [p, x{1:ncomp-1}, sW, wellvars{:}, so, w{1:ncomp-1}, sg] = initfn(...
     p, x{1:ncomp-1}, sW, wellvars{:}, so, w{1:ncomp-1}, sg);
    primaryVars = {'pressure', xnames{1:end-1}, 'satw', wellVarNames{:}, 'sato', ynames{1:end-1}, 'satg'};
else
    [p, x{1:ncomp-1}, wellvars{:}, so, w{1:ncomp-1}, sg] = initfn(...
     p, x{1:ncomp-1}, wellvars{:}, so, w{1:ncomp-1}, sg);
    primaryVars = {'pressure', xnames{1:end-1}, wellVarNames{:}, 'sato', ynames{1:end-1}, 'satg'};
    sW = zeros(model.G.cells.num, 1);
end
sample = p;

sO = model.AutoDiffBackend.convertToAD(sO, sample);
sO(freeOil) = so;
sG = model.AutoDiffBackend.convertToAD(sG, sample);
sG(freeGas) = sg;


[sO, sG] = setMinimums(model, state, sW, sO, sG, pureVapor, pureLiquid);
[sO0, sG0] = setMinimums(model, state0, sW0, sO0, sG0, pureVapor0, pureLiquid0);

if isempty(twoPhaseIx) || opt.resOnly
    reorder = [];
else
    n2ph = nnz(twoPhase);
    nVars = sum(p.getNumVars());
    reorder = 1:nVars;
    start = nc + twoPhaseIx;
    stop = nc*(ncomp+model.water) + nwellvar + (1:n2ph);
    reorder(start) = stop;
    reorder(stop) = start;
end


% Property pressure different from flow potential
p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end


cellJacMap = cell(numel(primaryVars), 1);

offset = ncomp + model.water + numel(wellvars);
if any(twoPhase) && ~all(twoPhase)
    for i = 1:(ncomp+1)
        cellJacMap{i + offset} = twoPhaseIx;
    end
end
x{end} = ones(model.G.cells.num, 1);
w{end} = ones(nnz(twoPhase), 1);

for i = 1:ncomp-1
    x{end} = x{end}-x{i};
    if any(twoPhase)
        w{end} = w{end}-w{i};
    end
end

for i = 1:ncomp
    y{i} = ~pureLiquid.*x{i} + value(x{i}).*pureLiquid;
    if any(twoPhase)
        if ~opt.resOnly
            assert(isa(y{i}, 'ADI'));
        end
        y{i}(twoPhase) = w{i};
    end
    x{i}(pureVapor) = value(x{i}(pureVapor));
end


% Compute properties and fugacity
[xM,  yM,  rhoO,  rhoG,  muO,  muG, f_L, f_V, xM0, yM0, rhoO0, rhoG0] = ...
                  model.getTimestepPropertiesEoS(state, state0, p_prop, temp, x, y, z, sO, sG, cellJacMap);
%%
if model.water
    sat = {sW, sO, sG};
else
    sat = {sO, sG};
end


if model.water
    [krW, krO, krG] = model.evaluateRelPerm(sat);
else
    [krO, krG] = model.evaluateRelPerm(sat);
end

% Compute transmissibility
T = s.T;
if isfield(model.G.rock.shaleMechanisms,'NNCDiffusion') || model.G.rock.shaleMechanisms.Diffusion
    Tdiff_gas = s.Tdiff_gas; Tdiff_oil = s.Tdiff_oil;
end
if isfield(model.G.rock.shaleMechanisms,'TraditionalDiffusion')
    Tdiff_gas = s.Tdiff_gas; Tdiff_oil = s.Tdiff_oil;
    Tdiff_gas(s.NNC_startID:end,:) = 0; Tdiff_oil(s.NNC_startID:end,:) = 0; %model traditional diffusion without M-F
end
% Gravity gradient per face
gdz = model.getGravityGradient();

% Oil flux
rhoOf  = s.faceAvg(sO.*rhoO)./max(s.faceAvg(sO), 1e-8);
mobO   = krO./muO;
dpO    = s.Grad(p) - rhoOf.*gdz;
upco  = (value(dpO)<=0);
vO = -s.faceUpstr(upco, mobO).*T.*dpO;
if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'nonDarcyinFracs')
    %HR-OMO edit: Non-Darcy oil flow
    turbFactor=3.2808*1.485e9./(model.G.rock.perm/(milli*darcy)).^1.021;
    turbFactor = (~s.isMatrix).*turbFactor;
    rhoO_X_mobO_X_k = rhoO.*mobO.*model.G.rock.perm.*turbFactor;
    
    rhoO_X_mobO_X_k_2face = s.faceUpstr(upco, rhoO_X_mobO_X_k);
%     vO = 2*vO./(1 + power((1 - 4*rhoO_X_mobO_X_k_2face.*vO),0.5)); %old
     vO = 2*vO./(1 + power((1 + 4*rhoO_X_mobO_X_k_2face.*vO),0.5)); %+
%     vO = 2*vO./(1 - power((1 + 4*rhoO_X_mobO_X_k_2face.*vO),0.5)); %-
end

if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'nonDarcyEverywhere')
    %HR-OMO edit: Non-Darcy oil flow
%     turbFactor = model.G.rock.turbFactor;
    turbFactor=3.2808*1.485e9./(model.G.rock.perm/(milli*darcy)).^1.021;
    rhoO_X_mobO_X_k = rhoO.*mobO.*model.G.rock.perm.*turbFactor;
    
    rhoO_X_mobO_X_k_2face = s.faceUpstr(upco, rhoO_X_mobO_X_k);
    vO = 2*vO./(1 + power((1 - 4*rhoO_X_mobO_X_k_2face.*vO),0.5));
end

%HR edit: Incorporate Gangi model by scaling oil velocity
if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms, 'Gangi')
    F_appStar = model.fluid.KGangiFn(p);
    F_appStar(~s.isMatrix) = 1;
    [KGangif, ~] = s.splitFaceCellValue(s, upco, F_appStar);
    vO = (KGangif.*vO);
end
% Gas flux
rhoGf  = s.faceAvg(sG.*rhoG)./max(s.faceAvg(sG), 1e-8);
mobG   = krG./muG;
if isfield(fluid, 'pcOG')
    pcOG  = fluid.pcOG(sG);
    pG = p + pcOG;
else
    pG = p;
end
dpG    = s.Grad(pG) - rhoGf.*gdz;

upcg  = (value(dpG)<=0);
vG = -s.faceUpstr(upcg, mobG).*T.*dpG;
if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'nonDarcyinFracs')
    %HR-OMO edit: Non-Darcy gas flow
    rhoG_X_mobG_X_k = rhoG.*mobG.*model.G.rock.perm.*turbFactor;
    rhoG_X_mobG_X_k_2face = s.faceUpstr(upco, rhoG_X_mobG_X_k);
    vG = 2*vG./(1 + power((1 - 4*rhoG_X_mobG_X_k_2face.*vG),0.5));
end
if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'nonDarcyEverywhere')
    %HR-OMO edit: Non-Darcy gas flow
    rhoG_X_mobG_X_k = rhoG.*mobG.*model.G.rock.perm.*turbFactor;
    rhoG_X_mobG_X_k_2face = s.faceUpstr(upco, rhoG_X_mobG_X_k);
    vG = 2*vG./(1 + power((1 - 4*rhoG_X_mobG_X_k_2face.*vG),0.5));
end
%HR edit: Incorporate Gangi model by scaling gas velocity
if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'Gangi')
    vG = KGangif.*vG;
end
rOvO = s.faceUpstr(upco, rhoO).*vO;
rGvG = s.faceUpstr(upcg, rhoG).*vG;

pv = model.operators.pv;
pv0 = pv;
if isfield(fluid, 'pvMultR')
    pv = pv.*fluid.pvMultR(p_prop);
    pv0 = pv0.*fluid.pvMultR(p0);
end


% EQUATIONS -----------------------------------------------------------
if model.water
    % Water flux
    if isfield(fluid, 'pcOW')
        pcOW  = fluid.pcOW(sW);
    else
        pcOW = 0;
    end
    pW = p_prop;
    pW0 = p0;
    muW = f.muW(pW);
    bW     = fluid.bW(pW);
    rhoW   = bW.*fluid.rhoWS;
    bW0 = fluid.bW(pW0);

    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./muW;

    dpW    = s.Grad(p - pcOW) - rhoWf.*gdz;
    upcw  = (value(dpW)<=0);
    vW = -s.faceUpstr(upcw, mobW).*T.*dpW;
    if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'nonDarcyinFracs')
        %HR-OMO edit: Non-Darcy water flow
        rhoW_X_mobW_X_k = rhoW.*mobW.*model.G.rock.perm.*turbFactor;
        rhoW_X_mobW_X_k_2face = s.faceUpstr(upco, rhoW_X_mobW_X_k);
        vW = 2*vW./(1 + power((1 - 4*rhoW_X_mobW_X_k_2face.*vW),0.5));
    end
    
    if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'nonDarcyEverywhere')
        %HR-OMO edit: Non-Darcy water flow
        rhoW_X_mobW_X_k = rhoW.*mobW.*model.G.rock.perm.*turbFactor;
        rhoW_X_mobW_X_k_2face = s.faceUpstr(upco, rhoW_X_mobW_X_k);
        vW = 2*vW./(1 + power((1 - 4*rhoW_X_mobW_X_k_2face.*vW),0.5));
    end
    %OMO edit: Incorporate Gangi model by scaling water velocity
    if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'Gangi')
        vW = KGangif.*vW; 
    end
    rWvW = s.faceUpstr(upcw, bW).*vW;
    water = (1/dt).*(pv.*bW.*sW - pv0.*bW0.*sW0);
else
    [vW, mobW, upcw, bW, rhoW] = deal([]);
end

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG);
end

if model.extraStateOutput
    bO = rhoO./fluid.rhoOS;
    bG = rhoG./fluid.rhoGS;
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeViscosities(state, muW,muO, muG); %Hassan's Edit: store phase viscosity
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
    state = model.storeDensities(state, rhoW, rhoO, rhoG);
end

state = model.storeDensities(state, rhoW, rhoO, rhoG);
% water equation + n component equations
[eqs, types, names] = deal(cell(1, ncomp + model.water));


if opt.reduceToPressure
    C = cell(2*ncomp + model.water, 1);
end

fluxes = cell(ncomp, 1);

%HR-OMO edit: modified the accumulation term in the mass balance equation
    %to include the contributions of the adsorption in shale formation.
if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'sorption')
    sumIso = y{1}./model.EOSModel.fluid.isotherm(1,1);
    sumIso0 = y0*(1./model.EOSModel.fluid.isotherm(1,:)');
    for ii=2:numel(y)
        sumIso = sumIso + y{ii}./model.EOSModel.fluid.isotherm(1,ii);
    end
end

compFlux = zeros(model.G.faces.num, ncomp);

for i = 1:ncomp
    names{i} = compFluid.names{i};
    types{i} = 'cell';
    
    if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'sorption')
        eqs{i} = (1/dt).*( ...
            rhoO.*pv.*sO.*xM{i} - rhoO0.*pv0.*sO0.*xM0{i} + ...
            rhoG.*pv.*sG.*yM{i} - rhoG0.*pv0.*sG0.*yM0{i}+...
            (s.gv.*model.EOSModel.fluid.isotherm(2,i)./model.EOSModel.fluid.isotherm(1,i)).*...
            ((y{i}.*p)./(1+p.*sumIso)- (y0(:,i).*p0)./(1+p0.*sumIso0) )  );

    else
        eqs{i} = (1/dt).*( ...
            pv.*rhoO.*sO.*xM{i} - pv0.*rhoO0.*sO0.*xM0{i} + ...
            pv.*rhoG.*sG.*yM{i} - pv0.*rhoG0.*sG0.*yM0{i});
    end
    vi = rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i});
    fluxes{i} = vi;
    compFlux(model.operators.internalConn,i) = value(vi);
    if opt.reduceToPressure
        C{i} = eqs{i};
    end
end
state.componentFluxes = compFlux;
if model.water
    mf = [model.fluid.rhoWS.*value(rWvW), value(rOvO), value(rGvG)];
else
    mf = [value(rOvO), value(rGvG)];
end
state.massFlux = zeros(model.G.faces.num, 2 + model.water);
state.massFlux(model.operators.internalConn, :) = mf;

if model.water
    eqs{ncomp+1} = water;
    names{ncomp+1} = 'water';
    types{ncomp+1} = 'cell';
    C{ncomp+1} = water;
    
    rho = {rhoW, rhoO, rhoG};
    mob = {mobW, mobO, mobG};
    pressures = {pW, p, pG};
else
    rho = {rhoO, rhoG};
    mob = {mobO, mobG};
    pressures = {p, pG};
end
comps = cellfun(@(x, y) {x, y}, xM, yM, 'UniformOutput', false);

[eqs, state, ~] = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                 pressures, sat, mob, rho, ...
                                                 {}, comps, ...
                                                 drivingForces);

% Finally, add in and setup well equations
if opt.staticWells
    compSrc = vertcat(wellSol.components);
    for i = 1:ncomp
        wc = vertcat(W.cells);
        eqs{i+model.water}(wc) = eqs{i+model.water}(wc) - compSrc(:, i);
    end
    nwelleqs = 0;
else
    wellSol = state.wellSol;
    [eqs, names, types, state.wellSol, src] = model.insertWellEquations(eqs, names, ...
                                                     types, wellSol0, wellSol, ...
                                                     wellvars, ...
                                                     wellMap, p, mob, rho, ...
                                                     {}, comps, ...
                                                     dt, opt); %#ok
end

eq_offset = nwelleqs + ncomp + model.water;

Da_V = zeros(1,ncomp); Da_L = zeros(1,ncomp);
%% debugging ECL Diffusion
% persistent n
% if isempty(n)
%     n = 0;
% end
% n = n+1
% test_oil = s.Grad(f_L);
% test_gas = s.Grad(f_V);
% if any(~isfinite(s.Grad(f_L)))
%     disp('Liquid-phase fugacity has a non-finite value(s)!')
% end
%%
for i = 1:ncomp
    eqs{i} = s.AccDiv(eqs{i}, fluxes{i});
    if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'ECLdiffusion')
        RT = 8.314*temp(1);

        saturation_Gas = s.faceUpstr(upcg, sG);
        saturation_Oil = s.faceUpstr(upco, sO);

        density_Gas = s.faceAvg(rhoG);
        density_Oil = s.faceAvg(rhoO);
        
        yM_Di = s.faceAvg(yM{i}.*model.G.rock.Di(i).*model.G.rock.poro./(model.G.rock.tau.*RT)); %use this to make flux in kg/(m^2.s)
        xM_Di = s.faceAvg(xM{i}.*model.G.rock.Di_o(i).*model.G.rock.poro./(model.G.rock.tau.*RT)); %use this to make flux in kg/(m^2.s)

        
        eqs{i} = eqs{i} - s.Div(yM_Di.*s.dd.*saturation_Gas.*density_Gas.*(s.Grad(RT.*log(f_V{i}))-compFluid.molarMass(i).*s.g_gradz))...
                        - s.Div(xM_Di.*s.dd.*saturation_Oil.*density_Oil.*(s.Grad(RT.*log(f_L{i}))-compFluid.molarMass(i).*s.g_gradz));....
%                             Combined chemical potential and gravity effects in kg/(m^2.s)

%         eqs{i} = eqs{i} - s.Div(yM_Di.*s.dd.*saturation_Gas.*density_Gas.*s.Grad(RT.*log(f_V{i})))...
%                         - s.Div(xM_Di.*s.dd.*saturation_Oil.*density_Oil.*s.Grad(RT.*log(f_L{i})));....
%                              chemical potential only in kg/(m^2.s)

    end
    if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'MS_ECL_Diffusion')
        %Hassan edit: implement MS diffusion using NNC diff for dusty-gas
        %model at static mode
        RT = 8.314*temp(1);
        DgXGrad_g = 0; DoXGrad_o = 0;
        
        saturation_Gas = s.faceUpstr(upcg, sG);
        saturation_Oil = s.faceUpstr(upco, sO);
        
        constant_g = s.faceAvg(rhoG.*yM{i}.*model.G.rock.poro./(model.G.rock.tau.*RT));
        constant_g = constant_g .* saturation_Gas;
        
        constant_o = s.faceAvg(rhoO.*xM{i}.*model.G.rock.poro./(model.G.rock.tau.*RT));
        constant_o = constant_o .* saturation_Oil;
        
        column_id = getColumnId_MS_Diff(i,ncomp); %column_id has the same size of ncomp
        for j = 1:ncomp
            DgXGrad_g = DgXGrad_g + Tdiff_gas(:,column_id(j)).*(s.Grad(RT.*log(f_V{j}))-compFluid.molarMass(j).*gdz); 
            DoXGrad_o = DoXGrad_o + Tdiff_oil(:,column_id(j)).*(s.Grad(RT.*log(f_L{j}))-compFluid.molarMass(j).*gdz); 
        end
        eqs{i} = eqs{i} - s.Div(constant_g.*DgXGrad_g) - s.Div(constant_o.*DoXGrad_o);   
    end
    if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'MS_Diffusion')
        %Hassan edit to implement diffusion modeling using Maxwell-Stefan
        %(MS) diffusion matrix and the basic diffusion flux equation at static mode: 
        %J = (-phi.S/tau).pho.D.Grad(X)
        RT = 8.314*temp(1);
        DgXGrad_g = 0; DoXGrad_o = 0;
        Dg = model.G.rock.Dg; Do = model.G.rock.Do;
        saturation_Gas = s.faceUpstr(upcg, sG);
        saturation_Oil = s.faceUpstr(upco, sO);
        
        constant_g = s.faceAvg(rhoG.*yM{i}.*model.G.rock.poro./(model.G.rock.tau.*RT));
        constant_g = constant_g .* saturation_Gas;
        
        constant_o = s.faceAvg(rhoO.*xM{i}.*model.G.rock.poro./(model.G.rock.tau.*RT));
        constant_o = constant_o .* saturation_Oil;
        
        column_id = getColumnId_MS_Diff(i,ncomp); %column_id has the same size of ncomp
        for j = 1:ncomp
            DgXGrad_g = DgXGrad_g + s.faceAvg(Dg(:,column_id(j))).*(s.Grad(RT.*log(f_V{j}))-compFluid.molarMass(j).*gdz); 
            DoXGrad_o = DoXGrad_o + s.faceAvg(Do(:,column_id(j))).*(s.Grad(RT.*log(f_L{j}))-compFluid.molarMass(j).*gdz); 
        end
%         eqs{i} = eqs{i} - s.Div(s.dd.*constant_g.*DgXGrad_g) - s.Div(s.dd.*constant_o.*DoXGrad_o); 
        eqs{i} = eqs{i} - s.Div(constant_g.*DgXGrad_g) - s.Div(constant_o.*DoXGrad_o); 
    end
    if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'dispersion')
        %OMO and Harun edit to add diffusion in gases only:
        %add isMatrix to the different variables
%         Sg_rhoG_poro = s.faceUpstr(upcg,sG.*rhoG.*model.G.rock.poro);
%         Sg_rhoG_poro = s.faceUpstr(upcg,sG.*model.G.rock.poro);
%         Sg_rhoG_poro = s.faceAvg(sG.*model.G.rock.poro);

%         Sg_poro = s.faceUpstr(upcg,sG).*s.faceAvg(model.G.rock.poro);
%         So_poro = s.faceUpstr(upco,sO).*s.faceAvg(model.G.rock.poro);
%         eqs{i} = eqs{i} - s.Div(model.G.rock.Di(i)./model.G.rock.tau.*Sg_poro.*s.Grad(yM{i}.*rhoG)) - s.Div(0.01*model.G.rock.Di(i)./model.G.rock.tau.*So_poro.*s.Grad(xM{i}.*rhoO));
        
%         Sg_rhoG_poro = s.faceUpstr(upcg,sG.*rhoG.*model.G.rock.poro);
        %Sg_rhoG_poro = s.faceUpstr(upcg,sG.*model.G.rock.poro);
%         Sg_rhoG_poro = s.faceUpstr(upcg,sG.*model.G.rock.poro.*s.isMatrix);%Sg_rhoG_poro = s.faceUpstr(upcg,sG.*model.G.rock.poro.*s.isMatrix);
        
        % % Dispersion 
%         Sg_rhoG_poro = s.faceUpstr(upcg,sG.*model.G.rock.poro);
% %         Sg_rhoG_poro = s.faceUpstr(upcg,sG.*model.G.rock.poro);
%         eqs{i} = eqs{i} - s.Div(model.G.rock.Di(i)./model.G.rock.tau.*Sg_rhoG_poro.*s.Grad(yM{i}.*rhoG));
        %calculate fugacity fi at reference pressure
%         [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = model.getMixtureFugacityCoefficients(P, T, x, y, model.fluid.acentricFactors);
% 
% 
%         f_L = EOSModel.computeFugacity(p, x, Z_L, A_L, B_L, Si_L, Bi);
%         f_V = EOSModel.computeFugacity(p, y, Z_V, A_V, B_V, Si_V, Bi);


        Sg_poro = s.faceUpstr(upcg,sG);
        So_poro = s.faceUpstr(upco,sO);
        eqs{i} = eqs{i} - s.Div(model.G.rock.Di(i).*s.faceAvg(model.G.rock.poro)./model.G.rock.tau.*Sg_poro.*s.Grad(yM{i}.*rhoG))- s.Div(model.G.rock.Di_o(i).*s.faceAvg(model.G.rock.poro)./model.G.rock.tau.*So_poro.*s.Grad(xM{i}.*rhoO));
    end
    
    if isfield(model.G.rock, 'shaleMechanisms') && isfield(model.G.rock.shaleMechanisms,'diffusion')
        %OMO and Harun edit to add diffusion in gases only:
%         Sg_rhoG_poro = s.faceUpstr(upcg,rhoG.*sG.*model.G.rock.poro);
        Sg_rhoG_poro = s.faceAvg(rhoG.*sG.*model.G.rock.poro);
        eqs{i} = eqs{i} - s.Div(model.G.rock.Di(i)./model.G.rock.tau.*Sg_rhoG_poro.*s.Grad(yM{i}));
        
    end
    ix = i + eq_offset;
    names{ix}= ['f_', compFluid.names{i}];
    types{ix} = 'fugacity';
    if isfield(model.fluid, 'alpha')
        a = model.fluid.alpha{i}(p, temp, z);
        a = a(twoPhase);
    else
        a = zeros(nnz(twoPhase), 1);
    end
    eqs{ix} = (f_L{i}(twoPhase) - f_V{i}(twoPhase) + a)/barsa;
    
    absent = state.components(twoPhase, i) <= 10*z_tol;
    if model.water
        absent = absent | pureWater(twoPhase);
    end
    if any(absent) && isa(eqs{ix}, 'ADI')
        eqs{ix}.val(absent) = 0;
    end
    
    if opt.reduceToPressure
        C{ix - nwelleqs} = eqs{ix};
    end
end

if model.water
    eqs{ncomp+1} = s.AccDiv(eqs{ncomp+1}, rWvW);
end

cloix = eq_offset + ncomp + 1;

if any(multiPhase)
    eqs{cloix} = sW(multiPhase) + sO(multiPhase) + sG(multiPhase) - 1;
else
    eqs{cloix} = [];
end

types{cloix} = 'saturation';
names{cloix} = 'volclosure';
    
if model.water
    eqs{ncomp+1} = eqs{ncomp+1}.*model.fluid.rhoWS;
end
if opt.reduceToPressure
    C{cloix - nwelleqs} = eqs{cloix};
    
    if model.water
        C{ncomp+1} = C{ncomp+1}.*model.fluid.rhoWS;
    end

    problem = PressureReducedLinearSystem(eqs, types, names, primaryVars, state, dt);
    problem.accumulationTerms = C;
    problem.model = model;
    problem.wellVarIndices = nc*(ncomp+model.water) + (1:nwellvar);
    problem.nwellvar = nwelleqs;
    problem.wellvars = wellvars;
    problem.wellvarNames = wellVarNames;
else
    if model.reduceLinearSystem
        problem = ReducedLinearizedSystem(eqs, types, names, primaryVars, state, dt);
    else
        problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    end
end

if isa(problem, 'ReducedLinearizedSystem')
    % problem.keepNum = model.G.cells.num*(ncomp+model.water);
    problem.keepNum = nc*(ncomp+model.water) + nwellvar;
    problem.reorder = reorder;
end

problem.iterationNo = opt.iteration;
end


function [sO, sG] = setMinimums(model, state, sW, sO, sG, pureVapor, pureLiquid)
    stol = 1e-8;
    if model.water
        sT = sum(state.s, 2);
        if any(pureVapor)
            sG(pureVapor) = sT(pureVapor) - sW(pureVapor);
            if isa(sG, 'ADI')
                sG.val(pureVapor) = max(sG.val(pureVapor), stol);
            else
                sG(pureVapor) = max(sG(pureVapor), stol);
            end
        end

        if any(pureLiquid)
            sO(pureLiquid) = sT(pureLiquid) - sW(pureLiquid);
            if isa(sO, 'ADI')
                sO.val(pureLiquid) = max(sO.val(pureLiquid), stol);
            else
                sO(pureLiquid) = max(sO(pureLiquid), stol);
            end
        end
    end

end

function [c_g, c_o] = molarDensity(compFluid,rhoG,rhoO,ncomp,y,x)
    M_gas = 0; M_oil=0;
    for i = 1:ncomp
        M_gas = M_gas + compFluid.molarMass(i).*y{i};
        M_oil = M_oil + compFluid.molarMass(i).*x{i};
    end
    c_g = rhoG ./ M_gas;
    c_o = rhoO ./ M_oil;
end

%Hassan Edit: When called, getColumnId_MS_Diff function returns column
%indices corresponding to matrix row in MS diffusion matrix. For example, a
%3-component system should have a 3*3 MS diffusion matrix (Assumin we added
%a fourth dummy component), and it will have a Tdiff of size (#faces,9).
%So, this function gets column index for Tdiff corresponing to i-th
%component. For 2nd component, column indices should be [4,5,6], and so on.
function id = getColumnId_MS_Diff(i,ncomp)
    start_id = (i-1)*ncomp + 1;
    end_id = start_id + ncomp-1;
    id = [start_id:end_id];
end
%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
