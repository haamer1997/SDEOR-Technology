function [error,cost,fpath,outFolder] = f_of_x(X,runNumber)
%F_OF_X Summary of this function goes here
%   Detailed explanation goes here
error = 0;
fpath = '/Users/folorode/Documents/MATLAB/Mac_Nash_Code/SingFrac_HM/';
fullPath = [fpath 'Input.dat']; 
fid = fopen(fullPath,'r');
outFolder = ['Output_', num2str(runNumber)];

tline = fgetl(fid);
while ischar(tline)
   eval(tline)
   tline = fgetl(fid);
end

fclose(fid);

% frac.x_f = Fracture_half_length;
frac.x_f = X(5);

if (is_History_Match == 1)
    %DeltaT	    q(Mscf/d)	Pwf(psi)	pwf (psi) Tubing	"Daily Water Rate(M3)"
    histFileName = 'Hist.csv';
    histFileName = fullfile(fpath,histFileName);
    hist = csvread(histFileName);
%     deltaT = hist(:,1);
    Input_Times = cumsum(hist(:,1));
    Input_Times = [0,Input_Times'];
%     qHist = [repmat(hist(1,2),6,1); hist(:,2)];
    qHist = hist(:,2);
    qHist = 0.00032774128.*qHist; %convert from Mscf/d to s.m^3/s
    qHist = qHist./numStages; %divide by number of fractures
end

% hist(:,3) = 600 + hist(:,3);

fracFiles = dir([fpath,'Frac*.mat']);
frac.numFracs = length(fracFiles);

Bound(2:3,1) = X(8);

% load([fpath,'Wells.mat'])
Wells = mean(Bound);
[numWells, ~] = size(Wells);

% load([fpath,'Bound.mat'])
Frac1 = Bound(1:2,1:2);
Frac1(:,1) = Wells(1);
Frac1(1,2) = Wells(2)-frac.x_f;
Frac1(2,2) = Wells(2)+frac.x_f;

numSegs = ones(1,frac.numFracs);
for i=1:frac.numFracs
    str = [fpath,'Frac'];
    str = sprintf('%s%d.mat',str,i);
%     load(str)
    currFrac = eval(sprintf('Frac%d',i));
    [nRows, ~] = size(currFrac);
    numSegs(i) = nRows - 1;
end

fracMeshSize = frac.x_f / 7;
boxMeshSize = fracMeshSize;

meshFileName = [fpath,'caseMesh'];

meshGeoFile = sprintf('%s.geo',meshFileName);
meshID = fopen(meshGeoFile,'w');
fprintf(meshID ,'boundMeshSize = %d;\n',boundMeshSize);
fprintf(meshID ,'boxMeshSize = %d;\n',boxMeshSize);
fprintf(meshID ,'fracMeshSize = %d;\n',fracMeshSize);
fprintf(meshID ,'\n');


[boundDim,~] = size(Bound);
z = 0; %2D
indx = 0;
% Print points corresponding to reservoir boundary
for i=1:boundDim
    indx = indx + 1;
    fprintf(meshID ,'Point(%d) = {%d, %d, %d, %d};\n',indx,Bound(i,1),Bound(i,2),z,boundMeshSize);
end
fprintf(meshID ,'\n');


% Print points corresponding to fractures
for i=1:frac.numFracs
    currFrac = eval(sprintf('Frac%d\n',i));
    currFrac = sortrows(currFrac,-2);
    [currDim,~] = size(currFrac); 
    for j=1:currDim
        if(j==1)
            if (currFrac(2,1)==currFrac(1,1))
                xTop = currFrac(2,1);
                yTop = currFrac(1,2) + 1;
                indx = indx + 1;
                if(yTop>max(Bound(:,2)))
                   fprintf('Frac top tip is outside reservoir boundary');
                   cost = 1e10;
                   return;
                end
                fprintf(meshID ,'Point(%d) = {%d, %d, %d, %d};\n',indx,xTop,yTop,z,1.0);                
            else
                m = (currFrac(2,2)-currFrac(1,2)) / (currFrac(2,1)-currFrac(1,1)); 
                if (m > 0)
                    xTop = currFrac(1,1) + sqrt(1.0/(m*m+1.0));  
                else
                    xTop = currFrac(1,1) - sqrt(1.0/(m*m+1.0));
                end
                yTop = currFrac(1,2) + m*(xTop-currFrac(1,1));    
                indx = indx + 1;
                if(yTop>max(Bound(:,2)))
                   fprintf('Frac top tip is outside reservoir boundary');
                   cost = 1e10;
                   return;
                end
                fprintf(meshID ,'Point(%d) = {%d, %d, %d, %d};\n',indx,xTop,yTop,z,1.0);
            end

        end

        indx = indx + 1;
        fprintf(meshID ,'Point(%d) = {%d, %d, %d, %d};\n',indx,currFrac(j,1),currFrac(j,2),z,fracMeshSize);
        
        if(j==currDim)
            if (currFrac(j,1)==currFrac(j-1,1))
                xBot = currFrac(j,1);
                yBot = currFrac(j,2) - 1;
                indx = indx + 1;
                if(yBot<min(Bound(:,2)))
                   fprintf('Frac bottom tip is outside reservoir boundary');
                   cost = 1e10;
                   return;
                end
                fprintf(meshID ,'Point(%d) = {%d, %d, %d, %d};\n',indx,xBot,yBot,z,1.0);                
            else 
                m = (currFrac(j,2)-currFrac(j-1,2)) / (currFrac(j,1)-currFrac(j-1,1));
                if (m > 0)
                    xBot = currFrac(j,1) - sqrt(1.0/(m*m+1.0));  
                else
                    xBot = currFrac(j,1) + sqrt(1.0/(m*m+1.0));
                end
                yBot = currFrac(j,2) + m*(xBot-currFrac(j,1));
                indx = indx + 1;
                if(yBot<min(Bound(:,2)))
                   fprintf('Frac bottom tip is outside reservoir boundary');
                   cost = 1e10;
                   return;
                end
                fprintf(meshID ,'Point(%d) = {%d, %d, %d, %d};\n',indx,xBot,yBot,z,1.0);
            end
        end
    end
    fprintf(meshID ,'\n');
end
fprintf(meshID ,'\n');

% Print lines delineating the reservoir boundary
indx = 0;  %reset indx to zero for lines
for i=1:(boundDim-1)
    indx = indx + 1;
    fprintf(meshID ,'Line(%d) = {%d, %d};\n',indx,indx,indx+1);
end
indx = indx + 1;
fprintf(meshID ,'Line(%d) = {%d, %d};\n',indx,indx,1);
fprintf(meshID ,'\n');

% Print lines corresponding to fractures
pointIndx = indx + 1;
for i=1:frac.numFracs
    pointIndx = pointIndx + 1;
    currFrac = eval(sprintf('Frac%d\n',i));
    [currDim,~] = size(currFrac); 
    currDim = currDim + 2;
    for j=1:(currDim-1)
        indx = indx + 1;
        fprintf(meshID ,'Line(%d) = {%d, %d};\n',indx,pointIndx-1,pointIndx);
        pointIndx = pointIndx + 1;
    end
    fprintf(meshID ,'\n');
end
fprintf(meshID ,'\n');

surfIndx = 12; %surface index
fprintf(meshID ,'Line Loop(%d) = {%d, %d, %d, %d};\n',10,1,2,3,4);
fprintf(meshID ,'Plane Surface(%d) = {%d};\n',surfIndx,10);
fprintf(meshID ,'\n');

% Print fractures in their respective container surfaces
indx = boundDim;
fracIndxStr = [sprintf('')];
for i=1:frac.numFracs
    currFrac = eval(sprintf('Frac%d\n',i));
    [currDim,~] = size(currFrac);
    currDim = currDim + 2;
    for j=1:(currDim-1)
        indx = indx + 1;
        fprintf(meshID ,'Line {%d} In Surface {%d};\n',indx,surfIndx);
        if (i==1 && j==1)
            fracIndxStr = [fracIndxStr, sprintf('%d',indx)];
        else
            fracIndxStr = [fracIndxStr, sprintf(',%d',indx)];
        end
    end
    fprintf(meshID ,'\n');
end
fprintf(meshID ,'\n');


boundIndxStr = [sprintf('')];
for i=1:boundDim
    if (i==1)
        boundIndxStr = [boundIndxStr, sprintf('%d',i)];
    else
        boundIndxStr = [boundIndxStr, sprintf(',%d',i)];
    end
end

% Print Well Physical Entity
% fprintf(meshID ,'Physical Point(1) = {1,2};\n\n');

% Print Frac Physical Entity
fprintf(meshID ,'BoundID = 200;\n');

str = [sprintf('Physical Line(BoundID) = {'), boundIndxStr, sprintf('}; \n')];
fprintf(meshID, str);

% Print Frac Physical Entity
fprintf(meshID ,'\nfracEntityID = 300;\n');

str = [sprintf('Physical Line(fracEntityID) = {'), fracIndxStr, sprintf('}; \n')];
fprintf(meshID, str);

% Print Surface Physical Entity
fprintf(meshID ,'\nPhysical Surface(400) = {%d};\n',surfIndx);

fclose(meshID);

fpath2 =  fpath(1:end-1); %remove last backslash
idcs   = strfind(fpath2,'\');
newdir = fpath2(1:idcs(end)-1); %directory up one folder
%Running Gmsh
gmshStr = [newdir,'\Gmsh\gmsh -2 -algo ', meshAlgo];
command2Run = sprintf('%s %s', gmshStr, meshGeoFile);
[status, results] = system(command2Run);

if(status)
    cost = 1e10;
    return;
end


% gmshStr = [fpath,'gmsh -2 -algo del2d'];
% gmshStr = [fpath,'gmsh -2 -algo frontal'];
% gmshStr = [fpath,'gmsh -2 -algo meshadapt'];
% Mesh.Algorithm = 5 %(1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad)
% Mesh.Algorithm3D = 1 %(1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)

%read mesh
[mesh.nodes, mesh.elements, mesh.numCVs, mesh.numElems,mesh.boundEls,mesh.fracEls] = readGmesh(sprintf('%s.msh',meshFileName));

%Store values read into input file
% matrix = Matrix(Matrix_porosity*ones(1,numCVs), Matrix_permeability*ones(1,numCVs), Young_Modulus, ...
%     Poisson_ratio, Biot_coefficient, Viscoelastic_Shear_Coefficient, Tortuosity); %call constructor

simManager.inputTimes = Input_Times;
simManager.lastTime = Last_Time;
simManager.isHM = is_History_Match;
simManager.tolConv = convergence_tolerance;
simManager.eps = epsilon;
simManager.numStages = numStages;
simManager.isOptim = is_Optim;

% matrix.poro = Matrix_porosity;
% matrix.perm = Matrix_permeability*ones(1,mesh.numCVs);
% matrix.E = Young_Modulus;
matrix.E = X(9);
matrix.nu = Poisson_ratio;
matrix.Biot = Biot_coefficient;
matrix.tau = Tortuosity;
matrix.tortuInv = 1/Tortuosity;
matrix.h = Reservoir_thickness;

matrix.confineStres = Total_Stress*6894.757293178/1e6;
% matrix.maxStress = Max_Confining_Stress*6894.757293178/1e6;
matrix.alpha = alpha_k;
% matrix.mExp = m_exponent;

% matrix.epsL = Large_Pore_Threshold;
% matrix.phiL = Large_Pore_Threshold*Matrix_porosity*ones(1,mesh.numCVs);

matrix.poro = X(1);
matrix.maxStress = X(3)*6894.757293178/1e6;
matrix.mExp = X(4);
matrix.epsL = X(6);
matrix.perm = X(7)*ones(1,mesh.numCVs);
matrix.phiL = Large_Pore_Threshold*Matrix_porosity*ones(1,mesh.numCVs); 

matrix.kappa_xx = matrix.E/(1.0-matrix.nu*matrix.nu);
matrix.kappa_xy = matrix.E/(2.0*(1.0+matrix.nu));


% if (simManager.isHM)
%     frac.xf = X(5);
% end


if (simManager.isHM)
    well.BHP = 0; 
else
    well.BHP = Bottomhole_pressure;
end

well.skin_m = matrix_skin;
well.skin_f = frac_skin;
well.radius = Well_radius;

mesh.fracMeshSize = fracMeshSize;



state.T = Temperature;
state.p_init = Initial_pressure/1e6; % MPa 
state.zi_init = Initial_mole_fractions;
state.U_init = Initial_displacements;
state.stress_init = Initial_stress;

% fluid.compList = Composition;
% [~,fluid.n_c] = size(fluid.compList); 
% fluid.ncMinus1 = fluid.n_c - 1;
% fluid.isotherm = Isotherm_Parameters';

Constants.RR = 8.3144621;       % in Pa m3 / mol-K
Constants.psi2Pa = 6894.757293178;

chemProps = ChemicalProps(Composition);
[~,fluid.n_c] = size(chemProps.compList); 
fluid.ncMinus1 = fluid.n_c - 1;
fluid.isotherm = Isotherm_Parameters';

[mesh.S, mesh.numNeighs] = createDataStruct(mesh.elements,mesh.numCVs,mesh.numElems);

% n_cMinus1 = fluid.n_c - 1;
simManager.numEquTotal = fluid.n_c+2;
simManager.numEquTotalMin1 = simManager.numEquTotal - 1;

[n, ~] = size(mesh.boundEls);
boundIndx = ones(1,boundDim);
idx = 1;
lenAllBounds = n + boundDim;

mesh.boundCVs = zeros(1,lenAllBounds);
for i=1:n-1
    if (mesh.boundEls(i,1)==mesh.boundEls(i+1,1))
        mesh.boundCVs(i+idx-1) = mesh.boundEls(i,2);
    else
        mesh.boundCVs(i+idx-1) = mesh.boundEls(i,2);
        idx = idx+1;
        mesh.boundCVs(i+idx-1) = mesh.boundEls(i,3);
        boundIndx(idx) = i+idx;
    end
end
mesh.boundCVs(n+idx-1) = mesh.boundEls(n,2);
mesh.boundCVs(n+idx) = mesh.boundEls(n,3);

[n, ~] = size(mesh.fracEls);
mesh.fracIndx = ones(1,frac.numFracs);
idx = 1;
lenAllFracs =  n + frac.numFracs;

fracCVs = zeros(1,lenAllFracs);
for i=1:n-1
%     if (mesh.fracEls(i,1)==mesh.fracEls(i+1,1) || (length(mesh.fracEls((mesh.fracEls(:,1)==mesh.fracEls(i,1)),1)) == 1) )
    if (mesh.fracEls(i,1)==mesh.fracEls(i+1,1) )
        fracCVs(i+idx-1) = mesh.fracEls(i,2);
    else
        fracCVs(i+idx-1) = mesh.fracEls(i,2);
        
        if(mesh.fracEls(i,3) ~= mesh.fracEls(i+1,2))
            idx = idx+1;
            fracCVs(i+idx-1) = mesh.fracEls(i,3);
            mesh.fracIndx(idx) = i+idx;
        end
    end
end

fracCVs(n+idx-1) = mesh.fracEls(n,2);
fracCVs(n+idx) = mesh.fracEls(n,3);
mesh.fracIndx = [mesh.fracIndx length(fracCVs)+1];

% idx = ones(1,frac.numFracs+1);
% for i=2:(frac.numFracs+1)
%     idx(i) = idx(i-1) + numSegs(i-1);
% end
% mesh.fracIndx = mesh.fracIndx(idx);

mesh.sortedFracNodesXYZ = [fracCVs' mesh.nodes(fracCVs,:)];
mesh.indxWinF = zeros(numWells,1);
for i=1:numWells
    dist = sqrt(sum(bsxfun(@minus,mesh.sortedFracNodesXYZ(:,2:3),Wells(i,:)).^2, 2));
    [minDist, mesh.indxWinF(i)] = min(dist);
end
mesh.wellCells = fracCVs(mesh.indxWinF);
mesh.wellIndex = simManager.numEquTotal*(mesh.wellCells-1)+1;  

frac.numCVs = length(fracCVs);
mesh.numMatrixPVs = mesh.numCVs*simManager.numEquTotal;
mesh.numPVs = mesh.numMatrixPVs + frac.numCVs*1;

% zi = Initial_mole_fractions;

frac.poro = Fracture_porosity*ones(frac.numCVs,1);
% frac.perm = Fracture_permeability*ones(frac.numCVs,1);
frac.perm = X(2)*ones(frac.numCVs,1);

frac.permDyn = frac.perm;
frac.E = Fracture_Young_Modulus;
frac.nu = Fracture_Poisson_ratio;
frac.eta = Viscoelastic_Shear_Coefficient; 
frac.width = Fracture_width;
frac.kappa = frac.E/(1.0-frac.nu*frac.nu);

% if (simManager.isHM)
%     matrix.poro = X(1);
%     frac.perm = X(2);
%     matrix.maxStress = X(3);
%     matrix.mExp = X(4);
%     matrix.epsL = X(6);
%     matrix.phiL = Large_Pore_Threshold*Matrix_porosity*ones(1,mesh.numCVs); 
% end

mesh.fracCVs = mesh.sortedFracNodesXYZ(:,1)';

%Initialize Solution vectors
% state.X_guess = zeros(mesh.numPVs,1);
state.X_prev = zeros(mesh.numPVs,1);
state.X_nri = zeros(mesh.numPVs,1);

% pi_prev = zeros(mesh.numCVs,1); 
% z1 = zeros(mesh.numCVs,1);
% z2 = zeros(mesh.numCVs,1);

% p = zeros(mesh.numCVs,1);
% Ux = zeros(mesh.numCVs,1);
% Uy = zeros(mesh.numCVs,1);


idx_f_pG = zeros(frac.numCVs,1);
indxI_m_pG = (1:simManager.numEquTotal:mesh.numMatrixPVs)';
%First initializing every control volume's p to p_init
state.X_prev(indxI_m_pG) = state.p_init;
pi_prev= state.X_prev(1:simManager.numEquTotal:mesh.numMatrixPVs); 

switch fluid.n_c
    case 2
        state.X_prev(indxI_m_pG+1) = state.zi_init(1);
%         z1 = state.X_prev(indxI_m_pG+1);
      
    case 3
        state.X_prev(indxI_m_pG+1) = state.zi_init(1);
        state.X_prev(indxI_m_pG+2) = state.zi_init(2);
%         z1 = state.X_prev(indxI_m_pG+1);
%         z2 = state.X_prev(indxI_m_pG+2);      
        
    case 4
        state.X_prev(indxI_m_pG+1) = state.zi_init(1);
        state.X_prev(indxI_m_pG+2) = state.zi_init(2);
        state.X_prev(indxI_m_pG+3) = state.zi_init(3);
%         z1 = state.X_prev(indxI_m_pG+1);
%         z2 = state.X_prev(indxI_m_pG+2);
%         z3 = state.X_prev(indxI_m_pG+3);
        
    case 5
        state.X_prev(indxI_m_pG+1) = state.zi_init(1);
        state.X_prev(indxI_m_pG+2) = state.zi_init(2);
        state.X_prev(indxI_m_pG+3) = state.zi_init(3);
        state.X_prev(indxI_m_pG+4) = state.zi_init(4); 
%         z1 = state.X_prev(indxI_m_pG+1);
%         z2 = state.X_prev(indxI_m_pG+2);
%         z3 = state.X_prev(indxI_m_pG+3);
%         z4 = state.X_prev(indxI_m_pG+4);
        
    case 6
        state.X_prev(indxI_m_pG+1) = state.zi_init(1);
        state.X_prev(indxI_m_pG+2) = state.zi_init(2);
        state.X_prev(indxI_m_pG+3) = state.zi_init(3);
        state.X_prev(indxI_m_pG+4) = state.zi_init(4); 
        state.X_prev(indxI_m_pG+5) = state.zi_init(5); 
%         z1 = state.X_prev(indxI_m_pG+1);
%         z2 = state.X_prev(indxI_m_pG+2);
%         z3 = state.X_prev(indxI_m_pG+3);
%         z4 = state.X_prev(indxI_m_pG+4);
%         z5 = state.X_prev(indxI_m_pG+5);
        
    otherwise
        disp('Number of components is not between 2 and 5 (inclusive)')
end


for count = 1:frac.numCVs
   idx_f = indxI_m_pG(end)+ simManager.numEquTotalMin1 + count;
   idx_f_pG(count) = idx_f;

   state.X_prev(idx_f)   = state.p_init;                                    
end
%Implementing well pressure in well cell
% state.X_prev(wellCells*numEqu-1) = 2000.0*6894.757293178;
pi_prev(fracCVs) = 2500.0*6894.757293178/1e6;
pi_prev(mesh.wellCells) = 1000.0*6894.757293178/1e6;

% well.wellIndx = simManager.numEquTotal*(mesh.wellCells-1)+1;
state.X_init = state.X_prev; 


frac.numAllFracs = length(mesh.sortedFracNodesXYZ); 
fracMid_x = mesh.sortedFracNodesXYZ(1,2)+0.001;  
matrix.phi_ut = (1.0 - matrix.epsL)*matrix.poro;
matrix.phiL = matrix.epsL*matrix.poro*ones(mesh.numCVs,1);

L_f = zeros(frac.numFracs,1);
frac.Vol = zeros(1,frac.numAllFracs);
frac.Vol(1:frac.numAllFracs) = frac.width*matrix.h*abs(mesh.sortedFracNodesXYZ(5,3)-mesh.sortedFracNodesXYZ(4,3));
frac.topTip = zeros(frac.numFracs,1);
frac.botTip = zeros(frac.numFracs,1);
for zz = 1:frac.numFracs  
    frac.topTip(zz) = mesh.fracIndx(zz);
    frac.botTip(zz) = mesh.fracIndx(zz+1)-1;
    L_f(zz) = abs(mesh.sortedFracNodesXYZ(mesh.fracIndx(zz)+2,3)-mesh.sortedFracNodesXYZ(mesh.fracIndx(zz)+1,3));
    frac.Vol([mesh.fracIndx(zz)+1,mesh.fracIndx(zz+1)-2]) = frac.width*matrix.h*((abs(mesh.sortedFracNodesXYZ(mesh.fracIndx(zz)+1,3)-mesh.sortedFracNodesXYZ(mesh.fracIndx(zz),3)) ...
                                + abs(mesh.sortedFracNodesXYZ(mesh.fracIndx(zz)+2,3)-mesh.sortedFracNodesXYZ(mesh.fracIndx(zz)+1,3)))/2.0);
    frac.Vol([frac.topTip(zz),frac.botTip(zz)])= frac.width*matrix.h*abs(mesh.sortedFracNodesXYZ(mesh.fracIndx(zz)+1,3)-mesh.sortedFracNodesXYZ(mesh.fracIndx(zz),3))/2.0;
end

fluid.Pci = chemProps.Pci./1e6;
fluid.Tci = chemProps.Tci;
fluid.accFact_i = chemProps.AccFacti;
fluid.MWi = chemProps.MWi;
fluid.Vci = chemProps.Vci;
% fluid.Parachor = chemProps.Parachori;  %Parachors are dimensionless constants
fluid.deltaij = chemProps.bic;

%From Taylor and Krishna pg 69
VofC = 15.9;
VofH = 2.31;
%         VofO = 6.11;
VofC1 = VofC + 4.0*VofH;
VofC2 = 2.0*VofC + 6.0*VofH;
%         VofCO2 = VofC + 2.0*VofO;
%         VofCO2 = 26.7;
VofC3 = 3.0*VofC + 8.0*VofH;
VofC4 = 4.0*VofC + 10.0*VofH;
VofC5 = 5.0*VofC + 12.0*VofH;
VofC6 = 6.0*VofC + 14.0*VofH;

switch fluid.n_c
    case 2
        fluid.D12xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(2),VofC1,VofC2);
    case 3
        fluid.D12xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(2),VofC1,VofC2);
        fluid.D13xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(3),VofC1,VofC3);
        fluid.D23xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(3),VofC2,VofC3);
        
        fluid.DInit = computeD(fluid.D12xP,fluid.D13xP,fluid.D23xP,state.p_init*1e6,...
                         fluid.ncMinus1,state.T,fluid.Pci*1e6,fluid.Tci,fluid.accFact_i, ...
                         state.zi_init, fluid.deltaij,Constants.RR);
    case 4
        fluid.D12xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(2),VofC1,VofC2);
        fluid.D13xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(3),VofC1,VofC3);
        fluid.D14xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(4),VofC1,VofC4);
        fluid.D23xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(3),VofC2,VofC3);
        fluid.D24xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(4),VofC2,VofC4);
        fluid.D34xP = computeGasDijxP(state.T,fluid.MWi(3),fluid.MWi(4),VofC3,VofC4);
 
        fluid.DInit = computeD4comps_nc(fluid.D12xP,fluid.D13xP,fluid.D14xP,...
                         fluid.D23xP,fluid.D24xP,fluid.D34xP,...
                         state.p_init*1e6,fluid.ncMinus1,state.T,fluid.Pci*1e6,...
                         fluid.Tci,fluid.accFact_i,state.zi_init,fluid.deltaij,Constants.RR);
    case 5
        fluid.D12xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(2),VofC1,VofC2);
        fluid.D13xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(3),VofC1,VofC3);
        fluid.D14xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(4),VofC1,VofC4);
        fluid.D15xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(5),VofC1,VofC5);
        fluid.D23xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(3),VofC2,VofC3);
        fluid.D24xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(4),VofC2,VofC4);
        fluid.D25xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(5),VofC2,VofC5);
        fluid.D34xP = computeGasDijxP(state.T,fluid.MWi(3),fluid.MWi(4),VofC3,VofC4);
        fluid.D35xP = computeGasDijxP(state.T,fluid.MWi(3),fluid.MWi(5),VofC3,VofC5);
        fluid.D45xP = computeGasDijxP(state.T,fluid.MWi(4),fluid.MWi(5),VofC4,VofC5);
 
        fluid.DInit = computeD5comps_nc(fluid.D12xP,fluid.D13xP,fluid.D14xP,fluid.D15xP,...
                         fluid.D23xP,fluid.D24xP,fluid.D25xP,fluid.D34xP,fluid.D35xP,...
                         fluid.D45xP,state.p_init*1e6,fluid.ncMinus1,state.T,fluid.Pci*1e6,...
                         fluid.Tci,fluid.accFact_i,state.zi_init,fluid.deltaij,Constants.RR);
    case 6
        fluid.D12xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(2),VofC1,VofC2);
        fluid.D13xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(3),VofC1,VofC3);
        fluid.D14xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(4),VofC1,VofC4);
        fluid.D15xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(5),VofC1,VofC5);
        fluid.D16xP = computeGasDijxP(state.T,fluid.MWi(1),fluid.MWi(6),VofC1,VofC6);
        fluid.D23xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(3),VofC2,VofC3);
        fluid.D24xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(4),VofC2,VofC4);
        fluid.D25xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(5),VofC2,VofC5);
        fluid.D26xP = computeGasDijxP(state.T,fluid.MWi(2),fluid.MWi(6),VofC2,VofC6);
        fluid.D34xP = computeGasDijxP(state.T,fluid.MWi(3),fluid.MWi(4),VofC3,VofC4);
        fluid.D35xP = computeGasDijxP(state.T,fluid.MWi(3),fluid.MWi(5),VofC3,VofC5);
        fluid.D36xP = computeGasDijxP(state.T,fluid.MWi(3),fluid.MWi(6),VofC3,VofC6);
        fluid.D45xP = computeGasDijxP(state.T,fluid.MWi(4),fluid.MWi(5),VofC4,VofC5);
        fluid.D46xP = computeGasDijxP(state.T,fluid.MWi(4),fluid.MWi(6),VofC4,VofC6);
        fluid.D56xP = computeGasDijxP(state.T,fluid.MWi(5),fluid.MWi(6),VofC5,VofC6);
 
        fluid.DInit = computeD6comps_nc(fluid,state,Constants);                     
                     
    otherwise
        disp('Number of components is not between 2 and 6 (inclusive)')
end


fluid.DInit = diag(diag(fluid.DInit));

Constants.nOverV_sc = 101325.0/(Constants.RR*288.706); %g-mol/s.m^3

% Z_prev  = zeros(mesh.numCVs,1);
% c_prev  = zeros(mesh.numCVs,1);       
mesh.volCV  = zeros(mesh.numCVs,1);
%Compute geometrical part of flux computation for all cells
maxNeighs = max(mesh.numNeighs); 
GeomProps.geomf1_4_a1 = zeros(mesh.numCVs,maxNeighs);
GeomProps.geomf2_4_a1 = zeros(mesh.numCVs,maxNeighs);
GeomProps.geomf1_4_a2 = zeros(mesh.numCVs,maxNeighs);
GeomProps.geomf2_4_a2 = zeros(mesh.numCVs,maxNeighs);
GeomProps.geomf1_4_a3 = zeros(mesh.numCVs,maxNeighs);
GeomProps.geomf2_4_a3 = zeros(mesh.numCVs,maxNeighs);

GeomProps.a1x  = zeros(mesh.numCVs,maxNeighs); GeomProps.a1y  = zeros(mesh.numCVs,maxNeighs); 
GeomProps.a2x  = zeros(mesh.numCVs,maxNeighs); GeomProps.a2y  = zeros(mesh.numCVs,maxNeighs); 
GeomProps.a3x  = zeros(mesh.numCVs,maxNeighs); GeomProps.a3y  = zeros(mesh.numCVs,maxNeighs); 
GeomProps.b1x  = zeros(mesh.numCVs,maxNeighs); GeomProps.b1y  = zeros(mesh.numCVs,maxNeighs); 
GeomProps.b2x  = zeros(mesh.numCVs,maxNeighs); GeomProps.b2y  = zeros(mesh.numCVs,maxNeighs); 
GeomProps.b3x  = zeros(mesh.numCVs,maxNeighs); GeomProps.b3y  = zeros(mesh.numCVs,maxNeighs); 
GeomProps.Dphi1Dx = zeros(mesh.numCVs,maxNeighs); GeomProps.Dphi1Dy  = zeros(mesh.numCVs,maxNeighs); 
GeomProps.Dphi2Dx = zeros(mesh.numCVs,maxNeighs); GeomProps.Dphi2Dy  = zeros(mesh.numCVs,maxNeighs); 
GeomProps.Dphi3Dx = zeros(mesh.numCVs,maxNeighs); GeomProps.Dphi3Dy  = zeros(mesh.numCVs,maxNeighs);
mesh.sCVs = zeros(mesh.numCVs,maxNeighs);

fIdx = 0;
c_subL = zeros(frac.numAllFracs,1);
c_subR = zeros(frac.numAllFracs,1);
sCVL = zeros(frac.numAllFracs,1);
sCVR = zeros(frac.numAllFracs,1);
iFrac = zeros(frac.numAllFracs,1);
for i=1:mesh.numCVs
    if any(mesh.fracCVs==i)
        fIdx = fIdx + 1;
        iFrac(fIdx) = i;
    end
    %This is a loop through all the neighbors of a cell, m
    for j=1:mesh.numNeighs(i)
        x1 = mesh.nodes(i,1);          y1 = mesh.nodes(i,2);
        x2 = mesh.nodes(mesh.S(i,j),1);     y2 = mesh.nodes(mesh.S(i,j),2);
        x3 = mesh.nodes(mesh.S(i,j+1),1);   y3 = mesh.nodes(mesh.S(i,j+1),2);

        mesh.sCVs(i,j) = (matrix.h/2.0*((x2*y3-x3*y2)-x1*(y3-y2)+y1*(x3-x2)))/3;

        area123_X_2 = 6.0*mesh.sCVs(i,j)/matrix.h;    % 2*3=6
        mesh.volCV(i) = mesh.volCV(i) + mesh.sCVs(i,j);               

        deltaXf1 = matrix.h*(2.0*x3-x2-x1)/6.0;
        deltaYf1 = matrix.h*(2.0*y3-y2-y1)/6.0;
        deltaXf2 = matrix.h*(x1-2.0*x2+x3)/6.0;
        deltaYf2 = matrix.h*(y1-2.0*y2+y3)/6.0;

        GeomProps.Dphi1Dx(i,j) = (y2-y3)/area123_X_2;
        GeomProps.Dphi1Dy(i,j) = (x3-x2)/area123_X_2;
        GeomProps.Dphi2Dx(i,j) = (y3-y1)/area123_X_2;
        GeomProps.Dphi2Dy(i,j) = (x1-x3)/area123_X_2;
        GeomProps.Dphi3Dx(i,j) = (y1-y2)/area123_X_2;
        GeomProps.Dphi3Dy(i,j) = (x2-x1)/area123_X_2;

        %Computing and storing the geometric terms for each face & "a" term
        GeomProps.geomf1_4_a1(i,j) = (GeomProps.Dphi1Dx(i,j)*deltaYf1 - GeomProps.Dphi1Dy(i,j)*deltaXf1);
        GeomProps.geomf2_4_a1(i,j) = (GeomProps.Dphi1Dx(i,j)*deltaYf2 - GeomProps.Dphi1Dy(i,j)*deltaXf2);
        GeomProps.geomf1_4_a2(i,j) = (GeomProps.Dphi2Dx(i,j)*deltaYf1 - GeomProps.Dphi2Dy(i,j)*deltaXf1);
        GeomProps.geomf2_4_a2(i,j) = (GeomProps.Dphi2Dx(i,j)*deltaYf2 - GeomProps.Dphi2Dy(i,j)*deltaXf2);
        GeomProps.geomf1_4_a3(i,j) = (GeomProps.Dphi3Dx(i,j)*deltaYf1 - GeomProps.Dphi3Dy(i,j)*deltaXf1);
        GeomProps.geomf2_4_a3(i,j) = (GeomProps.Dphi3Dx(i,j)*deltaYf2 - GeomProps.Dphi3Dy(i,j)*deltaXf2);

        %Constant parameters for the Ux equation
        GeomProps.a1x(i,j) = matrix.kappa_xx*GeomProps.Dphi1Dx(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.kappa_xy*GeomProps.Dphi1Dy(i,j)*(deltaXf1 + deltaXf2);
        GeomProps.a2x(i,j) = matrix.kappa_xx*GeomProps.Dphi2Dx(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.kappa_xy*GeomProps.Dphi2Dy(i,j)*(deltaXf1 + deltaXf2);
        GeomProps.a3x(i,j) = matrix.kappa_xx*GeomProps.Dphi3Dx(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.kappa_xy*GeomProps.Dphi3Dy(i,j)*(deltaXf1 + deltaXf2); 
        GeomProps.b1x(i,j) = matrix.nu*matrix.kappa_xx*GeomProps.Dphi1Dy(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.kappa_xy*GeomProps.Dphi1Dx(i,j)*(deltaXf1 + deltaXf2);
        GeomProps.b2x(i,j) = matrix.nu*matrix.kappa_xx*GeomProps.Dphi2Dy(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.kappa_xy*GeomProps.Dphi2Dx(i,j)*(deltaXf1 + deltaXf2);
        GeomProps.b3x(i,j) = matrix.nu*matrix.kappa_xx*GeomProps.Dphi3Dy(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.kappa_xy*GeomProps.Dphi3Dx(i,j)*(deltaXf1 + deltaXf2);    
        %Constant parameters for the Uy equation
        GeomProps.a1y(i,j) = matrix.kappa_xy*GeomProps.Dphi1Dx(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.kappa_xx*GeomProps.Dphi1Dy(i,j)*(deltaXf1 + deltaXf2);
        GeomProps.a2y(i,j) = matrix.kappa_xy*GeomProps.Dphi2Dx(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.kappa_xx*GeomProps.Dphi2Dy(i,j)*(deltaXf1 + deltaXf2);
        GeomProps.a3y(i,j) = matrix.kappa_xy*GeomProps.Dphi3Dx(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.kappa_xx*GeomProps.Dphi3Dy(i,j)*(deltaXf1 + deltaXf2); 
        GeomProps.b1y(i,j) = matrix.kappa_xy*GeomProps.Dphi1Dy(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.nu*matrix.kappa_xx*GeomProps.Dphi1Dx(i,j)*(deltaXf1 + deltaXf2);
        GeomProps.b2y(i,j) = matrix.kappa_xy*GeomProps.Dphi2Dy(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.nu*matrix.kappa_xx*GeomProps.Dphi2Dx(i,j)*(deltaXf1 + deltaXf2);
        GeomProps.b3y(i,j) = matrix.kappa_xy*GeomProps.Dphi3Dy(i,j)*(deltaYf1 + deltaYf2) - ...
                   matrix.nu*matrix.kappa_xx*GeomProps.Dphi3Dx(i,j)*(deltaXf1 + deltaXf2); 
        if any(mesh.fracCVs==i)                     
            c_i = (x1 + x2 + x3)/3.0;
            m1 =  (x1 + x2)/2.0;
            m2 =  (x1 + x3)/2.0;
            if any([x1,x2,x3]>fracMid_x)
                c_subR(fIdx) = c_subR(fIdx) + mesh.sCVs(i,j)*(x1+m1+c_i + x1+c_i+m2)/6.0;
                sCVR(fIdx) = sCVR(fIdx) + mesh.sCVs(i,j);
            else
                c_subL(fIdx) = c_subL(fIdx) + mesh.sCVs(i,j)*(x1+m1+c_i + x1+c_i+m2)/6.0;   
                sCVL(fIdx) = sCVL(fIdx) + mesh.sCVs(i,j);
            end
        end
    end
end 

c_subL = c_subL ./ sCVL; 
c_subR = c_subR ./ sCVR; 

well.ro_f = 0.14*sqrt(L_f.*L_f + matrix.h*matrix.h);
well.ro_m = sqrt((mesh.volCV(mesh.wellCells)/matrix.h)/pi);

[~,fracSortIdx]=ismember(mesh.fracCVs,iFrac);
c_subL = c_subL(fracSortIdx);
c_subR = c_subR(fracSortIdx);

frac.d_nnc = (sCVL(fracSortIdx).*abs(mesh.sortedFracNodesXYZ(:,2) - c_subL) + ...
    sCVR(fracSortIdx).*abs(c_subR - mesh.sortedFracNodesXYZ(:,2))) ./ mesh.volCV(mesh.fracCVs);
L_f =  mesh.fracMeshSize.*ones(frac.numAllFracs,1);
L_f([1 frac.numAllFracs]) = 0.5;
L_f([2 (frac.numAllFracs-1)]) = (1.0+mesh.fracMeshSize)/2;
frac.A_nnc = matrix.h.*L_f; %just temporary approx (not accounting for shorter frac segment at tip)

% mesh.numMatrixPVs = simManager.numEquTotal*mesh.numCVs;
















lastTime = Last_Time;
Input_Times = Input_Times*86400;

lengthUserTime = length(Input_Times);

p_results = zeros(mesh.numCVs,lengthUserTime);
pf_results = zeros(lenAllFracs,lengthUserTime);
% z1_results = zeros(mesh.numCVs,lengthUserTime);
% z2_results = zeros(mesh.numCVs,lengthUserTime);
switch fluid.n_c
    case 2
        z1_results = zeros(mesh.numCVs,lengthUserTime);
        
    case 3
        z1_results = zeros(mesh.numCVs,lengthUserTime);
        z2_results = zeros(mesh.numCVs,lengthUserTime);
        
    case 4
        z1_results = zeros(mesh.numCVs,lengthUserTime);
        z2_results = zeros(mesh.numCVs,lengthUserTime);
        z3_results = zeros(mesh.numCVs,lengthUserTime);
        
    case 5
        z1_results = zeros(mesh.numCVs,lengthUserTime);
        z2_results = zeros(mesh.numCVs,lengthUserTime);
        z3_results = zeros(mesh.numCVs,lengthUserTime);
        z4_results = zeros(mesh.numCVs,lengthUserTime);
        
    case 6
        z1_results = zeros(mesh.numCVs,lengthUserTime);
        z2_results = zeros(mesh.numCVs,lengthUserTime);
        z3_results = zeros(mesh.numCVs,lengthUserTime);
        z4_results = zeros(mesh.numCVs,lengthUserTime);  
        z5_results = zeros(mesh.numCVs,lengthUserTime);  
    otherwise
        disp('Number of components is not between 2 and 5 (inclusive)')
end

qg_results = zeros(1,lengthUserTime);
pwf_results = zeros(lengthUserTime,1);
Ux_results = zeros(mesh.numCVs,lengthUserTime);
Uy_results = zeros(mesh.numCVs,lengthUserTime);

sigma_xx_results = zeros(mesh.numCVs,lengthUserTime);
sigma_xxFrac_results = zeros(length(fracCVs),lengthUserTime); 
sigma_xxFracCum_results  = zeros(length(fracCVs),lengthUserTime); 
sigma_xxFracCum = zeros(length(fracCVs),1); 
wf_results = zeros(length(fracCVs),lengthUserTime); 
kfwf_results = zeros(length(fracCVs),lengthUserTime); 


simManager.maxNRI = 12;
simManager.maxNumTimestepCuts = 5;

simManager.maxTimeSteps = 40;
currTime = 0.0;
simManager.deltaT = Input_Times(2) - Input_Times(1);
if (simManager.isHM)
   qHistTS = qHist(1);
else
   qHistTS = 0;
end
 
timeStep = 1;
timeArray = zeros(100,1);
usertimeIndx = 1;

%%
% phiPrev = phiI';
sigma_xxFrac_prev = zeros(length(fracCVs),1);
simManager.isFirstTime = 1;
simManager.isX_prevUpdated = 1;
%timestep loop begins here
for nt = 2:simManager.maxTimeSteps   %timeStep = 2:length(time)
   prevTime = currTime;
      
   if ((prevTime >= Input_Times(usertimeIndx))&& usertimeIndx<lengthUserTime)
       simManager.deltaT = Input_Times(usertimeIndx+1) - prevTime;
       if (simManager.isHM)
           qHistTS = qHist(usertimeIndx); 
       else
           qHistTS = 0;
       end
       usertimeIndx = usertimeIndx+1;     
   end
   
   if  (currTime > lastTime*86400)      
       break;
   end
   
   %Just an initial guess
   state.X_guess = state.X_prev;
   
   simManager.isX_prevUpdated = 1;
   %Adaptive Timestep cutting loop
   for timeStepCut = 1:simManager.maxNumTimestepCuts
      state.X_nri = state.X_guess;     
      for nri = 1:simManager.maxNRI
         if (any(imag(state.X_nri)) || any(isnan(state.X_nri)) )  
             fprintf('One of the primary variables is complex or NaN\n');
             cost = 1e10;
             error = 5;
             return;
         end

         perm = matrix.perm.*(1.0-((matrix.confineStres-matrix.alpha*(state.X_nri(1:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs))')/matrix.maxStress).^matrix.mExp).^3;
         if (any(imag(perm)) )
             fprintf('Perm is complex or NaN\n');
             cost = 1e10;
             error = 3;
             return;
         end

 
         [error,DX,myNorm,rate_sc,sigma_xx,wf_curr,kfwf,BHP] = jacobian(state.X_nri, sigma_xxFrac_prev, ...
            sigma_xxFracCum,nt,qHistTS,matrix,frac,well,simManager,state,fluid,Constants,mesh,GeomProps);
         if(error)
             cost = 1e10;
             return;
         end        
                   
         if simManager.isFirstTime
              p_results(:,timeStep) = state.X_nri(1:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
              switch fluid.n_c
                case 2
                    z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    
                case 3
                    z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z2_results(:,timeStep) = state.X_nri(3:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
              
                case 4
                    z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z2_results(:,timeStep) = state.X_nri(3:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z3_results(:,timeStep) = state.X_nri(4:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    
                case 5
                    z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z2_results(:,timeStep) = state.X_nri(3:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z3_results(:,timeStep) = state.X_nri(4:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z4_results(:,timeStep) = state.X_nri(5:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    
                case 6
                    z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z2_results(:,timeStep) = state.X_nri(3:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z3_results(:,timeStep) = state.X_nri(4:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z4_results(:,timeStep) = state.X_nri(5:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                    z5_results(:,timeStep) = state.X_nri(6:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
                otherwise
                    disp('Number of components is not between 2 and 6 (inclusive)')
              end
              Ux_results(:,timeStep) =  state.X_nri(simManager.numEquTotalMin1:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
              Uy_results(:,timeStep) =  state.X_nri(simManager.numEquTotal:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);

              pf_results(:,timeStep) = state.X_nri(simManager.numEquTotal*mesh.numCVs+1:1:end);
              qg_results(timeStep) = rate_sc(1)*3.051187204736614e+03; %convert from s.m^3/s to Mscf/d
              pwf_results(timeStep) = BHP*1e6/Constants.psi2Pa;
         end          
                   
         simManager.isFirstTime = 0;   
         simManager.isX_prevUpdated = 0;
         state.X_nri = state.X_nri + DX;
         
%          %Forcing zi's to be btw 0 and 1 (just made it <1 for now)
         if(~isempty(find(state.X_nri>1 & state.X_nri<3, 1)))
             fprintf('Unphysical Solution..\n');             
             switch fluid.n_c
                case 2
                    state.X_nri(2:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(2:simManager.numEquTotal:mesh.numMatrixPVs);
                    
                case 3
                    state.X_nri(2:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(2:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(3:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(3:simManager.numEquTotal:mesh.numMatrixPVs);
             
                case 4
                    state.X_nri(2:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(2:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(3:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(3:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(4:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(4:simManager.numEquTotal:mesh.numMatrixPVs);
                    
                case 5
                    state.X_nri(2:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(2:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(3:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(3:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(4:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(4:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(5:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(5:simManager.numEquTotal:mesh.numMatrixPVs);
                case 6
                    state.X_nri(2:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(2:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(3:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(3:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(4:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(4:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(5:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(5:simManager.numEquTotal:mesh.numMatrixPVs);
                    state.X_nri(6:simManager.numEquTotal:mesh.numMatrixPVs) = state.X_init(6:simManager.numEquTotal:mesh.numMatrixPVs);
                    
                otherwise
                    disp('Number of components is not between 2 and 6 (inclusive)')
             end
             if(timeStep<7)
                cost = 1e10;
                return;
             end
         end

         if(myNorm <= simManager.tolConv )          
            break; 
         end 
         %check if solution is diverging instead of converging
         if(nri > 2 && myNorm/myNorm_prev > 0.85)
            fprintf('System is diverging! Cutting timestep..\n');
            if(timeStep<7)
               cost = 1e10;
               return;
            end
            break;
         end
         myNorm_prev = myNorm;
      end
      
%       if (any(imag(state.X_nri)) || any(isnan(state.X_nri)) )  
%          fprintf('One of the primary variables is complex or NaN\n');
%          cost = 1e10;
%          return;
%       end
      %Exit NRI if there is convergence
      if(myNorm <= simManager.tolConv)
         break; 
      end

      %Cut timestep size if there is no convergence
      if(myNorm > simManager.tolConv)
         simManager.deltaT = simManager.deltaT/4.0;
         fprintf('NRI did not converge! Cutting timestep to: %6.2f secs..\n',simManager.deltaT);
      end
      
      
   end
   
   if(myNorm > simManager.tolConv)
      fprintf('Execution Failed!!!\nNRI did not converge after %d timestep cut-backs!!!..\n',simManager.maxNumTimestepCuts);
      break; 
   end
   
   %storing results at user specified frequencies
   if (any(abs(currTime-Input_Times)<1e-10))    %(mod(currTime,30) == 0)
      timeStep = timeStep + 1;   
      
      p_results(:,timeStep) = state.X_nri(1:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
      switch fluid.n_c
         case 2
            z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);

         case 3
            z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z2_results(:,timeStep) = state.X_nri(3:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);

         case 4
            z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z2_results(:,timeStep) = state.X_nri(3:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z3_results(:,timeStep) = state.X_nri(4:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);

         case 5
            z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z2_results(:,timeStep) = state.X_nri(3:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z3_results(:,timeStep) = state.X_nri(4:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z4_results(:,timeStep) = state.X_nri(5:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);

         case 6
            z1_results(:,timeStep) = state.X_nri(2:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z2_results(:,timeStep) = state.X_nri(3:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z3_results(:,timeStep) = state.X_nri(4:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z4_results(:,timeStep) = state.X_nri(5:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
            z5_results(:,timeStep) = state.X_nri(6:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);

         otherwise
            disp('Number of components is not between 2 and 6 (inclusive)')
      end
      Ux_results(:,timeStep) =  state.X_nri(simManager.numEquTotalMin1:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
      Uy_results(:,timeStep) =  state.X_nri(simManager.numEquTotal:simManager.numEquTotal:simManager.numEquTotal*mesh.numCVs);
      
      pf_results(:,timeStep) = state.X_nri(simManager.numEquTotal*mesh.numCVs+1:1:end);
      
      qg_results(timeStep) = rate_sc(1)*86400.0/(0.3048^3)/1000.0; %Mscf/d
      pwf_results(timeStep) = BHP*1e6/Constants.psi2Pa; %psi
   end
   
   %update solution and advance in time
   state.X_prev = state.X_nri;
   
%    if(any(wf_curr<=0))
%        fprintf('fracture width is negative!!! Terminating execution...\n');
%        break;
%    end
   
   wf_results(:,timeStep) = wf_curr; 
   kfwf_results(:,timeStep) = kfwf; 
   
   sigma_xxFrac = sigma_xx(fracCVs);
   sigma_xxFrac_results(:,timeStep) = sigma_xxFrac;
   sigma_xx_results(:,timeStep) = sigma_xx;
   if(nt == 2)
       sigma_xxFracCum = sigma_xxFracCum + 0.5*simManager.deltaT*sigma_xxFrac;
   else
       sigma_xxFracCum = sigma_xxFracCum + 0.5*simManager.deltaT*(sigma_xxFrac+sigma_xxFrac_prev);
   end
   sigma_xxFracCum = state.stress_init*simManager.deltaT + sigma_xxFracCum;
   
   sigma_xxFrac_prev = sigma_xxFrac;
   sigma_xxFracCum_results(:,timeStep) = sigma_xxFracCum;
   
   
   fprintf('Current time = %6.2f days \n', currTime/86400);
   str = sprintf('Current time = %6.2f days', currTime/86400);
   
   currTime = currTime + simManager.deltaT;
   timeArray(nt) = currTime;
   
   simManager.isX_prevUpdated = 1;
   
%    %Double Timestep size if NRI converges too quickly
%    if(nri<4)
%       simManager.deltaT = simManager.deltaT * 2.0;
%       fprintf('Converged in %2.0f iterations... Doubling timestep to %6.2f secs \n',nri, simManager.deltaT);
%    end
   
end  

%%
  %Compute total concentrations from component concentrations
 
% for count = 1:mesh.numCVs
%   indxI = indxI_m_pG(count);
%   p(count)   = state.X_nri(indxI);
%   z1(count)  = state.X_nri(indxI+1);
%   z2(count)  = state.X_nri(indxI+2);
%   Ux(count)  = state.X_nri(indxI+3);
%   Uy(count)  = state.X_nri(indxI+4);
% end

tinDays = Input_Times(1:timeStep-1)' ./ 86400;   
qg_results = qg_results(1:timeStep-1);
qg_results = qg_results';
prod = 0.5*diff(tinDays).*(2*qg_results(2:end) - diff(qg_results));
cumProd = [0; cumsum(prod)]./1000;

production = [tinDays, qg_results, cumProd];
fprintf('Cumulative production is %2.0f MMscf \n',cumProd(end));
p_results = p_results(:,1:timeStep);
pf_results = pf_results(:,1:timeStep);
z1_results = z1_results(:,1:timeStep);
z2_results = z2_results(:,1:timeStep);
Ux_results = Ux_results(:,1:timeStep);
Uy_results = Uy_results(:,1:timeStep);

if(~simManager.isOptim)
    mkdir(fpath, outFolder)
    % save('Output/CaseSummary.txt','CaseTitle','CaseDescription', '-ascii') 
    save([fpath,outFolder,'\p_results.txt'],'p_results','-ascii') 
    save([fpath,outFolder,'\z1_results.txt'],'z1_results','-ascii') 
    save([fpath,outFolder,'\z2_results.txt'],'z2_results','-ascii') 
    save([fpath,outFolder,'\pf_results.txt'],'pf_results','-ascii') 
    save([fpath,outFolder,'\Ux_results.txt'],'Ux_results','-ascii') 
    save([fpath,outFolder,'\Uy_results.txt'],'Uy_results','-ascii') 
    save([fpath,outFolder,'\production.txt'],'production','-ascii')

    csvwrite([fpath,'nodes.csv'],mesh.nodes) 
    csvwrite([fpath,outFolder,'\production.csv'],production) 
    csvwrite([fpath,outFolder,'\pwf_results.csv'],pwf_results) 
end


if(simManager.isHM)
    if(currTime/86400<simManager.lastTime)
        cost = 1.5*norm(pwf_results(2:timeStep-1)-hist(1:timeStep-2,3))/(timeStep-2) 
%       cost = 1e10;
        return
    else
        cost = norm(pwf_results(2:timeStep-1)-hist(1:timeStep-2,3))/(timeStep-2) 
    end
end

if(~simManager.isOptim)
    figure(10)
        plot(tinDays(2:end),pwf_results(2:timeStep-1),tinDays(2:end),hist(1:timeStep-2,3),'*')
        axis([0 max(tinDays) 0 6000])
end

end

