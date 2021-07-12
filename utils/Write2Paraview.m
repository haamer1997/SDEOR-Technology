function Write2Paraview(fname,G,states,varargin)
opt = struct('NumHFs', -1,...
             'NumNFs',  -1,...
             'temp',-1);
opt = merge_options(opt, varargin{:});

[filepath,name,ext] = fileparts(fname);

%Folder Name
mkdir('ParaviewResults_BaseCase')

%test

%Write FracPlane to Paraview
if isfield(opt, 'fracplanes') 
    fname=fullfile('ParaviewResults',[name '_Frac.vtk'])
    Frac2vtk(fname,fracplanes,opt.NumHFs);
end
%Write Wells to Paraview
if isfield(opt, 'wells')
    fname=fullfile('ParaviewResults',[name '_Well.vtk'])
    well2vtk(fname,wells);
end
%Write Solution into Paraview
for i=1:size(states,1)
    fname=fullfile('ParaviewResults',[name '_Field_Step' num2str(i) '.vtk']);
    grid2vtk(fname,G,...
        {'Pressure','Saturation','Components'},...
        {states{i}.pressure,states{i}.s,states{i}.components});
end
end


function Frac2vtk(fname,fracplanes,NumHFs)
%GRID2VTK Summary of this function goes here
%Detailed explanation goes here

NumFracs=size(fracplanes,2);

%Extract all points
pts=[];
for fi=1:NumFracs
    pts=[pts; fracplanes(fi).points];
end

NumPts=size(pts,1);

%Find the poly node idx
FracPolyId={};
id=0;
for fi=1:NumFracs
    NumFracPts=size(fracplanes(fi).points,1);
    id_ST=id;
    id_END=id+NumFracPts-1;
    FracPolyId{end+1}=[NumFracPts id_ST:id_END];
    id=id_END+1;
end

%Frac Plane props
frac_aperture=[fracplanes(1:NumFracs).aperture];
frac_perm=[fracplanes(1:NumFracs).perm];
frac_poro=[fracplanes(1:NumFracs).poro];

%Hydraulic Fractures tag=0, Natural Fracs=1
if(NumHFs==-1)
    FracTags=1:NumFracs;
else
    FracTags=ones(NumFracs,1);
    FracTags(1:NumHFs)=0;
end

fid = fopen(fname, 'w'); 
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Fracture Geometry\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');
fprintf(fid, 'POINTS %d float\n', NumPts);
fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', pts');fprintf(fid, '\n');
fprintf(fid, 'POLYGONS %d %d\n', NumFracs,NumFracs+size(pts,1));
for fi=1:NumFracs
    format=getPrintRowFormat(numel(FracPolyId{fi}),'%d');
    fprintf(fid, format, FracPolyId{fi}');
end
fprintf(fid, '\n');

fprintf(fid, 'CELL_DATA %d\n', NumFracs);
fprintf(fid, '\n');

fprintf(fid, ['SCALARS ', 'FractureTag', ' int\n']);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid, '%d %d %d %d %d %d\n', FracTags');
fprintf(fid, '\n');

fprintf(fid, ['SCALARS ', 'Aperture', ' float\n']);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', frac_aperture');
fprintf(fid, '\n');

fprintf(fid, ['SCALARS ', 'Poro', ' float\n']);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', frac_poro');
fprintf(fid, '\n');

fprintf(fid, ['SCALARS ', 'Perm', ' float\n']);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid, '%0.15f %0.15f %0.15f %0.15f %0.15f %0.15f\n', frac_perm');

fclose(fid);

end

function well2vtk(fname,wells)
%Write well into vtk

NumWells=size(wells,2);

%Extract all points
pts=[];
for wi=1:NumWells
    pts=[pts; wells(wi).points];
end

NumPts=size(pts,1);

%Find the poly node idx
WellPolyId={};
id=0;
for wi=1:NumWells
    NumWellPts=size(wells(wi).points,1);
    id_ST=id;
    id_END=id+NumWellPts-1;
    WellPolyId{end+1}=[NumWellPts id_ST:id_END];
    id=id_END+1;
end

WellTags=1:NumWells;

fid = fopen(fname, 'w'); 

fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Well Geometry\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');
fprintf(fid, 'POINTS %d float\n', NumPts);
fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', pts');fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, 'LINES %d %d\n', NumWells,NumWells+size(pts,1));
for wi=1:NumWells
    format=getPrintRowFormat(numel(WellPolyId{wi}),'%d');
    fprintf(fid, format, WellPolyId{wi}');
end

fprintf(fid, '\n');
fprintf(fid, 'CELL_DATA %d\n', NumWells);
fprintf(fid, '\n');

fprintf(fid, ['SCALARS ', 'WellTag', ' int\n']);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid, '%d ', WellTags');
fprintf(fid, '\n');

fclose(fid);


end

function grid2vtk(fname,G,dataNames,dataArrays)
%GRID2VTK Summary of this function goes here
%Detailed explanation goes here

NX=G.cartDims(1);
NY=G.cartDims(2);
if(numel(G.cartDims)==3), NZ=G.cartDims(3);end

%Row order in tensor grid
G_coordX=G.nodes.coords(1:NX+1,1);
G_coordY=G.nodes.coords(1:NX+1:(NX+1)*(NY+1),2);
if(numel(G.cartDims)==3)
    G_coordZ=G.nodes.coords(1:(NX+1)*(NY+1)+1:(NX+1)*(NY+1)*(NZ+1),3);
else
    G_coordZ=[0 G_coordX(end)/10.0];
end

fid = fopen(fname, 'w'); 
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Rectilinear grid\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET RECTILINEAR_GRID\n');
fprintf(fid, 'DIMENSIONS %d %d %d\n',NX+1,NY+1,NZ+1);
fprintf(fid, 'X_COORDINATES %d float\n', NX+1);
fprintf(fid, '%0.10f ', G_coordX(:)');fprintf(fid, '\n');
fprintf(fid, 'Y_COORDINATES %d float\n', NY+1);
fprintf(fid, '%0.10f ', G_coordY(:)'); fprintf(fid, '\n');
fprintf(fid, 'Z_COORDINATES %d float\n', NZ+1);
fprintf(fid, '%0.10f ', G_coordZ(:)');fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, 'CELL_DATA %d\n', NX*NY*NZ);

for i=1:length(dataNames)
    name=dataNames{i};
    data=dataArrays{i};
    if(size(data,2)==1)
        %Write SCALAR
        fprintf(fid, ['SCALARS ', name, ' float\n']);
        fprintf(fid,'LOOKUP_TABLE default\n');
        fprintf(fid, '%0.15f %0.15f %0.15f %0.15f %0.15f %0.15f\n', data(1:NX*NY*NZ)');
        fprintf(fid, '\n');
    end
    if(size(data,2)>1)
        %MultiComp Filed
        NumComp=size(data,2);
        fprintf(fid, 'FIELD FieldData 1\n');
        fprintf(fid, [name,' ', num2str(NumComp),' ', num2str(NX*NY*NZ), ' float\n']);
        %format=['%0.10f'];
        %for fi=1:NumComp-1, format=[format ' %0.10f'];end;
        %format=[format '\n'];
        format=getPrintRowFormat(NumComp,'%0.10f');
        fprintf(fid, format, data(1:NX*NY*NZ,:)');
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
end

fclose(fid);


end

function [format] = getPrintRowFormat(NumEntries,data_type)
    format=[data_type];
    for fi=1:NumEntries-1, format=[format ' ' data_type];end;
    format=[format '\n'];
end