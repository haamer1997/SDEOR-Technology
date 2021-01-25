%% This function creates a single fracture based on the type curve problem geometry for fractional flow
% make sure to copy this file to mrst-2019b\modules\pEDFM\utils
% Edit: Hassan
function [G,frac_ids,well_ids]=explicitStencilSingleFrac(physdim, varargin)
    opt = struct('aperture', 0.02,...
                 'numStages'        ,  1,                      ...
                 'fracSpacing'      , 100,                      ...
                 'fracHalfLength'   , 100,                      ...
                 'fracHeight', 60,                             ...
                 'nyRefine',10,...
                 'tipNX',5,...
                 'physdim1',200,...
                 'nx',11, ...
                 'nz',11, 'gridType','geomspacing');
    opt = merge_options(opt, varargin{:});
    
    switch(lower(opt.gridType)) %only linspacing works for now
        case 'geomspacing'
            ptY = geomspace(opt.aperture/2, opt.fracSpacing/2, opt.nyRefine,opt.aperture);
            wellCellcentroidY = opt.fracSpacing/2;
            facesXcoords = wellCellcentroidY + [-fliplr(ptY) ,ptY];
        case 'logspacing'
            ptX = logspace(log10(opt.aperture/2),log10(opt.fracSpacing/2),opt.nxRefine);
            wellCellcentroidX = opt.fracSpacing/2;
            facesXcoords = wellCellcentroidX + [-fliplr(ptX) ,ptX];
        case 'linspacing'
            ptX = linspace(opt.aperture/2, opt.fracSpacing/2, opt.nx);
            wellCellcentroidX = opt.fracSpacing/2;
            facesXcoords = wellCellcentroidX + [-fliplr(ptX) ,ptX];
        case 'cartesian'
            facesXcoords = linspace(0, opt.fracSpacing, 2*opt.nxRefine);             
        otherwise
            error('Unknown case')
    end
    ptY = geomspace(opt.aperture/2, physdim(2)/2, opt.nyRefine,opt.aperture);
    wellCellcentroidY = physdim(2)/2;
    facesYcoords = wellCellcentroidY + [-fliplr(ptY) ,ptY];      
    facesZcoords = linspace(0,physdim(3),opt.nz);
    
    G = tensorGrid(facesXcoords, facesYcoords, facesZcoords);
    %plotGrid(G)
    Frac_I = 1:G.cartDims(1);
    Frac_J = opt.nyRefine;
    Frac_K = 1:opt.nz-1;
    z_idx = ceil(numel(Frac_K)/2);
    [II,JJ,KK]=meshgrid(Frac_I,Frac_J,Frac_K);
    frac_ids = sub2ind(G.cartDims, II,JJ,KK);
    frac_ids = reshape(frac_ids,1,[]);
    well_ids=[median(frac_ids)]; %place well cell to be in the middle of the frac plane
%     centerId=ceil(numel(frac_ids(:))/2);
%     well_ids=[frac_ids(centerId)];
end
