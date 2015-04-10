function [ h_ ] = renderSurf(s, data)
%renderSurf - use patch to render a surface struct
%
%      usage: [ h_ ] = renderSurf( s, data )
%         by: ds1
%       date: Mar 09, 2015
%     inputs: s
%    outputs: h_
%
%    purpose: shortcut for rendering a surface stored in a struct with some
%    nice default settings 
%
%        e.g: 
%             sr = loadSurfVTK('surf/rh.white.vtk')
%             figure, renderSurf( sr )
%
%       -or-
%             s = loadSurfVTK('surf/rh.sphere.vtk')
%             c = loadSurfVTK('surf/rh.curv.vtk')
%             figure, renderSurf( s, sign(c.data) )

if nargin < 2
    data = [];
end

if nargin < 1 || ~isstruct(s)
    help renderSurf
    return
end

if ~all(isfield(s, {'vtcs', 'tris', 'Nvtcs', 'Ntris', 'Nedges', 'filename'}))
    disp('(renderSurf) something wrong with the struct you passed in')
    return
end

if exist('data') && ~isempty('data') && numel(data) == s.Nvtcs
    % user wants to render data provides
    facevertexcdata = data(:);
elseif (isfield(s, 'data') && ~isempty(s.data)) && not(numel(s.data)==1 && all(s.data == 0))
   disp('(renderSurf) data present')
    facevertexcdata = s.data;
else
    % simply constant
    facevertexcdata = ones(size(s.vtcs));
end

if s.Ntris == 0
    disp('(renderSurf) no faces - cannot render - using point cloud instead')  
    usePointCloud = true;
else
    usePointCloud = false;
end

if usePointCloud
    h_ = scatter3(s.vtcs(:,1), s.vtcs(:,2), s.vtcs(:,3), 10, facevertexcdata); 
    light
else
    h_ = patch('vertices', s.vtcs, 'faces', s.tris, ...
        'facevertexcdata', facevertexcdata );
    material dull
    shading interp
    colormap gray
    light
end
axis equal
axis vis3d
axis off

end