function [ h_ ] = renderSurf(s)
%renderSurf - use patch to render a surface struct
%
%      usage: [ h_ ] = renderSurf( s )
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
%             firgure, renderSurf( sr )

if nargin < 1 || ~isstruct(s)
    help renderSurf
    return
end

if ~all(isfield(s, {'vtcs', 'tris', 'Nvtcs', 'Ntris', 'Nedges', 'filename'}))
    disp('(renderSurf) something wrong with the struct you passed in')
    return
end

h_ = patch('vertices', s.vtcs, 'faces', s.tris, ...
    'facevertexcdata', ones(size(s.vtcs)) );
shading interp
colormap gray
light
axis equal
axis vis3d

end