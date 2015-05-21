function [  ] = showV1Patch(spCoords, curv, v1label, thr)
%showV1Patch - render the patch of v1
%
%      usage: [  ] = showV1Patch( spCoords, curv, v1label, thr )
%         by: lpzds1
%       date: May 20, 2015
%        $Id$
%     inputs: spCoords, curv, v1label, thr
%    outputs:
%
%    purpose: display a set of labeled coordinates with an overlay
%
%   see also: fitV1ellipse, mapV1

oldStyle = false; % true -> scatter plot; false -> patch

% threshold the V1 label probabilites according to a value
thresholdedV1label = v1label(v1label(:,5)>thr , :);

% and pick the corresponding curvature values
fvertexcdata = curv(thresholdedV1label(:,1));

if oldStyle
    % used to do something like this...
    v1label = v1label(v1label(:,5)>thr , :);
    scatter(spCoords(:,1), spCoords(:,2),35, fvertexcdata , 'filled'), colormap(gray);   
else
    % pick out the triangulation of the points in V1label, then show them with
    % patch
    dt = delaunayTriangulation(spCoords(:,[1 2]));
    patch('vertices', dt.Points , 'faces', dt.ConnectivityList, 'facevertexcdata', fvertexcdata)
    shading flat 
    axis equal
end

end