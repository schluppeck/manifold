function [  ] = showV1Patch(spCoords, curv, v1label, thr, binarized)
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

% flag that determines plotting style
oldStyle = false; % true -> scatter plot; false -> patch

if ieNotDefined('binarized'), binarized = false; end

% threshold the V1 label probabilites according to a value
% return all columns ID, sphericalCoords [1 2 3], label prob
thresholdedV1label = v1label(v1label(:,5)>thr , :);

% and pick the corresponding curvature values
fvertexcdata = curv(thresholdedV1label(:,1));

if binarized, fvertexcdata = -sign(fvertexcdata); end

if oldStyle
    % used to do something like this...
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