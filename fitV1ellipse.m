function [ xform, H, p ] = fitV1ellipse(spCoords,  plotEllipse, useShear, debugPlot)
%fitV1ellipse - fit ellipse around coordinates of labels that define V1
%
%      usage: [ xform ] = fitV1ellipse( coords,[plotEllipse], [useShear], [debugPlot] )
%         by: lpzds1
%       date: May 20, 2015
%        $Id$
%     inputs: coords [N x 3], usually theta/phi/r
%             
%    outputs: xform, H (coordinates of the convex hull)
%             plotEllipse [true] - draw resulting ellipse?
%             useShear [false] - apply a shear transform from Benson et al
%             debugPlot [false] - show some additional information
%             
%    purpose: fit an ellipse around a collection of points in 2d (cartesian
%             or spherical coordinates). this is done by finding the convex
%             hull and then using a robust ellipse fitting method
%             (Fitzgibbon et al, 1999).
%
%             the function plots the ellipse if the apporpriate flag is set
%             and returns an xform that allows un-shifting / unrotating
%             data points to be in a common frame of reference.
%
%   see also: fitellipse, drawellipse, convhull
% 
%        e.g: 
%             x = 5 .* randn(100,1) + 3;
%             y = 0.15 .* x + randn(100,1) + 0.5;
%             figure, scatter(x,y,'r.'), hold on
%             fitV1ellipse([x y], true)
%  

% set up some flags
if ieNotDefined('plotEllipse'), plotEllipse = 1; end
if ieNotDefined('useShear'), useShear = 0; end
if ieNotDefined('debugPlot'), debugPlot = 0; end

% a selection of points in 2d / e.g. where V1 is labelled:
x = spCoords(:,1);
y = spCoords(:,2);

if debugPlot
    % for debugging, plot the position of the points
    scatter(x, y, [], 'r', 'o'); % markertype
    axis equal
end

% benson et al. apply a shear transform to flatten out the 2d plots of
% spherical coordinates. see e.g. the PLoS Comp Biol 2014 paper for a short
% description of this.

if useShear
    % there is also a shear transformation, according to Benson et al  {{1, 0.65}, {0, 1}};
    
    % we usually deal with this in homogeneous coordinates, although there
    % it's not strictly necessary (as there is no translation)
    sMatrix = [1 0.65 0;
        0 1 0;
        0 0 1];

    xySheared = sMatrix*toHomogeneous([x(:), y(:)]);
    
    if debugPlot
        hold on
        scatter(xySheared(1,:), xySheared(2,:), [], 'g', 'o'); 
    end
    
    % and re-assign the transformed coordinates.
    x = xySheared(1,:)';
    y = xySheared(2,:)';
end

% calculate the convex hull of the points. this will allow us to robustly
% fit an ellipse. this is different from Benson et al's approach, but
% seemed less complicated to me and apears to work just as well

k = convhull([x,y]);

% the following nicely does the job of fitting an ellipse - see comments in
% function for a reference to the Fitzgibbon et al paper
p = fitellipse(x(k), y(k));

if debugPlot || plotEllipse
    [ellipseP_, xyEllipse] = drawellipse(p, [] ,'color', 'b', 'linewidth',3);
end

% plot the symbols on top - the fitted ellipse underneath
if debugPlot || plotEllipse
    hold on
    plot(x(k), y(k), 'ko', 'markerfacecolor', 'w', 'linewidth',2, 'markeredgecolor', 'k', 'markersize',15);
end


% the params in |p| are: (Cx, Cy, Rx, Ry, theta_radians)
% ultimately, to undo the rotation / shift... can convert to homogeneous coords and
% then apply the inverse transform

% rotation could be Rot-pi/2 and Rx Ry swapped... as per comments in
% Fitzgibbon code. make this unambiguous here:

if p(5) < -pi/4
    p(5) = p(5) + pi/2;
    p([3 4]) = p([4 3]); % swap radii
    fprintf('(fitV1ellipse) theta < -pi/4');
elseif p(5) > +pi/4
    p(5) = p(5) - pi/2;
    p([3 4]) = p([4 3]); % swap radii
    fprintf('(fitV1ellipse) theta > +pi/4');
end

% xform contains the parameters of theta, centre shift...
xform = [cos(p(5)) -sin(p(5)) p(1)
         sin(p(5)) cos(p(5))  p(2)
         0          0           1 ];

% also, optionally return the coordinates of the convex hull:
H = [x(k), y(k)];
     
end