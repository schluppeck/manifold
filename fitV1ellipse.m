function [ xform, H, p ] = fitV1ellipse(spCoords,  plotEllipse, useShear, debugPlot, spCoordsExpanded)
%fitV1ellipse - fit ellipse around coordinates of labels that define V1
%
%      usage: [ xform ] = fitV1ellipse( coords,[plotEllipse], [useShear], [debugPlot], [spCoordsExpanded] )
%         by: lpzds1
%       date: May 20, 2015
%        $Id$
%     inputs: coords [N x 3], usually theta/phi/r
%             plotEllipse [true] - draw resulting ellipse?
%             useShear [false] - apply a shear transform from Benson et al
%             debugPlot [false] - show some additional information
%             spCoordsExpanded - needed for new rendering (and benson style ellipse fitting)
%             
%             
%    outputs: xform, H (coordinates of the [convex] hull)
%             
%             
%    purpose: fit an ellipse around a collection of points in 2d (cartesian
%             or spherical coordinates). this is done by finding the convex
%             hull and then using a robust ellipse fitting method
%             (Fitzgibbon et al, 1999).
%
%             ALSO - implemented the Benson et al way of fitting V1 for
%             consistency... SET internal flag useConvexHull=false
%
%             the function plots the ellipse if the apporpriate flag is set
%             and returns an xform that allows un-shifting / unrotating
%             data points to be in a common frame of reference.
%
%   see also: fitellipse, drawellipse, convhull, plateau
% 
%        e.g: 
%             x = 5 .* randn(100,1) + 3;
%             y = 0.15 .* x + randn(150,1) + 0.5;
%             xE = 7 .* randn(150,1) + 3;
%             yE = 0.5 .*  randn(150,1) + 0.5;
%             figure, scatter(x,y,'r.'), hold on
%             fitV1ellipse([x y], true) % old style
%             fitV1ellipse([x y], true, [], true, [xE yE])  %  benson ellipse fit
%  

% set up some flags
if ieNotDefined('plotEllipse'), plotEllipse = 1; end
if ieNotDefined('useShear'), useShear = 0; end
if ieNotDefined('debugPlot'), debugPlot = 0; end
if ieNotDefined('spCoordsExpanded'), spCoordsExpanded = []; end

% - the benson-style fitting: set to false
% - the uon-style fitting: set to true (old-style)
if ieNotDefined('useConvexHull'), useConvexHull = false; end

% if useConvexHull is false, we also need to have expanded V1 coords for
% fitting! check here.
if ~useConvexHull 
    if isempty(spCoordsExpanded)
        disp('(uhoh) to use the BENSON et al method for fitting ellipse')
        disp('       need to have expanded V1');
        keyboard % for now!
    else
        % have access to expande V1 definition:
        xE = spCoordsExpanded(:,1);
        yE = spCoordsExpanded(:,2);
    end
    
end

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

if useConvexHull

    % this method is very quick, but it's not very robust to the shape of
    % the V1 patch on the cortical surface. in particular we noticed that
    % the calcarine sulcus is not very well aligned with the V1 model.
    %
    % the alternative method useConvexHull == false implements the Benson
    % et al version from the Curr Biol paper
    
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

else % don't use convex hull for fitting.
   
    % the alternative method useConvexHull == false implements the Benson
    % et al version from the Curr Biol paper
    
    % 1 - 1/(1 +
    %   Exp[-20*
    %     (Sqrt[\[Sigma]x*((pt[[1]] - x0)*
    %              Cos[\[Theta]] + (pt[[2]] - y0)*Sin[-\[Theta]])^2 +
    %         \[Sigma]y*((pt[[1]] - x0)*Sin[\[Theta]] + (pt[[2]] - y0)*
    %              Cos[\[Theta]])^2] - 2)])]},

    if debugPlot && false
        pSimulated = [pi/4, 0, 0, 25, 106];
        [xSim, ySim] = meshgrid(-0.5:0.01:0.5, -0.5:0.01:0.5); 
        v = plateau(pSimulated, xSim, ySim);
        figure, surf(xSim, ySim, v)
        shading interp
        material metal
        lighting phong
        light
        camlight HEADLIGHT
        axis off
        title('plateau function')
    end
    
    % determine which points of xE/yE are inside convex hull / V1
    
    k = convhull([x,y]);
    insideV1 = inpolygon(xE,yE,x(k),y(k));
    
    
    % make the V1 valueLandscape (1s inside, 0s outside)
    valueLandscape = zeros(size(xE));
    % now label stuff inside as 1
    valueLandscape(insideV1) = 1;
    
    
    % now do the minimization:
    % initial values for params:
    % [theta, x0, y0, sigmaX, sigmaY]
    p0 = [pi/8, mean(x), mean(y), 5, 20]; % benson et al
    pFinal = fmincon(@(params) plateauFit(params, [xE, yE], valueLandscape ), p0, ...
        [],[],[],[],[-pi, -5 -5 0 0], [pi 5 5 inf inf ]);
    
    insidePlateau = plateau(pFinal, xE, yE) > 0.5;
    xEinside = xE(insidePlateau);
    yEinside = yE(insidePlateau);
    
    kInside = convhull(xEinside, yEinside);
    
    p = fitellipse(xEinside(kInside), yEinside(kInside));

    

    

    if debugPlot
        % figure
        % subplot(1,2,1)
        % plot(xE, yE, 'ro', xE(~insideV1), yE(~insideV1), 'k+',...
        % xE(insidePlateau), yE(insidePlateau), 'gs', 'markerfacecolor', 'g' )
        % hold on
        % 
        % [ellipseP_, xyEllipse] = drawellipse(p, [] ,'color', 'b', 'linewidth',3);
        % 
        % subplot(1,2,2)
        % scatter(xE, yE, 15, valueLandscape); % 
        % 
        % figure, scatter(xE, yE,[], valueLandscape)
        % hold on
        % plot(xEinside(kInside), yEinside(kInside), 'r-')
        % keyboard

       
    end

    
    % plot the data used for alignment on top - the fitted ellipse underneath
    if debugPlot || plotEllipse
        [ellipseP_, xyEllipse] = drawellipse(p, [] ,'color', 'b', 'linewidth',3);
        %hold on
        %plot(x(k), y(k), 'ko', 'markerfacecolor', 'w', 'linewidth',2, 'markeredgecolor', 'k', 'markersize',15);
    end
    
end


    

% the params in |p| are: (Cx, Cy, Rx, Ry, theta_radians)
% ultimately, to undo the rotation / shift... can convert to homogeneous coords and
% then apply the inverse transform

% rotation could be Rot-pi/2 and Rx Ry swapped... as per comments in
% Fitzgibbon code. make this unambiguous here:

if p(5) < -pi/4
    p(5) = p(5) + pi/2;
    p([3 4]) = p([4 3]); % swap radii
    fprintf('(fitV1ellipse) theta < -pi/4\n');
elseif p(5) > +pi/4
    p(5) = p(5) - pi/2;
    p([3 4]) = p([4 3]); % swap radii
    fprintf('(fitV1ellipse) theta > +pi/4\n');
end

% at this point also check that the ordering of major minor axes is as we
% expect
assert(p(3) >= p(4), '(uhoh) major / minor axes swapped?' )

% xform contains the parameters of theta, centre shift...
xform = [cos(p(5)) -sin(p(5)) p(1)
         sin(p(5)) cos(p(5))  p(2)
         0          0           1 ];

% also, optionally return the coordinates of the convex hull:
H = [x(k), y(k)];
     
end


function val = plateau(p, x, y)
% plateau - for benson et al function fitting
%
% x's will be a list of V1inside and V1expanded
% y's will be ..........    1             0
%
% p [theta, x0, y0, sigmaX, sigmaY]
%
%      explanation: this function returns an elliptical looking plateau
%                   for fitting we will use a V1 region containing 1's and
%                   a surround filled with 0's - then we want to adjust the
%                   parameters in such a way that the squared errors are
%                   minimized...

theta = p(1);
x0 = p(2);
y0 = p(3);
sigmaX = p(4); 
sigmaY = p(5);
k = 20;

val = 1- 1./ ... 
    (1 + exp(-k .* (sqrt( ...
        sigmaX .* ((x - x0) .* cos(theta) + (y - y0) .* sin(-theta)).^2 ...
      + sigmaY .* ((x - x0) .* sin(theta) + (y - y0) .* cos(+theta)).^2) - 2)));

end

function e = plateauFit(p, xy, valueLandscape)
% plateauFit - function to pass to fminsearch...
%

val = plateau(p, xy(:,1), xy(:,2));

% SSE
e = sum((val - valueLandscape).^2);

end