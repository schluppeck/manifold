function [ h, xy ] = drawellipse(params, nPoints)
%drawellipse - draw an ellipse
%
%      usage: [ h, xy ] = drawellipse( params, nPoints )
%         by: lpzds1
%       date: Feb 18, 2015
%        $Id$
%     inputs: params, nPoints
%             []
%    outputs: h, xy
%
%    purpose: draw an ellipse given parameters returned from the function
%    fitellipse from fitzgibbon and others.
%
%        e.g: %        [rx,  ry,  cx, cy, rot]
%             params = [2.5, 1.2, 2,  5, pi/8];
%             figure, drawellipse(params, 360)
%             axis equal
%
%  see also: fitellipse
%

if nargin < 2
    nPoints = 100;
end

if nargin < 1
    help drawellipse
    return
end

if numel(params) ~= 5
    error('params needs to 5 elements long')
end

% params as returned by
% params = fitellipse(nx,ny);

% Draw the returned ellipse
t = linspace(0,pi*2, nPoints);
x = params(3) * cos(t); % scaled x
y = params(4) * sin(t); % scaled y

% rotated and shifted ellipse
nx = x*cos(params(5))-y*sin(params(5)) + params(1);
ny = x*sin(params(5))+y*cos(params(5)) + params(2);

% xy to return:
xy = [nx(:), ny(:)];

hold on
h = plot(nx,ny,'r-');


end