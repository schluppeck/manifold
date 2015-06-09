function [ h, xy ] = drawellipse(params, nPoints, varargin)
%drawellipse - draw an ellipse
%
%      usage: [ h, xy ] = drawellipse( params, nPoints, [linspec] )
%         by: lpzds1
%       date: Feb 18, 2015
%        $Id$
%     inputs: a) params, nPoints, [linespec]
%             b) params, nPoints, 'fillColor', [1x3, rgb]
%             
%    outputs: h, xy
%
%    purpose: draw an ellipse given parameters returned from the function
%    fitellipse from fitzgibbon and others.
%
%        e.g: %        [rx,  ry,  cx, cy, rot]
%             params = [2.5, 1.2, 2,  5, pi/8];
%             figure, drawellipse(params, 360)
%             hold on, drawellipse(params, 360, 'color', 'b', 'linewidth', 5)
%             axis equal
%
%             % - or - filled
%             drawellipse(params, 360, 'fillColor', [1 1 0])
%
%  see also: fitellipse
%

if nargin < 2 || isempty(nPoints)
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
if isempty(varargin)
    h = plot(nx,ny,'r-');
elseif strcmp(varargin{1}, 'fillColor')
    disp('(!) drawing a filled ellipse')
    h_ = patch(nx, ny, varargin{2});
else
    h = plot(nx,ny,varargin{:});
end

end