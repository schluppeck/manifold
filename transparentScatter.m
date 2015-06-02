function [ scatterPoints ] = transparentScatter(x, y, sizeOfCircle, color, opacity)
%transparentScatter - plot data scatter with transparency
%
%      usage: [ scatterPoints ] = transparentScatter( x, y, sizeOfCircle, color, opacity )
%         by: lpzds1
%       date: Jun 02, 2015
%
%     inputs: x, y, sizeOfCircle, color, opacity
%    outputs: scatterPoints
%
%    purpose: make transparent scatter plot
%            [arguments are vectorized, so should work ok] with inputs
%
%            inspired by http://stackoverflow.com/questions/6366404/semi-transparent-markers-in-matlab-figures
%
%            NB! could add option to deal with different pbaspect ratios,
%            for now, assuming that we want to show things with axis equal
%            Otherwise assumption about little circles is wrong and need to
%            multiply with appropriate factor
%
%        e.g:
%             n = 150;
%             cmap = hot(64);
%             colors = cmap( randi(64, n,1),:);
%             scatterPoints = transparentScatter(randn(n,1),randn(n,1),0.1,colors,randn(n,1));
%             axis equal % !to get circles rather than ellipses

if ieNotDefined('sizeOfCircle')
    sizeOfCircle = range(x) ./ 10; % silly heuristic?
end

if ieNotDefined('color')
    defaultColors = get(0,'DefaultAxesColorOrder');
    color = defaultColors(1,:);
else
    % check that either 1x3 of nx3
    szCheck = ( all(size(color) == [1 3]) ) ||  ( (size(color,1) == size(x,1)) && size(color,2) ==3);
    assert(szCheck, 'color needs to have 1 or n rows, 3 columns!');
    % could add option for string based colors, too...
end

if ieNotDefined('opacity'), opacity = 0.2; end


defaultColors = get(0,'DefaultAxesColorOrder');
assert(size(x,2)  == 1 && size(y,2)  == 1 , 'x and y should be column vectors');

% size of disks
t= 0:pi/10:2*pi;

rep_x = repmat(x',[size(t,2),1]);
rep_y = repmat(y',[size(t,2),1]);
rep_t = repmat(t',[ 1, size(x,1)]);

% reshape color so it can be passed to patch
rsColor = reshape(color, [1 size(color)]);

scatterPoints = patch((sizeOfCircle*sin(rep_t)+ rep_x),(sizeOfCircle*cos(rep_t)+rep_y), rsColor,'edgecolor','none');
alpha(scatterPoints,opacity);

end