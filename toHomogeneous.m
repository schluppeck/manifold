function coordsH = toHomogeneous(coords)
% toHomogeneous - convert from simple to homogeneous coords (2D)
%
% add a row of 1's to take coords (x,y,1 ) from homogeneous to (x,y)
% coords should be [x,y] - row vectors
%
%  ds - 2015/02

if size(coords,1) ~=2
    % assume that they are flipped)
    coords = coords';
end

coordsH = [coords; ones(1,size(coords,2))];

end