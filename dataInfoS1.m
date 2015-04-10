%dataInfoS1 - quick walkthrough of how data are organized
%
%      usage: [  ] = dataInfoS1(  )
%         by: ds1
%       date: Mar 05, 2015
%     inputs: 
%    outputs: 
%
%    purpose: function that illustrates how data are organized and how to
%    access info in the different data structures 
%
%        e.g: 
%             dataInfoS1

%% 4d functional imaging data is stored in this mat file
% we usually store data in a different format (NIFTI), but for convenience,
% here I saved out the same info into a MAT file
fname = 'somato-fMRI-periodic';  % big data file. > 60mb
load(fname)

% provides "data" and "hdr" (some info)
% dims 1-3 are space (x,y,z), dim 4 is time

load('S1-definition-ds20100728')
% provides S1 (some coordinates and also a matrix/transform for getting
% between the scan (data) space and the anatomy space in which the cortical
% surfaces are stores.
disp('4x4 transform that takes ANATOMY coordinates to SCAN space')
disp(S1.base2scan)

% for this data set, separation between timepoints is 2.4s
t = 2.4 .* [0:(size(data,4)-1)] ;

% so plotting time series at x,y,z: 37,61,17:
voxCoord = [37,61,17];
figure, plot(t, squeeze(data(voxCoord(1), voxCoord(2), voxCoord(3),:)), 'r-')
title(sprintf('fMRI signal over time at [%i,%i,%i]',voxCoord ))
xlabel('Time (s)')
ylabel('fMRI reponse (image intensity)')

% or 2d image at z = 17, t=52
whichZ = 17; whichT = 52;
figure
imagesc(squeeze(data(:,:, whichZ, whichT))); 
axis image
colormap(gray), colorbar
title(sprintf('Slice through data at z=%i, t=%i', whichZ, whichT))

%% onto the surface
%
% how does this map onto the 2d surface?
% fMRI data were obtained with a stimulus on the left hand 
% therefore look on right hemisphere

% if we load in a surface and try to convert the coords into scan space:
% note that the VTK files are read in using an additional transform 
%   - the 1 sets transformToSurfRelax on, which shifts the origin
% to take into account differences in how data formats are used in our analysis tools 
% (geomview [.off] and freesurfer [-> to .VTK])
%
% we also have to deal with the fact that in matlab indeces go 1 .. n
%                                          and in c indices go 0 .. n-1
 
s = loadSurfVTK('surf/rh.white.vtk', 1); 

% now xformSurfaceWorld2Array to shift into frame of ref of data
s = xformSurfaceWorld2Array(s, base.hdr);

figure, renderSurf(s)
alpha(0.9)

% find where that voxel would have been:
voxCoord = [37,61,17];
% need inv(xform) because we are going the other way
voxCoordInAnatomy = inv(S1.base2scan) * [ voxCoord, 1]';

% and plot the point
hold on
p_ = plot3(voxCoordInAnatomy(1), voxCoordInAnatomy(2), voxCoordInAnatomy(3), 'ro');
set(p_, 'markersize',15, 'markerfacecolor','r')

% plot some of the vertices corresponding to S1, say every 250th, to make
% rendering abit quicker
skip = 250;
s1_ = plot3(S1.volumeCoords(1:skip:end,1), ...
    S1.volumeCoords(1:skip:end,2), ...
    S1.volumeCoords(1:skip:end,3), 'b.');

%% onto the sphere
%
% now load in sphere representation of the hemisphere
sSphere = loadSurfVTK('surf/rh.sphere.vtk', 1); 

% and the curvature info.
sCurve = loadSurfVTK('surf/rh.curv.vtk', 1); 

% across all those surfaces for this one subject, the vertex number (its
% ID) remains the same, only their locations in 3D change...
%
% renderSurf is just a small wrapper function around matlab
% |patch| to make the code a bit cleaner.

figure, subplot(2,1,1)
renderSurf(sSphere, -sign(sCurve.data))
title('sign(curvature) rendered on spherical representation')

subplot(2,1,2)
renderSurf(sCurve)
title('curvature rendered on normal/white matter representation')
colormap(gray)

%% to find the points of interest on the sphere 
%
% discussion c. 2015-04-10
%
% find vertices of interest in e.g. the "white" surface. this corresponds 
% to the surface at the interface between white matter and gray matter.
%
% then, knowing the vertex IDs, we can look up spatial location from the
% .sphere surfaces for the same cortical hemisphere & subject

% coordinates of voxels in 3d space (homogeneous)
S1.volumeCoords;

% vertices of the "white matter surface"
s.vtcs;

% vertices of the "gray matter surface"
g = loadSurfVTK('surf/rh.pial.vtk', 1); 

% mid cortex is the halfway point between white and gray matter
% the surfaces have corresponding vertices, so we can use the following
% trick

m.vtcs = s.vtcs + 0.5 .* (g.vtcs - s.vtcs); % half way between them
m.tris = s.tris; % triangulation stays the same
m.Ntris = s.Ntris;
m.Nvtcs = s.Nvtcs;
m.filename = 'calculated from wm and gm surfaces';

% there is a c/mex implementation of this, courtesy of wade/wandell
% (VISTA soft) to speed things up (assignToNearest.m/mexmaci)
% 
% but matlab's knnsearch should do the right thing:
%  first input: point-cloud, second input: query points (making sure not to use homogeneous coords!).

% can try to find out which vertices are closest in WM, mid cortex and GM
% in our software, something slightly different happens, but this should
% work.
idxS = knnsearch(s.vtcs, S1.volumeCoords(:,1:3));
idxM = knnsearch(m.vtcs, S1.volumeCoords(:,1:3));
idxG = knnsearch(g.vtcs, S1.volumeCoords(:,1:3));

% some voxels will point to the same 3d point, so could take "unique", but
% set operation will take care of duplicates, too

idx = union(union(idxS, idxM), idxG);


% now render, taking the INDICES of the original surface, but the 3d
% locations of those vertices mapped into the sphere by the freesurfer
% program

figure
renderSurf(sSphere, -sign(sCurve.data));

title('points of interest in 3d')

hold on

% here we take the 3d coords on a sphere (using indices from the closest
% points in the folded surface)
spherePoints_ = plot3(sSphere.vtcs(idx,1), ...
    sSphere.vtcs(idx,2), ...
    sSphere.vtcs(idx,3), 'b.', 'markersize',5);

axis vis3d

% so the appropriate 3d coords (on the sphere) are just
sSphere.vtcs(idx,:);

% show some:
sSphere.vtcs(idx(1:10),:)

% and on to spherical coords via
% [TH,PHI,R] = cart2sph(X,Y,Z)


