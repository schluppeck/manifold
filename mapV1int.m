function [  ] = mapV1int(v, overlayNum, scan, x, y, s, roi)
%mapV1int - mapV1 interrogator function to be run from mrTools
%
%      usage: [  ] = mapV1int( v, overlayNum, scan, x, y, s, roi )
%         by: lpzds1
%       date: Jun 01, 2015
%        $Id$
%     inputs: v, overlayNum, scan, x, y, s, roi
%    outputs: 
%
%    purpose: load fMRI data with freesurfer / benson et al ellipses, etc
%
%        e.g: 

% figure out which segmentation to go to

% flat patch, label, etc

% display in extra figure window: flat map and also overlay data:

subject = lower(viewGet(v, 'subject'));
base = viewGet(v ,'base');

% strip away the "surfRelax" and "subect part, so we end up in 
subjectRoot = strrep(base.coordMap.path, sprintf('%s/surfRelax',subject), '');

% the name of the base contains a hint whether it should be lh or rh
% e.g. : 'bm_right_WM_Flat_109_53_99_Rad65'
rightMatch = strfind(base.name,'_right_');
leftMatch = strfind(base.name,'_left_');

if ~isempty(rightMatch) && rightMatch(1)
    hemi = 'rh'
elseif ~isempty(leftMatch) && leftMatch(1)
    hemi = 'lh'
else
    disp('(mapV1int) big uhoh - no matching side...')
    return
end

% make a new graph window or select one if it's open
g_ = selectGraphWin;

showSphere = true;
fsData = mapV1(subject, hemi, subjectRoot, showSphere );
colormap(gray)
caxis([-2 4])

% then add (as a scatter plot for now) data from overlay!
overlays = viewGet(v, 'overlays');
curOverlay = overlays(overlayNum);
alphaOverlay = overlays(1); % that's how it is...
hold on

% get vertices
inner = loadSurfOFF( fullfile(base.coordMap.path, base.coordMap.innerCoordsFileName) );
outer = loadSurfOFF( fullfile(base.coordMap.path, base.coordMap.outerCoordsFileName) );
inner = xformSurfaceWorld2Array(inner,base.hdr);
outer = xformSurfaceWorld2Array(outer,base.hdr);

midCortex = inner; % initialize, then fix.
midCortex.fname = 'midcortex';
midCortex.vtcs = inner.vtcs + 0.5.*(outer.vtcs - inner.vtcs);
midCortex.tris = inner.tris;

base2scan = viewGet(v, 'base2scan');
scan2base = inv(base2scan);

% indices where V1 is labelled:
idx = fsData.patch.v1label(:,1);
vtcs = base2scan * toHomogeneous( midCortex.vtcs(idx,:) );
% [vtcs, IA, IC] = unique(round(vtcs(1:3,:)'), 'rows');
% the corresponding points on the FS flat patch are in W

if ~isempty(vtcs) 

  % check scan dimensions
  scanDims = viewGet(v,'scandims');

  % make sure we are inside scan dimensions
  xCheck = (vtcs(1,:) >= 1) & (vtcs(1,:) <= scanDims(1));
  yCheck = (vtcs(2,:) >= 1) & (vtcs(2,:) <= scanDims(2));
  sCheck = (vtcs(3,:) >= 1) & (vtcs(3,:) <= scanDims(3));

  % only return ones that are in bounds
  vtcs = vtcs(:,find(xCheck & yCheck & sCheck));
  % now convert to columns 
  vtcs = round(vtcs(1:3,:)');

  % that means we also need to chop out the values from flat patch
  inbounds = find(xCheck & yCheck & sCheck);
  fsData.patch.W = fsData.patch.W(inbounds,:);
  
else
    disp('something is very wrong - no vtcs!')
    return
end


% now get values from the whole data array
overlayData = curOverlay.data{scan};
alphaData = alphaOverlay.data{scan};

overlay_values = overlayData(sub2ind(size(overlayData),vtcs(:,1), vtcs(:,2), vtcs(:,3)));
alpha_values = alphaData(sub2ind(size(overlayData),vtcs(:,1), vtcs(:,2), vtcs(:,3)));

% interp from 0...1 * colormap range, then use (4) for alpha?! Undocumented

sizeOfCircle = 0.01;
%get colormap info from mrTools
cmapRange = curOverlay.range;
cmap = curOverlay.colormap;
nColors = size(cmap,1);

colorIdx = interp1(linspace(cmapRange(1), cmapRange(2), nColors), 1:nColors, overlay_values, 'nearest'); 

badIdx = isnan(colorIdx);
colorIdx(badIdx) = 1; % lowest value in colormap - but also change alpha to 0
opacity = alpha_values;
opacity(badIdx) = 0;


[ scatterPoints ] = transparentScatter( fsData.patch.W(:,1), fsData.patch.W(:,2), sizeOfCircle, cmap(colorIdx,:), opacity )
% s = scatter(fsData.patch.W(:,1), fsData.patch.W(:,2), 100, overlay_values, 'filled') ;

% make sure we can go from 3D coords back to vertices in the Freesurfer
% world!
keyboard


end


