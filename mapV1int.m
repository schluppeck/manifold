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

% flat map only. figure placement is chosen for this at the moment
showSphere = false; 

% viewGet information about the current session.
subject = lower(viewGet(v, 'subject'));
base = viewGet(v ,'base');

% strip away the "surfRelax" and "subect part, so we end up in
subjectRoot = strrep(base.coordMap.path, sprintf('%s/surfRelax',subject), '');

if ~(exist('subjectRoot','dir') == 7)
    % try local
    subjectRoot = getenv('SUBJECTS_DIR');
    if ~exist('subjectRoot','dir') == 7
        error(sprintf('(uhoh) really cannot find %s\n ', subjectRoot))
    end
end
    


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

% the left column of plots: PA, the right column: ECC
sp1_ = axes('Position',[0 0.5 0.5 0.5])
fsData = mapV1(subject, hemi, subjectRoot, showSphere );
colormap(gray)
caxis([-2 4])

sp2_ = axes('Position',[0.5 0.5 0.5 0.5])
mapV1(subject, hemi, subjectRoot, showSphere );
colormap(gray)
caxis([-2 4])

% then add (as a scatter plot for now) data from overlay!
overlays = viewGet(v, 'overlays');
curOverlay = overlays(overlayNum);
alphaOverlay = overlays(1); % that's how things are organized: r2, polaAngle, eccentricity, ... 
paOverlay = overlays(2);
eccOverlay = overlays(3);

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
    % vtcs(:,~(xCheck & yCheck & sCheck)) = nan ; % this leaves sizes the same, but sets bad indeces to nan
    
    % now convert to columns
    vtcs = round(vtcs(1:3,:)');
    
    % that means we also need to chop out the values from flat patch
    inbounds = find(xCheck & yCheck & sCheck);
    fsData.patch.W = fsData.patch.W(inbounds,:); % this causes problems with differeing sizes
    % fsData.patch.W(~inbounds,:) = nan; % keep but don't display?!
    
else
    disp('something is very wrong - no vtcs!')
    return
end


% now get values from the whole data array
overlayData = curOverlay.data{scan};
alphaData = alphaOverlay.data{scan};
paData = paOverlay.data{scan};
eccData = eccOverlay.data{scan};


% transform alpha data
% hard threshold on r2 or "gamma correct" with an exponent
r2threshold = 0.01;
alphaExponent = 0.5;
alphaData = real(alphaData.^alphaExponent); % avoid complex numbers...

% get all the appropriate overlays, not just the one that's shown in mrTools
overlay_values = overlayData(sub2ind(size(overlayData),vtcs(:,1), vtcs(:,2), vtcs(:,3)));
alpha_values = alphaData(sub2ind(size(overlayData),vtcs(:,1), vtcs(:,2), vtcs(:,3)));
pa_overlay_values = paData(sub2ind(size(overlayData),vtcs(:,1), vtcs(:,2), vtcs(:,3)));
ecc_overlay_values = eccData(sub2ind(size(overlayData),vtcs(:,1), vtcs(:,2), vtcs(:,3)));

% transparentScatter needs a circle size
sizeOfCircle = 0.01;

% get colormap info from mrTools
cmapRange = curOverlay.range;
cmap = curOverlay.colormap;
nColors = size(cmap,1);

% plotBensonModel functions - eccentricity
bensonECC = @(q, xCoord) 90.*exp(q .* (xCoord-1));
% the benson et al paper uses 0...180deg, we use radians
% bensonPA = @(q, yCoord) d2r(90 + 90 .*sign(yCoord) .* (abs(yCoord)).^q);
bensonPA = @(q, yCoord) d2r(90 .*sign(yCoord) .* (abs(yCoord)).^q);

% deal with data in each of two hemispheres
pafit.coords = fsData.patch.elCoords(inbounds,1); 
if strcmp(hemi, 'lh')
    % PA
    pafit.cmap = paOverlay.colormap;% (use the color map that's going)  
    pafit.cmapRange = [pi -pi]; 
elseif strcmp(hemi, 'rh')
    % PA
    pafit.cmap = circshift(paOverlay.colormap, round(nColors./2)); % shift color map
    pafit.cmapRange = [-pi pi]; 
    pa_overlay_values = mod(overlay_values+2*pi, 2*pi)-pi; % and make values in the same convenient range around 0
end

% ECC is the same for both hemispheres...
eccfit.coords = fsData.patch.elCoords(inbounds,2);
eccfit.cmap = eccOverlay.colormap% (use the color map that's going)  
eccfit.cmapRange = [0 8];

% display the data on top of the prepped flat patches

axes(sp1_)
% column 1 - PA
hold on
getColorsAndDisplay(fsData, pa_overlay_values, pafit.cmap, paOverlay.range, alpha_values);

axes(sp2_)
% column 2 - ECC
hold on
getColorsAndDisplay(fsData, ecc_overlay_values, eccOverlay.colormap, eccOverlay.range, alpha_values);

% and now fit and display the fits...

% only try to fit points for which there are data (therefore the
% restriction to inbounds)

initQ = 1;

% do the fitting - estimate least squares solution on each individually.
eccfit.q = fminsearch(@(val) nansum((bensonECC(val, eccfit.coords) - ecc_overlay_values).^2), ...
    initQ);

pafit.q = fminsearch(@(val) nansum((bensonPA(val, pafit.coords) - pa_overlay_values).^2), ...
    initQ);

% model fit
eccfit.model = bensonECC(eccfit.q, fsData.patch.elCoords(inbounds,2)); % PA first, then ecc
pafit.model = bensonPA(pafit.q, fsData.patch.elCoords(inbounds,1));


% in a separate subplot, show the fit(s)
sp_3 = axes('Position',[0 0.2 0.5 0.5]);
fit_alpha_values = alpha_values; % or % ones(size(pafit.model));
getColorsAndDisplay(fsData, pafit.model, pafit.cmap, pafit.cmapRange, fit_alpha_values);

hold on
e_ = drawellipse([0 0 fsData.ellipse.p(3:4) 0]);
set(e_, 'linewidth',2, 'color', [0 0 1]);

% and the fit value...
text(0,-0.25,sprintf('fit value: %.2f',pafit.q),'fontsize',14)

sp4_ = axes('position',[0.5 0.2 0.5 0.5]);
getColorsAndDisplay(fsData, eccfit.model, eccOverlay.colormap, eccfit.cmapRange, fit_alpha_values);

hold on
e_ = drawellipse([0 0 fsData.ellipse.p(3:4) 0]);
set(e_, 'linewidth',2, 'color', [0 0 1]);

text(0,-0.25,sprintf('fit value: %.2f',eccfit.q),'fontsize',14)



% make the next snippet of code contingent on what is being fitted...
% little diagnostic plot:
nBins = 50;
axes('position', [0.1 0.1 0.1 0.1])
hist(pa_overlay_values,nBins), xlabel('pa (overlay\_values)')
axis([-pi pi 0 inf])

axes('position', [0.25 0.1 0.1 0.1])
hist(pafit.model,nBins), xlabel('fit')
axis([-pi pi 0 inf])

% necc
axes('position', [0.6 0.1 0.1 0.1])
hist(ecc_overlay_values,nBins), xlabel('ecc (overlay\_values)')
axis([0 inf 0 inf])

axes('position', [0.75 0.1 0.1 0.1])
hist(eccfit.model,nBins), xlabel('fit')
axis([0 inf 0 inf])


end


function [colorVals, opacityVals] = mapValuesToColormap(overlay_values, cmap, cmapRange, alpha_values)
% mapValuesToColormap - map overlay_values into a triplet from a given colormap 
%
%    helper function for making color mainpulations a bit easier 
%
%    see also: mapV1int, mapV1, showV1patch, transparentScatter
%

% need to have whole numbers to index into color map, so no other method
% appropriate!
interpMethod = 'nearest'; 

% if no alpha_values are passed in, assume we want to show all!
if ieNotDefined('alpha_values'), alpha_values = ones(size(overlay_values)); end

nColors = size(cmap,1); % # of rows in cmap

colorIdx = interp1(linspace(cmapRange(1), cmapRange(2), nColors), 1:nColors, overlay_values, interpMethod);

badIdx = isnan(colorIdx);
colorIdx(badIdx) = 1; % lowest value in colormap - but also change alpha to 0

% package for return
opacityVals = alpha_values;
opacityVals(badIdx) = 0;

% and color vals are indexed into the map
colorVals = cmap(colorIdx,:);

end

function scatterPoints = getColorsAndDisplay(data, ovals, cmap, cmapRange, alphaVals, sizeOfCircle)
% getColorsAndDisplay - helper function for displaying maps transparently

if ieNotDefined('sizeOfCircle'), sizeOfCircle = 0.01; end

% get colors
[colorVals, opacityVals] = mapValuesToColormap(ovals, cmap, cmapRange, alphaVals);

% and plot them
[ scatterPoints ] = transparentScatter( data.patch.W(:,1), data.patch.W(:,2), sizeOfCircle, colorVals, opacityVals );

axis equal
axis off

end
