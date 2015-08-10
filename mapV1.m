function [ returnData ] = mapV1(subject, hemi, subjectRoot, showSphere, showPlot, firstIteration)
%mapV1 - get V1 definition automatically in a subject
%
%      usage: [  ] = mapV1(subject, hemi, subjectRoot, showSphere )
%         by: lpzds1
%       date: May 18, 2015
%     inputs: subject, hemi, subjectRoot, showSphere (true/false)
%    outputs:
%
%    purpose: using freesurfer functions and segmentations.
%
%   see also: bensonCoordinates, fitBensonCoords, mapV1int
%
%        e.g:
%             subjectRoot = '/Volumes/research/VNG/data/anatomy/freesurfer/subjects-7T'
%             % - or -
%             subjectRoot = getenv('SUBJECTS_DIR');
%             % - or -
%             subjectRoot = '/data/anatomy/freesurfer/subjects';
%             subject = 'ab';
%             hemi = 'rh';
%             figure, mapV1(subject,hemi,subjectRoot)


if ieNotDefined('subjectRoot'), subjectRoot = '/data/anatomy/freesurfer/subjects-7T/'; end
if ieNotDefined('subject'), subject = 'fsaverage'; end
if ieNotDefined('hemi'), hemi = 'rh'; end
if ieNotDefined('showSphere'), showSphere = false; end
if ieNotDefined('showPatchForDebug'), showPatchForDebug = false; end

if ieNotDefined('showPlot'), showPlot = false; end
if ieNotDefined('firstIteration'), firstIteration = false; end

% probability threshold for deciding what is INSIDE V1
if ieNotDefined('thr'), thr = 0.8; end % 0.8

% should be undo the residual rotation in the V1 ellipse?
% if set to TRUE then an affine transform is applied to center the ellipse
% on 0,0 with the major axis aligned to the xaxis
% if set to FALSE, then the points are left in place.
if ieNotDefined('unrotate'), unrotate = true; end

% show curvature binarized true/false?
binarized = true;

% check that code is available
% modified version of read_surf
if ~(exist('read_surfFS','file') == 2)
    disp('(uhoh) some important functions are not on the path!?')
    disp('       addpath(''~/matlab/manifold/'')' )
    return
end

% which depends on having /Applications/freesurfer/matlab on the path!
if ~(exist('read_surf','file') == 2)
    disp('(uhoh) some important functions are not on the path!?')
    disp('       addpath(''/Applications/freesurfer/matlab/'')' )
    return
end

% check if the SUBJECTS_DIR is consistent with what user asks for
% FS read_label relies on this to be set correctly...
%
% !TODO refactor this at some point to remove this assumption
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
if ~strcmpi(SUBJECTS_DIR, subjectRoot)
    oldSUBJECTS_DIR = SUBJECTS_DIR;
    setenv('SUBJECTS_DIR', subjectRoot);
    didResetSubjectsDir = true;
    disp('(!!) had to reset env variable for freesurfer')
else
    didResetSubjectsDir = false;
end

pname = fullfile(subjectRoot, subject);

% set some relative locations
surffolder = 'surf';
% special case: if RH and fsaverage_sym:
if strcmp(hemi,'rhx')
    surffolder = 'xhemi/surf';
    hemi = 'lh'; % and reset
elseif strcmp(hemi,'lhx')
    surffolder = 'surf';
    hemi = 'lh'; % and reset
end

labelfolder = 'label';

% ... and naming conventions
% for making this function more general (S1? Brodman areas) we could
% consider passing these in as option arguments
surfname = sprintf('%s.white',hemi);
curvname = sprintf('%s.curv',hemi);
labelname = sprintf('%s.v1.prob',hemi);

% load fsaverage sphere...
fsaverage_surfname = sprintf('%s.sphere.reg',hemi);
full_surfname = fullfile(pname,surffolder,fsaverage_surfname);
full_curvname = fullfile(pname,surffolder,curvname);

% define some helper functions view this in spherical coords?
mycart2sph = @(X, ts) cart2sph(ts.vtcs(X(:,1),1), ts.vtcs(X(:,1),2),ts.vtcs(X(:,1),3));
if any(strcmp(subject,{'fsaverage','fsaverage_sym'}))
    if strcmp(hemi, 'lh')  % 0.65 and -0.65
       myshear = @(X) [1 0.65 0; 0 1 0; 0 0 1] * toHomogeneous([X(:,1), X(:,2)]); % only for FSAVERAGE sphere?
    elseif strcmp(hemi, 'rh')
       myshear = @(X) [1 -0.65 0; 0 1 0; 0 0 1] * toHomogeneous([X(:,1), X(:,2)]); % only for FSAVERAGE sphere?
    end
       disp('using special SHEAR matrix!!')
else
    myshear = @(X) [1 0 0; 0 1 0; 0 0 1] * toHomogeneous([X(:,1), X(:,2)]);
end
fromHomogeneous = @(X) transpose(X([1 2],:));
myUnrotate = @(X, xform) fromHomogeneous(xform * toHomogeneous([X(:,1), X(:,2)]));
myScatter = @(X,ts) scatter(  ts.vtcs(X(:,1),1), ts.vtcs(X(:,1),2), 5, 'r.'  );
myScatter3 = @(X,ts) scatter3(  ts.vtcs(X(:,1),1), ts.vtcs(X(:,1),2),ts.vtcs(X(:,1),3), 5, 'r.'  );

% load sphere...
try
    sSphere = read_surfFS(full_surfname);
catch
    fprintf('(ugh) could not read %s\n', full_surfname)
    return
end

% corresponding curvature
try
    c = read_curv(full_curvname);
catch
    fprintf('(ugh) could not read %s\n', full_curvfname)
    return
end

% read the label containing V1
try
    v1labelRaw = read_label(subject,labelname); % uses GETENV, etc
catch
    fprintf('(ugh) could not read data for subject: %s, label: %s\n', subject, labelname);
    return
end

% make sure vertices are 1-offset
v1labelRaw = oneoffset_label(v1labelRaw);
% and threshold the label
[ v1label ] = threshold_label(v1labelRaw, thr);

% convert v1 label and coordinates into SPHERICAL coordinates
S = zeros(size(v1label,1), 3);
[S(:,1), S(:,2), S(:,3) ] = mycart2sph(v1label(:,1), sSphere);
W = transpose(myshear(S(:,[1 2]))); % long rows, so need to transpose

Sraw = zeros(size(v1labelRaw,1), 3);
[Sraw(:,1), Sraw(:,2), Sraw(:,3) ] = mycart2sph(v1labelRaw(:,1), sSphere);
Wraw = transpose(myshear(Sraw(:,[1 2]))); % long rows, so need to transpose

if showPlot
    figHandle = gcf;
end

if showSphere && showPlot
    % plot into a subplot
    subplot(1,3,1)
    renderSurf(sSphere,sign(c))
    hold on,
    myScatter3(v1label, sSphere)
    
    % set view to home in on V1 location on RH and LH, respectively...
    if strcmp(hemi,'rh')
        view([-45 -30]);
    else
        view([30 -30]);
    end
    camlight('left')
end


if unrotate
    % fit an ellipse around the points and don't show it
    % - now uses the expanded V1 stored in Wraw as well...
    
    [xform, h,p] = fitV1ellipse(W, false,[],[], Wraw);
    
    % this overwrites the points stored in W, which are then used further
    % along
    
    % keep a copy of rotated W / raw (for debug plot)
    Wrotated = W;
    WrawRotated = Wraw;
    W = myUnrotate(W, inv(xform));
    Wraw = myUnrotate(Wraw, inv(xform));
    
else
    Wrotated = []; % no need to keep an extra copy. W is still rotated
end

if showPlot && ~firstIteration
    if showSphere
        % and make the 2nd patch a subplot figure, too
        subplot(1,3,[2 3])
    end
    
    % show curvature pattern in 2d scatter plot or another funky version, which
    % is implemented in showV1Patch
    showV1Patch(Wraw, c, v1labelRaw, thr, binarized);
    alpha(0.5)
    showV1Patch(W, c, v1label, thr, binarized);
    
    hold on
end

% fit and show ellipse to "de/un-rotated points". Also get the parameters
% again. The center should be 0,0, the radii the same as before and the
% rotation 0 rad!
[xform_after, hA,pA] = fitV1ellipse(W, showPlot,[],[],Wraw); % either show or don't

if showPlot && ~firstIteration
    colormap(gray)
    caxis([-5 5])
    daspect([2 1 1]) % axis equal
    axis off
end
% package up some data to return
returnData.ellipse.xform = xform;
returnData.ellipse.p = p;
returnData.ellipse.pA = pA;

returnData.ellipse.h = h;

returnData.patch.W = W;
returnData.patch.c = c;
returnData.patch.sSphere = sSphere;
returnData.patch.v1label = v1label;
returnData.patch.v1LabelRaw = v1labelRaw;

if showPatchForDebug
    fd_ = figure; % debug figure, remember to switch back to main figure later
    subplot(2,1,1)
    scatter(W(:,1), W(:,2), 'ro')
    daspect([2 1 1]) % axis equal
    title(sprintf('ellipse; sub: %s, hemi: %s', subject, hemi))
    if ~isempty(Wrotated)
        subplot(2,1,2)
        scatter(WrawRotated(:,1), WrawRotated(:,2), 'ro')
        hold on
        scatter(Wrotated(:,1), Wrotated(:,2), 'b+')
        daspect([2 1 1]) % axis equal
        xNeg = sum(Wrotated(:,1) <= -pi/2);
        xPos = sum(Wrotated(:,1) > -pi/2);
        title(sprintf('orig ellipse; sub: %s, hemi: %s [< / > -pi/2] %i/%i', subject, hemi, xNeg, xPos))
    end
end

if showPlot
    figure(figHandle); % switch back to main figure
end

% spherical coords are stored in W for this patch... ellipse params in pA
% (after un-rotation step)
returnData.patch.elCoords  = fitBensonCoords( W, pA, hemi);



% and reset UNIX environment variable to what we found...
if didResetSubjectsDir
    setenv('SUBJECTS_DIR', oldSUBJECTS_DIR);
    disp('(!!) returned env variable for freesurfer')
end

end
