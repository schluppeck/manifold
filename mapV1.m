function [ returnData ] = mapV1(subject, hemi, subjectRoot, showSphere, data)
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
%        e.g:
%             subjectRoot = '/Volumes/research/VNG/data/anatomy/freesurfer/subjects-7T'
%             % - or -
%             subjectRoot = getenv('SUBJECTS_DIR');
%             % - or -
%             subjectRoot = '/data/anatomy/freesurfer/subjects';
%             subject = 'ab';
%             hemi = 'rh';
%             figure, mapV1(subject,hemi,subjectRoot)


if ieNotDefined('subjectRoot'), subjectRoot = '/data/anatomy/freesurfer/subjects'; end
if ieNotDefined('subject'), subject = 'ab'; end
if ieNotDefined('hemi'), hemi = 'lh'; end
if ieNotDefined('showSphere'), showSphere = true(); end
if ieNotDefined('data'), data = []; end
if ieNotDefined('thr'), thr = 0.8; end % 0.8

% should be undo the residual rotation?
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
labelfolder = 'label';

% ... and naming conventions
surfname = sprintf('%s.white',hemi);
curvname = sprintf('%s.curv',hemi);
labelname = sprintf('%s.v1.prob',hemi);

% load fsaverage sphere...
fsaverage_surfname = sprintf('%s.sphere.reg',hemi);
full_surfname = fullfile(pname,surffolder,fsaverage_surfname);
full_curvname = fullfile(pname,surffolder,curvname);

% define some helper functions view this in spherical coords?
mycart2sph = @(X, ts) cart2sph(ts.vtcs(X(:,1),1), ts.vtcs(X(:,1),2),ts.vtcs(X(:,1),3));
% myshear = @(X) [1 0.65 0; 0 1 0; 0 0 1] * toHomogeneous([X(:,1), X(:,2)]); only for FSAVERAGE sphere?     
myshear = @(X) [1 0 0; 0 1 0; 0 0 1] * toHomogeneous([X(:,1), X(:,2)]);
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

if showSphere
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
    
if showSphere
    % and make the 2nd patch a subplot figure, too
    subplot(1,3,[2 3])
end


if unrotate
    % fit an ellipse around the points and don't show it
    [xform, h,p] = fitV1ellipse(W, false);
    W = myUnrotate(W, inv(xform));
end
[xform_after, h,p] = fitV1ellipse(W, true);

% show curvature pattern in 2d scatter plot or another funky version, which
% is implemented in showV1Patch
showV1Patch(W, c, v1labelRaw, thr, binarized);

axis equal
axis off

% show elliptical coordinate frame
% quadratic: (x/a)^2 + (y/b)^2 = 1
% geometric: x = a * sin(t);
%            y = b * cos(t);


% package up some data to return
returnData.ellipse.xform = xform;
returnData.ellipse.p = p;
returnData.ellipse.h = h;

returnData.patch.W = W;
returnData.patch.c = c;
returnData.patch.sSphere = sSphere;
returnData.patch.v1label = v1label;
returnData.patch.v1LabelRaw = v1labelRaw;


% and reset UNIX environment variable to what we found...
if didResetSubjectsDir
    setenv('SUBJECTS_DIR', oldSUBJECTS_DIR);
    disp('(!!) returned env variable for freesurfer')
end

end
