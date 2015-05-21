function [ s ] = read_surfFS(fname)
%read_surfFS - read freesurfer surface and return surface struct
%
%      usage: [ s ] = read_surfFS( fname )
%         by: lpzds1
%       date: May 20, 2015
%        $Id$
%     inputs: fname
%    outputs: s
%
%    purpose: read freesurfer surfaces in standard format and return in
%             struct format as used for OFF files in mrTools 
%
%   see also: read_surf, read_label, renderSurf
%
%        e.g: 
%            fname = '/data/anatomy/freesurfer/subjects/ab/surf/lh.white'
%            s = read_surfFS(fname)
%            renderSurf(s)
%   

if nargin < 1
    help read_surfFS
    return
end

if ~(exist(fname)==2)
    % then the file does not exist
    disp('(read_surfFS) cannot find file!')
    return
end

% otherwise good to go
[v, f] = read_surf(fname);

% renderSurf wants:
% {'vtcs', 'tris', 'Nvtcs', 'Ntris', 'Nedges', 'filename'})

s.vtcs = v;
s.tris = f + 1; % 0-offset to 1 all the triangles are made up of vertices +1
s.Nvtcs = size(v,1);
s.Ntris = size(f,1);
s.Nedges = 0;
s.filename = fname;

end