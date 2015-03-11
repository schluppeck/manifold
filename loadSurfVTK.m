function [ s ] = loadSurfVTK(filename, transformToSurfRelax, volumeSize, verbose)
%loadVTK - load VTK file into structure
%
%      usage: [ s ] = loadSurfVTK( filename, [verbose] )
%         by: lpzds1
%       date: Feb 10, 2015
%     inputs: filename
%             transformToSurfRelax [default 0]
%               set this to 1, to transform to coordinate frame used in mrTools
%             volumeSize = [175 256 256]  (only used with previous option)
%             verbose
%
%    outputs: s - a struct w/ the following fields
%                 vtcs
%                 tris
%                 Nvtcs
%                 Ntris
%                 Nedges [empty here]
%
%                 not implemented (nParent,nPatch,patch2parent)
%
%    purpose: load VTK data into a structure that follows the same
%             conventions as loadSurfOFF from mrTools distribution
%
%        e.g: 
%              % fname = 'surf/rh.pial.vtk'
%              fname = 'surf/rh.inflated.vtk';
%              cname = 'surf/rh.curv.vtk'; 
%              s = loadSurfVTK(fname);
%              c = loadSurfVTK(cname);
%              figure
%              patch('vertices',s.vtcs, 'faces',s.tris, ...
%                    'facevertexcdata', sign(c.data(:)) )
%              caxis([-1 1])
%              axis equal, axis off, axis vis3d 
%              shading interp
%              colormap(gray), material dull, light
%
%   see also: loadSurfOFF, loadSurfCaret


% error checking - does file exist, etc.
if nargin < 4
    verbose = 0;
end

% default volume size in RAS coordinates
if nargin < 3
  volumeSize = [176 256 256];
end

if nargin < 2
    % by default, read data as it is in file...
    % is this is set to 1, then subtract 1 from tris and vtcs
    % and shift by 0.5 x, y and z dim respectively to "uncenter" coords in
    % cube of "volumeSize"
    transformToSurfRelax = 0;
end

try
    [vertex,faces,data] = read_vtkData(filename, verbose);
catch
    disp('(loadSurfVTK) problems reading vtk data')
    s = [];
    return
end

if transformToSurfRelax
    % for OFF compatibility need to subtract on, but not here...
    % faces   = faces  -1; % not tris, as these should be read in 1-index
    % vertex  = vertex -1;

    % center image
    vertex(:,1) = vertex(:,1) + volumeSize(1)/2;   % higher, more right
    vertex(:,2) = vertex(:,2) + volumeSize(2)/2;   % higher, more anterior
    vertex(:,3) = vertex(:,3) + volumeSize(3)/2;   % higher, more superior
end

s.vtcs = vertex;
s.tris = faces;
s.Nvtcs = size(s.vtcs,1);
s.Ntris = size(s.tris,1);
s.Nedges = []; % to do / need to check how to get this from VTK
s.data = data(:);

s.filename = filename; % but could be relative

end

function [vertex,face,data] = read_vtkData(filename, verbose)

%% read_vtkData - read data from VTK file.
%
%   [vertex,face, pointData] = read_vtkData(filename, verbose);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%   'pointData' is a 'nb.vert x 1' array specifying data at each vertex [ds]
%
%   2015-02-01 / mods by ds
%
%       e.g:
%             filename = 'metadata/V1-predict.mh.vtk'
%             [vertex,face, pointData] = read_vtkData(filename, 1);
%
%  based on read_vtk.m by
%   Copyright (c) Mario Richtsfeld

if nargin < 2
    verbose = 1;
end

fid = fopen(filename,'r');
if( fid == -1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(3:5), 'vtk')
    error('The file is not a valid VTK one.');
end

%% read header %%
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
nvert = sscanf(str,'%*s %d %*s', 1);

%% read vertices
[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A';

%% read polygons
str = fgets(fid);
str = fgets(fid);

info = sscanf(str,'%c %*s %*s', 1);

if((info ~= 'P') && (info ~= 'V'))
    str = fgets(fid);
    info = sscanf(str,'%c %*s %*s', 1);
end

if(info == 'P')
    
    nface = sscanf(str,'%*s %d %*s', 1);
    
    [A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
    if cnt~=4*nface
        warning('Problem in reading faces.');
    end
    
    A = reshape(A, 4, cnt/4);
    face = transpose( A(2:4,:)+1 );
    
end

if(info ~= 'P')
    face = 0;
end


%% read vertex ids
if(info == 'V')
    
    nv = sscanf(str,'%*s %d %*s', 1);
    
    [A,cnt] = fscanf(fid,'%d %d \n', 2*nv);
    if cnt~=2*nv
        warning('Problem in reading faces.');
    end
    
    A = reshape(A, 2, cnt/2);
    face = transpose( repmat(A(2,:)+1, 3, 1) );
end

if((info ~= 'P') && (info ~= 'V'))
    face = 0;
end

%% and now try to read point dat
str = fgets(fid);

if str == -1
    % no data present beyond this point
    data = 0;
else
    % try reading it out
    pat = '\w+';
    m = regexpi(str, pat, 'match');
    
    % info = sscanf(str,'%c %*s %*s', 1);
    isPointData = strcmpi('POINT_DATA',m);
    if any(isPointData)
        np = str2num(m{find(isPointData)+1});
        
        % two additional lines describing data
        str = fgets(fid);
        str = fgets(fid);
        
        % then the data
        [data,cnt] = fscanf(fid,'%f\n', np);
        %info = sscanf(str,'%c %*s %*s', 1);
    else
        data = 0;
    end
    % and close file
    fclose(fid);
end

end

