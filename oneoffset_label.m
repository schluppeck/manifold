function [ label ] = oneoffset_label(label)
%oneoffset_label - make label array 1-offset (from 0-offset)
%
%      usage: [ label ] = oneoffset_label( label )
%         by: lpzds1
%       date: May 21, 2015
%        $Id$
%     inputs: label
%    outputs: label
%
%    purpose: read_label and assorted functions return vertices with
%             0-offset, we need 1-offset in Matlab
%

label(:,1) = label(:,1) + 1; % 1-offset

end