function [ thresholdedLabel ] = threshold_label(label, thr)
%threshold_label - threshold a freesurfer label array
%
%      usage: [ thresholdedLabel ] = threshold_label( label, thr )
%         by: lpzds1
%       date: May 21, 2015
%        $Id$
%     inputs: label, thr
%    outputs: thresholdedLabel
%
%    purpose: take an array loaded with read_label and return tresholded
%             version of the array
%
%             [ id, x,y,z, stat_value ]



% threshold V1 probability to some level that was determined by Hind's
% study
if ieNotDefined('thr'), thr = 0.8; end

thresholdedLabel = label(label(:,5) > thr , :);

end