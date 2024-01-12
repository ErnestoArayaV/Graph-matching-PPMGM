function [matches ] = merge_ground_truth( gt_M, gt_N )
%MERGE_GROUND_TRUTH Creates GT from kidi to kidj with information from
%kidi_ref and kidj_ref

maximum = max(max(gt_M(:,2)), max(gt_N(:,2)));

bins_M = [ (1:maximum)', zeros(maximum, 1) ];
bins_N = [ (1:maximum)', zeros(maximum, 1) ];

% sort in all matches
for i=1:size(gt_M,1)
   bins_M(gt_M(i,2),2) = gt_M(i,1); 
end

for i=1:size(gt_N,1)
   bins_N(gt_N(i,2),2) = gt_N(i,1); 
end

% find transitive matches

ind_M = bins_M(bins_M(:,2) ~= 0, 1);
ind_N = bins_N(bins_N(:,2) ~= 0, 1);

joint_ind = intersect(ind_M, ind_N);

matches = [bins_M(joint_ind, 2), bins_N(joint_ind, 2) ];

end

