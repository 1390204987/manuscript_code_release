function [p_angle, shuffle_angle, shuffle_angle_mean, shuffle_angle_CI, real_angle] = ...
    angle_shuffle_test(allvec, nShuffle)
% Shuffle test for pairwise angles between coding axes.
%
% This version keeps the angle range from 0 to 180 deg.
% It does NOT use abs(cos_angle).
%
% INPUT:
%   allvec    : nNeuron x nAxis matrix
%               each column is one coding axis
%
%   nShuffle : number of shuffles, e.g. 1000
%
% OUTPUT:
%   p_angle            : nAxis x nAxis p-value matrix
%                        p = probability that shuffled angle <= real angle
%
%   shuffle_angle      : nAxis x nAxis x nShuffle shuffled angle distribution
%
%   shuffle_angle_mean : nAxis x nAxis mean shuffled angle
%
%   shuffle_angle_CI   : nAxis x nAxis x 2, 95% shuffle range
%                        (:,:,1) = lower bound
%                        (:,:,2) = upper bound
%
%   real_angle         : nAxis x nAxis real pairwise angle, 0-180 deg

if nargin < 2 || isempty(nShuffle)
    nShuffle = 1000;
end

[nNeuron, nAxis] = size(allvec);

real_angle = nan(nAxis, nAxis);
shuffle_angle = nan(nAxis, nAxis, nShuffle);
p_angle = nan(nAxis, nAxis);

%% Calculate real pairwise angles, 0-180 deg
for i = 1:nAxis
    for j = 1:nAxis

        x = allvec(:, i);
        y = allvec(:, j);

        valid_idx = ~isnan(x) & ~isnan(y);
        x = x(valid_idx);
        y = y(valid_idx);

        if length(x) > 2
            real_angle(i,j) = local_axis_angle_0to180(x, y);
        end
    end
end

%% Shuffle test
for s = 1:nShuffle

    for i = 1:nAxis
        for j = 1:nAxis

            x = allvec(:, i);
            y = allvec(:, j);

            valid_idx = ~isnan(x) & ~isnan(y);
            x = x(valid_idx);
            y = y(valid_idx);

            if length(x) > 2

                % Shuffle neuron identity of the second axis
                y_shuff = y(randperm(length(y)));

                shuffle_angle(i,j,s) = local_axis_angle_0to180(x, y_shuff);
            end
        end
    end
end

%% p value
% Smaller angle means stronger positive alignment.
% p_angle tests whether the real angle is smaller than expected by shuffle.
for i = 1:nAxis
    for j = 1:nAxis

        null_dist = squeeze(shuffle_angle(i,j,:));
        null_dist = null_dist(~isnan(null_dist));

        if ~isempty(null_dist) && ~isnan(real_angle(i,j))
%             p_angle(i,j) = mean(null_dist <= real_angle(i,j));
           p_two_tail = 2 * min( ...
                mean(null_dist <= real_angle(i,j)), ...
                mean(null_dist >= real_angle(i,j)) ...
                );
           p_two_tail = min(p_two_tail, 1);
           p_angle(i,j) = p_two_tail ;
        end
    end
end

%% Summary of shuffle distribution
shuffle_angle_mean = mean(shuffle_angle, 3, 'omitnan');

shuffle_angle_CI = nan(nAxis, nAxis, 2);
shuffle_angle_CI(:,:,1) = prctile(shuffle_angle, 2.5, 3);
shuffle_angle_CI(:,:,2) = prctile(shuffle_angle, 97.5, 3);

end


function angle_deg = local_axis_angle_0to180(x, y)
% Calculate angle between two vectors, keeping the range 0-180 deg.

dot_xy = nansum(x .* y);
norm_xy = sqrt(nansum(x.^2) .* nansum(y.^2));

cos_angle = dot_xy ./ norm_xy;

% Avoid numerical error outside [-1, 1]
cos_angle = min(max(cos_angle, -1), 1);

angle_deg = acosd(cos_angle);

end