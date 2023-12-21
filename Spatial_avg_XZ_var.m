function [var_avg] = Spatial_avg_XZ_var(var,num_points_x,num_points_y,num_points_z, fi, varargin)

if ~isempty(varargin)
    idx = varargin{1}; % idx = 17;    % To avoid near-wall general effects (Larsson 2015) at y/delta_BL = 0.2
else
    idx = fi + 1;      % To avoid near-wall filtering effects
end

num_points_xz     = (num_points_x-(idx-1)*2)*(num_points_z-(idx-1)*2);

var_avg       = zeros(1,num_points_y);




for jj = idx:num_points_y-(idx-1)
    for ii = idx:num_points_x-(idx-1)
        for kk = idx:num_points_z-(idx-1)
            % Compute across entire channel
            aux_j = jj;
            if isnan(var(ii,jj,kk))
                var(ii,jj,kk) = 0;
            end
            var_avg(aux_j)       = var_avg(aux_j)  + var(ii,jj,kk)/num_points_xz;

        end

    end
end


end