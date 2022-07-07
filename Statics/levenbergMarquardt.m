function [Y, wrench_internal_base_guess] = levenbergMarquardt(p0, q0, n0, m0, cosserat_rod_list, number_of_discretizations, positions_magnets_cosserat_array)

% Make an array of ring magnet classes
global array_ring_magnets
array_ring_magnets = {RingMagnet(p0, q0)}; % Base magnet

% Guessed internal wrench at base
wrench_internal_base_guess = [n0; m0];

% Distal boundary condition
boundary_condition_distal = [0; 0; 0; 0; 0; 0]; % Unloaded

% Number of rods
cosserat_rods_total = size(cosserat_rod_list,2);


Y = shoot([p0; q0; wrench_internal_base_guess]);
residuals = getResiduals(Y);
residuals_norm_minus = norm(residuals);
residuals_norm_plus = 1e5;
residuals_norm_change = residuals_norm_plus - residuals_norm_minus;
lambda = 1e-4;

while (abs(residuals_norm_change) > 1e-8)
    jacobian = getJacobian(residuals);
    lhs = jacobian' * jacobian;
    rhs = -jacobian' * residuals;
    lambda = lambda / 4;
    lhs = lhs + lambda * eye(size(jacobian,2));
    step = lhs \ rhs;

    wrench_internal_base_guess = wrench_internal_base_guess + step;
    Y = shoot([p0; q0; wrench_internal_base_guess]);
    residuals = getResiduals(Y);
    residuals_norm_plus = norm(residuals);
    wrench_internal_base_guess = wrench_internal_base_guess - step;

    if (residuals_norm_plus < residuals_norm_minus)
        wrench_internal_base_guess = wrench_internal_base_guess + step;
        residuals_norm_change = residuals_norm_plus - residuals_norm_minus;
        residuals_norm_minus = residuals_norm_plus;
    else
        while (residuals_norm_plus > residuals_norm_minus)
            lambda = lambda * 4;
            lhs = lhs + lambda * eye(size(jacobian,2));
            step = lhs \ rhs;
            wrench_internal_base_guess = wrench_internal_base_guess + step;
            Y = shoot([p0; q0; wrench_internal_base_guess]);
            residuals = getResiduals(Y);
            residuals_norm_plus = norm(residuals);
            wrench_internal_base_guess = wrench_internal_base_guess - step;
        end
        wrench_internal_base_guess = wrench_internal_base_guess + step;
        residuals_norm_change = residuals_norm_plus - residuals_norm_minus;
        residuals_norm_minus = residuals_norm_plus;
    end

end



% Shoot forward through cosserat rods
    function Y = shoot(y0)
        
        % Array of ring magnets
        array_ring_magnets = {RingMagnet(p0, q0)}; % Base magnet
        
        % Initialize
        y_base = y0;
        Y = [];
        
        % Forward integrate
        for i = 1:cosserat_rods_total
            cosserat_rod = cosserat_rod_list{i};
            Y = [ Y, cosserat_rod.forwardIntegrate(y_base, number_of_discretizations(i)) ];
            y_base = Y(:,end);
            
            % If we just passed a ring magnet, add to ring magnet list
            if ismember(i, positions_magnets_cosserat_array)
                % Get ring magnet base pose
                p_base = Y(1:3,end);    q_base = Y(4:7,end);
                % Add to cell
                array_ring_magnets{size(array_ring_magnets,2)+1} = RingMagnet(p_base, q_base);
            end
        end
    end


% Get boundary residuals
    function residuals_distal = getResiduals(Y)
        residuals_distal = Y(8:13,end) - boundary_condition_distal;
    end


% Compute Jacobian matrix
    function jacobian = getJacobian(residuals)
        jacobian = NaN(size(residuals,1), size(wrench_internal_base_guess,1));
        new_guess = wrench_internal_base_guess;
        increment = 1e-8;
        for i = 1:size(new_guess,1)
            new_guess(i) = new_guess(i) + increment;
            new_y0 = [p0; q0; new_guess];
            new_Y = shoot(new_y0);
            new_residuals = getResiduals(new_Y);
            jacobian_column = (new_residuals - residuals) / increment;
            jacobian(:,i) = jacobian_column;
            new_guess(i) = new_guess(i) - increment;
        end
    end

end