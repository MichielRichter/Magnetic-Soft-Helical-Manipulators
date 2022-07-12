classdef RingMagnet
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        transform_center
    end
    
    methods
        function obj = RingMagnet(position_surface_global, orientation_surface_global)
            % Assign pose of magnet
            obj.transform_center = [quat2rotm(orientation_surface_global') position_surface_global;
                zeros(1,3) 1];
        end
        
        function [field, field_gradient] = getRingGradientField(obj,position_global)
            
            % Transform global position to local frame
            position_local = obj.transform_center \ [position_global; 1];
            x_local = position_local(1);
            y_local = position_local(2);
            z_local = position_local(3);
            
            % Check if the z-position is further than some threshold and
            % compute local field and gradient
            if z_local < 0.5e-3 % Only consider positions that are in front of the magnet. But not too close
                field_local = zeros(3,1);
                field_gradient_local = zeros(3,3);
            else 
                field_local = 1.8*get_ring_field(x_local, y_local, z_local); % Factor 1.8 scales the model to an N45 magnet
                field_gradient_local = 1.8*get_ring_field_gradient(x_local, y_local, z_local);
            end
            
            % Transform field to global coordinates
            field = obj.transform_center(1:3,1:3) * field_local;
            field_gradient = obj.transform_center(1:3,1:3) * field_gradient_local * obj.transform_center(1:3,1:3)';
            
            
        end
    end
end

