function [field, field_gradient] = getExternalGradientField(p)

field_magnitude = 25e-3; 

field_direction = [0; 1; 0];

field = field_magnitude * field_direction;
field_gradient = zeros(3,3);

end

