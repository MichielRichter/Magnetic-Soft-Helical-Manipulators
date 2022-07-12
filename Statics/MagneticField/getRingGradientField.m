function [field, field_gradient] = getRingGradientField(p)

global array_ring_magnets

% Initialize
field = zeros(3,1);
field_gradient = zeros(3,3);


% Loop through ring magnets and add the field
for i = 1:size(array_ring_magnets, 2)
    magnet = array_ring_magnets{i};
    [field_magnet, field_gradient_magnet] = magnet.getRingGradientField(p);
    field = field + field_magnet;
    field_gradient = field_gradient + field_gradient_magnet;
end

end

