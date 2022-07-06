clear all;
addpath('./MagneticFields');

%% Simulation variables
number_of_magnets = 3;
angle_bending_range = deg2rad( 1:40 );
length_magnet = 1e-3;
length_helix = 4e-3;


%% Storage variables
Torque = zeros(3, length(angle_bending_range)); % Sum of torques on tip magnet due to field and gradient
Torque_magnet_field = zeros(3, length(angle_bending_range)); % Torque on tip magnet due to field
Torque_magnet_gradient = zeros(3, length(angle_bending_range)); % Torque on tip magnet due to field gradient
Force = zeros(3, length(angle_bending_range));


%% Simulate
iteration = 1;
for angle_bending_segment = angle_bending_range

    disp([num2str(iteration) ' / ' num2str(size(angle_bending_range,2))]);
    
    % Compute axial directions of magnets
    r = length_helix / angle_bending_segment; % Helix arc radius
    angles = (0:number_of_magnets-1) * angle_bending_segment; % Angles (absolute) at which magnets are located
    directions = [zeros(1,number_of_magnets); r*cos(angles); -r*sin(angles)]; % Axial direction vector of magnets

    % Compute center positions of magnets
    positions = zeros(3, number_of_magnets); % initialize
    matrix_transform_magnet = [eye(3) [0;length_magnet/2;0]; zeros(1,3) 1]; % Transformation matrix for magnets
    matrix_transform_arc = [ ...
        [1 0 0; 0 cos(-angle_bending_segment) -sin(-angle_bending_segment); 0 sin(-angle_bending_segment) cos(-angle_bending_segment)], ...
        [0; r*sin(angle_bending_segment); r*cos(angle_bending_segment)-r];
        0 0 0 1]; % Transformation matrix for helix
    matrix_transform_initial = [eye(3) zeros(3,1); zeros(1,3) 1]; % Initial center pose of first magnet
    for i = 2:number_of_magnets
        matrix_transform_initial = matrix_transform_initial * matrix_transform_magnet * matrix_transform_arc * matrix_transform_magnet; % Transform to next magnet center position
        positions(:,i) = matrix_transform_initial(1:3,4);
    end

    % Make cell array of ring magnets
    magnets_array = cell(1, number_of_magnets);
    for i = 1:number_of_magnets
        magnets_array{i} = RingMagnet(positions(:,i), directions(:,i));
    end

    % Compute torque and force
    magnet_target = magnets_array{number_of_magnets}; % Get target magnet
    for i = 1 : number_of_magnets-1 % Loop over all preceding magnets
        % Compute
        magnet_source = magnets_array{i}; % Get source magnet
        [torque, force, torque_magnet_field, torque_magnet_gradient] = magnet_target.get_torques_and_forces(magnet_source); % Compute torque and force on target magnet from source magnet
        % Store
        Torque(:,iteration) = Torque(:,iteration) + torque;
        Force(:,iteration) = Force(:,iteration) + force;
        Torque_magnet_field(:,iteration) = Torque_magnet_field(:,iteration) + torque_magnet_field;
        Torque_magnet_gradient(:,iteration) = Torque_magnet_gradient(:,iteration) + torque_magnet_gradient;
    end

    iteration = iteration + 1;
end


%% Plot bending torques from field, gradient, and combination

% Variables
lw = 6; % line width
fs = 40; % font size
tscale = 1e3; % Scaling factor for torque: 1e3 -> Nm to Nmm

% Convert relative angular deflections from radians to degrees
angles_degrees = rad2deg(angle_bending_range);

% Plot
figure

yyaxis left
plot(angles_degrees, Torque_magnet_field(1,:)*tscale,'r:','LineWidth',lw); % Bending torques due to field
hold on; grid on;
plot(angles_degrees, Torque_magnet_gradient(1,:)*tscale,'r--','LineWidth',lw); % Bending torques due to field gradient
ylim([-50e-3 50e-3]);
ylabel({'Torque (Nmm)'}, 'Interpreter', 'Latex')


yyaxis right
plot(angles_degrees, Torque(1,:)*tscale,'b','LineWidth',lw); % Bending torque due to both field and gradient
ylim([-8.5e-3 8.5e-3]);
ylabel({'Torque (Nmm)'}, 'Interpreter', 'Latex')
xlabel({'Segment deflection (deg)'}, 'Interpreter', 'Latex')

% Legend and title
legend({'$\tau_{x,B}$', '$\tau_{x,B_\nabla}$', '$\tau_{x,B}+\tau_{x,B_\nabla}$'},'Interpreter','Latex','Location','northwest')
title(['Magnet ' num2str(number_of_magnets)],'Interpreter','Latex')

% Color axes
ax = gca;   ax.YAxis(1).Color = 'r';    ax.YAxis(2).Color = 'b';
set(gca, 'FontSize', fs)



%% Plot field from ring magnets
% Make 2D coordinate mesh
y = -5e-3:0.05e-3:20e-3;
z = -20e-3:0.05e-3:5e-3;
[Y,Z] = meshgrid(y,z);

% Initialize mesh for storing field directions and norm
field_norm = zeros(size(Y));
field_Bx = zeros(size(Y));
field_By = zeros(size(Y));
field_Bz = zeros(size(Y));

% Loop through 2D coordinate mesh
for i = 1:size(Y,1)
    for j = 1:size(Z,2)
        % Position in yz-plane
        p = [0; Y(i,j); Z(i,j)];

        % Loop over all magnets
        for k = 1:number_of_magnets

            % Compute field
            magnet_source = magnets_array{k}; % Get source magnet
            [field, ~] = magnet_source.get_field_and_gradient(p); % Get field
            
            % Store
            field_norm(i,j) = field_norm(i,j) + norm(field);
            field_Bx(i,j) = field_Bx(i,j) + field(1);
            field_By(i,j) = field_By(i,j) + field(2);
            field_Bz(i,j) = field_Bz(i,j) + field(3);
        end

    end
end

% Convert 2D mesh from meters to millimeters
Y = Y*1e3;
Z = Z*1e3;

% Plot surface of field magnitude
figure
surf(Y,Z,field_norm, 'EdgeAlpha',0); colormap('jet'); 
clim([0 0.3]) % Colorbar limit
hold on

% Plot contour of each magnet
for k = 1:number_of_magnets; magnets_array{k}.draw_magnet; end

% Plot quiver of field
step = 10; % dilute number of data points
quiver3(Y(1:step:end,1:step:end), ...
    Z(1:step:end,1:step:end), ...
    3*ones(size(Y(1:step:end,1:step:end))), ...
    field_By(1:step:end,1:step:end), ...
    field_Bz(1:step:end,1:step:end),...
    field_Bx(1:step:end,1:step:end),'r');

axis equal;
xlabel('$y$ (mm)', 'Interpreter', 'Latex'); ylabel('$z$ (mm)', 'Interpreter', 'Latex'); colorbar;
xlim([min(min(Y)) max(max(Y))]);
ylim([min(min(Z)) max(max(Z))]);
set(gca, 'FontSize', fs)

% Place camera on top
campos([0 0 20])