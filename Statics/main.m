clear all; clc;
addpath('./MagneticField', './CosseratRods')

%% External magnetic field
global B_global Bgrad_global

% Magnetic field
field_range = (0:1:20) * 1e-3;


%% Internal permanent magnet values
% Length
length_magnet = 1e-3;

% Radii
radius_magnet_outer = 2e-3;
radius_magnet_inner = 0.75e-3;

% Magnetic properties
residual_flux_density_magnet = 1.35; 

% Density
density_magnet = 7000;

% Elastic and shear modulus
elastic_modulus_magnet = 152e9;
shear_modulus_magnet = elastic_modulus_magnet / (2*(1+0.3));


%% Backbone values
% Radii
radius_backbone_outer = 0.75e-3;
radius_backbone_inner = 0.25e-3;

% Magnetic properties
residual_flux_density_backbone = 0;

% Density
density_backbone = 1000;

% Elastic and shear modulus
shore_silicone = 50;
elastic_modulus_backbone = 0*(0.0981*(56 + 7.66*shore_silicone) / (0.137505*(254-2.54*shore_silicone))) * 1e6; 
shear_modulus_backbone = elastic_modulus_backbone / (2*(1+0.4));


%% Helix values
% Helix length and boundary deflection
length_helix = ( 4 ) * 1e-3;
bending_angle_design = deg2rad( 40 );

% Number of windings
number_of_windings = 2;

% Radii
radius_helix_outer = radius_magnet_outer;
radius_helix_inner = 0.95e-3;
radius_core_outer = radius_helix_inner;
radius_core_inner = radius_backbone_outer;

% Helix pitch and width
pitch = radius_helix_outer * (bending_angle_design / number_of_windings);
height = (length_helix - number_of_windings * pitch) / (2 + number_of_windings);

% Magnetic properties
volume_fraction_magnetic = 0.25;
residual_flux_density_helix = volume_fraction_magnetic * 1;

% Density
density_helix = volume_fraction_magnetic*7000 + (1-volume_fraction_magnetic)*1000;

% Elastic and shear modulus
shore_pdms = 33;
E_pdms = (0.0981*(56 + 7.66*shore_pdms) / (0.137505*(254-2.54*shore_pdms))) * 1e6;
shear_modulus_pdms = E_pdms/(2*(1+0.4));
shear_modulus_helix = (shear_modulus_pdms * ...
    exp(2.5*volume_fraction_magnetic / (1-1.35*volume_fraction_magnetic)));
elastic_modulus_helix = 2*(1+0.4) * shear_modulus_helix;


%% Solve
% ---------------------------------------------------------
% MAKE COSSERAT RODS AND SEGMENT
% ---------------------------------------------------------
Magnet = CosseratCylinder(length_magnet, radius_magnet_outer, ...
    radius_magnet_inner, radius_backbone_outer, radius_backbone_inner, ...
    density_magnet, density_backbone, residual_flux_density_magnet, ...
    residual_flux_density_backbone, elastic_modulus_magnet, elastic_modulus_backbone, ...
    shear_modulus_magnet, shear_modulus_backbone);

Disc = CosseratCylinder(height, radius_magnet_outer, ...
    radius_magnet_inner, radius_backbone_outer, radius_backbone_inner, ...
    density_helix, density_backbone, residual_flux_density_helix, ...
    residual_flux_density_backbone, elastic_modulus_helix, elastic_modulus_backbone, ...
    shear_modulus_helix, shear_modulus_backbone);

Helix = CosseratHelix(length_helix-2*height, radius_helix_outer, ...
    radius_helix_inner, radius_core_outer, radius_core_inner, ...
    radius_backbone_outer, radius_backbone_inner, ...
    density_helix, density_backbone, residual_flux_density_helix, ...
    residual_flux_density_backbone, elastic_modulus_helix, ...
    elastic_modulus_backbone, shear_modulus_helix, shear_modulus_backbone, ...
    bending_angle_design, number_of_windings);

Segment = {Magnet,Disc,Helix,Disc,Magnet};%, Disc,Helix,Disc,Magnet};
number_of_discretizations = [10 10 180 10 10];% 10 180 10 10];
positions_magnets_cosserat_array = [1 5];% 8];

% Storage variables
Deflections = NaN(1,size(field_range,2));

% Guess proximal boundary variables
force_internal_base = [0; 0; 0];
moment_internal_base = [0; 0; 0];

% Solve
cycle = 1;
for i = 1:size(field_range,2)
    
    % Keep track of iteration
    disp(['Iteration: ' num2str(i) ' / ' num2str(size(field_range,2))]);
    
    % Set external field and field gradient
    B_global = field_range(i) * (rotx(-90) * [0; 0; 1]);
    Bgrad_global = zeros(3,3);
    
    % Solve statics
    tic
    [Y, solution] = levenbergMarquardt([0;0;0], [1;0;0;0], force_internal_base, ...
        moment_internal_base, Segment, number_of_discretizations, positions_magnets_cosserat_array);
    toc
    % Update guess
    force_internal_base = solution(1:3);
    moment_internal_base = solution(4:6);
    
    % Compute deflection angle
    orientation_tip = quat2rotm(Y(4:7,end)');
    orientation_z_tip = orientation_tip(:,3);
    deflection = rad2deg(acos(dot(orientation_z_tip/norm(orientation_z_tip),[0;0;1])));
    
    % Store deflection angle
    Deflections(i) = deflection;
    
    cycle = cycle + 1;
end

