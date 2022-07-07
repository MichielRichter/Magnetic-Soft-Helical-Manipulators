classdef CosseratCylinder
    % Cosserat rod class for cylinders with an additional central backbone
    
    properties
        % Assigned
        length;
        
        radius_outer;
        radius_inner;
        
        radius_backbone_outer;
        radius_backbone_inner;
        
        density;
        density_backbone;
        
        residual_flux_density;
        residual_flux_density_backbone;
        
        elastic_modulus;
        elastic_modulus_backbone;
        
        shear_modulus;
        shear_modulus_backbone;
        
        % Derived
        Kb;
        Ks;
    end
    
    methods
        
        %% Constructor
        function obj = CosseratCylinder(length, radius_outer, ...
                radius_inner, radius_backbone_outer, ...
                radius_backbone_inner, density, density_backbone, ...
                residual_flux_density, residual_flux_density_backbone, ...
                elastic_modulus, elastic_modulus_backbone, ...
                shear_modulus, shear_modulus_backbone)
            
            % Assign properties
            obj.length = length;
            
            obj.radius_outer = radius_outer;
            obj.radius_inner = radius_inner;
            
            obj.radius_backbone_outer = radius_backbone_outer;
            obj.radius_backbone_inner = radius_backbone_inner;
            
            obj.density = density;
            obj.density_backbone = density_backbone;
            
            obj.residual_flux_density = residual_flux_density;
            obj.residual_flux_density_backbone = residual_flux_density_backbone;
            
            obj.elastic_modulus = elastic_modulus;
            obj.elastic_modulus_backbone = elastic_modulus_backbone;
            
            obj.shear_modulus = shear_modulus;
            obj.shear_modulus_backbone = shear_modulus_backbone;
            
            % Derived properties
            % area moments of inertia for composite core and backbone
            moment_area_x = pi * (radius_outer^4 - radius_backbone_inner^4) / 4;
            moment_area_y = pi * (radius_outer^4 - radius_backbone_inner^4) / 4;
            moment_area_z = moment_area_x + moment_area_y;
            
            % Elastic modulus for composite
            volume_core = pi * (radius_outer^2 - radius_inner^2) * length;
            volume_backbone = pi * (radius_backbone_outer^2 - radius_backbone_inner^2) * length;
            volume_fraction_core = volume_core / (volume_core + volume_backbone);
            volume_fraction_backbone = volume_backbone / (volume_core + volume_backbone);
            elastic_modulus_composite = volume_fraction_core*elastic_modulus + volume_fraction_backbone*elastic_modulus_backbone;
            shear_modulus_composite = volume_fraction_core*shear_modulus + volume_fraction_backbone*shear_modulus_backbone;
            
            % Stiffness matrices
            stiffness_bending = diag([elastic_modulus_composite, elastic_modulus_composite, shear_modulus_composite]) * ...
                diag([moment_area_x, moment_area_y, moment_area_z]);
            stiffness_shear = diag([shear_modulus_composite, shear_modulus_composite, elastic_modulus_composite]) * ...
                (pi * (radius_outer^2-radius_backbone_inner^2));

            % Bending and shear stiffness matrices
            obj.Kb = stiffness_bending;
            obj.Ks = stiffness_shear;

        end
        
        function Y = forwardIntegrate(obj, y0, number_of_discretizations)
            % Input base states and number of discretizations
            % Output shape solution matrix
            
            % Initialize
            N = number_of_discretizations;
            delta_s = obj.length / N;
            Y = NaN(size(y0,1), N+1);
            Y(:,1) = y0;
            
            % Integrate
            for i = 1:N
                y1 = Y(:,i);
                k1 = delta_s * obj.ODE(y1, delta_s);
                y2 = y1 + 0.5*k1;
                k2 = delta_s * obj.ODE(y2, delta_s);
                y3 = y1 + 0.5*k2;
                k3 = delta_s * obj.ODE(y3, delta_s);
                y4 = y1 + k3;
                k4 = delta_s * obj.ODE(y4, delta_s);
                Y(:,i+1) = y1 + (k1 + 2*k2 + 2*k3 + k4)/6;
            end
            
        end
        
        function dy = ODE(obj, y, delta_s)
            
            global B_global Bgrad_global
            
            % Extract
            p = y(1:3);
            q = y(4:7);     R = quat2rotm(q');
            n = y(8:10);
            m = y(11:13);
            
            % Compute shear variables
            v = obj.Ks\R'*n + [0;0;1];
            u = obj.Kb\R'*m + [0;0;0];
            
            % Discretize cross-section of cylinder (backbone non-magnetic)
            Nangle = 36;    dangle = 2*pi/Nangle;
            Nr = 12;        dr = (obj.radius_outer - obj.radius_inner) / Nr;
            
            % Initialize magnetic force and torque vector
            force_magnetic = 0;
            torque_magnetic = 0;
            
            % Compute forces and torques
            for theta = 0 : dangle : 2*pi-dangle % Loop over angle
                for r = obj.radius_inner : dr : obj.radius_outer-dr % Loop over radius
                    % Position in local frame
                    px_local = (r+0.5*dr) * cos(theta+0.5*dangle);
                    py_local = (r+0.5*dr) * sin(theta+0.5*dangle);
                    pz_local = 0;
                    
                    % Position in global frame
                    p_global = p + R * [px_local; py_local; pz_local];
                    
                    % Sector area and volume
                    area = 0.5 * dangle * ((r+dr)^2 - r^2);
                    volume = area * delta_s;
                    
                    % Magnetic moment
                    moment_magnetic = ((obj.residual_flux_density * volume)/(pi*4e-7)) * R(:,3);
                    
                    % Get magnetic field
                    [field_ring, field_gradient_ring] = getRingGradientField(p_global);
                    [field_external, field_gradient_external] = getExternalGradientField(p_global);
                    field_global = field_ring + B_global;
                    
                    % Get magnetic field gradients
                    field_gradient_matrix_global = field_gradient_ring + Bgrad_global;
                    field_gradient_global = [field_gradient_matrix_global(1,1);
                        field_gradient_matrix_global(1,2); field_gradient_matrix_global(1,3);
                        field_gradient_matrix_global(2,2); field_gradient_matrix_global(2,3)];
                    
                    % Force and torque
                    force = [moment_magnetic(1) moment_magnetic(2) moment_magnetic(3) 0 0;
                        0 moment_magnetic(1) 0 moment_magnetic(2) moment_magnetic(3);
                        -moment_magnetic(3) 0 moment_magnetic(1) -moment_magnetic(3) moment_magnetic(2)] * ...
                        field_gradient_global;
                    torque = cross(moment_magnetic, field_global) + cross(p_global-p, force);
                    
                    % Add to total force and torque
                    force_magnetic = force_magnetic + force;
                    torque_magnetic = torque_magnetic + torque;
                end
            end
            
            % Partial derivatives
            p_prime = R * v;
            q_prime = 0.5 * [-q(2:4)'; q(1)*eye(3)-obj.skew(q(2:4))] * R * u;
            n_prime = -force_magnetic / delta_s;
            m_prime = cross(-p_prime, n) - (torque_magnetic / delta_s);
            dy = [p_prime; q_prime; n_prime; m_prime];
        end
        
        function S = skew(obj, x)
           S = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0]; 
        end
    end
end

