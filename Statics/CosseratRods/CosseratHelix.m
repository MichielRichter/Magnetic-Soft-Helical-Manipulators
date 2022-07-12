classdef CosseratHelix
    % Cosserat class for helices with a core cylinder and backbone,
    % starting at onset of the helix (no bottom and top disks,
    % use CosseratCylinder for that). So length should be:
    % total segment length-2*width
    
    % The helix starts on the local x-axis
    
    % IMPORTANT: For now we assume that the material in the helix windings
    % does NOT contribute to the bending stiffness matrix.
    
    properties
        % Assigned
        length;
        
        radius_outer;
        radius_inner;
        
        radius_core_outer;
        radius_core_inner;
        
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
        
        angle_bending_design;
        
        % Derived
        phi;
        number_of_windings;
        pitch;
        height;
        width;
        
        Kb;
        Ks;
        
        J0;
    end
    
    methods
        %% Constructor
        function obj = CosseratHelix(length, radius_outer, radius_inner, ...
                radius_core_outer, radius_core_inner, ...
                radius_backbone_outer, radius_backbone_inner, ...
                density, density_backbone, ...
                residual_flux_density, residual_flux_density_backbone, ...
                elastic_modulus, elastic_modulus_backbone, ...
                shear_modulus, shear_modulus_backbone, angle_bending_design, ...
                number_of_windings)
            
            % Assign properties
            obj.length = length;
            
            obj.radius_outer = radius_outer;
            obj.radius_inner = radius_inner;
            
            obj.radius_core_outer = radius_core_outer;
            obj.radius_core_inner = radius_core_inner;
            
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
            
            obj.angle_bending_design = angle_bending_design;
            obj.number_of_windings = number_of_windings;
            
            % Derived
            obj.width = radius_outer - radius_inner;
            obj.pitch = radius_outer * (angle_bending_design / number_of_windings);
            obj.height = (length - number_of_windings * obj.pitch) / (2 + number_of_windings);
            obj.phi = (obj.height / (obj.height+obj.pitch)) * 2*pi;
            
            % Derived properties
            % area moments of inertia for core 
            moment_area_x_core = pi * (radius_core_outer^4 - radius_core_inner^4) / 4;
            moment_area_y_core = pi * (radius_core_outer^4 - radius_core_inner^4) / 4;
            moment_area_z_core = moment_area_x_core + moment_area_y_core;
            
            % Derived properties
            % area moments of inertia for core 
            moment_area_x_backbone = pi * (radius_backbone_outer^4 - radius_backbone_inner^4) / 4;
            moment_area_y_backbone = pi * (radius_backbone_outer^4 - radius_backbone_inner^4) / 4;
            moment_area_z_backbone = moment_area_x_backbone + moment_area_y_backbone;
            
            % Bending stiffness matrices
            stiffness_bending_core = diag([elastic_modulus, elastic_modulus, shear_modulus]) * ...
                diag([moment_area_x_core, moment_area_y_core, moment_area_z_core]);
            stiffness_bending_backbone = diag([elastic_modulus_backbone, elastic_modulus_backbone, shear_modulus_backbone]) * ...
                diag([moment_area_x_backbone, moment_area_y_backbone, moment_area_z_backbone]);
            
            % Shear stiffness matrices
            stiffness_shear_core = diag([shear_modulus, shear_modulus, elastic_modulus]) * ...
                (pi * (radius_core_outer^2-radius_core_inner^2));
            stiffness_shear_backbone = diag([shear_modulus_backbone, shear_modulus_backbone, elastic_modulus_backbone]) * ...
                (pi * (radius_backbone_outer^2-radius_backbone_inner^2));
            
            % Bending and shear stiffness matrices
            obj.Kb = stiffness_bending_core + stiffness_bending_backbone;
            obj.Ks = stiffness_shear_core + stiffness_shear_backbone;
            
        end
        
        function Y = forwardIntegrate(obj, y0, number_of_discretizations)
            % Input base states and number of discretizations
            % Output shape solution matrix
            
            % Initialize
            N = number_of_discretizations;
            delta_s = obj.length / N;
            delta_beta = (2*pi*obj.number_of_windings*delta_s) / (obj.length); %\beta(L) = 2*pi*W, \beta(s) = s/L * 2*pi*W, \dbeta = ds/L * 2*pi*W
            Y = NaN(size(y0,1), N+1);
            Y(:,1) = y0;
            
            % Integrate
            beta = 0;
            for i = 1:N
                y1 = Y(:,i);
                k1 = delta_s * obj.ODE(y1, delta_s, beta);
                y2 = y1 + 0.5*k1;
                k2 = delta_s * obj.ODE(y2, delta_s, beta);
                y3 = y1 + 0.5*k2;
                k3 = delta_s * obj.ODE(y3, delta_s, beta);
                y4 = y1 + k3;
                k4 = delta_s * obj.ODE(y4, delta_s, beta);
                Y(:,i+1) = y1 + (k1 + 2*k2 + 2*k3 + k4)/6;
                
                beta = beta + delta_beta;
            end
        end
        
        function dy = ODE(obj, y, delta_s, beta)
            
            global B_global Bgrad_global
            
            % Extract
            p = y(1:3);
            q = y(4:7);     R = quat2rotm(q');
            n = y(8:10);
            m = y(11:13);
            
            % Update stiffness matrix
            stiffness_shear = obj.Ks;
            stiffness_bending = obj.Kb;
            
            %stiffness_bending = stiffness_bending + ...
            %    diag([obj.elastic_modulus, obj.elastic_modulus, obj.shear_modulus])*(rotz(beta)*obj.J0*rotz(beta)');
            
            % Compute shear variables
            v = stiffness_shear\R'*n + [0;0;1];
            u = stiffness_bending\R'*m + [0;0;0];
            
            % Initialize magnetic force and torque vector
            force_magnetic = 0;
            torque_magnetic = 0;
            
            %------------------------------------------------------------
            % MAGNETIC FORCES AND TORQUES ON CORE CYLINDER
            %------------------------------------------------------------
            Nangle = 36;    dangle = 2*pi / Nangle;
            Nr = 1;         dr = (obj.radius_core_outer - obj.radius_core_inner) / Nr;
            for theta = 0 : dangle : 2*pi-dangle % Loop over angle
                for r = obj.radius_core_inner : dr : obj.radius_core_outer-dr % Loop over radius
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
            
            %------------------------------------------------------------
            % MAGNETIC FORCES AND TORQUES ON HELIX SECTOR
            %------------------------------------------------------------
            Nangle = 36;    dangle = obj.phi / Nangle;
            Nr = 12;        dr = (obj.radius_outer - obj.radius_inner) / Nr;
            for theta = beta : dangle : beta+obj.phi-dangle % Loop over angle
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
            
            %------------------------------------------------------------
            % STATE GRADIENTS
            %------------------------------------------------------------
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

