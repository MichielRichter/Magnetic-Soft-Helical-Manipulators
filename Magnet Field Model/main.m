clear all; clc;
%{
    This program requires running COMSOL Multiphysics with Matlab. Field estimates from a COMSOL MODEL
    are obtained. Then the program aims to fit an axially symmetric field model to the estimates. For the
    model, see Richter et al. (2021) "Multi-point orientation control of discretely magnetized continuum
    manipulators". The magnet axis is assumed the z-axis. 
%}

%% Make model and coordinates grid
model = Ring_Magnet_COMSOL;

% Coordinates
z1 = (0.55 : 0.5 : 60) * 1e-3;
y1 = (0 : 0.5 : 50) * 1e-3;
z2 = (0 : 0.5 : 20) * 1e-3;
y2 = (2.25 : 0.5 : 20) * 1e-3;

Z = [z1];
Y = [y1];

coord_3d = [];
for z = Z
    for y = Y
        coord_3d = [coord_3d; [0 y z]];
    end
end


%% Field data
Bx = mphinterp(model, 'mfnc.Bx', 'coord', [coord_3d(:,1), coord_3d(:,2), coord_3d(:,3)]', 'dataset', 'dset1');
By = mphinterp(model, 'mfnc.By', 'coord', [coord_3d(:,1), coord_3d(:,2), coord_3d(:,3)]', 'dataset', 'dset1');
Bz = mphinterp(model, 'mfnc.Bz', 'coord', [coord_3d(:,1), coord_3d(:,2), coord_3d(:,3)]', 'dataset', 'dset1');
field_3d = [Bx; By; Bz]';


%% Fitting
fitOrder = 6;
dispVars = 2;
[field_3d_mod, ~, Bsym, R2Y, R2Z] = start_fit(coord_3d, field_3d, 0, 1, fitOrder, dispVars);
syms x y z
dxBx = diff(Bsym(1),x);
dxBy = diff(Bsym(2),x);
dxBz = diff(Bsym(3),x);
dyBy = diff(Bsym(2),y);
dyBz = diff(Bsym(3),y);
dzBz = diff(Bsym(3),z);
dBsym = [dxBx dxBy dxBz;
    dxBy dyBy dyBz;
    dxBz dyBz dzBz];
field = matlabFunction(Bsym);
fieldGrad = matlabFunction(dBsym);

%matlabFunction(Bsym,'File','get_ring_field_z0_to_6_y3_to_6','Vars',{x,y,z});
%matlabFunction(dBsym,'File','get_ring_gradient_z0_to_6_y3_to_6','Vars',{x,y,z});
%% Plotting
N = 120;
ymin = -3e-3;
ymax = 3e-3;
ystep = (ymax-ymin)/N;
zmin = 3e-3;
zmax = 5e-3;
zstep = (zmax-zmin)/N;
y = [ymin:ystep:ymax];
z = [zmin:zstep:zmax];
[Y,Z] = meshgrid(y,z);

Bnorm = NaN(size(Y));
By = NaN(size(Y));
Bz = NaN(size(Y));
dxBx = NaN(size(Y));
dyBx = NaN(size(Y));
dzBx = NaN(size(Y));
dyBy = NaN(size(Y));
dzBy = NaN(size(Y));
dzBz = NaN(size(Y));

% Get data
f = waitbar(0,'Getting field data');
count = 1; maxCount = size(Y,1)*size(Y,2);
for i = 1:size(Y,1)
    for j = 1:size(Y,2)
        waitbar(count/maxCount,f,['Getting field data (' num2str(count) '/' num2str(maxCount) ')']);
        
        % Coordinate
        y = Y(i,j);
        z = Z(i,j);
        
        % Field
        B = field(0,y,z);
        dB = fieldGrad(0,y,z);
        
        % Assign
        Bnorm(i,j) = norm(B);
        By(i,j) = B(2);
        Bz(i,j) = B(3);
        dxBx(i,j) = dB(1,1);
        dyBx(i,j) = dB(1,2);
        dzBx(i,j) = dB(1,3);
        dyBy(i,j) = dB(2,2);
        dzBy(i,j) = dB(2,3);
        dzBz(i,j) = dB(3,3);
        
    end
end
close(f)


%%
% Plot
figure(1)
surf(Y*1e3, Z*1e3, Bnorm*1e3,'EdgeAlpha',0); hold on;
colorbar; colormap('jet')
c = 16;
h1 = quiver3(Y(1:c:end,1:c:end)*1e3, Z(1:c:end,1:c:end)*1e3, ...
    100*ones(size(Z(1:c:end,1:c:end))), 0.6*By(1:c:end,1:c:end), ...
    0.6*Bz(1:c:end,1:c:end), zeros(size(By(1:c:end,1:c:end))), ...
    'LineWidth',2.25, 'Color', [0,0,0]);
set(h1,'AutoScale','on', 'AutoScaleFactor', 0.2)
xlabel('$y\mathrm{[mm]}$','Interpreter','latex'); ylabel('$z\mathrm{[mm]}$','Interpreter','Latex'); title('Field norm','Interpreter','Latex');
axis([min(min(Y*1e3)) max(max(Y*1e3)) min(min(Z*1e3)) max(max(Z*1e3)) -100 100])    
view(90,90); axis square;
set(gca,'FontSize',25)

figure(2)
surf(Y*1e3, Z*1e3, dyBy,'EdgeAlpha',0)
colorbar; colormap('jet')
xlabel('$y\mathrm{[mm]}$','Interpreter','latex'); ylabel('$z\mathrm{[mm]}$','Interpreter','Latex'); title('dyBy');
axis([min(min(Y*1e3)) max(max(Y*1e3)) min(min(Z*1e3)) max(max(Z*1e3)) -100 100])    
view(90,90); axis square;
set(gca,'FontSize',25)

figure(3)
surf(Y*1e3, Z*1e3, dzBy,'EdgeAlpha',0)
colorbar; colormap('jet')
xlabel('$y\mathrm{[mm]}$','Interpreter','latex'); ylabel('$z\mathrm{[mm]}$','Interpreter','Latex'); title('dzBy');
axis([min(min(Y*1e3)) max(max(Y*1e3)) min(min(Z*1e3)) max(max(Z*1e3)) -100 100])    
view(90,90); axis square;
set(gca,'FontSize',25)

figure(4)
surf(Y*1e3, Z*1e3, dzBz,'EdgeAlpha',0)
colorbar; colormap('jet')
xlabel('$y\mathrm{[mm]}$','Interpreter','latex'); ylabel('$z\mathrm{[mm]}$','Interpreter','Latex'); title('dzBz');
axis([min(min(Y*1e3)) max(max(Y*1e3)) min(min(Z*1e3)) max(max(Z*1e3)) -100 100])    
view(90,90); axis square;
set(gca,'FontSize',25)
