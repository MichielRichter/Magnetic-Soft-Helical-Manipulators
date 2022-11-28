function [Rs] = objective(disp, B, fitOrder, coilMeas, activA, activB)

global a
%% Get source displacements, fields, positions
if activA && activB
    dA = disp(1:2:end);
    dB = disp(2:2:end);
elseif activA
    dA = disp;
elseif activB
    dB = disp;
end

Bx = coilMeas(:,1);
By = coilMeas(:,2);
Bz = coilMeas(:,3);
Px = coilMeas(:,4);
Py = coilMeas(:,5);
Pz = coilMeas(:,6);


%% Fit
A = NaN(3*size(coilMeas,1),fitOrder*size(disp,2));
b = NaN(size(A,1),1);

for i = 1:size(Px,1)
    % Position
    px = Px(i);
    py = Py(i);
    pz = Pz(i);
    
    if activA && activB
        for j = 1:size(dA,2)
            d_A = dA(j);
            d_B = dB(j);
            A((i-1)*3+1:i*3, 2*(j-1)*fitOrder+1:2*j*fitOrder) = B(d_A, ...
                d_B, px, py, pz);
        end
        
    elseif activA
        for j = 1:size(dA,2)
            d_A = dA(j);
            A((i-1)*3+1:i*3, (j-1)*fitOrder+1:j*fitOrder) = B(d_A, ...
                px, py, pz);
        end
        
    elseif activB
        for j = 1:size(dB,2)
            d_B = dB(j);
            A((i-1)*3+1:i*3, (j-1)*fitOrder+1:j*fitOrder) = B(d_B, ...
                px, py, pz);
        end
    end

    
    b((i-1)*3+1:i*3) = [Bx(i); By(i); Bz(i)];

end

%%
a = pinv(A)*b; %coefficients
Bmod = A*a; %Model field

Bxmod = Bmod(1:3:end);
Bymod = Bmod(2:3:end);
Bzmod = Bmod(3:3:end);

%% Calculate R squared values
y = [Bx By Bz];
f = [Bxmod Bymod Bzmod];


SSresX = sum((y(:,1)-f(:,1)).^2);
SSresZ = sum((y(:,3)-f(:,3)).^2);

SStotX = sum((y(:,1)-mean(y(:,1))).^2);
SStotZ = sum((y(:,3)-mean(y(:,3))).^2);

R2X = 1-(SSresX/SStotX);
R2Z = 1-(SSresZ/SStotZ);


%Ox = sum( (y(:,1) - f(:,1)).^2./((y(:,1)-mean(y(:,1))).^2) + ((y(:,1) - f(:,1))./y(:,1)).^2 );
%Oz = sum( (y(:,3) - f(:,3)).^2./((y(:,3)-mean(y(:,3))).^2) + ((y(:,3) - f(:,1))./y(:,3)).^2 );
%%
Rs = -R2X + -R2Z;
%Rs = ( mean(abs(y(:,1)-f(:,1))) + mean(abs(y(:,3)-f(:,3))) );
%Rs = -Ox + -Oz;

end

