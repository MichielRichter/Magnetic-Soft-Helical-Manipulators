function [field_3d_mod, B, Bsymbolic, R2Y, R2Z] = start_fit(coords_3d, field_3d, ...
    activA, activB, fitOrder, dispVars)

global a

%% Start fitting
coilMeas = [field_3d coords_3d];
syms x y z d_A d_B;
p = [x;y;z];

Rs = NaN(1,dispVars);
coeffCells = cell(1,dispVars);
dispVecs = cell(1,dispVars);
bestOrder = NaN(1,dispVars);

for i = 1:dispVars
    fVal_min = 1e10;
    % Make displacement vector
    if i == 1
        if activA && activB
            dispVec0 = [0e-3, 0e-3]; %dA1 dB1 dA2 dB2 ...
            UB = [10, 10]; %mm,mm,deg,deg,deg,mm
            LB = -UB;
        else
            dispVec0 = [0e-3];
            UB = [10]; %mm,mm,deg,deg,deg,mm
            LB = -UB;
        end
    else
        if activA && activB
            dispVec0 = [dispVec (dispVec(end-1)+5e-3) (dispVec(end)+5e-3)];
            UB = [UB UB(1)*ones(1,2)];
            LB = [LB -UB(end-1:end)*0];
        else
            dispVec0 = [dispVec (dispVec(end)+0.5e-3)];
            UB = [UB UB(1)*ones(1,1)];
            LB = [LB -UB(end)*0];
        end
    end
    
    for order = 1:fitOrder
        % Make symbolic matlab function for A and B coefficients
        BsymA = sym('BA',[3, order]);
        BsymB = sym('BB',[3, order]);
        psiAB = 1/((x)^2 + (y)^2 + (z+d_A)^2)^(1/2);
        psiB = 1/((x)^2 + (y)^2 + (z+d_B)^2)^(1/2);
        for j = 1:order
            psiB = diff(psiB,z);
            psiAB = diff(psiAB,z);
            psiA = (((x)^2 + (y)^2 + (z+d_A)^2)^(1/2))^(2*j+1) * psiAB;
            
            BsymA(:,j) = gradient(psiA,p);
            BsymB(:,j) = gradient(psiB,p);
        end
        
        % Make combined symbolic matlab function
        if activA && activB
            Bsym = [BsymA BsymB];
            vars = {d_A d_B, x, y, z};
        elseif activA
            Bsym = BsymA;
            vars = {d_A, x, y, z};
        elseif activB
            Bsym = BsymB;
            vars = {d_B, x, y, z};
        end
        B = matlabFunction(Bsym,'Vars',vars);
        
        %{
        Obtained so far:
        BsymA = Symbolic field matrix 3xN for scalar potentials with A
        BsymB = Symbolic field matrix 3xN for scalar potentials with B
        Bsym  = Symbolic field matrix 3xN(activA + activB)
        B     = matlab function of Bsym
    
        Todo:
        find optimal coefficients and optimal displacement variables
        %}
        
        % Fitting
        options = optimset('Largescale','off','Display','iter');
        [dispVec,fVal] = fmincon(@(disp)objective(disp, B, order, ...
            coilMeas, activA, activB), dispVec0, [], [],[], [], LB, UB, ...
            [], options);
        if fVal < fVal_min
            fVal_min = fVal;
            Rs(i) = fVal;
            dispVecs{i} = dispVec;
            coeffCells{i} = a;
            bestOrder(i) = order;
        end
    end
    dispVec = dispVecs{i};
    a = coeffCells{i};
end

%% Choose best fit
idx_bestFit = find(min(Rs) == Rs);
dispVec = cell2mat(dispVecs(idx_bestFit));
coeffs = cell2mat(coeffCells(idx_bestFit));
order = bestOrder(idx_bestFit);

% Make symbolic matlab function for A and B coefficients
if activA && activB
    dA = dispVec(1:2:end);
    dB = dispVec(2:2:end);
    
    B = [];
    for i = 1:size(dA,2)
        BsymA = sym('BA',[3, order]);
        BsymB = sym('BB',[3, order]);
        d_A = dA(i);
        d_B = dB(i);
        psiAB = 1/((x)^2 + (y)^2 + (z+d_A)^2)^(1/2);
        psiB = 1/((x)^2 + (y)^2 + (z+d_B)^2)^(1/2);
        for j = 1:order
            psiB = diff(psiB,z);
            psiAB = diff(psiAB,z);
            psiA = ((((x)^2 + (y)^2 + (z+d_A)^2))^(1/2))^(2*j+1) * psiAB;
            
            BsymA(:,j) = gradient(psiA,p);
            BsymB(:,j) = gradient(psiB,p);
        end
        B = [B BsymA BsymB];
    end
    B = B*coeffs;
    
elseif activA
    dA = dispVec;
    dB = [];
    B = [];
    for i = 1:size(dA,2)
        BsymA = sym('BA',[3, order]);
        d_A = dA(i);
        psiAB = 1/((x)^2 + (y)^2 + (z+d_A)^2)^(1/2);
        for j = 1:order
            psiAB = diff(psiAB,z);
            psiA = ((((x)^2 + (y)^2 + (z+d_A)^2))^(1/2))^(2*j+1) * psiAB;
            
            BsymA(:,j) = gradient(psiA,p);
        end
        B = [B BsymA];
    end
    B = B*coeffs;
elseif activB
    dA = [];
    dB = dispVec;
    B = [];
    for i = 1:size(dB,2)
        BsymB = sym('BB',[3, order]);
        d_B = dB(i);
        psiB = 1/((x)^2 + (y)^2 + (z+d_B)^2)^(1/2);
        for j = 1:order
            psiB = diff(psiB,z);
            BsymB(:,j) = gradient(psiB,p);
        end
        B = [B BsymB];
    end
    B = B*coeffs;
end

% Make matlab function
Bsymbolic = B;
B = matlabFunction(B);

%% Get field at 3d coordinates
field_3d_mod = NaN(size(field_3d));

for i = 1:size(coords_3d,1)
    px = coords_3d(i,1);
    py = coords_3d(i,2);
    pz = coords_3d(i,3);
    
    field_mod = B(px,py,pz);
    field_3d_mod(i,:) = field_mod';
end

%% R squared value
y = field_3d;
f = field_3d_mod;

SSresY = sum((y(:,2)-f(:,2)).^2);
SSresZ = sum((y(:,3)-f(:,3)).^2);

SStotY = sum((y(:,2)-mean(y(:,2))).^2);
SStotZ = sum((y(:,3)-mean(y(:,3))).^2);

R2Y = 1-(SSresY/SStotY);
R2Z = 1-(SSresZ/SStotZ);


end

