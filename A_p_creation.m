function [A_p, b_p, Ap_u, Ap_v] = A_p_creation(u_star, v_star, p_guess, A_u, A_v, ...
    dx1, dy1,dx2,dy2, p_type, A_p)

[Ny_u, Nx_u] = size(u_star);
[Ny_v, Nx_v] = size(v_star);
[Ny_p, Nx_p] = size(p_guess);

 %% Get the velocity at the midpoints of the Pressure cell

uw_star = u_star(:,1:end-1);
ue_star = u_star(:,2:end);

vn_star = v_star(2:end,:);
vs_star = v_star(1:end-1,:);

%% Set the value of aP for the u and v velocities

Ap_u = reshape(diag(A_u),Ny_u,Nx_u);    

Ap_v = reshape(diag(A_v),Ny_v,Nx_v);

%% Calculate the pressure correction matrix coefficients


Cw1 = 1./Ap_u(:,1:end-1)/dx1^2;
Ce1 = 1./Ap_u(:,2:end)/dx1^2;
Cw2 = 1./Ap_u(:,1:end-1)/dx2^2;
Ce2 = 1./Ap_u(:,2:end)/dx2^2;
Cn = 1./Ap_v(2:end,:)/dy1^2;
Cs = 1./Ap_v(1:end-1,:)/dy1^2;

%% Create the b vector

b_p = zeros(length(p_guess(:)),1);

%% Set the pressure correction coefficients

for i = 1:Ny_p*Nx_p      

    % Boundaries on the west and south face
    if p_type(i) == 1
        A_p(i,i) = Ce1(i) + Cn(i);
        A_p(i,i+Ny_p) = -Ce1(i);
        A_p(i,i+1)    = -Cn(i);

    % Boundaries on the west and north face
    elseif p_type(i) == 2
        A_p(i,i) = Ce1(i) + Cs(i);
        A_p(i,i+Ny_p) = -Ce1(i);
        A_p(i,i-1)    = -Cs(i);

    % Boundaries on the west face
    elseif p_type(i) == 3            
        A_p(i,i) = Ce1(i) + Cn(i) + Cs(i);
        A_p(i,i+Ny_p) = -Ce1(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = -Cs(i);

    % Boundaries on the east and south face
    elseif p_type(i) == 4
        A_p(i,i) = Cw1(i) + Cn(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-Ny_p) = -Cw1(i);

    % Boundaries on the east and north face
    elseif p_type(i) == 5
        A_p(i,i) = Cw1(i) + Cs(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = -Cw1(i);

    % Boundaries on the east face
    elseif p_type(i) == 6            
        A_p(i,i) = Cw1(i) + Cn(i) + Cs(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = -Cw1(i);

    % Boundaries on the south face
    elseif p_type(i) == 7
        if i<=0.55/dx1*Ny_p || i>=((0.55/dx1)+(0.9/dx2))*Ny_p
        A_p(i,i) = Cw1(i) + Ce1(i) + Cn(i);
        A_p(i,i+Ny_p) = -Ce1(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-Ny_p) = -Cw1(i);      
        else
        A_p(i,i) = Cw2(i) + Ce2(i) + Cn(i);
        A_p(i,i+Ny_p) = -Ce2(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-Ny_p) = -Cw2(i);
        end


    % Boundaries on the north face
    elseif p_type(i) == 8
        if i<=(0.55/dx1) *Ny_p|| i>=((0.55/dx1)+(0.9/dx2))*Ny_p
            A_p(i,i) = Cw1(i) + Ce1(i) + Cs(i);
            A_p(i,i+Ny_p) = -Ce1(i);
            A_p(i,i-1) = -Cs(i);
            A_p(i,i-Ny_p) = -Cw1(i);
        else
            A_p(i,i) = Cw2(i) + Ce2(i) + Cs(i);
            A_p(i,i+Ny_p) = -Ce2(i);
            A_p(i,i-1) = -Cs(i);
            A_p(i,i-Ny_p) = -Cw2(i);
        end


    elseif p_type(i) == -1
        A_p(i,i) = 1;        
        
    % Interior cells
    else
        if i<=(0.55/dx1)*Ny_p || i>=((0.55/dx1)+(0.9/dx2))*Ny_p
        A_p(i,i) = Cw1(i) + Ce1(i) + Cn(i) + Cs(i);
        A_p(i,i+Ny_p) = -Ce1(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = -Cw1(i);      
        else
        A_p(i,i) = Cw2(i) + Ce2(i) + Cn(i) + Cs(i);
        A_p(i,i+Ny_p) = -Ce2(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = -Cw2(i);
        end


    end            

    if p_type(i) == -1
        b_p(i) = 0;
    elseif i<=(0.55/dx1)*Ny_p || i>=((0.55/dx1)+(0.9/dx2))*Ny_p
        b_p(i) = -(ue_star(i)-uw_star(i))/dx1 -(vn_star(i)-vs_star(i))/dy1;
    else
        b_p(i) = -(ue_star(i)-uw_star(i))/dx2 -(vn_star(i)-vs_star(i))/dy2;
    end

    
    
end

%% Set the correction of the pressure at the middle of the domain to zero

ind = Ny_p*Nx_p-100; 

A_p(ind,:) = 0;
A_p(ind,ind) = 1;
b_p(ind) = 0;


