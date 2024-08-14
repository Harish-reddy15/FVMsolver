function [A_u, b_u] = A_u_creation(u_guess, v_guess, p_guess, dx1, dy1,dx2,dy2 ,dt,...
    Re, u_type, u_bot, u_prevTime, u_lef,u_top,A_u,x_u)

[Ny_u, Nx_u] = size(u_guess);

u_P = u_guess;
u_E = [u_guess(:,2:end) zeros(Ny_u,1)];
u_W = [zeros(Ny_u,1) u_guess(:,1:end-1)];

v_ne = [v_guess(2:end,:) zeros(Ny_u,1)];
v_nw = [zeros(Ny_u,1) v_guess(2:end,:)];
v_se = [v_guess(1:end-1,:) zeros(Ny_u,1)];
v_sw = [zeros(Ny_u,1) v_guess(1:end-1,:)];

%% Calculation of the velocity at the midpoints (nonlinear terms)
    x_e = [0 x_u(:,2:end)-x_u(:,1:end-1)];
    x_w = [x_u(:,2:end)-x_u(:,1:end-1) 0];
    x_ew=x_e+x_w;
    v_nwe=v_nw.*x_e;
    v_new =v_ne.*x_w;
    v_swe=v_sw.*x_e;
    v_sew =v_se.*x_w;

ue = (u_P + u_E)/2;
uw = (u_P + u_W)/2;
vn = (v_nwe+v_new)./x_ew;
vs = (v_swe+v_sew)./x_ew;

%% Calculate the pressures at the east and west midpoint

p_e = [p_guess zeros(Ny_u,1)];
p_w = [zeros(Ny_u,1) p_guess];

%% Coefficients for the A matrix

    for s=1:Nx_u
        if s<=0.55/dx1 || s>=(0.55/dx1)+(0.9/dx2)
            Fe_u(:,s) = ue(:,s)/2/dx1  ; 
        else
            Fe_u(:,s) = ue(:,s)/2/dx2 ; 
        end
    end

    for o=1:Nx_u
        if o<=0.55/dx1 || o>=(0.55/dx1)+(0.9/dx2)
            Fw_u(:,o) = uw(:,o)/2/dx1  ; 
        else
            Fw_u(:,o) = uw(:,o)/2/dx2 ; 
        end
    end
Fn_u = vn/2/dy1; Fs_u = vs/2/dy1;
De1 = 1/dx1^2/Re;
De2 = 1/dx2^2/Re;
Dw1 = 1/dx1^2/Re;
Dw2 = 1/dx2^2/Re;
Dn = 1/dy1^2/Re; 
Ds = 1/dy1^2/Re;        

%% Create the b vector

b_u = zeros(length(u_guess(:)),1);

for i = 1:length(u_type(:))
    
    %% If the cell is ON the right boundary
    if u_type(i) == 1            

        % Outlet boundary condition
        A_u(i,i) = 1;  
        A_u(i,i-Ny_u) = -1;
        b_u(i) = 0;

    %% If the cell is ON the left boundary at the bottom
    elseif u_type(i) == 2
        
        % Inlet boundary condition
        A_u(i,i) = 1;
        b_u(i) = u_lef;           

    %% If the cell is Neighboring the bottom boundary 
    elseif u_type(i) == 3
        if i<=0.55/dx1*Ny_u || i>=(0.55/dx1+0.9/dx2)*Ny_u
        A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-2*Fs_u(i)+De1+Dw1+Dn+2*Ds ;      
        A_u(i,i-Ny_u) = -Fw_u(i)-Dw1;   
        A_u(i,i+Ny_u) =  Fe_u(i)-De1;    
        A_u(i,i+1)    =  Fn_u(i)-Dn;            
        b_u(i)        = -(p_e(i)-p_w(i))/dx1 + 2*u_bot*Ds ;
        else
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-2*Fs_u(i)+De2+Dw2+Dn+2*Ds ;      
        A_u(i,i-Ny_u) = -Fw_u(i)-Dw2;   
        A_u(i,i+Ny_u) =  Fe_u(i)-De2;    
        A_u(i,i+1)    =  Fn_u(i)-Dn;            
        b_u(i)        = -(p_e(i)-p_w(i))/dx2 + 2*u_bot*Ds ;
        end


    %% If the cell is Neighboring the top boundary
    elseif u_type(i) == 4
        if i<=0.55/dx1*Ny_u || i>=(0.55/dx1)+(0.9/dx2)*Ny_u
        A_u(i,i)      =  Fe_u(i)-Fw_u(i)+2*Fn_u(i)-Fs_u(i)+De1+Dw1+Ds+2*Dn ;      
        A_u(i,i-Ny_u) = -Fw_u(i)-Dw1;   
        A_u(i,i+Ny_u) =  Fe_u(i)-De1;                 
        A_u(i,i-1)    = -Fs_u(i)-Ds;      
        b_u(i)        = -(p_e(i)-p_w(i))/dx1+ 2*u_top*Dn;
        else 
        A_u(i,i)      =  Fe_u(i)-Fw_u(i)+2*Fn_u(i)-Fs_u(i)+De2+Dw2+Ds+2*Dn ;      
        A_u(i,i-Ny_u) = -Fw_u(i)-Dw2;   
        A_u(i,i+Ny_u) =  Fe_u(i)-De2;                 
        A_u(i,i-1)    = -Fs_u(i)-Ds;      
        b_u(i)        = -(p_e(i)-p_w(i))/dx2 +2*u_top*Dn;
        end


    %% If the cell is in or on the obstruction
    elseif u_type(i) == -1
        A_u(i,i) = 1;
        b_u(i) = 0;    
        
    %% If the cell is IN the domain
    else

        if i<=0.55/dx1*Ny_u || i>=(0.55/dx1)+(0.9/dx2)*Ny_u
        A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De1+Dw1+Dn+Ds ;        
        A_u(i,i-Ny_u) = -Fw_u(i)-Dw1;   
        A_u(i,i+Ny_u) =  Fe_u(i)-De1;   
        A_u(i,i+1)    =  Fn_u(i)-Dn;                
        A_u(i,i-1)    = -Fs_u(i)-Ds; 
        b_u(i) = -(p_e(i)-p_w(i))/dx1 ;
        else
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De2+Dw2+Dn+Ds ;        
        A_u(i,i-Ny_u) = -Fw_u(i)-Dw2;   
        A_u(i,i+Ny_u) =  Fe_u(i)-De2;   
        A_u(i,i+1)    =  Fn_u(i)-Dn;                
        A_u(i,i-1)    = -Fs_u(i)-Ds; 
        b_u(i) = -(p_e(i)-p_w(i))/dx2 ;
        end


    end
end
