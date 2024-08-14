function [A_v, b_v] = A_v_creation(u_guess, v_guess, p_guess, dx1, dy1,dx2,dy2, dt,...
            Re, v_type, v_bot, v_top, v_prevTime,v_lef,A_v)
       
[Ny_v, Nx_v] = size(v_guess);

%% Calculation of the velocity at the midpoints (nonlinear terms)

u_ne = [u_guess(:,2:end); zeros(1,Nx_v)];
u_se = [zeros(1,Nx_v); u_guess(:,2:end)];
u_nw = [u_guess(:,1:end-1); zeros(1,Nx_v)];
u_sw = [zeros(1,Nx_v); u_guess(:,1:end-1)];

v_P = v_guess;
v_N = [v_guess(2:end,:); zeros(1,Nx_v)];
v_S = [zeros(1,Nx_v); v_guess(1:end-1,:)];

ue = (u_ne+u_se)/2;
uw = (u_nw+u_sw)/2;
vn = (v_P+v_N)/2;
vs = (v_P+v_S)/2;

%% Calculate of pressure at the north and south midpoints

p_n = [zeros(1,Nx_v); p_guess(2:end,:); zeros(1,Nx_v)];
p_s = [zeros(1,Nx_v); p_guess(1:end-1,:); zeros(1,Nx_v)];

%% Calculation of the coefficients for the matrix

Fe_v = ue/2; Fw_v = uw/2; Fn_v = vn/2/dy1; Fs_v = vs/2/dy1;
De1 = 1/dx1^2/Re;
De2 = 1/dx2^2/Re;
Dw1 = 1/dx1^2/Re;
Dw2 = 1/dx2^2/Re;
Dn = 1/dy1^2/Re; 
Ds = 1/dy1^2/Re; 

%% Create the b vector for the vertical velocity calculation

b_v = zeros(length(v_guess(:)),1);

%% Set the values of the matrix

for i = 1:Ny_v*Nx_v

    %% If the cell is ON the bottom boundary    
    if v_type(i) == 1
        A_v(i,i) = 1;    
        b_v(i) = v_bot;

    %% If the cell is ON the top boundary
    elseif v_type(i) == 2
        A_v(i,i) = 1; 
        b_v(i) = v_top;

    %% If the cell is Neighboring the left boundary
    elseif v_type(i) == 3
        A_v(i,i)      =  Fe_v(i)/dx1 +Fn_v(i)/dx1-Fs_v(i)+De1+2*Dw1+Dn+Ds ;
        A_v(i,i+1)    =  Fn_v(i)-Dn;  
        A_v(i,i-1)    = -Fs_v(i)-Ds;            
        A_v(i,i+Ny_v) =  Fe_v(i)/dx1-De1;
        b_v(i) = -(p_n(i)-p_s(i))/dy1  + 2*(Fw_v(i)/dx1+Dw1)*v_lef;




    %% If the cell is Neighboring the right boundary
    elseif v_type(i) == 4
        A_v(i,i)      =  2*Fe_v(i)/dx1-Fw_v(i)/dx1+Fn_v(i)-Fs_v(i)+Dw1+Dn+Ds ;
        A_v(i,i+1)    =  Fn_v(i)-Dn;  
        A_v(i,i-1)    = -Fs_v(i)-Ds; 
        A_v(i,i-Ny_v) = -Fw_v(i)/dx1-Dw1;            
        b_v(i) = -(p_n(i)-p_s(i))/dy1 ;  

        

    %% If the cell has an East face on a No-slip boundary
    elseif v_type(i) == 5
        if i<=0.55/dx1*Ny_v || i>=(0.55/dx1)+(0.9/dx2)*Ny_v
        A_v(i,i)      =   -Fw_v(i)/dx1+Fn_v(i)-Fs_v(i)+2*De1+Dw1+Dn+Ds ;
        A_v(i,i+1)    =  Fn_v(i)-Dn;  
        A_v(i,i-1)    = -Fs_v(i)-Ds; 
        A_v(i,i-Ny_v) = -Fw_v(i)/dx2-Dw1; 
        b_v(i) = -(p_n(i)-p_s(i))/dy1 ;  
        else
        A_v(i,i)      =   -Fw_v(i)/dx2+Fn_v(i)-Fs_v(i)+2*De2+Dw2+Dn+Ds ;
        A_v(i,i+1)    =  Fn_v(i)-Dn;  
        A_v(i,i-1)    = -Fs_v(i)-Ds; 
        A_v(i,i-Ny_v) = -Fw_v(i)/dx2-Dw2; 
        b_v(i) = -(p_n(i)-p_s(i))/dy2 ; 
        end


    
    %% If the cell is in or on the obstruction
    elseif v_type(i) == -1
        A_v(i,i) = 1;
        b_v(i) = 0;
        
    %% If the cell is In the domain
    else
        if i<=0.55/dx1*Ny_v|| i>=(0.55/dx1)+(0.9/dx2)*Ny_v
        A_v(i,i)      =  Fe_v(i)/dx1-Fw_v(i)/dx1+Fn_v(i)-Fs_v(i)+De1+Dw1+Dn+Ds ;
        A_v(i,i+1)    =  Fn_v(i)-Dn;  
        A_v(i,i-1)    = -Fs_v(i)-Ds; 
        A_v(i,i-Ny_v) = -Fw_v(i)/dx1-Dw1;             
        A_v(i,i+Ny_v) =  Fe_v(i)/dx1-De1;
        b_v(i) = -(p_n(i)-p_s(i))/dy1 ;
        else
        A_v(i,i)      =  Fe_v(i)/dx2-Fw_v(i)/dx2+Fn_v(i)-Fs_v(i)+De2+Dw2+Dn+Ds ;
        A_v(i,i+1)    =  Fn_v(i)-Dn;  
        A_v(i,i-1)    = -Fs_v(i)-Ds; 
        A_v(i,i-Ny_v) = -Fw_v(i)/dx2-Dw2;             
        A_v(i,i+Ny_v) =  Fe_v(i)/dx2-De2;
        b_v(i) = -(p_n(i)-p_s(i))/dy2 ;
        end

    end
end