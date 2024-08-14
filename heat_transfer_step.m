%% Perform a step of the heat transfer analysis

A_phi = zeros(length(phi_n(:)));
b_phi = zeros(length(phi_n(:)),1);

%% Set the heat transfer properties

phi_in = 0;
phi_wall = 1;
kappa = 1;

 %% Get the velocity at the midpoints of the Pressure/Temperature cell

uw = u_star(:,1:end-1);
ue = u_star(:,2:end);

vn = v_star(2:end,:);
vs = v_star(1:end-1,:);

%% Coefficients for the A matrix

Fe1 = ue/2/dx1;
Fe2 = ue/2/dx2;
Fw2 = ue/2/dx2;
Fw1 = uw/2/dx1;   Fn = vn/2/dy1;   Fs = vs/2/dy1;
De1 = kappa/dx1^2; Dw1 = kappa/dx1^2;
De2 = kappa/dx2^2; Dw2 = kappa/dx2^2;

Dn = kappa/dy1^2; Ds = kappa/dy1^2;

%% Create the A matrix for the phi system

for i = 1:Ny_p*Nx_p      

    % Boundaries on the west and south face
    if p_type(i) == 1
        A_phi(i,i) = Fe1(i)+Fn(i)-2*Fs(i) + De1+2*Dw1+Dn +2*Ds;
        A_phi(i,i+Ny_p) = Fe1(i)-De1;
        A_phi(i,i+1)    = Fn(i)-Dn;
        b_phi(i) =   2*(Fw1(i)+Dw1)*phi_in+ 2*(Fs(i)+Ds)*(phi_wall);

    % Boundaries on the west and north face
    elseif p_type(i) == 2
        A_phi(i,i) = Fe1(i)-Fs(i) + De1+2*Dw1+2*Dn+Ds ;
        A_phi(i,i+Ny_p) = Fe1(i)-De1;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        b_phi(i) =  2*(Fw1(i)+Dw1)*phi_in + 2*(-Fn(i)+Dn)*phi_wall;

    % Boundaries on the west face
    elseif p_type(i) == 3        
        A_phi(i,i) = Fe1(i)+Fn(i)-Fs(i) + De1+2*Dw1+Dn+Ds ;
        A_phi(i,i+Ny_p) = Fe1(i)-De1;
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        b_phi(i) =  2*(Fw1(i)+Dw1)*phi_in;

    % Boundaries on the east and south face
    elseif p_type(i) == 4
        A_phi(i,i) = 2*Fe1(i)-Fw1(i)+Fn(i)-2*Fs(i) +Dw1+Dn +2*Ds;        
        A_phi(i,i+1)    = Fn(i)-Dn;        
        A_phi(i,i-Ny_p) =-Fw1(i)-Dw1;
        b_phi(i) = 2*(Fs(i)+Ds)*phi_wall;

    % Boundaries on the east and north face
    elseif p_type(i) == 5
        A_phi(i,i) = 2*Fe1(i)-Fw1(i)-Fs(i) + Dw1+Ds+2*Dn ;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-Ny_p) =-Fw1(i)-Dw1;
        b_phi(i) =  2*(-Fn(i)+Dn)*phi_wall;

    % Boundaries on the east face
    elseif p_type(i) == 6         
        A_phi(i,i) = 2*Fe1(i)-Fw1(i)+Fn(i)-Fs(i) + +Dw1+Dn+Ds ;        
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-Ny_p) =-Fw1(i)-Dw1;
        b_phi(i) = 0;

    % Boundaries on the south face
    elseif p_type(i) == 7
        
        if i<=(0.55/dx1)*Ny_u || i>=(0.55/dx1+0.9/dx2)*Ny_u
        
            A_phi(i,i) = Fe1(i)-Fw1(i)+Fn(i)-2*Fs(i) + De1+Dw1+Dn+2*Ds ;
        A_phi(i,i+Ny_p) = Fe1(i)-De1;
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-Ny_p) =-Fw1(i)-Dw1;
        b_phi(i) = 0+ 2*(Fs(i)+Ds)*phi_wall; 

        else
         A_phi(i,i) = Fe2(i)-Fw2(i)+Fn(i)-2*Fs(i) + De2+Dw2+Dn+2*Ds ;
        A_phi(i,i+Ny_p) = Fe2(i)-De2;
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-Ny_p) =-Fw2(i)-Dw2;
        b_phi(i) = 0+2*(Fs(i)+Ds)*phi_wall;
        end


    % Boundaries on the north face
    elseif p_type(i) == 8
        if i<=(0.55/dx1)*Ny_u || i>=((0.55/dx1)+(0.9/dx2))*Ny_u
        A_phi(i,i) = Fe1(i)-Fw1(i)   -Fs(i) + De1+Dw1+2*Dn+Ds ;
        A_phi(i,i+Ny_p) = Fe1(i)-De1;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-Ny_p) =-Fw1(i)-Dw1;
        b_phi(i) = 0+ 2*(-Fn(i)+Dn)*phi_wall;
        else
            A_phi(i,i) = Fe2(i)-Fw2(i)   -Fs(i) + De2+Dw2+2*Dn+Ds ;
        A_phi(i,i+Ny_p) = Fe2(i)-De2;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-Ny_p) =-Fw2(i)-Dw2;
        b_phi(i) = 0+ 2*(-Fn(i)+Dn)*phi_wall;
        end


    elseif p_type(i) == -1
        A_phi(i,i) = 1;  
        b_phi(i) = 0;      
        
    % Interior cells
    else
        if i<=(0.55/dx1)*Ny_u || i>=((0.55/dx1)+(0.9/dx2))*Ny_u
        A_phi(i,i) = Fe1(i)-Fw1(i)+Fn(i)-Fs(i) + De1+Dw1+Dn+Ds ;
        A_phi(i,i+Ny_p) = Fe1(i)-De1;
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-Ny_p) =-Fw1(i)-Dw1;
        b_phi(i) = 0;
        else
        A_phi(i,i) = Fe2(i)-Fw2(i)+Fn(i)-Fs(i) + De2+Dw2+Dn+Ds ;
        A_phi(i,i+Ny_p) = Fe2(i)-De2;
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-Ny_p) =-Fw2(i)-Dw2;
        b_phi(i) = 0;
        end
    end            

    
end

phi_np1 = A_phi\b_phi;

figure(6)
imagesc(x_p,y_p,reshape(phi_np1,length(y_p),length(x_p)));
set(gca,'YDir','normal')
colormap jet
