%% Mesh Generation

% This script should generate the three types of cell centers


%% Create the domain
x_min=0;
x_max=2;
y_min=0;
y_max=1;


obj_start=0.55;
obj_end =(obj_start+0.9);


dx1 = 1/40;
dy1 = 1/80;
dx2 = 1/80;
dy2 = 1/80;




[x_u, y_u, x_v, y_v, x_p, y_p] = grid_generation(dx1,dy1,dx2,dy2,x_max,x_min,y_min,y_max,obj_start,obj_end);

%% Calculate the size of each of the meshs

Ny_u = length(y_u); 
Nx_u = length(x_u);
Ny_v = length(y_v); 
Nx_v = length(x_v);
Ny_p = length(y_p); 
Nx_p = length(x_p);

%% Create cell types

[u_type, v_type, p_type] = set_cell_type(Ny_u, Nx_u, Ny_v, Nx_v, Ny_p, Nx_p,dx1,dy1,dx2,dy2,x_u, y_u, x_v, y_v, x_p, y_p,x_max,x_min,y_min,y_max,obj_start,obj_end);
