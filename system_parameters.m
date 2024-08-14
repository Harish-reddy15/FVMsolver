%% System parameters

%% Set the Reynolds number

Re = 1;

%% Set the relaxation parameters

alpha_u = .2;
alpha_v = .2;
alpha_p = .0025;

%% Set the boundary speeds

u_bot =  0; v_bot = 0;
u_top =  0; v_top = 0;
u_lef =  1; v_lef = 0; 
u_rig =  0; v_rig = 0; 

%% Set the temporal components

dt = 0.02;
t = 0:dt:4;
Nt = length(t);

%% Set the stop conditions

II_max = 500;
u_tol = 10^-6;
II = 1;
u_change = 1;