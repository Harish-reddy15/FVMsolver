%% Time step initialization


%% Set the the u^n, and v^n values    
n=1;
u_prevTime = u_time(:,:,n);
v_prevTime = v_time(:,:,n);
p_prevTime = p_time(:,:,n);

%% Set the initial guess as the previous time step velocity

u_guess = u_prevTime;
v_guess = v_prevTime;
p_guess = p_prevTime;

%% Reset the iteration count and convergence criteria

II = 1;
u_change = 1;

%% Initializae the sparse A matricies

A_u = sparse([],length(u_guess(:)),length(u_guess(:)));
A_v = sparse([],length(v_guess(:)),length(v_guess(:)));
A_p = sparse([],length(p_guess(:)),length(p_guess(:)));
 