%% Guess initialization


%% Create the guessed values

u_guess = ones(length(y_u),length(x_u)); 
% Set u-velocity to zero in and on the obstruction
u_guess(u_type==-1) = 0;

v_guess = zeros(length(y_v),length(x_v)); 
p_guess = zeros(length(y_p),length(x_p)); 

%% Store the values of the guess for a time history

u_time = u_guess;
v_time = v_guess;
p_time = p_guess;
phi_n = zeros(length(y_p),length(x_p));