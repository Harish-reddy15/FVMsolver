clear; close all; clc;

%% Create the domain

mesh_generation


 

%% Set system and solver parameters

system_parameters

%% Generate the initial guess for the system

guess_initialization

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%% TIME MARCHING %%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    %% Set the variables to begin a new time step
    
    time_step_initialization
    
    %% Begin recursive calculation of u^n+1 and v^n+1

    simple_algorithm
    
    %% Perform a heat transfer step
    
    heat_transfer_step
    
    %% Plot results and store the u^n+1, v^n+1, and p^n+1 fields
    
    plot_store_results

    %%% mass fluxes 

    m_start_obj=sum((u_guess(:,0.55/dx1)*dy1));
    m_start_2step=sum((u_guess(:,0.55/dx1+0.3/dx2)*dy1));
    m_start_3step=sum((u_guess(:,0.55/dx1+0.6/dx2)*dy1));
    m_end_obj=sum((u_guess(:,0.55/dx1+0.9/dx2)*dy1));

    for i=1:length(u_guess)
    avg_vel(i)=mean(nonzeros(u_guess(:,i)));
    end




    
