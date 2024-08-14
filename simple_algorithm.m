%% Perform simple algorithm to calculate u^n+1, v^n+1, and p^n+1

while II <= II_max && abs(u_change(end)) > u_tol

    u_guess_old = u_guess;
    v_guess_old = v_guess;
    p_guess_old = p_guess;    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of u-star
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [A_u, b_u] = A_u_creation(u_guess, v_guess, p_guess, dx1, dy1,dx2,dy2 ,dt,...
    Re, u_type, u_bot, u_prevTime, u_lef,u_top,A_u,x_u);
    %% Calculate the u star value    
    u_star = reshape(A_u\b_u,Ny_u,Nx_u);

    %% Update the u_star value with relaxation
    u_star = u_guess + alpha_u*(u_star-u_guess);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of v-star
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   [A_v, b_v] = A_v_creation(u_guess, v_guess, p_guess, dx1, dy1,dx2,dy2, dt,...
            Re, v_type, v_bot, v_top, v_prevTime,v_lef,A_v);

    %% Calculate the v-star value

    v_star = reshape(A_v\b_v,Ny_v,Nx_v);

    %% Use relaxation to set the v-star value

    v_star = v_guess + alpha_v*(v_star-v_guess);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of the pressure correction
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A_p, b_p, Ap_u, Ap_v] = A_p_creation(u_star, v_star, p_guess, A_u, A_v, ...
    dx1, dy1,dx2,dy2, p_type, A_p);

%% Calculate the pressure correction and apply it with under-relaxation

    p_correction = reshape(A_p\b_p,Ny_p,Nx_p);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of the velocity correction
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the pressure correction on the east and west faces of u-cells

    pc_E = p_correction(:,2:end);
    pc_W = [p_correction(:,1:end-1)];

%% Calculate the velocity correction for all cells IN the domain

    u_correction = 0*u_guess;
     u_correction(:,2:end-1) = (pc_W-pc_E)/dx1./Ap_u(:,2:end-1);

       u_correction(u_type==-1) = 0;
    
%% Calculate the pressure correction on the north and south faces of v-cells

    pc_N = [p_correction(2:end,:)];
    pc_S = [p_correction(1:end-1,:)];

%% Calculate the velocity correction for all cells IN the domain

    v_correction = 0*v_guess;
    v_correction(2:end-1,:) = (pc_S-pc_N)/dy2./Ap_v(2:end-1,:);

    % Set the correction for cells with a center ON the boundary
    v_correction(v_type==1) = 0;
    v_correction(v_type==2) = 0;
    v_correction(v_type==-1) = 0;
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Corect the velocity directly and pressure with under-relaxation
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u_guess = u_star + u_correction;

    v_guess = v_star + v_correction;

    p_guess = p_guess + alpha_p*p_correction;


%% Impose flux conservation

    m_in = sum(u_guess(:,1)*dy1);
    m_out = sum(u_guess(:,end-1)*dy1);

    u_guess(:,end) = m_in/m_out*u_guess(:,end-1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the continuity residule and changes in u,v, and p
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ue = u_guess(:,2:end);
    uw = u_guess(:,1:end-1);
    vn = v_guess(2:end,1:end);
    vs = v_guess(1:end-1,1:end);

    continuity_residule(II) = sqrt(sum(((ue(:)-uw(:))/dx2+(vn(:)-vs(:))/dy2).^2));

    u_change(II) = sqrt(sum((u_guess(:)-u_guess_old(:)).^2))/Nx_u/Ny_u;
    v_change(II) = sqrt(sum((v_guess(:)-v_guess_old(:)).^2))/Nx_v/Ny_v;
    p_change(II) = sqrt(sum((p_guess(:)-p_guess_old(:)).^2))/Nx_p/Ny_p;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the continuity residule
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     close all;

    if mod(II,5) == 0
        
        figure(4)
        subplot(311)
        imagesc(u_guess); colorbar; caxis([-1 1]*max(abs(u_guess(:))));
        set(gca,'YDir','normal')
        subplot(312)
        imagesc(v_guess); colorbar; caxis([-1 1]*max(abs(v_guess(:))));
        set(gca,'YDir','normal')
        subplot(313)
        imagesc(p_guess); colorbar; caxis([-1 1]*max(abs(p_guess(:))));
        set(gca,'YDir','normal')
        colormap jet
        drawnow

        figure(5)
        plot(u_change,'.-r');
        hold on;
        plot(v_change,'.-g');
        plot(p_change,'.-k');
        plot(continuity_residule,'.-b');
        set(gca,'YScale','log')

    end

    II = II + 1;

end