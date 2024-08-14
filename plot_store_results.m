

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the results
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col = jet(201);
if mod(n,2) == 0 
    
    figure(2);
    hold off; 
    u_temp = [u_bot+zeros(1,Nx_u); u_guess; u_top+zeros(1,Nx_u)];
    v_temp = [v_lef+zeros(Ny_v,1)  v_guess  v_rig+zeros(Ny_v,1)];

    x_vt = [min(x_u(:)) x_v max(x_u(:))]; [X_v, Y_v] = meshgrid(x_vt,y_v);
    y_ut = [min(y_v(:)) y_u max(y_v(:))]; [X_u, Y_u] = meshgrid(x_u,y_ut);

    [X,Y] = meshgrid(linspace(min(x_u(:)),max(x_u(:)),100),linspace(min(y_v(:)),max(y_v(:)),40));

    U = griddata(X_u(:),Y_u(:),u_temp(:),X,Y);
    V = griddata(X_v(:),Y_v(:),v_temp(:),X,Y);

    imagesc(x_p,y_p,reshape(p_guess,Ny_p,Nx_p)); hold on;

    colormap jet
    quiver(X,Y,U,V,'k')
%     plot(u_guess,y_u,'.-','MarkerSize',30,'LineWidth',4,'Color',col(n,:))
    axis([min(x_u(:)) max(x_u(:)) min(y_v(:)) max(y_v(:))])
    set(gca,'YDir','normal')
    colorbar
    drawnow

end
    
    
%% Store the results
 
u_time(:,:,n+1) = u_guess;
v_time(:,:,n+1) = v_guess;
p_time(:,:,n+1) = p_guess;
