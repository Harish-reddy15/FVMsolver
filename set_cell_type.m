function [u_type, v_type, p_type] = set_cell_type(Ny_u, Nx_u, Ny_v, Nx_v, Ny_p, Nx_p,dx1,dy1,dx2,dy2,x_u, y_u, x_v, y_v, x_p, y_p,x_max,x_min,y_min,y_max,obj_start,obj_end);


step_1_front=obj_start;
step_1_top=0.1;
step_2_top=0.2;
step_2_front=0.3;
step_3_top=0.3;
step_3_front=0.6;
step_end_front=0.9;

t=round(step_1_front/dx1);
t1=round(step_end_front/dx2)+t;
t2=round(step_1_top/dy2);
t3=round(step_2_top/dy2);
t4=round(step_3_top/dy2);
t5=(step_2_front/dx2)+t;
t6=(step_3_front/dx2)+t;
t6=round(t6);

t7=round((y_max-step_1_top)/dy2);
t8=round((y_max-step_2_top)/dy2);
t9=round((y_max-step_3_top)/dy2);
t0=(step_2_front/dy2)+t;
t11=(step_3_front/dy2)+t;
t11=round(t11);
%% Create cell type variables
u_type =zeros(Ny_u,Nx_u);
v_type =zeros(Ny_v,Nx_v);
p_type =zeros(Ny_p,Nx_p);
%% Set the u_type

u_type(:,1)=2;% Cells on the left get the in let

u_type(:,end)=1;% Cells on the right get the outlet

u_type(1,2:t) =3;
u_type(1,t1:end-1)=3;% Cells along the bottom of the domain that are not 
% on the left/right wall get a no slip to the south

u_type(end,2:t) =4;
u_type(end,t1:end-1) =4;% Cells along the top of the domain that are not 
% on the left/right wall get a no slip to the
% north

u_type(1:t2,t+1)=-1;%%%first vertical wall on obstraction on the bottom
u_type(t2+1,t+1:t5+2)=3;%%%first horizontal wall on obstraction on the bottom

u_type(t2+1:t3,t5+1)=-1; %%%second vertical wall on obstraction on the bottom
u_type(t3+1,t5+1:t1+1)=3;%%%second horizontal wall on obstraction on the bottom

u_type(t3+1:t4,t6+1)=-1;%%%third vertical wall on obstraction on the bottom
u_type(t4+1,t6+1:t1+1)=3;%%%third vertical wall on obstraction on the bottom

u_type(1:t4,t1+1)=-1; %%% final vertical wall on the bottom


u_type(t7+1:end,t+1)=-1;%%%first vertical wall on obstraction on the top
u_type(t7,t+1:t0+1)=4;%%%first horizontal wall on obstraction on the top

u_type(t8+1:t7,t0+1)=-1;%%%second vertical wall on obstraction on the top
u_type(t8,t0+1:t11+1)=4;%%%second horizontal wall on obstraction on the top

u_type(t9+1:t8,t11+1)=-1;%%%third vertical wall on obstraction on the top
u_type(t9,t11+1:t1+1)=4;%%%third vertical wall on obstraction on the top

u_type(t9+1:end,t1+1)=-1;  %%% final vertical wall on the top





u_type(1:t2,t+2:t1)=-1;%%%%bottom
u_type(t2:t3,t5+2:t1)=-1;%%%%bottom
u_type(t3:t4,t6+2:t1)=-1;%%%%bottom

u_type(t7+1:end,t+2:t1)=-1;%%%top
u_type(t8+1:t7+1,t0+2:t1)=-1;%%%top
u_type(t9+1:t8+1,t11+2:t1)=-1;%%%top



[X_u, Y_u] = meshgrid(x_u, y_u);
figure(1)
plot([0 2 2 0 0],[0 0 1 1 0],'k','LineWidth',2);
hold on
plot([0.55 0.55 0.85 0.85 1.15 1.15 1.45 1.45 ],[0 0.1 0.1 0.2 0.2 0.3 0.3 0],'k','LineWidth',2);
s = scatter(X_u(:),Y_u(:),1000,u_type(:),'Marker','.');
colormap jet;

%% Set the v_type
                         
v_type(1,1:t)=1;
v_type(1,t1:end)=1;% Cells on the bottom get the No slip

v_type(end,1:t)=2;
v_type(end,t1:end)=2;% Cells on the top get the Dirichlet equation v=v_top

v_type(2:end-1,1)=3;% Cells along the left of the domain that are not 
% on the top/bottom wall get a inlet
                         
v_type(2:end-1,end)=4;% Cells along the right of the domain that are not 
% on the top/bottom wall get a outlet
% north

v_type(1:t2+1,t)=5;%%%first vertical wall on obstraction on the bottom
v_type(t2+1,t+1:t5+1)=1;%%%first horizontal wall on obstraction on the top

v_type(t2+2:t3+1,t5)=5; %%%second vertical wall on obstraction on the bottom
v_type(t3+1,t5+1:t1)=1;%%%second horizontal wall on obstraction on the bottom

v_type(t3+2:t4+1,t6)=5;%%%third vertical wall on obstraction on the bottom
v_type(t4+1,t6+1:t1)=1;%%%third vertical wall on obstraction on the bottom

v_type(1:t4+1,t1+1)=5; %%% final vertical wall on the bottom


v_type(t7+1:end,t)=5;%%%first vertical wall on obstraction on the top
v_type(t7+1,t+1:t0+1)=2;%%%first horizontal wall on obstraction on the top

v_type(t8+1:t7,t0)=5;%%%second vertical wall on obstraction on the top
v_type(t8+1,t0+1:t11+1)=2;%%%second horizontal wall on obstraction on the top

v_type(t9+1:t8,t11)=5;%%%third vertical wall on obstraction on the top
v_type(t9+1,t11+1:t1)=2;%%%third vertical wall on obstraction on the top

v_type(t9+1:end,t1+1)=5;  %%% final vertical wall on the top


v_type(1:t2,t+1:t1)=-1;%%%%bottom
v_type(t2:t3,t5+1:t1)=-1;%%%%bottom
v_type(t3:t4,t6+1:t1)=-1;%%%%bottom

v_type(t7+2:end,t+1:t1)=-1;%%%top
v_type(t8+2:t7+1,t0+1:t1)=-1;%%%top
v_type(t9+2:t8+1,t11+1:t1)=-1;%%%top


[X_v, Y_v] = meshgrid(x_v, y_v);
figure(2)
plot([0 2 2 0 0],[0 0 1 1 0],'k','LineWidth',2);
hold on
plot([0.55 0.55 0.85 0.85 1.15 1.15 1.45 1.45 ],[0 0.1 0.1 0.2 0.2 0.3 0.3 0],'k','LineWidth',2);
s = scatter(X_v(:),Y_v(:),1000,v_type(:),'Marker','.');
colormap jet;


%% Set the p_type                         
                         
p_type (1,1:t)=7;
p_type (1,t1:end)=7;% Cells at the bottom boundary but not in the corners so only v_s' = 0;
p_type (end,1:t)=8;
p_type (end,t1:end)=8;% Cells at the top boundary but not in the corners so only v_n' = 0;
p_type (2:end-1,1)=3;% Cells at the left boundary but not in the corners so only u_w' = 0;
p_type (2:end-1,end)=6;% Cells at the right boundary but not in the corners so only u_e' = 0;

p_type (1,1)=1;% Cell in the bottom left corner has a u_w and v_s that do not need correcting
p_type (1,end)=4;% Cell in the bottom right corner has a u_e and v_s that do not need correcting
p_type(1,t)=4;
p_type(t2+1,t5)=4;
p_type(t3+1,t6)=4;
p_type (end,1)=2;% Cell in the top left corner has a u_w and v_n that do not need correcting
p_type (end,end)=5;% Cell in the top right corner has a u_e and v_n that do not need correcting


p_type(2:t2,t)=6;%%%first vertical wall on obstraction on the bottom
p_type(t2+1,t+1:t5-1)=7;%%%first horizontal wall on obstraction on the top

p_type(t2+2:t3,t5)=6; %%%second vertical wall on obstraction on the bottom
p_type(t3+1,t5+1:t11)=7;%%second horizontal wall on obstraction on the bottom

p_type(t3+2:t4,t6)=6;%%%third vertical wall on obstraction on the bottom
p_type(t4+1,t6+1:t1)=7;%%%third vertical wall on obstraction on the bottom

p_type(1:t4,t1+1)=3; %%% final vertical wall on the bottom


p_type(t7+1:end,t)=6;%%%first vertical wall on obstraction on the top
p_type(t7,t+1:t0+2)=8;%%%first horizontal wall on obstraction on the top

p_type(t8+1:t7,t0)=6;%%%second vertical wall on obstraction on the top
p_type(t8,t0+1:t11+2)=8;%%%second horizontal wall on obstraction on the top

p_type(t9+1:t8,t11)=6;%%%third vertical wall on obstraction on the top
p_type(t9,t11+1:t1)=8;%%%third vertical wall on obstraction on the top

p_type(t9+1:end,t1+1)=3;  %%% final vertical wall on the top


p_type(1:t2,t+1:t1)=-1;%%%%bottom
p_type(t2:t3,t5+1:t1)=-1;%%%%bottom
p_type(t3:t4,t6+1:t1)=-1;%%%%bottom

p_type(t7+1:end,t+1:t1)=-1;%%%top
p_type(t8+1:t7+1,t0+1:t1)=-1;%%%top
p_type(t9+1:t8+1,t11+1:t1)=-1;%%%top

p_type(t3+1,t6)=4;

p_type(end,t)=5;
p_type(t7,t0)=5;
p_type(t8,t11)=5;
p_type(1,t1+1)=1;
p_type(end,t1+1)=2;


% Plot the p-type results
[X_p, Y_p] = meshgrid(x_p, y_p);
figure(3)
plot([0 2 2 0 0],[0 0 1 1 0],'k','LineWidth',2);
hold on
plot([0.55 0.55 0.85 0.85 1.15 1.15 1.45 1.45 ],[0 0.1 0.1 0.2 0.2 0.3 0.3 0],'k','LineWidth',2);

s = scatter(X_p(:),Y_p(:),1000,p_type(:),'Marker','.');

colormap jet;
return


