function [x_u, y_u, x_v, y_v, x_p, y_p] = grid_generation(dx1,dy1,dx2,dy2,x_max,x_min,y_min,y_max,obj_start,obj_end)


% Because we know the value of u at the left and right boundary
% We don't need ghost cells


x_u_1 = x_min:dx1:obj_start;
x_u_2 = obj_start+dx2:dx2:obj_end -dx2;
x_u_3 = obj_end:dx1:x_max;

x_u =[x_u_1 x_u_2 x_u_3];


y_u =y_min+dy2/2:dy2:y_max-dy2/2;

x_v_1 = x_min+dx1/2:dx1:obj_start-dx1/2;
x_v_2 = obj_start+dx2/2:dx2:obj_end-dx2/2;
x_v_3 = obj_end+dx1/2:dx1:x_max-dx1/2;

x_v=[x_v_1 x_v_2 x_v_3];


y_v =y_min:dy2:y_max;

x_p_1 = x_min+dx1/2:dx1:obj_start-dx1/2;
x_p_2 = obj_start+dx2/2:dx2:obj_end-dx2/2;
x_p_3 = obj_end+dx1/2:dx1:x_max-dx1/2;

x_p=[x_p_1 x_p_2 x_p_3];

y_p=y_min+dy2/2:dy2:y_max-dy2/2;

end

