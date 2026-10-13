function [d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y] = delta_X(x_f, x_l)


a_f     = x_f(1);
e_f     = x_f(2);
i_f     = x_f(3);
omega_f = x_f(4);
Omega_f = x_f(5);
M_f     = x_f(6);

a_l     = x_l(1);
e_l     = x_l(2);
i_l     = x_l(3);
omega_l = x_l(4);
Omega_l = x_l(5);
M_l     = x_l(6);



u_l    = omega_l + M_l;
u_f    = omega_f + M_f;
u_diff =  ssa(M_f -M_l) + ssa(omega_f - omega_l); 

d_a          = (a_f - a_l) / a_l;
d_lambda     = u_diff + (Omega_f - Omega_l)*cos(i_l);
d_e_x        = e_f*cos(omega_f) - e_l*cos(omega_l);
d_e_y        = e_f*sin(omega_f) - e_l*sin(omega_l);
d_i_x        = i_f - i_l;
d_i_y        = (Omega_f -Omega_l)*sin(i_l);


end

function angle = ssa(angle)
    angle = mod( angle + pi, 2 * pi ) - pi; 
end