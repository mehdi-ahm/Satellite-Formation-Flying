function y = nmpc_orbit_output(x,u,p)
% NMPC_ORBIT_OUTPUT
%
% x =
% [delta_a
%  delta_lambda
%  delta_ex
%  delta_ey
%  delta_ix
%  delta_iy]
%
% Outputs controlled by NMPC:
% y = [delta_a; delta_lambda]

    %#ok<INUSD>

    y = [x(1);
         x(2)];

end