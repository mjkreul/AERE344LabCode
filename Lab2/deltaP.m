function ans = deltaP(rho, h_i, h_i0, h_ref, h_ref0)


g = 9.8; %m/s^2

% -rho*g*sin(theta)*((h_i - h_i0) - (h_ref - h_ref0)); but theta is 90 so
% it's equal to 1.  Maybe we need to change it in later labs?
ans = -rho*g*((h_i - h_i0) - (h_ref - h_ref0));
end