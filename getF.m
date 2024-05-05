function F = getF(rv)
%{
This function finds the capital F matrix which is in the dynamics equation.
It takes in the vector of the position state rv(ECI, km) which is a column vector
%}
arguments
    rv (3,1) double
end
mu = 398600;
r = norm(rv);
A = ((3*mu*(rv*rv'))/r^5)-mu*eye(3)/(r^3);
F = [zeros(3), eye(3);
A, zeros(3)];
end