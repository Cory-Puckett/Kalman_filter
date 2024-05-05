function [x2,P2] = propogate(t1, t2, x, P)
%{
This function does the propogate step
Inputs:
t1: time of given state
t2: time of state you want to find
x: state at t1
P: covariance matrix at t1

Outputs:
x2: state at t2
P2: covariance matrix at t2
%}
arguments
    t1 (1,1) double
    t2 (1,1) double
    x (6,1) double
    P (6,6) double
end

%Givens
mu = 398600;
G = [zeros(3);eye(3)];
sigw = 0.05; %Meritt told me to use this
Q = sigw^2*eye(3);

%Finding new P
F = getF(x(1:3));
Pdot = F*P + P*F' + G*Q*G';
P2 = P+(t2-t1)*Pdot;

%Finding new x
r = x(1:3);
v = x(4:6);
a = (-mu*r)/(norm(r)^3);
xdot = [v;a];
x2 = x+xdot*(t2-t1);

%Finding new x using ODE 45
% [~, xlist] = ode45(@odefun,[t1,t2],x);
% x2 = (xlist(end,:))';
end

%Probably don't need
% function xdot = odefun(~,x)
% %{
% Function used in ODE 45 to find new x
% %}
% mu = 398600;
% xdot = zeros(size(x));
% r = x(1:3);
% v = x(4:6);
% xdot(1:3) = v;
% xdot(4:6) = -mu*r/(norm(r)^3);
% end