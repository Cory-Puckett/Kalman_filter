function [x2,P2]= rangeRateUpdate(t, x, P, rhoDot)
%{
This function does the update step using a given range rate measurement
Inputs:
t: time of measurement
x: state at before update
P: covariance matrix before update
rhoDot: range rate measurement

Outputs:
x2: state after update
P2: covariance matrix after update
%}
arguments 
    t (1,1) double
    x (6,1) double
    P (6,6) double
    rhoDot (1,1) double
end

rotRate = 2*pi/(86164.0905); %Earth rotation rate rad/s
omega = [0,0,rotRate]';
R = getR(t); %Vector to station
rho = norm(x(1:3)-R);%Estimate of range based on preupdate state
bearing = (x(1:3)-R)/rho;

%Finding derivative of b with respect to r
dbdr = eye(3)/rho - ((x(1:3)-R)*(x(1:3)-R)')/(rho^3);

%Finding derivative of rho_dot with respect to r
dpddr = dbdr'*(x(4:6)-cross(omega,R));

%Defining H matrix
H = [dpddr',bearing'];

rhoDot_est = bearing'*(x(4:6)-cross(omega,R));

%Coppied from rangeUpdate
sig_rhodot = .1; %Uncertainty in range measure
R_kal = sig_rhodot^2; %The r variable in the kalman gain equation.
%Calculating Kalman gain
K = P*H'*(H*P*H' + R_kal)^-1;
%Calculating updated state
x2 = x+K*(rhoDot - rhoDot_est);
%Calculating updated covariance matrix
P2 = (eye(6)-K*H)*P*(eye(6)-K*H)'+K*R_kal*K';

end