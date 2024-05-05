function [x2,P2]= bearingUpdate(t, x, P, bearingE)
%{
This function does the update step using a given bearing measurement
Inputs:
t: time of measurement
x: state at before update
P: covariance matrix before update
bearingE: bearing measurement in ENU

Outputs:
x2: state after update
P2: covariance matrix after update
%}
%Pretty sure this is problem
arguments 
    t (1,1) double
    x (6,1) double
    P (6,6) double
    bearingE (3,1) double
end

R = getR(t); %Vector to station
rho = norm(x(1:3)-R);%Estimate of range based on preupdate state

%Finding derivative of b with respect to r
dbdr = eye(3)/rho - ((x(1:3)-R)*(x(1:3)-R)')/(rho^3);

%Finding bearing measurement and estimate in inertial frame
bearing = bE2bI(t,bearingE);
bearing_est = (x(1:3)-R)/rho;

H = [dbdr, zeros(3)];

sig_bear = (1/60)*(pi/180); %Uncertainty in range measure
R_kal = sig_bear^2*(eye(3)-bearing*bearing'); %The R variable in the kalman gain equation.
%Calculating Kalman gain
K = P*H'*(H*P*H' + R_kal)^-1;
%Calculating updated state
x2 = x+K*(bearing - bearing_est);
%Calculating updated covariance matrix
P2 = (eye(6)-K*H)*P*(eye(6)-K*H)'+K*R_kal*K';

end