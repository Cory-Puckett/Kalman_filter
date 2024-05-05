function [x2,P2]= rangeUpdate(t, x, P, rho)
%{
This function does the update step using a given range measurement
Inputs:
t: time of measurement
x: state at before update
P: covariance matrix before update
rho: range measurement

Outputs:
x2: state after update
P2: covariance matrix after update
%}
arguments 
    t (1,1) double
    x (6,1) double
    P (6,6) double
    rho (1,1) double
end

sig_rho = 10; %Uncertainty in range measure
R = getR(t); %Vector to station
rho_est = norm(x(1:3)-R);%Estimate of range based on preupdate state
bearing = (x(1:3)-R)/rho_est;
H_rho = [bearing' , 0, 0, 0];
R_kal = sig_rho^2; %The r variable in the kalman gain equation.
%Calculating Kalman gain
K = P*H_rho'*(H_rho*P*H_rho' + R_kal)^-1;
%Calculating updated state
x2 = x+K*(rho - rho_est);
%Calculating updated covariance matrix
P2 = (eye(6)-K*H_rho)*P*(eye(6)-K*H_rho)'+K*R_kal*K';

end