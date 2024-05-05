function R = getR(t)
%{
This function finds the capital R (ECI, km) vector which is the vector from the
center of the earth to the deep space network station in Goldstone, CA.

It takes in the time in seconds since the start of the last siderial day 
that you want to find it at.
%}
arguments
    t (1,1) double
end
%To get RF once
% Re = 6378; %Radius of earth in km
% lat = 35.426667;
% long = -116.89;
% z = Re*sind(lat);
% xymag = Re*cosd(lat);
% x = Re*cosd(lat)*cosd(long);
% y = Re*cosd(lat)*sind(long);
% RF = [x,y,z]';%R vector in ECEF cord in kms

%Fixing RF that way it's not running as much trig
RF = [-2.350568838222485e+03;-4.635229164173814e+03;3.697074618177178e+03];
TF2I = getTF2I(t);
R = TF2I*RF;
end