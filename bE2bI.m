function bI = bE2bI(t,bE)
%{
This function converts a bearing measurement in the ENU frame to the ECI
frame.
It takes in t which is time in seconds since the start of the last siderial
day. b is then the bearing measurement taken from Golstone, CA
Frames abbreviations ECI = I, ECEF = F, ENU = E
%}
arguments
    t (1,1) double
    bE (3,1) double
end
%Check that bearing is unit vector
if abs(1-norm(bE))>=10^-6
    warning("Bearing measurement isn't unit vector. It is off by " + abs(1-norm(bE)));
end

%Position of observatory
Re = 6378; %Radius of earth in km
lat = 35.426667;
long = -116.89;
z = Re*sind(lat);
x = Re*cosd(lat)*cosd(long);
y = Re*cosd(lat)*sind(long);
RF = [x,y,z]';%R vector in ECEF cord in kms

%Constructing E, N, and U vectors
u = RF/norm(RF);
e = cross(-u,[0;0;1]);
e = e/norm(e);
n = cross(u,e);
n = n/norm(n);
%Transition matrices
TF2E = [e';n';u'];
TE2F = TF2E';
TF2I = getTF2I(t);
bI = TF2I*TE2F*bE;
end