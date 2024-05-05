function T = getTF2I(t)
%{
This function finds the state transition matrix from ECEF to ECI

It takes in the time in seconds since the start of the last siderial day 
that you want to find it at.
%}
arguments
    t (1,1) double
end
rotRate = 2*pi/(86164.0905); %Earth rotation rate rad/s
ang = rotRate*t;
T = [cos(ang),sin(ang),0;
    -sin(ang),cos(ang),0;
    0,0,1]'; %Attitude transformation matrix from fixed to inertial frames
end