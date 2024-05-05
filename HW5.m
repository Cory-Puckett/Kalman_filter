%{
Homework 5
This homework assignment creates a Extended Calman filter to solve for the
state of the Artemis Spacecraft based on range, range rate, and bearing
measurements.
%}
clear all
close all

t0 = 54000;
s0 = [500,6500,3500,-10,2,3]';
P0 = ([100,100,100,0.1,0.1,0.1].*eye(6)).^2;
measures = readtable("MANE6964_HW6_meas.xlsx");
traj = readtable("MANE6964_HW6_traj.xlsx");
tf = max(measures.t);
dt = 60;%Time between each propogate step
states = zeros(361,13); %Matrix to store all states in same format as trajectory file
%States format: t,x,y,z,vx,vy,vz,varx,vary,varz,varvx,varvy,varvz


%% Looping through all measurements
%{
This section will create all of the states using the propogate and update
steps.
%}
%Update for first entry only
[x0,P0plus] = rangeUpdate(t0,s0,P0,measures.ran(1));
states(1,:) = [t0, x0',diag(P0plus)'];

%P and X to be updated and propogated
lastP = P0plus;
lastx = x0;
% lastP = P0;
% lastx = s0;

warncount = 0;

%Update for every other entry
for j = 1:(length(states)-1)
    t = t0+dt*j;

    %Propogate step
    [lastx,lastP] = propogate(t-dt,t,lastx,lastP);

    %If there is a measurement then update
    if rem((t0-t)/300,1)==0
    % if false
        i = find(measures.t==t,1);
        %First not empty measure value
        measure_num = find(table2array(measures(i,2:end)),1);
        if isempty(measure_num)
            error("No measure at time " + t + " seconds");
        elseif measure_num == 1
            %Update for range
            rho = measures.ran(i);
            [lastx,lastP] = rangeUpdate(t,lastx,lastP,rho);
        elseif measure_num == 2
            %Update for range rate
            rho_dot = measures.randot(i);
            [lastx,lastP] = rangeRateUpdate(t,lastx,lastP,rho_dot);
        elseif measure_num == 3
            %Update for bearing
            bearingE = table2array(measures(i,4:6))';
            [lastx,lastP] = bearingUpdate(t,lastx,lastP,bearingE);
        end
    end
    %Store new values
    states(j+1,:) = [t, lastx',diag(lastP)'];
    posdif = lastx(1:3)-table2array(traj(j+1,2:4))';
    veldif = lastx(4:6)-table2array(traj(j+1,5:7))';
    % if norm(posdif)>=300 || norm(veldif)>=0.3
    %     warning("State diverged at time " + t);
    %     warncount = warncount+1;
    %     % if warncount >=4
    %     %     break
    %     % end
    % end
end

%% 3D Plotting
%X value
close all;
figure
hold on;
grid on;
plot3(states(:,2),states(:,3),states(:,4),DisplayName="EKF results");
plot3(table2array(traj(:,2)),table2array(traj(:,3)),table2array(traj(:,4)), ...
    DisplayName="True Trajectory");
ellipsoid(0,0,0,6378,6378,6378,20);%Adding in earth
legend("EKF Results","True Trajectory","Earth",Location="northwest");
xlabel("X axis (km)");
ylabel("Y axis (km)");
zlabel("Z axis (km)");

view(3);
axis equal;
exportgraphics(gca,"3Dplot.png",Resolution=600);

%% Tiles for x,y,z cordinates
figure("Name","X,Y,Z Position");
numtiles = 3;
xyzplot = tiledlayout(numtiles,1);

for i = 1:numtiles
    nexttile;
    hold on;
    plot(states(:,1),states(:,i+1),"-",DisplayName="Kalman",Color=[0.8500 0.3250 0.0980]);
    std = sqrt(states(:,i+7));%Standard deviation of current measure
    plot(states(:,1),states(:,i+1)+3*std,"-b",DisplayName="+3*sig");
    plot(states(:,1),states(:,i+1)-3*std,"-m",DisplayName="-3*sig");
    plot(table2array(traj(:,1)),table2array(traj(:,i+1)),DisplayName="Trajectory");
    legend(Location="eastoutside");
    grid on;
    xlabel("Time (s)");
    if i == 1
        ylabel("X position (km)");
    elseif i==2
        ylabel("Y position (km)");
    elseif i==3
        ylabel("Z position (km)");
    end
end
exportgraphics(xyzplot,"XYZPos.png",Resolution=600);

%% Tiles for x,y,z velocity
figure("Name","X,Y,Z Velocity");
numtiles = 3;
vxvyvzplot = tiledlayout(numtiles,1);

for i = 4:numtiles+3
    nexttile;
    hold on;
    plot(states(:,1),states(:,i+1),"-",DisplayName="Kalman",Color=[0.8500 0.3250 0.0980]);
    std = sqrt(states(:,i+7));%Standard deviation of current measure
    plot(states(:,1),states(:,i+1)+3*std,"-b",DisplayName="+3*sig");
    plot(states(:,1),states(:,i+1)-3*std,"-m",DisplayName="-3*sig");
    plot(table2array(traj(:,1)),table2array(traj(:,i+1)),DisplayName="Trajectory");
    legend(Location="eastoutside");
    grid on;
    xlabel("Time (s)");
    if i == 4
        ylabel("Velocity in X (km/s)");
    elseif i==5
        ylabel("Velocity in Y (km/s)");
    elseif i==6
        ylabel("Velocity in Z (km/s)");
    end
end
exportgraphics(vxvyvzplot,"XYZVel.png",Resolution=600);

%% Tiles for All error
figure("Name","Error of state");
numtiles = 2;
errplot = tiledlayout(numtiles,1);

err = states(:,1:7);
err(:,2:7) = err(:,2:7)-table2array(traj(:,2:7));

for j = 1:numtiles
    nexttile;
    hold on;
    ax = ["X","Y","Z"];
    specs = ["-r","--b","-.m"];
    for i = 1:3
        plot(err(:,1),err(:,i+1+(j-1)*3),specs(i),DisplayName=ax(i));
    end
    grid on;
    legend;
    xlabel("Time (s)");
    if j == 1
        ylabel("Error in pos(km)");
    elseif j==2
        ylabel("Error in vel(km/s)");
    end
end
exportgraphics(errplot,"StateErr.png",Resolution=600);

%% Tiles for All error with STD
figure("Name","Error of state with Deviation");
numtiles = 2;
errplot = tiledlayout(numtiles,1);

err = states(:,1:7);
err(:,2:7) = err(:,2:7)-table2array(traj(:,2:7));

for j = 1:numtiles
    nexttile;
    hold on;
    ax = ["X","Y","Z"];
    axSTD = ["3 sig X","3 sig Y","3 sig Z"];
    specs = ["-r","--b","-.m"];
    specsSTD = [":r",":b",":m"];
    for i = 1:3
        plot(err(:,1),err(:,i+1+(j-1)*3),specs(i),DisplayName=ax(i));
        plot(states(:,1),states(:,i+7+(j-1)*3),specsSTD(i),DisplayName=axSTD(i));
        plot(states(:,1),-states(:,i+7+(j-1)*3),specsSTD(i),DisplayName=axSTD(i));
    end
    grid on;
    legend(Location="eastoutside");
    xlabel("Time (s)");
    if j == 1
        ylabel("Error in pos(km)");
        ylim([-400,400]);
    elseif j==2
        ylabel("Error in vel(km/s)");
        ylim([-0.75,0.75]);
    end
end
exportgraphics(errplot,"StateErrStd.png",Resolution=600);


%% 3D Plotting with sig
% Lecture 11 to get this to work. Need all P matrices
figure
hold on;
grid on;
plot3(states(:,2),states(:,3),states(:,4),DisplayName="Kal");
plot3(table2array(traj(:,2)),table2array(traj(:,3)),table2array(traj(:,4)), ...
    DisplayName="Traj X");
xlabel("X axis (km)");
ylabel("Y axis (km)");
zlabel("Z axis (km)");


%Adding elipses
std = sqrt(states(:,8:10))*3;
for i =1:length(std)
    xc = states(i,2);yc = states(i,3);zc = states(i,4);
    xr = std(i,1);yr = std(i,2);zr = std(i,3);
    ellipsoid(xc,yc,zc,xr,yr,zr,10);
end
view(3);
axis equal;
exportgraphics(gca,"3DplotSig.png",Resolution=600);






