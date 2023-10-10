%% Reading data
clc
close all
clear

load 'DATA/NominalUERE/dataset_6_20180328T121701'
addpath('Utilities')

T = 3600;

%% GPS visibility check

GPS_rho = RHO.GPS; %
n = zeros(T,1);
N_sat = length(SAT_POS_ECEF.GPS);

for t = 1 : T
    for k = 1 : N_sat
        if GPS_rho(k,t) > 10
            n(t) = n(t) + 1;
        end
    end
end

figure
plot(n)
title('GPS visible satellites')
xlabel('Time','FontWeight','bold')
ylabel('Number of satellites','FontWeight','bold')

figure
title('GPS measured pseudorange')
xlabel('Time','FontWeight','bold')
ylabel('Pseudorange','FontWeight','bold')
hold on
for k = 1:N_sat
    plot(GPS_rho(k,:))
end
hold off

%% GLONASS visibility check
GLO_rho = RHO.GLO;
n = zeros(T,1);
N_sat = length(SAT_POS_ECEF.GLO);

for t = 1 : T
    for k = 1 : N_sat
        if GLO_rho(k,t) > 10
            n(t) = n(t) + 1;
        end
    end
end

figure
plot(n)
title('GLONASS visible satellites')
xlabel('Time','FontWeight','bold')
ylabel('Number of satellites','FontWeight','bold')

figure
title('GLONASS measured pseudorange')
xlabel('Time','FontWeight','bold')
ylabel('Pseudorange','FontWeight','bold')
hold on
for k = 1:N_sat
    plot(GLO_rho(k,:))
end
hold off

%% BEIDOU
BEI_rho = RHO.BEI;
n = zeros(T,1);
N_sat = length(SAT_POS_ECEF.BEI);

for t = 1 : T
    for k = 1 : N_sat
        if BEI_rho(k,t) > 10
            n(t) = n(t) + 1;
        end
    end
end

figure
plot(n)
title('BEIDOU visible satellites')
xlabel('Time','FontWeight','bold')
ylabel('Number of satellites','FontWeight','bold')

figure
title('BEIDOU measured pseudorange')
xlabel('Time','FontWeight','bold')
ylabel('Pseudorange','FontWeight','bold')
hold on
for k = 1:N_sat
    plot(BEI_rho(k,:))
end
hold off

%% GALILEO
GAL_rho = RHO.GAL;
n = zeros(T,1);
N_sat = length(SAT_POS_ECEF.GAL);

for t = 1 : T
    for k = 1 : N_sat
        if GAL_rho(k,t) > 10
            n(t) = n(t) + 1;
        end
    end
end

figure
plot(n)
title('GALILEO visible satellites')
xlabel('Time','FontWeight','bold')
ylabel('Number of satellites','FontWeight','bold')

figure
title('GALILEO measured pseudorange')
xlabel('Time','FontWeight','bold')
ylabel('Pseudorange','FontWeight','bold')
hold on
for k = 1:N_sat
    plot(GAL_rho(k,:))
end
hold off

%% RLS soluion GPS

x = [];

for t = 1:T % for loop for each time instant
    x_hat = zeros(1,4)'; % each time I start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.GPS(:,t))); % active satellites
        xyz_sats = []; 
        rho = [];

        % take satellites measurements
        for k = 1 : length(Sat_ind)


            xyz_sat = SAT_POS_ECEF.GPS(Sat_ind(k)).pos(t,:); % take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat];                   % store them in a vector
            rho = [rho; RHO.GPS(Sat_ind(k),t)];              % take rho measurements
        end

        % geometrical distance between the linearization ponit and the satellite

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);

        % coefficients computation

        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;

        % Geometrical matrix

        H = [ax ay az ones(length(Sat_ind),1)];

        % computation of delta rho

        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;

        % Definition of delta x 

        x_delta = H\rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt = mean(x);

figure
pcshow(x(:,1:3))
hold on
pcshow(xt(1:3),'o')
xlabel('x')
ylabel('y')
zlabel('z')
hold off
title('3D plot GPS')

figure
plot(x(:,1),x(:,2),'o')
hold on
plot(xt(1),xt(2),'o')
hold off
title('xy GPS')

x_lla_GPS = ecef2lla(xt(1:3))'
'Cape Town coordinates'

xt_ls = mean(x);  

%% GPS Position error for each time instant

er = zeros(T,4);
er_dist = zeros(T,1);
for t = 1:T
    er(t,:) = xt_ls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end

cov_GPS_ls = cov(er(:,1:3))
std_GPS_ls = sqrt(trace(cov_GPS_ls.^2))

%mean error over 450 sec

for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t GPS')

%% RLS soluion GALILEO

x = [];

for t = 1:T % for loop for each time instant
    x_hat = zeros(1,4)'; % each time I start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.GAL(:,t))); % active satellites
        xyz_sats = []; 
        rho = [];

        % take satellites measurements
        for k = 1 : length(Sat_ind)


            xyz_sat = SAT_POS_ECEF.GAL(Sat_ind(k)).pos(t,:); % take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat];                   % store them in a vector
            rho = [rho; RHO.GAL(Sat_ind(k),t)];              % take rho measurements
        end

        % geometrical distance between the linearization ponit and the satellite

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);

        % coefficients computation

        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;

        % Geometrical matrix

        H = [ax ay az ones(length(Sat_ind),1)];

        % computation of delta rho

        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;

        % Definition of delta x 

        x_delta = H\rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt = mean(x);

% figure
% plot(x(:,1))
% title('x over time')
% figure
% plot(x(:,2))
% title('y over time')
% figure
% plot(x(:,3))
% title('z over time')

figure
pcshow(x(:,1:3))
hold on
pcshow(xt(1:3),'o')
xlabel('x')
ylabel('y')
zlabel('z')
hold off
title('3D plot GALILEO')

figure
plot(x(:,1),x(:,2),'o')
hold on
plot(xt(1),xt(2),'o')
hold off
title('xy GALILEO')

x_lla_GAL = ecef2lla(xt(1:3))'
'Cape Town coordinates'

%% GALILEO position error for each time instant

er = zeros(T,4);
er_dist = zeros(T,1);
for t = 1:T
    er(t,:) = xt_ls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end

cov_GAL_ls = cov(er(:,1:3))
std_GAL_ls = sqrt(trace(cov_GAL_ls.^2))

%mean error over 450 sec

for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t GALILEO')

%% RLS soluion GLO

x = [];

for t = 1:T % for loop for each time instant
    x_hat = zeros(1,4)'; % each time I start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.GLO(:,t))); % active satellites
        xyz_sats = []; 
        rho = [];

        % take satellites measurements
        for k = 1 : length(Sat_ind)


            xyz_sat = SAT_POS_ECEF.GLO(Sat_ind(k)).pos(t,:); % take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat];                   % store them in a vector
            rho = [rho; RHO.GLO(Sat_ind(k),t)];              % take rho measurements
        end

        % geometrical distance between the linearization ponit and the satellite

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);

        % coefficients computation

        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;

        % Geometrical matrix

        H = [ax ay az ones(length(Sat_ind),1)];

        % computation of delta rho

        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;

        % Definition of delta x 

        x_delta = H\rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt = mean(x);

% figure
% plot(x(:,1))
% title('x over time')
% figure
% plot(x(:,2))
% title('y over time')
% figure
% plot(x(:,3))
% title('z over time')

figure
pcshow(x(:,1:3))
hold on
pcshow(xt(1:3),'o')
xlabel('x')
ylabel('y')
zlabel('z')
hold off
title('3D plot GLONAS')

figure
plot(x(:,1),x(:,2),'o')
hold on
plot(xt(1),xt(2),'o')
hold off
title('xy GLONAS')

x_lla_GLO = ecef2lla(xt(1:3))'
'Cape Town coordinates'

%% GLONASS position error for each time instant

er = zeros(T,4);
er_dist = zeros(T,1);
for t = 1:T
    er(t,:) = xt_ls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end

cov_GLO_ls = cov(er(:,1:3))
std_GLO_ls = sqrt(trace(cov_GLO_ls.^2))

%mean error over 450 sec

for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t GLONASS')

%% RLS soluion BEI

x = [];

for t = 1:T % for loop for each time instant
    x_hat = zeros(1,4)'; % each time I start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.BEI(:,t))); % active satellites
        xyz_sats = []; 
        rho = [];

        % take satellites measurements
        for k = 1 : length(Sat_ind)


            xyz_sat = SAT_POS_ECEF.BEI(Sat_ind(k)).pos(t,:); % take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat];                   % store them in a vector
            rho = [rho; RHO.BEI(Sat_ind(k),t)];              % take rho measurements
        end

        % geometrical distance between the linearization ponit and the satellite

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);

        % coefficients computation

        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;

        % Geometrical matrix

        H = [ax ay az ones(length(Sat_ind),1)];

        % computation of delta rho

        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;

        % Definition of delta x 

        x_delta = H\rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt = mean(x);

% figure
% plot(x(:,1))
% title('x over time')
% figure
% plot(x(:,2))
% title('y over time')
% figure
% plot(x(:,3))
% title('z over time')

figure
pcshow(x(:,1:3))
hold on
pcshow(xt(1:3),'o')
xlabel('x')
ylabel('y')
zlabel('z')
hold off
title('3D plot BEIDOU')

figure
plot(x(:,1),x(:,2),'o')
hold on
plot(xt(1),xt(2),'o')
hold off
title('xy BEIDOU')

x_lla_BEI = ecef2lla(xt(1:3))'
'Cape Town coordinates'

%% BEIDOU position error for each time instant

er = zeros(T,4);
er_dist = zeros(T,1);
for t = 1:T
    er(t,:) = xt_ls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end

cov_BEI_ls = cov(er(:,1:3))
std_BEI_ls = sqrt(trace(cov_BEI_ls.^2))

%mean error over 450 sec

for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t BEIDOU')

%% Standard deviation comparison

 Constellations = {'GPS'; 'GALILEO'; 'GLONASS'; 'BEIDOU'};
 Standard_deviations = [std_GPS_ls std_GAL_ls std_GLO_ls std_BEI_ls];
 figure, bar(1:4,Standard_deviations)
 title('Standard deviation comparison')
 set(gca,'xticklabel',Constellations) 







