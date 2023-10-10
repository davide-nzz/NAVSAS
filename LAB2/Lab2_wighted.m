%------LS and WLS on realdataset------%
%% Reading data
clc
close all
clear

load 'DATA/RealisticUERE/dataset_1_20180329T160947'
addpath('Utilities')

T = 3600;

%% GPS
GPS_rho = RHO.GPS;
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
plot(n), title('GPS satellites')

figure, title('GPS rho'), hold on
for k = 1:N_sat
    plot(GPS_rho(k,:))
end
hold off

x = [];

for t = 1:T % for each time instant
    x_hat = zeros(1,4)'; % each time i start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.GPS(:,t))); % active satellites
        xyz_sats = []; 
        rho = [];

        %take satellites measurements
        for k = 1 : length(Sat_ind)


            xyz_sat = SAT_POS_ECEF.GPS(Sat_ind(k)).pos(t,:); %take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat]; % store them in a vector
            rho = [rho; RHO.GPS(Sat_ind(k),t)]; % take rho measurements
        end

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);
        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;
        H = [ax ay az ones(length(Sat_ind),1)];
        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;
        x_delta = H\rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt_ls = mean(x);

% GPS error for each time instant
er = zeros(T,4);
er_dist = zeros(T,1);
for t = 1:T
    er(t,:) = xt_ls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end
cov_GPS_ls = cov(er(:,1:3))
std_GPS_ls = sqrt(trace(cov_GPS_ls.^2))
cov_GPS_ls_2d = cov(er(:,1:2));
 %mean er over 450 sec
for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t GPS')




%--------WEIGHTED-----------%
x = [];

for k = 1:N_sat
    Sat_ind = find(~isnan(RHO.GPS(k,:))); % active period satellite
    if length(Sat_ind) > 1
        rho = RHO.GPS(k,Sat_ind)';
        rho_diff = diff(rho,2);
        UERE(k,1) = sqrt((rho_diff-mean(rho_diff))'*(rho_diff-mean(rho_diff))/length(rho_diff));
    else
        UERE(k,1) = 0;
    end
end

for t = 1:T % for each time instant
    x_hat = zeros(1,4)'; % each time i start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.GPS(:,t))); % active satellites
        R = diag(UERE(Sat_ind));
        W = R^(-1);
        xyz_sats = []; 
        rho = [];

        %take satellites measurements
        for k = 1 : length(Sat_ind)
            xyz_sat = SAT_POS_ECEF.GPS(Sat_ind(k)).pos(t,:); %take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat]; % store them in a vector
            rho = [rho; RHO.GPS(Sat_ind(k),t)]; % take rho measurements
        end

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);
        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;
        H = [ax ay az ones(length(Sat_ind),1)];
        H_w = (H'*W*H)\H'*W;
        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;
        x_delta = H_w*rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt_wls = mean(x);

figure
pcshow(x(:,1:3))
hold on
pcshow(xt_ls(1:3),'o')
pcshow(xt_wls(1:3),'o')
xlabel('x')
xlabel('y')
xlabel('z')
hold off
title('3D plot GPS')

figure
plot(x(:,1),x(:,2),'o')
hold on
plot(xt_ls(1),xt_ls(2),'o')
plot(xt_wls(1),xt_wls(2),'o')

title('xy GPS weighted estimation and regular')

%er for each time instant
er = zeros(T,4);
er_dist = zeros(T,1);
er_squared = zeros(T,4);
for t = 1:T
    er(t,:) = xt_wls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end
cov_GPS_wls = cov(er(:,1:3))
std_GPS_wls = sqrt(trace(cov_GPS_wls.^2))
cov_GPS_wls_2d = cov(er(:,1:2));
 %mean er over 450 sec
for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t GPS weighted')

[eigvec,eigval] = eig(cov_GPS_ls_2d);
ellipse_ls = ellipse(xt_ls(1:2),eigvec,eigval);
[eigvec,eigval] = eig(cov_GPS_wls_2d);
ellipse_wls = ellipse(xt_wls(1:2),eigvec,eigval);

figure
plot(x(:,1),x(:,2),'o')
hold on
plot(xt_ls(1),xt_ls(2),'o')
plot(xt_wls(1),xt_wls(2),'o')
plot(ellipse_ls(1,:),ellipse_ls(2,:),'linewidth',2);
plot(ellipse_wls(1,:),ellipse_wls(2,:),'linewidth',2);
legend('','ls','wls','ls','wls')
hold off
title('xy GPS wigthed and regular least squares')
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
plot(n), title('GALILEO satellites')

figure, title('GALILEO rho'), hold on
for k = 1:N_sat
    plot(GAL_rho(k,:))
end
hold off

x = [];

for t = 1:T % for each time instant
    x_hat = zeros(1,4)'; % each time i start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.GAL(:,t))); % active satellites
        xyz_sats = []; 
        rho = [];

        %take satellites measurements
        for k = 1 : length(Sat_ind)


            xyz_sat = SAT_POS_ECEF.GAL(Sat_ind(k)).pos(t,:); %take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat]; % store them in a vector
            rho = [rho; RHO.GAL(Sat_ind(k),t)]; % take rho measurements
        end

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);
        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;
        H = [ax ay az ones(length(Sat_ind),1)];
        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;
        x_delta = H\rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt_ls = mean(x);

% GALILEO error for each time instant
er = zeros(T,4);
er_dist = zeros(T,1);
for t = 1:T
    er(t,:) = xt_ls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end
cov_GAL_ls = cov(er(:,1:3))
std_GAL_ls = sqrt(trace(cov_GAL_ls.^2))

 %mean er over 450 sec
for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t GALILEO')

%--------WEIGHTED--------%
x = [];

for k = 1:N_sat
    Sat_ind = find(~isnan(RHO.GAL(k,:))); % active period satellite
    if length(Sat_ind) > 1
        rho = RHO.GAL(k,Sat_ind)';
        rho_diff = diff(rho,2);
        UERE(k,1) = sqrt((rho_diff-mean(rho_diff))'*(rho_diff-mean(rho_diff))/length(rho_diff));
    else
        UERE(k,1) = 0;
    end
end

for t = 1:T % for each time instant
    x_hat = zeros(1,4)'; % each time i start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.GAL(:,t))); % active satellites
        R = diag(UERE(Sat_ind));
        W = R^(-1);
        xyz_sats = []; 
        rho = [];

        %take satellites measurements
        for k = 1 : length(Sat_ind)
            xyz_sat = SAT_POS_ECEF.GAL(Sat_ind(k)).pos(t,:); %take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat]; % store them in a vector
            rho = [rho; RHO.GAL(Sat_ind(k),t)]; % take rho measurements
        end

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);
        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;
        H = [ax ay az ones(length(Sat_ind),1)];
        H_w = (H'*W*H)\H'*W;
        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;
        x_delta = H_w*rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt_wls = mean(x);

figure
pcshow(x(:,1:3))
hold on
pcshow(xt_ls(1:3),'o')
pcshow(xt_wls(1:3),'o')
xlabel('x')
xlabel('y')
xlabel('z')
hold off
title('3D plot GALILEO')

figure
plot(x(:,1),x(:,2),'o')
hold on
plot(xt_ls(1),xt_ls(2),'o')
plot(xt_wls(1),xt_wls(2),'o')
hold off
title('xy GALILEO weighted estimation and regular')

%er for each time instant
er = zeros(T,4);
er_dist = zeros(T,1);
er_squared = zeros(T,4);
for t = 1:T
    er(t,:) = xt_wls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end
cov_GAL_wls = cov(er(:,1:3))
std_GAL_wls = sqrt(trace(cov_GAL_wls.^2))

 %mean er over 450 sec
for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t GALILEO weighted')





%% GLONASS
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
plot(n), title('GLONASS satellites')

figure, title('GLONASS rho'), hold on
for k = 1:N_sat
    plot(GLO_rho(k,:))
end
hold off

x = [];

for t = 1:T % for each time instant
    x_hat = zeros(1,4)'; % each time i start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.GLO(:,t))); % active satellites
        xyz_sats = []; 
        rho = [];

        %take satellites measurements
        for k = 1 : length(Sat_ind)


            xyz_sat = SAT_POS_ECEF.GLO(Sat_ind(k)).pos(t,:); %take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat]; % store them in a vector
            rho = [rho; RHO.GLO(Sat_ind(k),t)]; % take rho measurements
        end

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);
        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;
        H = [ax ay az ones(length(Sat_ind),1)];
        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;
        x_delta = H\rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt_ls = mean(x);

% GPS error for each time instant
er = zeros(T,4);
er_dist = zeros(T,1);
for t = 1:T
    er(t,:) = xt_ls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end
cov_GLO_ls = cov(er(:,1:3))
std_GLO_ls = sqrt(trace(cov_GLO_ls.^2))

 %mean er over 450 sec
for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t GLONASS')

%--------WEIGHTED-----------%
x = [];

for k = 1:N_sat
    Sat_ind = find(~isnan(RHO.GLO(k,:))); % active period satellite
    if length(Sat_ind) > 1
        rho = RHO.GLO(k,Sat_ind)';
        rho_diff = diff(rho,2);
        UERE(k,1) = sqrt((rho_diff-mean(rho_diff))'*(rho_diff-mean(rho_diff))/length(rho_diff));
    else
        UERE(k,1) = 0;
    end
end

for t = 1:T % for each time instant
    x_hat = zeros(1,4)'; % each time i start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.GLO(:,t))); % active satellites
        R = diag(UERE(Sat_ind));
        W = R^(-1);
        xyz_sats = []; 
        rho = [];

        %take satellites measurements
        for k = 1 : length(Sat_ind)
            xyz_sat = SAT_POS_ECEF.GLO(Sat_ind(k)).pos(t,:); %take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat]; % store them in a vector
            rho = [rho; RHO.GLO(Sat_ind(k),t)]; % take rho measurements
        end

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);
        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;
        H = [ax ay az ones(length(Sat_ind),1)];
        H_w = (H'*W*H)\H'*W;
        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;
        x_delta = H_w*rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt_wls = mean(x);

figure
pcshow(x(:,1:3))
hold on
pcshow(xt_ls(1:3),'o')
pcshow(xt_wls(1:3),'o')
xlabel('x')
ylabel('y')
zlabel('z')
hold off
title('3D plot GLONASS')

figure
plot(x(:,1),x(:,2),'o')
hold on
plot(xt_ls(1),xt_ls(2),'o')
plot(xt_wls(1),xt_wls(2),'o')
hold off
title('xy GLONASS weighted estimation and regular')

%er for each time instant
er = zeros(T,4);
er_dist = zeros(T,1);
er_squared = zeros(T,4);
for t = 1:T
    er(t,:) = xt_wls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end
cov_GLO_wls = cov(er(:,1:3))
std_GLO_wls = sqrt(trace(cov_GLO_wls.^2))

 %mean er over 450 sec
for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t GLONASS weighted')




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
plot(n), title('BEIDOU satellites')

figure, title('BEIDOU rho'), hold on
for k = 1:N_sat
    plot(BEI_rho(k,:))
end
hold off

x = [];

for t = 1:T % for each time instant
    x_hat = zeros(1,4)'; % each time i start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.BEI(:,t))); % active satellites
        xyz_sats = []; 
        rho = [];

        %take satellites measurements
        for k = 1 : length(Sat_ind)


            xyz_sat = SAT_POS_ECEF.BEI(Sat_ind(k)).pos(t,:); %take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat]; % store them in a vector
            rho = [rho; RHO.BEI(Sat_ind(k),t)]; % take rho measurements
        end

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);
        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;
        H = [ax ay az ones(length(Sat_ind),1)];
        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;
        x_delta = H\rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt_ls = mean(x);

% GPS error for each time instant
er = zeros(T,4);
er_dist = zeros(T,1);
for t = 1:T
    er(t,:) = xt_ls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end
cov_BEI_ls = cov(er(:,1:3))
std_BEI_ls = sqrt(trace(cov_BEI_ls.^2))

 %mean er over 450 sec
for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t BEIDOU')

%--------WEIGHTED-----------%
x = [];

for k = 1:N_sat
    Sat_ind = find(~isnan(RHO.BEI(k,:))); % active period satellite
    if length(Sat_ind) > 1
        rho = RHO.BEI(k,Sat_ind)';
        rho_diff = diff(rho,2);
        UERE(k,1) = sqrt((rho_diff-mean(rho_diff))'*(rho_diff-mean(rho_diff))/length(rho_diff));
    else
        UERE(k,1) = 0;
    end
end

for t = 1:T % for each time instant
    x_hat = zeros(1,4)'; % each time i start from earth's center

    for s = 1:7 % 7 iterations for each time instant, since we start from far away

        Sat_ind = find(~isnan(RHO.BEI(:,t))); % active satellites
        R = diag(UERE(Sat_ind));
        W = R^(-1);
        xyz_sats = []; 
        rho = [];

        %take satellites measurements
        for k = 1 : length(Sat_ind)
            xyz_sat = SAT_POS_ECEF.BEI(Sat_ind(k)).pos(t,:); %take xyz of each active satellite
            xyz_sats = [xyz_sats;xyz_sat]; % store them in a vector
            rho = [rho; RHO.BEI(Sat_ind(k),t)]; % take rho measurements
        end

        r = sqrt((xyz_sats(:,1)-x_hat(1)).^2+(xyz_sats(:,2)-x_hat(2)).^2+(xyz_sats(:,3)-x_hat(3)).^2);
        ax = (xyz_sats(:,1) - x_hat(1))./r;
        ay = (xyz_sats(:,2) - x_hat(2))./r;
        az = (xyz_sats(:,3) - x_hat(3))./r;
        H = [ax ay az ones(length(Sat_ind),1)];
        H_w = (H'*W*H)\H'*W;
        rho_hat_x = xyz_sats(:,1) - x_hat(1);
        rho_hat_y = xyz_sats(:,2) - x_hat(2);
        rho_hat_z = xyz_sats(:,3) - x_hat(3);
        rho_hat = sqrt(rho_hat_x.^2+rho_hat_y.^2+rho_hat_z.^2);
        rho_delta = rho_hat - rho;
        x_delta = H_w*rho_delta;
        
        x_hat = x_hat + x_delta;
    end
    x = [x;x_hat'];
end

xt_wls = mean(x);

figure
pcshow(x(:,1:3))
hold on
pcshow(xt_ls(1:3),'o')
pcshow(xt_wls(1:3),'o')
xlabel('x')
ylabel('y')
zlabel('z')
hold off
title('3D plot BEIDOU')

figure
plot(x(:,1),x(:,2),'o')
hold on
plot(xt_ls(1),xt_ls(2),'o')
plot(xt_wls(1),xt_wls(2),'o')
hold off
title('xy BEIDOU weighted estimation and regular')

%er for each time instant
er = zeros(T,4);
er_dist = zeros(T,1);
er_squared = zeros(T,4);
for t = 1:T
    er(t,:) = xt_wls-x(t,:);
    er_dist(t) = norm(er(t,1:3));
end
cov_BEI_wls = cov(er(:,1:3))
std_BEI_wls = sqrt(trace(cov_BEI_wls.^2))

 %mean er over 450 sec
for k = 1:8
    er_mean(k) = mean(er_dist(1+(k-1)*450:k*450));
end

er_mean_hold = [ones(450,1)*er_mean(1);ones(450,1)*er_mean(2);ones(450,1)*er_mean(3);ones(450,1)*er_mean(4);ones(450,1)*er_mean(5);ones(450,1)*er_mean(6);ones(450,1)*er_mean(7);ones(450,1)*er_mean(8)];
figure, hold on
plot(er_dist),plot(er_mean_hold,'c')
hold off
title('error over t BEIDOU weighted')
