clear
clc

load('processed_data.mat')
% Contains landmarks, l_time, l_depth, l_bearing, odom_t, pos_odom_se2, vel_odom
load('truth_data.mat')
% Contains t_truth, x_truth, y_truth, and th_truth
% fix wrapping on th_truth
th_truth = rad_wrap_pi(th_truth);

% consistencize time vector names
lm_t = l_time;

% rng(0); % Fix the random number generator for debugging

% dt = 0.1;
Tf = odom_t(end);
% time = 0:dt:Tf;
time = odom_t;
NN = length(odom_t)+length(l_time);

% Field of operation is 20m by 20m, presumably 10 m from 0 in each
% direction.

% Robot initial conditions
x0 = 2.55; %m
y0 = -1.7; %m
th0 = 10 * pi/180; % pi/2; %rad, 90 degrees

% Create linear and angular velocity commands
uu = vel_odom;

% Create noisy output velocities
alph = [0.1, 0.1];

% Initialize truth data vectors, as well as observations
xx = [x0, zeros(1, NN-1)];
yy = [y0, zeros(1, NN-1)];
th = [th0, zeros(1, NN-1)];

% % Fill in truth data based on commanded velocities

XX = pos_odom_se2;

% Problem-defined range and bearing noise
% TODO: tune these values
sig_r = 0.1; % m
sig_th = 0.05; % rad
sig_rth = [sig_r; sig_th];

% Define landmark positions

% Lm = zeros(2,NL);
Lm = landmarks';
NLm = size(Lm,2); % Number of landmarks

% measurements = [ ]

% Initialize mu and sigma, then the extended Kalman filter
mu0 = [x0;y0;th0];
sig0 = [1 0 0; 0 1 0; 0 0 0.1];
ekf = EKF(mu0,sig0,sig_r,sig_th,alph,Lm,NLm,NN);

% Run the states through the EKF
l_idx = 1; % index of current l_time
o_idx = 1; % index of current odom_t
prev_t = 0;
for t_idx = 1:NN
    % if t_idx ~= NN || odom_t(o_idx) < l_time(l_idx)
    if o_idx <= length(odom_t) && (l_idx > length(l_time) || odom_t(o_idx) < l_time(l_idx))
    % if prev_t < odom_t(end) && (prev_t >= l_time(end) || odom_t(o_idx) < l_time(l_idx))

        tt = odom_t(o_idx);
        dt = tt - prev_t;

        ekf.predict(t_idx, dt, uu(:,o_idx));

        o_idx = o_idx + 1;
    % elseif t_idx == NN || l_time(l_idx) < odom_t(o_idx)
    elseif l_idx <= length(l_time) && (o_idx > length(odom_t) || l_time(l_idx) < odom_t(o_idx))
    % elseif prev_t < lm_t(end) && (prev_t >= odom_t(end) || odom_t(o_idx) < l_time(l_idx))

        tt = l_time(l_idx);
        dt = tt - prev_t;

        zz = [l_depth(:,l_idx)'; l_bearing(:,l_idx)'];
        ekf.correct(t_idx, dt, zz);

        l_idx = l_idx + 1;
    end

    prev_t = tt;
end


%% Plotting
[mu_h, sig_h, gains] = ekf.get_estimates();

labels = {'X Position', 'Y Position', 'Heading', 'X Error', 'Y Error',...
    'Heading Error'};

% Plot the robot animation
robotPlot = mobileRobotVis(XX(:,1),Lm,length(odom_t));
for tt = 2:length(odom_t)
    robotPlot.updatePlot(XX(:,tt));
    pause(0.001)
end

estimate_time = sort(unique([odom_t l_time']));
truth_data = zeros(3,length(estimate_time));
est_idx = 1;
for tru_idx = 1:length(t_truth)
    if t_truth(tru_idx) > estimate_time(est_idx)
        truth_data(:,est_idx) = [x_truth(tru_idx);y_truth(tru_idx);th_truth(tru_idx)];
        est_idx = est_idx + 1;
    end
end

% Plot true and estimated states versus time
figure(2)
for ii = 1:3
    subplot(3,1,ii); hold on;
    plot(odom_t,XX(ii,:),'color',[0 0.6588 0.8039]);
    plot(estimate_time,mu_h(ii,:),'color',[0 0.3098 0.4196],'LineStyle','--');
    xlabel('Time'); ylabel(labels(ii))
    legend('Actual','Measured');
end

figure(3)
% Plot error and covariance over time
for ii = 1:3
    subplot(3,1,ii); hold on;
    plot(odom_t,XX(ii,:)-mu_h(ii,:),'color',[0 0.6588 0.8039]);
    plot(estimate_time,2*sqrt(sig_h(ii,:)),'color',[0 0.3098 0.4196])
    plot(estimate_time,-2*sqrt(sig_h(ii,:)),'color',[0 0.3098 0.4196])
    xlabel('Time'); ylabel(labels(3+ii))
    legend('Error','95% Confidence');
end

var_labels = {'X ', 'Y ', 'Theta '};
meas_labels = {'Range Gain', 'Bearing Gain'};

% Plot gains over time
% figure(4)
% for ii = 1:3
%     for jj = 1:2
%         subplot(3,2,2*(ii-1)+jj)
%         plot(time,reshape(gains(ii,jj,:),[1,NN]));
%         xlabel('Time'); ylabel(strcat(var_labels(ii),meas_labels(jj)));
%     end
% end

% Plot actual range/bearing with measurments
% figure(3)
% for ii = 1:3
%     subplot(2,3,2*ii-1)
%     plot(time,reshape(z_true(1,ii,:),[1,NN]),'LineWidth',2);
%     hold on;
%     plot(time,reshape(zz(1,ii,:),[1,NN]),'LineWidth',0.5);
%     subplot(2,3,2*ii)
%     plot(time,reshape(z_true(2,ii,:),[1,NN]),'LineWidth',2);
%     hold on;
%     plot(time,reshape(zz(2,ii,:),[1,NN]),'LineWidth',0.5);
% end
