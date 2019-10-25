clear
clc

load('processed_data.mat')
% Contains landmarks, l_time, l_depth, l_bearing, odom_t, pos_odom_se2, vel_odom
load('truth_data.mat')
% Contains t_truth, x_truth, y_truth, and th_truth
% Fix wrapping on th_truth
th_truth = rad_wrap_pi(th_truth);
meas_t = l_time;

% Make time vector names consistent
% lm_t = meas_t;

% Define the total number of time steps
NN = length(odom_t)+length(meas_t);

% Robot initial conditions
x0 = 2.55; %m
y0 = -1.7; %m
th0 = 10 * pi/180; % pi/2; %rad, 90 degrees

% Create linear and angular velocity commands
uu = vel_odom;

% Create noisy output velocities
% alph = [0.1, 0.01, 0.01, 0.1];
alph = [0.2, 0.2];

% Initialize truth data vectors, as well as observations
xx = [x0, zeros(1, NN-1)];
yy = [y0, zeros(1, NN-1)];
th = [th0, zeros(1, NN-1)];

% Truth data from the odemetry data
XX = pos_odom_se2;

% Problem-defined range and bearing noise
% TODO: tune these values
sig_r = 0.6; % m
sig_th = 0.3; % rad
sig_rth = [sig_r; sig_th];

% Define landmark positions
Lm = landmarks';
NLm = size(Lm,2); % Number of landmarks

% Initialize mu and sigma, then the extended Kalman filter
mu0 = [x0;y0;th0];
sig0 = [0.1 0 0; 0 0.1 0; 0 0 0.1];
ekf = EKF(mu0,sig0,sig_r,sig_th,alph,Lm,NN);


%
[~, order] = sort([odom_t meas_t']);
all_times = [odom_t,                meas_t';
            ones(1,length(odom_t)), 2*ones(1,length(meas_t))];
% estimate_time = all_times(1,order);
estimate_type = all_times(2,order);

% Run the states through the EKF
all_t = [odom_t';meas_t];
all_t_sort = sort(all_t);
l_idx = 0; % index of current meas_t
% o_idx = 1; % index of current odom_t
command_idx = 0;
prev_t = 0;
uu_command = uu(:,1);
for t_idx = 1:NN
    % if o_idx <= length(odom_t) && (l_idx > length(meas_t) || odom_t(o_idx) < meas_t(l_idx))
    %     % o_idx = o_idx + 1;
    % end
    % use odom_t

    % tt = odom_t(min(o_idx, length(o_idx)));
    tt = all_t_sort(t_idx);
    dt = tt - prev_t;

    if estimate_type(t_idx) == 1  % odom
        command_idx = command_idx + 1;
        uu_command = uu(:, command_idx);
    end
    ekf.predict(t_idx, dt, uu_command); %, odom_dt, trdt);

    if estimate_type(t_idx) == 2  % measurement
        % update l_idx
        l_idx = l_idx + 1;
        zz = [l_depth(:,l_idx)'; l_bearing(:,l_idx)'];
        ekf.correct(t_idx, zz);

    end





    % if o_idx <= length(odom_t) && (l_idx > length(meas_t) || odom_t(o_idx) < meas_t(l_idx))
    %     o_idx = o_idx + 1;
    % end

    % if l_idx <= length(meas_t) && (o_idx > length(odom_t) || meas_t(l_idx) < odom_t(o_idx))
        % use meas_t

        % tt = meas_t(l_idx);
        % dt = tt - prev_t;




    % end

    prev_t = tt;
end


%% Plotting
[mu_h, sig_h, gains] = ekf.get_estimates();

labels = {'X Position', 'Y Position', 'Heading', 'X Error', 'Y Error',...
    'Heading Error'};

% % Plot the robot animation
% robotPlot = mobileRobotVis(XX(:,1),Lm,length(odom_t));
% for tt = 2:length(odom_t)
%     robotPlot.updatePlot(XX(:,tt));
%     pause(0.001)
% end

estimate_time = sort(unique([odom_t meas_t']));
truth_data = zeros(3,length(estimate_time));
odom_data = zeros(3,length(estimate_time));
for est_idx = 1:length(estimate_time)
    [~,tru_idx] = min(abs(t_truth-estimate_time(est_idx)));
    truth_data(:,est_idx) = [x_truth(tru_idx);
                            y_truth(tru_idx);
                            th_truth(tru_idx)];
    [~,o_idx] = min(abs(odom_t-estimate_time(est_idx)));
    odom_data(:,est_idx) = pos_odom_se2(:,o_idx);
end

% Plot true and estimated states versus time
f2 = figure(2);
clf(f2);
for ii = 1:3
    subplot(3,1,ii); hold on;
    plot(estimate_time,truth_data(ii,:),'color',[0 0.6588 0.8039]);
    plot(estimate_time,mu_h(ii,:),'color',[0 0.3098 0.4196],'LineStyle','--');
    xlabel('Time'); ylabel(labels(ii))
    legend('Actual','Predicted');
end

f3 = figure(3);
clf(f3);
% Plot error and covariance over time
for ii = 1:3
    subplot(3,1,ii); hold on;
    if ii < 3
        plot(estimate_time,truth_data(ii,:)-mu_h(ii,:),'color',[0 0.6588 0.8039]);
    elseif ii == 3
        heading_err = rad_wrap_pi(truth_data(ii,:)-mu_h(ii,:));
        plot(estimate_time, heading_err,'color',[0, 0.6588, 0.8039]);
    end
    plot(estimate_time,2*sqrt(sig_h(ii,:)),'color',[0, 0.3098, 0.4196])
    plot(estimate_time,-2*sqrt(sig_h(ii,:)),'color',[0, 0.3098, 0.4196])
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
