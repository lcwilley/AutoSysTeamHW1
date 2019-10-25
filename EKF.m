classdef EKF < handle
    properties
        % Estimate and Covariance
        mu
        sig
        % Histories
        mu_h
        sig_h
        K_h
        % Time Parameters
        tt
        % Landmark Positions
        Lm
        % Noise variables
        QQ
        alph
    end
    methods
        % NN is total number of time pts
        function self = EKF(mu0, sig0, sig_r, sig_th, alph, Lm, NN)
            self.mu = mu0;
            self.sig = sig0;
            self.mu_h = [self.mu, zeros(3,NN-1)];
            self.sig_h = [diag(self.sig), zeros(3,NN-1)];
            self.K_h = zeros(3,2,NN);
            self.Lm = Lm;
            self.QQ = [sig_r^2, 0; 0, sig_th^2];
            self.alph = alph;
        end

        function self = predict(self, t_idx, dt, uu) %, odom_dt, trdt)
            % Unicycle model: [ x+vt*cos(th)*dt;
            %                   y+vt*sin(th)*dt;
            %                   th+wt*dt]
            % Unpack parameters
            th = self.mu(3);
            vt = uu(1);
            wt = uu(2);

            % Jacobian wrt state
            GG = eye(3);
            GG(1,3) = -vt*sin(th)*dt;
            GG(2,3) = vt*cos(th)*dt;
            % Jacobian wrt control
            VV = [cos(th)*dt, 0;
                  sin(th)*dt, 0;
                  0, dt];
            % Motion noise
            % MM  = [self.alph(1)*vt^2 + self.alph(2)*wt^2, 0;
            %         0, self.alph(3)*vt^2 + self.alph(2)*wt^2];
            % %
            MM  = [self.alph(1)*vt^2, 0;
                    0, self.alph(2)*wt^2];

            % Update estimates
            mu_update = [vt*cos(th)*dt;
                        vt*sin(th)*dt;
                        rad_wrap_pi(wt*dt)];
            %
            % test
            % update_err = mu_update - odom_dt

            self.mu = self.mu + mu_update;
            % self.mu = self.mu + odom_dt;
            self.mu(3) = rad_wrap_pi(self.mu(3));
            self.sig = GG*self.sig*GG' + VV*MM*VV';

            % Update histories
            self.mu_h(:,t_idx) = self.mu;
            self.sig_h(:,t_idx) = diag(self.sig);
        end

        function self = correct(self, t_idx, zz)
            % ezpz = 1;
            KK = NaN;
            for ii = 1:size(self.Lm,2)
                if ~any(isnan(zz(:,ii)),1)
                    % Unpack landmark and estimate positions
                    mx = self.Lm(1,ii);
                    my = self.Lm(2,ii);
                    mbx = self.mu(1);
                    mby = self.mu(2);
                    mbth = self.mu(3);

                    % Calculate predicted obsevation and Kalman gain
                    qq = (mx - mbx)^2 + (my - mby)^2;
                    zhat = [sqrt(qq);
                            rad_wrap_pi(atan2((my-mby),(mx-mbx))-mbth)];
                    HH = [-(mx-mbx)/sqrt(qq), -(my-mby)/sqrt(qq), 0;
                        (my-mby)/qq, -(mx-mbx)/qq, -1];
                    SS = HH*self.sig*HH' + self.QQ;
                    KK = self.sig*HH'/SS;

                    % Update estimate
                    K_innovation = KK*(zz(:,ii)-zhat);
                    self.mu = self.mu + K_innovation;
                    self.mu(3) = rad_wrap_pi(self.mu(3));
                    self.sig = (eye(length(self.sig))-KK*HH)*self.sig;

                end
            end

            % Update histories
            self.mu_h(:,t_idx) = self.mu;
            self.sig_h(:,t_idx) = diag(self.sig);
            self.K_h(:,:,t_idx) = KK;
        end

        function [mh,sh,kh] = get_estimates(self)
            mh = self.mu_h;
            sh = self.sig_h;
            kh = self.K_h;
        end
    end
end
