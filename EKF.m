classdef EKF < handle
    properties
        % Estimate and Covariance
        mu_bar
        sig_bar
        mu
        sig
        % Histories
        mu_h
        sig_h
        mu_b_h
        sig_b_h
        K_h
        % Time Parameters
        dt
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
            self.mu_b_h = self.mu_h;
            self.sig_b_h = self.sig_h;
            self.K_h = zeros(3,2,NN);
            % self.dt = dt;
            self.Lm = Lm;
            self.QQ = [sig_r^2, 0; 0, sig_th^2];
            self.alph = alph;
        end
        function self = update(self, t_idx, dt, uu, zz)
            self.t_idx = t_idx;
            self.dt = dt;
            self.predict(uu);
            self.correct(zz);
        end
        function self = predict(self,uu)
            % Unicycle model: [x+vt*cos(th)*dt; y+vt*sin(th)*dt; th+wt*dt]
            % Unpack parameters
            th = self.mu(3);
            vt = uu(1);
            wt = uu(2);

            % Jacobian wrt state
            GG = eye(3);
            GG(1,3) = vt*sin(th)*self.dt;
            GG(2,3) = vt*cos(th)*self.dt;
            % Jacobian wrt control
            VV = [cos(th)*self.dt, 0;
                  sin(th)*self.dt, 0;
                  0, self.dt];
            % Motion noise
            MM = [self.alph(1)*vt^2, 0;
                  0, self.alph(2)*wt^2];

            % Update estimates
            self.mu_bar = self.mu + [vt*cos(th)*self.dt;...
                vt*sin(th)*self.dt; wt*self.dt];
            self.sig_bar = GG*self.sig*GG' + VV*MM*VV';

            % Update histories
            self.mu_b_h(:,self.t_idx) = self.mu_bar;
            self.sig_b_h(:,self.t_idx) = diag(self.sig_bar);
        end
        function self = correct(self,zz)
            pz = 1;
            for ii = 1:length(self.Lm)
                % Unpack landmark and estimate positions
                mx = self.Lm(1,ii);
                my = self.Lm(2,ii);
                mbx = self.mu_bar(1);
                mby = self.mu_bar(2);
                mbth = self.mu_bar(3);

                % Calculate predicted obsevation and Kalman gain
                qq = (mx - mbx)^2 + (my - mby)^2;
                zhat = [sqrt(qq);
                        atan2((my-mby),(mx-mbx))-mbth];
                HH = [-(mx-mbx)/sqrt(qq), -(my-mby)/sqrt(qq), 0;
                     (my-mby)/qq, -(mx-mbx)/qq, -1];
                SS = HH*self.sig_bar*HH' + self.QQ;
                KK = self.sig_bar*HH'/SS;

                % Update estimate
                self.mu_bar = self.mu_bar + KK*(zz(:,ii)-zhat);
                self.sig_bar = (eye(length(self.sig_bar))-KK*HH)*self.sig_bar;
                pz = pz * sqrt(det(2*pi*SS))*...
                    exp(-1/2*(zz(:,ii)-zhat)'/SS*(zz(:,ii)-zhat));
            end
            self.mu = self.mu_bar;
            self.sig = self.sig_bar;

            % Update histories
            self.mu_h(:,self.t_idx) = self.mu;
            self.sig_h(:,self.t_idx) = diag(self.sig);
            self.K_h(:,:,self.t_idx) = KK;
        end
        function [mh,sh,mbh,sbh,kh] = get_estimates(self)
            mbh = self.mu_b_h;
            sbh = self.sig_b_h;
            mh = self.mu_h;
            sh = self.sig_h;
            kh = self.K_h;
        end
    end
end
