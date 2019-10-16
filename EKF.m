classdef EKF < handle
    properties
        mu
        sig
        mu_h
        sig_h
        mu_bar
        sig_bar
        mu_b_h
        sig_b_h
        K_h
        dt
        Lm
        Q
        alph
        t
    end
    methods
        function self = EKF(mu0, sig0, dt, sig_r, sig_th, alph, Lm, N)
            self.mu = mu0;
            self.sig = sig0;
            self.mu_h = [self.mu, zeros(3,N-1)];
            self.sig_h = [diag(self.sig), zeros(3,N-1)];
            self.mu_b_h = self.mu_h;
            self.sig_b_h = self.sig_h;
            self.K_h = zeros(3,2,N);
            self.dt = dt;
            self.Lm = Lm;
            self.Q = [sig_r^2, 0; 0, sig_th^2];
            self.alph = alph;
        end
        function self = update(self,t,u,z)
            self.t = t;
            self.predict(u);
            self.correct(z);
        end
        function self = predict(self,u)
            % Unicycle model: [vt*cos(th); vt*sin(th); wt]
            th = self.mu(3);
            vt = u(1);
            wt = u(2);
            G = eye(3);
            G(1,3) = vt*sin(th);
            G(2,3) = vt*cos(th);
            V = [cos(th), 0;
                 sin(th), 0;
                 0, self.dt];
            M = [self.alph(1)*vt^2, 0;
                  0, self.alph(2)*wt^2];
              
            % Update estimates
            self.mu_bar = self.mu + [vt*cos(th); vt*sin(th); wt*self.dt];
            self.sig_bar = G*self.sig*G' + V*M*V';
            
            % Update histories
            self.mu_b_h(:,self.t) = self.mu_bar;
            self.sig_b_h(:,self.t) = diag(self.sig_bar);
        end
        function self = correct(self,z)
            pz = 1;
            for i = 1:length(self.Lm)
                mx = self.Lm(1,i);
                my = self.Lm(2,i);
                mbx = self.mu_bar(1);
                mby = self.mu_bar(2);
                mbth = self.mu_bar(3);
                q = (mx - mbx)^2 + (my - mby)^2;
                zhat = [sqrt(q);
                        atan2((my-mby),(mx-mbx))-mbth];
                H = [-(mx-mbx)/sqrt(q), -(my-mby)/sqrt(q), 0;
                     (my-mby)/q, -(mx-mbx)/q, -1];
                S = H*self.sig_bar*H' + self.Q;
                K = self.sig_bar*H'/S;
                self.mu_bar = self.mu_bar + K*(z(:,i)-zhat);
                self.sig_bar = (eye(length(self.sig_bar))-K*H)*self.sig_bar;
                pz = pz * sqrt(det(2*pi*S))*...
                    exp(-1/2*(z(:,i)-zhat)'/S*(z(:,i)-zhat));
            end
            self.mu = self.mu_bar;
            self.sig = self.sig_bar;
            self.mu_h(:,self.t) = self.mu;
            self.sig_h(:,self.t) = diag(self.sig);
            self.K_h(:,:,self.t) = K;
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