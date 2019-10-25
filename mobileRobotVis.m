classdef mobileRobotVis < handle
    properties
        radius
        robot_handle
        bearing_handle
        state_history
        history_handle
        tt
        robot_circ_pts
    end
    methods
        function self = mobileRobotVis(X0,Lm,NN)
            % Unpack state
            xx = X0(1);
            yy = X0(2);
            th = X0(3);

            self.radius = 0.5; % Robot radius

            f1 = figure(1)
            clf(f1)
            hold on;

            % Record and plot robot position history
            self.tt = 1;
            self.state_history = zeros(3,NN);
            self.state_history(:,self.tt) = X0;
            self.history_handle = plot(self.state_history(1,1),...
                self.state_history(2,1),'b');
            self.tt = self.tt + 1;

            % Define the robot coordinates
            % robot_cors = circle([X0(1),X0(2)],self.radius,'n',100);
            self.robot_circ_pts = tbot_circle(self.radius);
            robot_cors = self.robot_circ_pts + [xx; yy];
            self.robot_handle = fill(robot_cors(1,:),...
            	robot_cors(2,:),[0.4,0.6,1]);
            self.bearing_handle = plot([xx,xx+self.radius*cos(th)],...
                [yy,yy+self.radius*sin(th)],'k','LineWidth',2);

            % Plot the landmarks (assumes each landmark is a column vector)
            for ii = 1:length(Lm(1,:))
                scatter(Lm(1,ii),Lm(2,ii),100,'ro','filled')
            end

            % Define the arena size
            axis([-10,10,-10,10])
            axis equal
        end
        function self = updatePlot(self,XX)
            % Unpack state
            xx = XX(1);
            yy = XX(2);
            th = XX(3);

            %Update the robot trail
            self.state_history(:,self.tt) = XX;
            self.history_handle.XData = self.state_history(1,1:self.tt);
            self.history_handle.YData = self.state_history(2,1:self.tt);

            % Update robot position
            % robot_cors = circle([XX(1),XX(2)],self.radius,'n',100);
            robot_cors = self.robot_circ_pts + [xx; yy];
            self.robot_handle.XData = robot_cors(1,:);
            self.robot_handle.YData = robot_cors(2,:);
            self.bearing_handle.XData = [xx,xx+self.radius*cos(th)];
            self.bearing_handle.YData = [yy,yy+self.radius*sin(th)];

            % Update time step
            self.tt = self.tt + 1;
        end
    end
end
