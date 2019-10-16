classdef mobileRobotVis < handle
    properties
        radius
        robot_handle
        bearing_handle
        state_history
        history_handle
        t
    end
    methods
        function self = mobileRobotVis(X0,L,dt,Tf)
            % Unpack state
            x = X0(1);
            y = X0(2);
            th = X0(3);
            
            self.radius = 0.5; % Robot radius
            
            figure(1)
            hold on;
            
            % Record and plot robot position history
            self.t = 1;
            self.state_history = zeros(3,Tf/dt + 1);
            self.state_history(:,self.t) = X0;
            self.history_handle = plot(self.state_history(1,1),...
                self.state_history(2,1),'b');
            self.t = self.t + 1;
            
            % Define the robot coordinates
            robot_cors = circle([X0(1),X0(2)],self.radius,'n',100);
            self.robot_handle = fill(robot_cors(1,:),...
            	robot_cors(2,:),[0.4,0.6,1]);
            self.bearing_handle = plot([x,x+self.radius*cos(th)],...
                [y,y+self.radius*sin(th)],'k','LineWidth',2);
            
            % Plot the landmarks (assumes each landmark is a column vector)
            for i = 1:length(L(1,:))
                scatter(L(1,i),L(2,i),100,'ro','filled')
            end
            
            % Define the arena size
            axis([-10,10,-10,10])
            axis equal
        end
        function self = updatePlot(self,X)
            % Unpack state
            x = X(1);
            y = X(2);
            th = X(3);
            
            %Update the robot trail
            self.state_history(:,self.t) = X;
            self.history_handle.XData = self.state_history(1,1:self.t);
            self.history_handle.YData = self.state_history(2,1:self.t);
            
            % Update robot position
            robot_cors = circle([X(1),X(2)],self.radius,'n',100);
            self.robot_handle.XData = robot_cors(1,:);
            self.robot_handle.YData = robot_cors(2,:);
            self.bearing_handle.XData = [x,x+self.radius*cos(th)];
            self.bearing_handle.YData = [y,y+self.radius*sin(th)];
            
            % Update time step
            self.t = self.t + 1;
        end
    end
end