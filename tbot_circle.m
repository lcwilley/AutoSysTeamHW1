function circle_pts = tbot_circle(radius)

% receives radius, returns circle at origin

% ======================================

th = (0 : 0.02*pi : 2*pi)';
th(end) = th(1);

circ_pts = zeros(2,size(th,1));

r_vec = [radius; 0];

iter = 0;
for ii = th'
    iter = iter + 1;
    rotz = Rz_2d_I2B(ii);
    circle_pts(:,iter) = rotz * r_vec;
end

end

%
