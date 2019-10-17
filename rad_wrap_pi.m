function angle = rad_wrap_pi(angle)

% wrap an angle in rads, -pi <= theta < pi

% ======================================

angle = angle - 2*pi * floor((angle + pi) * (0.5 / pi));

end

%
