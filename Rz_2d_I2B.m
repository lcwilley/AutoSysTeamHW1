function rz_2d_i2b = Rz_2d_I2B(psi)

% Rotation matrix about z-axis, in 2-D.
% This is the inertial_2_body rotation
% For batch processing, the Transpose of the output
% will give you a set of b2i matrices. Transpose each
% in turn to get i2b matrices.

% ======================================

cps = cos( psi );
sps = sin( psi );

rz_2d_i2b = [cps, sps;-sps, cps];

%
