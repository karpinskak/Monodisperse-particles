function [ur,ufi_przez_r,uz] = velocity_field(A,r,z)
%VELOCITY_FIELD Gives values of dimensionless Burgers vortex velocity field
%elements
ur=-A*r;
ufi_przez_r=(1-exp(-r.^2/2))./(2*pi*r.^2);
uz=2*A*z;
end

