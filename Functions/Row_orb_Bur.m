function [ r0 ] = Row_orb_bur( par )
%ROWNANIE_ORBITY wyznacza promień orbity stacjonarnej wiru Burgersa bez
%wpływu grawitacji biorąc parametr par=St/A.

syms r
r0_sym=solve((1/par)^(1/2)*r^2-(1-exp(-r^2/2))/(2*pi)==0);

r0=abs(double(r0_sym(2)));
end