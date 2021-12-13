% Tara Urner
% November 2019

% The purpose of this code is to calculate the image distance behind the
% GRIN lenses in the V2 endoscope design. 


% https://www.grintech.de/fileadmin/Kunden/Datenblaetter/rev.07-2019/Introduction_Gradient_Index_Imaging_Optics.pdf

%n0 = 1.635;
%g  = pi/(2*2.23);
%l  = 2.46;
%fg = 1/n0 * 1/(g*sin(g*l))
%s  = 1/n0 * 1/(g*tan(g*l))

function [di,g,f,s,M] = ImageDistanceGRIN(n0,l,P,z)
% center refractive index = n0
% length = l
% pitch = P
% z is a range given as [far, near] to the GRIN lenses

%% GRIN Lens Parameters
g = pi/(2*2.23) % geometrical gradient constant
f = 1/(n0*g*sin(g*l)) % focal distance
s = 1/n0 * 1/(g*tan(g*l))

%% Image distance
do = z;

di = (1/f) - (1./do);

M = -di./do;

end