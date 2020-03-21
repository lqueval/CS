function [ K_4h_c,K_4h_s ] = fct_get_K4( Na,i_a,i_b,i_c,h_min,h_max,theta_1a,theta_2a,wa )
%---------------------------------------------------
%  NAME:      fct_get_K4.m
%  WHAT:      Get current sheet K4
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHORS:   20200321, J. Guo, L. Quéval (loic.queval@gmail.com)
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%
%  USE:
%    [ K_4h_c,K_4h_s ] = fct_get_K4( Na,i_a,i_b,i_c,h_min,h_max,theta_1a,theta_2a,wa );
%
%  INPUTS:
%    Na            = nb of turns of the armature windings
%    i_a,i_b,i_c   = instantaneous 3-ph currents
%    h_min,h_max   = min and max harmonic number to be calculated
%    theta_1a      = coil width electrical angle
%    theta_2a      = coil aperture electrical angle
%    wa            = coil width
%
%  OUTPUTS:
%    K_4h_c,K_4h_s = Fourier coefficients of the K3 current sheet
%---------------------------------------------------

for h=h_min:1:h_max
    K_temp = 4*Na/pi/wa/h*sin((theta_1a+theta_2a)*h/2)*sin(theta_1a*h/2); 
    K_4h_s(h)= K_temp*(i_a + i_b*cos(h*2*pi/3) + i_c*cos(h*4*pi/3));
    K_4h_c(h)= -K_temp*(i_b*sin(h*2*pi/3)+ i_c*sin(h*4*pi/3));
end

end

