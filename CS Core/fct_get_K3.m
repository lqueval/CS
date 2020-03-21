function [ K_3h_c,K_3h_s ] = fct_get_K3( h_min,h_max,Nf,i_f,theta_1f,theta_2f,alfa,P,wf )
%---------------------------------------------------
%  NAME:      fct_get_K3.m
%  WHAT:      Get current sheet K3
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHORS:   20200321, J. Guo, L. Quéval (loic.queval@gmail.com)
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%
%  USE:
%    [ K_3h_c,K_3h_s ] = fct_get_K3( h_min,h_max,Nf,i_f,theta_1f,theta_2f,alfa,P,wf );
%
%  INPUTS:
%    h_min,h_max = min and max harmonic number to be calculated
%    Nf         = nb of turns of the rotor winding
%    i_f        = instantaneous field coil current
%    theta_1f   = coil width electrical angle
%    theta_2f   = coil aperture electrical angle
%    alfa       = rotor mechanical angle
%    P          = nb of pairs of pole
%    wf         = coil width
%
%  OUTPUTS:
%    K_3h_c,K_3h_s = Fourier coefficients of the K3 current sheet
%---------------------------------------------------

for h=h_min:2:h_max %h is odd
    K_temp = 8*Nf*i_f/pi/wf/h*sin((theta_1f+theta_2f)*h/2)*sin(theta_1f*h/2);  
    K_3h_s(h) = K_temp*cos(h*P*alfa);
    K_3h_c(h) = -K_temp*sin(h*P*alfa);
end

end

