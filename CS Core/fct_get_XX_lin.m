function [a_lh,b_lh,c_lh,d_lh] = fct_get_XX_lin(h_min,h_max,P,r_l,mu_l,K_lh_s,K_lh_c)
%---------------------------------------------------
%  NAME:      fct_get_XX_lin.m
%  WHAT:      Build and solve linear system (Eq.6)
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHORS:   20200321, F. Trillaud, J. Guo, L. Quéval (loic.queval@gmail.com) 
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%
%  USE:
%    fct_get_XX_lin(h_min,h_max,P,r_l,mu_l,K_lh_s,K_lh_c)
%
%  INPUTS:
%    h_min,h_max         = min and max harmonic number to be calculated
%    P                   = nb of pairs of pole
%    r_l                 = vector of current sheet radius
%    mu_l                = vector of domain permeability
%    K_lh_s,K_lh_c       = Fourier coefficients of the current sheets
%        
%  OUTPUTS:
%    a_lh,b_lh,c_lh,d_lh = coefficients of the x-matrix (Eq.6)
%----------------------------------------------------

% initialization
N = length(mu_l); %nb of domains
RR = zeros(4*N,4*N);
BB = zeros(4*N,1);

RR(1,2) = 1;
RR(2,4) = 1;

BB(1)=0;
BB(2)=0;

RR(4*N-1,4*N-3) = 1;
RR(4*N,4*N-1) = 1;

BB(4*N-1) = 0;
BB(4*N) = 0;

for h=h_min:1:h_max
    % Build matrix
    for l=1:N-1
        i = 4*(l-1);
        r_l_phP = r_l(l)^(h*P);
        r_l_mhP = r_l(l)^(-h*P);
        
        RR(i+3,i+1) = -r_l_phP;
        RR(i+3,i+2) = -r_l_mhP;
        RR(i+3,i+5) = r_l_phP;
        RR(i+3,i+6) = r_l_mhP;
        
        RR(i+4,i+3) = -r_l_phP;
        RR(i+4,i+4) = -r_l_mhP;
        RR(i+4,i+7) = r_l_phP;
        RR(i+4,i+8) = r_l_mhP;
        
        RR(i+5,i+1) = r_l_phP/mu_l(l);
        RR(i+5,i+2) = -r_l_mhP/mu_l(l);
        RR(i+5,i+5) = -r_l_phP/mu_l(l+1);
        RR(i+5,i+6) = r_l_mhP/mu_l(l+1);
        
        RR(i+6,i+3) = r_l_phP/mu_l(l);
        RR(i+6,i+4) = -r_l_mhP/mu_l(l);
        RR(i+6,i+7) = -r_l_phP/mu_l(l+1);
        RR(i+6,i+8) = r_l_mhP/mu_l(l+1);
        
        BB(i+5) = K_lh_s(l,h)*r_l(l)/h/P;
        BB(i+6) = K_lh_c(l,h)*r_l(l)/h/P;
    end
    
    % Invert system
    XX = lsqlin(RR,BB);
    
    % get coeff from AA
    for l = 1:N
        i = 4*(l-1);
        a_lh(l,h)= XX(i+1);
        b_lh(l,h)= XX(i+2);
        c_lh(l,h)= XX(i+3);
        d_lh(l,h)= XX(i+4);
    end
end
