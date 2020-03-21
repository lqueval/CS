function [a_lh,b_lh,c_lh,d_lh,mur_l_new] = fct_get_XX_nonlin(h_min,h_max,P,r_l,mur_l,mu0,K_lh_s,K_lh_c,isnonlin)
%---------------------------------------------------
%  NAME:      fct_get_XX_nonlin.m
%  WHAT:      iterative procedure for the nonlinear case
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHORS:   20200321, J. Guo, L. Quéval (loic.queval@gmail.com), F. Trillaud, B. Roucaries
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%
%  USE:
%    [a_lh,b_lh,c_lh,d_lh,mur_l_new] = fct_get_XX_nonlin(h_min,h_max,P,r_l,mur_l,mu0,K_lh_s,K_lh_c,isnonlin)
%
%  INPUTS:
%    h_min,h_max         = min and max harmonic number to be calculated
%    P                   = nb of pairs of pole
%    r_l                 = vector of current sheet radius
%    mur_l               = vector of (initial) domain relative permeability
%    mu0                 = permeability of vacuum
%    K_lh_s,K_lh_c       = Fourier coefficients of the current sheets
%    isnonlin            = 0 for linear, 1 for nonlinear case
%        
%  OUTPUTS:
%    a_lh,b_lh,c_lh,d_lh = coefficients of the x-matrix (Eq.6)
%    mur_l_new           = updated material of the machine
%----------------------------------------------------

N = length(mur_l);

%Initialize loop
if sum(isnonlin)>0
    loop_counter = 1; %nonlinear case
else
    loop_counter = 99; %linear case
end

mur_l_old = mur_l.^2; %nonlinear subdomains are initialized at something big
mur_l_new = mur_l;
e_mur = 1; %target error   

while ( max(abs(mur_l_old-mur_l_new))>=e_mur && loop_counter<100 )
    
    %get new value of mu_l
    mu_l = mur_l_new*mu0; %vector of current sheet permeability
    
    %get coefficients of XX
    [a_lh,b_lh,c_lh,d_lh] = fct_get_XX_lin(h_min,h_max,P,r_l,mu_l,K_lh_s,K_lh_c);
  
    %update
    mur_l_old = mur_l_new;
    
    % Get new mur in nonlinear domains
    for l = 1:N
        if isnonlin(l)>0 % loop on nonlinear domains only
            [THETA,RHO] = ndgrid((0:0.5:360)*pi/180,(r_l(l-1)+r_l(l))/2); %define field points in the middle of domain l
            [~,~,~,~,NORMXXX_j] = fct_get_B(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA); %get normB in domain l
            mur_l_new(l) = fct_mur_B(max(max(NORMXXX_j))); %get new value of mur
            mur_l_new(l) = mur_l_old(l) + (mur_l_new(l)-mur_l_old(l))/(1+loop_counter/10); %damping
        end
    end
    
    fprintf('== Counter = %d ==  \n', loop_counter); fprintf('mur_j = %f \n', mur_l_new);
    loop_counter = loop_counter+1;
    
end
end
