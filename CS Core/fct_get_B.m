function [RHO,THETA,BR,BTH,NORMB] = fct_get_B(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA)
%---------------------------------------------------
%  NAME:      fct_get_B.m
%  WHAT:      Calculates B at field points.
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHORS:   20200321, J. Guo, L. Quéval (loic.queval@gmail.com)
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%
%  USE:
%    [RHO,THETA,BR,BTH,NORMB] = fct_get_B(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA);
%
%  INPUTS:
%    a_lh,b_lh,c_lh,d_lh   = coefficients of the x-matrix (Eq.6)
%    h_min,h_max           = min and max harmonic number to be calculated
%    P          = nb of pairs of pole
%    r_l        = geometry of the machine
%    RHO        = Field points rho-coordinate vector or matrix
%    THETA      = Field points theta-coordinate vector or matrix
%
%  OUTPUTS:
%    RHO        = Field points rho-coordinate vector or matrix
%    THETA      = Field points theta-coordinate vector or matrix
%    BR         = Field points B rho-component vector or matrix
%    BTH        = Field points B theta-component vector or matrix
%    NORMB      = Field points B norm vector or matrix
%---------------------------------------------------
            
% initialize (BR, BTH)
BR = zeros(size(RHO,1),size(RHO,2));
BTH = zeros(size(RHO,1),size(RHO,2));

for h=h_min:1:h_max %loop on harmonic number h
    for m=1:size(RHO,1) %loop on field points
        for q=1:size(RHO,2)
            rho_loop = RHO(m,q);
            theta_loop = THETA(m,q);
            sector = sum(rho_loop>r_l)+1; %find the sector number for a given radius
            
                BR(m,q) = (a_lh(sector,h)*rho_loop^(h*P-1) + b_lh(sector,h)*rho_loop^(-h*P-1))*h*P*cos(h*P*theta_loop) - (c_lh(sector,h)*rho_loop^(h*P-1) + d_lh(sector,h)*rho_loop^(-h*P-1))*h*P*sin(h*P*theta_loop)+BR(m,q);
                BTH(m,q) =(-a_lh(sector,h)*rho_loop^(h*P-1) + b_lh(sector,h)*rho_loop^(-h*P-1))*h*P*sin(h*P*theta_loop) + (-c_lh(sector,h)*rho_loop^(h*P-1) + d_lh(sector,h)*rho_loop^(-h*P-1))*h*P*cos(h*P*theta_loop)+BTH(m,q);
            
        end
    end
end
NORMB = sqrt(BR.^2+BTH.^2);
