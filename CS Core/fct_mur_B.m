function [mur]=fct_mur_B(normB)
%---------------------------------------------------
%  NAME:      fct_mur_B.m
%  WHAT:      Calculate the relative permeability
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHORS:   20200321, J. Guo, L. Quéval (loic.queval@gmail.com)
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%
%  USE:
%    [mur]=fct_mur_B(normB);
%    Call the function with normB=[] to plot the mur(B) curve
%
%  INPUTS:
%    normB   = norm of the magnetic flux density
%
%  OUTPUTS:
%    mur     = relative permeability
%---------------------------------------------------

mu0 = 4*pi*1e-7;

% BH curve soft iron from Comsol (added few points above 2T up to 5T)
H_index =[663.146,1067.5,1705.23,2463.11,3841.67,5425.74,7957.75,12298.3,20462.8,32169.6,61213.4,111408,0.5e+06,1.5e+06,3.9789e+06]; %norm of H
B_index =[1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.3,2.6,5];  %norm of B
mur_index = B_index./H_index/mu0; %relative permeability mur

if ~isempty(normB)
    mur = interp1([0,B_index,10], [1200,mur_index,1], normB, 'linear','extrap'); %add 2 points to avoid nonphysical values for B=0 and B=inf
else
    mur = interp1([0,B_index,10], [1200,mur_index,1], B_index, 'linear','extrap'); %add 2 points to avoid nonphysical values for B=0 and B=inf
    
    %plot
    fig_iron_murB = figure(); box on, grid on, hold on
        plot([0,B_index],fct_mur_B([0,B_index]),'ob','MarkerSize',3) %data
        plot([0:0.01:B_index(end)],fct_mur_B([0:0.01:B_index(end)]),'-r') %fit
    xlabel('|B| [T]'), ylabel('\mu_r [ ]');
    legend('data','fit'), legend boxoff
    xlim([0,5])
    ylim([-fct_mur_B(0)*0.05,fct_mur_B(0)*1.05])
end

end