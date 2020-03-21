function [] = fct_plot_geom_generic(r_l,is_currentsheet)
%---------------------------------------------------
%  NAME:      fct_plot_geom_generic.m
%  WHAT:      Plot current sheet model geometry
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHORS:   20200321, J. Guo, L. Quéval (loic.queval@gmail.com)
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%
%  USE:
%    [] = fct_plot_geom_generic(r_j,is_currentsheet)
%
%  INPUTS:
%    r_j              = vector of current sheet radius
%    is_currentsheet  = 0 for nothing, 1 for current sheet
%
%  OUTPUTS:
%    
%---------------------------------------------------

for j=1:length(r_l)
    
    if is_currentsheet(j) == 1
        [x,y,~] = cylinder(r_l(j),200); plot(x(1,:),y(1,:),'r'); %plot red circle at r_l
    else
        [x,y,~] = cylinder(r_l(j),200); plot(x(1,:),y(1,:),'k'); %plot black circle at r_l
    end
    
end

axis square