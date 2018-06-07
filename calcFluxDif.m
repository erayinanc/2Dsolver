%% Caculate the Diffusive fluxes
%  fluxDifX, fluxDifY over the cell surfaces as a function of
%  a conserved scalar field Phi, the diffusivity D and a 
%  grid-spacing of deltaX
function [fluxDifX,fluxDifY]=calcFluxDif(Phi,D)
    global dx Ifim Ifi Ila Ilap Jfim Jfi Jla Jlap;

    % Calculation
    fluxDifX = D/dx.*(Phi(Ifi:Ilap,:) - Phi(Ifim:Ila,:));
    fluxDifY = D/dx.*(Phi(:,Jfi:Jlap) - Phi(:,Jfim:Jla));
end
