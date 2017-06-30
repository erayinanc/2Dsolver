%% calcDivergenceFlux
%  calculates the divergence of the velocity field from velocity flux, div = du/dx + dv/dy

function [div]=calcDivergenceFlux(fluxU,fluxV)
    % Initialisation
    [Ima,~] = size(fluxU);
    [~,Jma] = size(fluxV);
    div = zeros(Ima+1,Jma+1);
    div(2:Ima-1,2:Jma-1) = ( fluxU(2:Ima-1,2:Jma-1)-fluxU(1:Ima-2,2:Jma-1)...
                           + fluxV(2:Ima-1,2:Jma-1)-fluxV(2:Ima-1,1:Jma-2));
end