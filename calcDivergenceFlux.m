%% calcDivergenceFlux
%  calculates the divergence of the velocity field from velocity flux, div = du/dx + dv/dy
function [div]=calcDivergenceFlux(fluxU,fluxV,rho)
    global Ima Jma dx nG Ifim Ifi Ilam Ila Jfim Jfi Jlam Jla ;

    % Calculation
    div = zeros(Ima+2*nG,Jma+2*nG);
    div(Ifi:Ila,Jfi:Jla) = 0.5*dx*rho*( fluxU(Ifi:Ila,Jfi:Jla) - fluxU(Ifim:Ilam,Jfi:Jla)...
                                      + fluxV(Ifi:Ila,Jfi:Jla) - fluxV(Ifi:Ila,Jfim:Jlam) );
end
