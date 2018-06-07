%% calcDivergence
%  calculates the divergence of the velocity field, div = du/dx + dv/dy
function [div]=calcDivergence(rhoU,rhoV)
    global Ima Jma nG Ifim Ifi Ifip Ilam Ila Ilap Jfim Jfi Jfip Jlam Jla Jlap;

    % Calculation
    div = zeros(Ima+2*nG,Jma+2*nG);
    div(Ifi:Ila,Jfi:Jla) = .5 * ( rhoU(Ifip:Ilap,Jfi:Jla) - rhoU(Ifim:Ilam,Jfi:Jla)...
                                + rhoV(Ifi:Ila,Jfip:Jlap) - rhoV(Ifi:Ila,Jfim:Jlam) );
end
