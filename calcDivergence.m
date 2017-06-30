%% calcDivergence
%  calculates the divergence of the velocity field, div = du/dx + dv/dy

function [div]=calcDivergence(U,V)
    % Initialisation
    [Ima,~] = size(U);
    [~,Jma] = size(V);
    div = zeros(Ima,Jma);
    div(2:Ima-1,2:Jma-1) = ( U(3:Ima,2:Jma-1)-U(1:Ima-2,2:Jma-1)...
                           + V(2:Ima-1,3:Jma)-V(2:Ima-1,1:Jma-2) );
end