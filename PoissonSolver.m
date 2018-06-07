%% PoissonSolver
%  Solve (approximately) the Poisson equation for the pressure correction.
%  Find a field p such that div(grad(p)) = div
function [p,nIt,eps]=PoissonSolver(div,dt,omega,nItMax,epsMax) 
    global Ima Jma dx nG Ifim Ifi Ifip Ilam Ila Ilap Jfim Jfi Jfip Jlam Jla Jlap;

    % Initialisation
    p = zeros(Ima+2*nG,Jma+2*nG);

    % Loop over iteration steps. 
    eps = 1.0E+9; nIt = 0;
    while (nIt < nItMax) && (eps > epsMax)
          % Increment iteration counter
            nIt = nIt + 1;
            
            % Store pressure field of last iteration step
            pold = p;

            % Compute new pressure with Jacobi point iteration
            p(Ifi:Ila,Jfi:Jla) = 0.25 * (pold(Ifip:Ilap,Jfi:Jla) + pold(Ifim:Ilam,Jfi:Jla) ...
                                       + pold(Ifi:Ila,Jfip:Jlap) + pold(Ifi:Ila,Jfim:Jlam) ...
                                       - omega*(div(Ifi:Ila,Jfi:Jla)*dx/dt) );
                                  
            % Set zero gradient BC at the inlet and fixed value of zero at all            
            p(Ifim,:) = p(Ifi,:); p(Ilap,:) = 0;
            
            % Compute maximum error
            eps = max (max (abs (p-pold)));  
    end
end   
