function [R, dRdx, T, U, Ud, PHI] = SHOOTINGRES(m, Uw, Fvec, harms, Ntcyc, Ncyc)
%SHOOTINGRES is the residue function one can use for shooting
%
%   USAGE:
%       [R, dRdx] = m.SHOOTINGRES(Uw, Fvec, harms)
%   INPUTS:
%       Uw      : (2*m.Ndofs+1, 1) [Udisp; Uvel; freq] solution vector
%       Fvec    : (m.Ndofs, Nh) Each column is the complex forcing vector 
%                   for each harmonic such that
%                   Fvec(:, 1)*exp(-1j*w*t)+conj(Fvec(:, 1)*exp(-1j*w*t))
%                   is the harms(1)^th harmonic forcing
%       harms   : (Nh, 1) List of excitation harmonics
%       Ntcyc   : Number of time points (excluding final point) in each
%                   period
%       Ncyc    : Number of periods
%   OUTPUTS:
%       R       : (2*m.Ndofs, 1) Residue (final state-initial state)
%       dRdx    : (2*m.Ndofs, 2*m.Ndofs) gradient 

    Wfrc = Uw(end);
    Tmax = 2*pi/Wfrc*Ncyc;
    
    Fex = @(t) 2*(real(Fvec)*cos(harms(:)'*t) - imag(Fvec)*sin(harms(:)'*t));
    
    opts = struct('Display', 'none');
    [T, U, Ud, ~, ~, PHI] = m.HHTAMARCH(0, Tmax, Tmax/(Ntcyc*Ncyc), ...
        Uw(1:m.Ndofs), Uw(m.Ndofs+(1:m.Ndofs)), Fex, opts);

    R = [U(:, end); Ud(:, end)]-Uw(1:end-1);
    dRdx = PHI-eye(m.Ndofs*2);
end
