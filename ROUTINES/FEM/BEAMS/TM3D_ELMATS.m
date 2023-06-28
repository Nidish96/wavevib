function [Me, Ke] = TM3D_ELMATS(Le, pars)
%TM3D_ELMATS returns the element mass and stiffness matrices for a
%linear 3D Timoshenko Beam element with DoF ordering:
%  [u, v, w, thx, thy, thz]
%
%  USAGE :
%    [Me, Ke] = TM3D_ELMATS(Le, pars);
%  INPUTS :
%    Le 	: Element Length
%    pars	: Structure containing,
%  	E, G, rho, A, Iy, Iz, k1, k2, k3
%  OUTPUTS :
%    Me 	: 12x12 mass matrix
%    Ke 	: 12x12 stiffness matrix

  rhoA  = pars.rho*pars.A;
  rhoIz = pars.rho*pars.Iz;
  rhoIy = pars.rho*pars.Iy;

  EA    = pars.E*pars.A;
  EIz   = pars.E*pars.Iz;
  EIy   = pars.E*pars.Iy;

  k1    = pars.k1;
  k2    = pars.k2;
  GA    = pars.G*pars.A;
  GIy   = pars.k3*pars.G*pars.Iy;
  GIz   = pars.k3*pars.G*pars.Iz;

  %% Mass Matrix
  Me = Le./6*[rhoA*2, 0, 0, 0, 0, 0, rhoA, 0, 0, 0, 0, 0;
              0, rhoA*2, 0, 0, 0, 0, 0, rhoA, 0, 0, 0, 0;
              0, 0, rhoA*2, 0, 0, 0, 0, 0, rhoA, 0, 0, 0;
              0, 0, 0, (rhoIy+rhoIz)*2, 0, 0, 0, 0, 0, (rhoIy+rhoIz), 0, 0;
              0, 0, 0, 0, rhoIy*2, 0, 0, 0, 0, 0, rhoIy, 0;
              0, 0, 0, 0, 0, rhoIz*2, 0, 0, 0, 0, 0, rhoIz;
              rhoA, 0, 0, 0, 0, 0, rhoA*2, 0, 0, 0, 0, 0;
              0, rhoA, 0, 0, 0, 0, 0, rhoA*2, 0, 0, 0, 0;
              0, 0, rhoA, 0, 0, 0, 0, 0, rhoA*2, 0, 0, 0;
              0, 0, 0, (rhoIy+rhoIz), 0, 0, 0, 0, 0, (rhoIy+rhoIz)*2, 0, 0;
              0, 0, 0, 0, rhoIy, 0, 0, 0, 0, 0, rhoIy*2, 0;
              0, 0, 0, 0, 0, rhoIz, 0, 0, 0, 0, 0, rhoIz*2];

  %% Stiffness Matrix
  % Ke = [EA/Le, 0, 0, 0, 0, 0, -EA/Le, 0, 0, 0, 0, 0;
  %       0, k1*GA/Le, 0, 0, 0, -k1*GA/2, 0, -k1*GA/Le, 0, 0, 0, -k1*GA/2;
  %       0, 0, k2*GA/Le, 0, k2*GA/2, 0, 0, 0, -k2*GA/Le, 0, k2*GA/2, 0;
  %       0, 0, 0, (GIz+GIy)/Le, 0, 0, 0, 0, 0, -(GIz+GIy)/Le, 0, 0;
  %       0, 0, k2*GA/2, 0, (k2*GA*Le^2+4*EIy)/(4*Le), 0, 0, 0, -k2*GA/2, 0, (k2*GA*Le^2-4*EIy)/(4*Le), 0;
  %       0, -k1*GA/2, 0, 0, 0, (k1*GA*Le^2+4*EIz)/(4*Le), 0, k1*GA/2, 0, 0, 0, (k1*GA*Le^2-4*EIz)/(4*Le);
  %       -EA/Le, 0, 0, 0, 0, 0, EA/Le, 0, 0, 0, 0, 0;
  %       0, -k1*GA/Le, 0, 0, 0, k1*GA/2, 0, k1*GA/Le, 0, 0, 0, k1*GA/2;
  %       0, 0, -k2*GA/Le, 0, -k2*GA/2, 0, 0, 0, k2*GA/Le, 0, -k2*GA/2, 0;
  %       0, 0, 0, -(GIz+GIy)/Le, 0, 0, 0, 0, 0, (GIz+GIy)/Le, 0, 0;
  %       0, 0, k2*GA/2, 0, (k2*GA*Le^2-4*EIy)/(4*Le), 0, 0, 0, -k2*GA/2, 0, (k2*GA*Le^2+4*EIy)/(4*Le), 0;
  %       0, -k1*GA/2, 0, 0, 0, (k1*GA*Le^2-4*EIz)/(4*Le), 0, k1*GA/2, ...
  %       0, 0, 0, (k1*GA*Le^2+4*EIz)/(4*Le)];
  
  %% New
  Ke = [EA/Le, 0, 0, 0, 0, 0, -EA/Le, 0, 0, 0, 0, 0;
       0, k1*GA/Le, 0, 0, 0, k1*GA/2, 0, -k1*GA/Le, 0, 0, 0, k1*GA/2;
       0, 0, k2*GA/Le, 0, -k2*GA/2, 0, 0, 0, -k2*GA/Le, 0, -k2*GA/2, 0;
       0, 0, 0, (GIz+GIy)/Le, 0, 0, 0, 0, 0, -(GIz+GIy)/Le, 0, 0;
       0, 0, -k2*GA/2, 0, (k2*GA*Le^2+4*EIy)/(4*Le), 0, 0, 0, k2*GA/2, 0, (k2*GA*Le^2-4*EIy)/(4*Le), 0;
       0, k1*GA/2, 0, 0, 0, (k1*GA*Le^2+4*EIz)/(4*Le), 0, -k1*GA/2, 0, 0, 0, (k1*GA*Le^2-4*EIz)/(4*Le);
       -EA/Le, 0, 0, 0, 0, 0, EA/Le, 0, 0, 0, 0, 0;
       0, -k1*GA/Le, 0, 0, 0, -k1*GA/2, 0, k1*GA/Le, 0, 0, 0, -k1*GA/2;
       0, 0, -k2*GA/Le, 0, k2*GA/2, 0, 0, 0, k2*GA/Le, 0, k2*GA/2, 0;
       0, 0, 0, -(GIz+GIy)/Le, 0, 0, 0, 0, 0, (GIz+GIy)/Le, 0, 0;
       0, 0, -k2*GA/2, 0, (k2*GA*Le^2-4*EIy)/(4*Le), 0, 0, 0, k2*GA/2, 0, (k2*GA*Le^2+4*EIy)/(4*Le), 0;
       0, k1*GA/2, 0, 0, 0, (k1*GA*Le^2-4*EIz)/(4*Le), 0, -k1*GA/2, 0, 0, 0, (k1*GA*Le^2+4*EIz)/(4*Le)];
end
