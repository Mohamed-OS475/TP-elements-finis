function val = f(x, y, omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x, y, omega)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%         omega : la pulsation
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choix de l'expression de la source
% A CHANGER EVENTUELLEMENT POUR LES SIMULATIONS

%choix_source = 'sin-sin' 
choix_source = 'Radiale-Exo1-et-2' 
%choix_source = 'Radiale-Exo3' 
%choix_source = 'validation' 
%choix_source = 'constante' 



if strcmp(choix_source,'validation')
  % A COMPLETER POUR LA VALIDATION :
  
  val = sin(3*pi*x).*sin(4*pi*y);
end

if strcmp(choix_source,'sin-sin')
  % Polynome trigonometrique
  val = sin(3*pi*x).*sin(4*pi*y);
end

if strcmp(choix_source,'constante')
  % Terme source constant
  % A COMPLETER 
  val = ones(size(x));
end

if strcmp(choix_source,'Radiale-Exo1-et-2')
  % Source radiale et sinusoidale
  A = 5;           % Amplitude de la source
  x0 = 0.6;        % Abscisse du centre 
  y0 = 2.5;        % Ordonnee du centre 
  rr = sqrt((x - x0).^2 + (y - y0).^2);   % Variable radiale
  val = A * sin(omega*rr);                % Expression de la source
end


if strcmp(choix_source,'Radiale-Exo3')
A= 5;
% Pour Pacman : 
x0=-0.5;
y0=0.;
rr=sqrt((x-x0).^2+(y-y0).^2);
val = A*cos(omega*rr);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%25
