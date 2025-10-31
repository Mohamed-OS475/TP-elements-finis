function val = mu_1(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MU_1 :
% Evaluation de la fonction mu_1.
%
% SYNOPSIS val = mu_1(x,y)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choix de l'expression du coefficient mu_1
% A CHANGER EVENTUELLEMENT POUR LES SIMULATIONS
%choix_mu1 = 'constante';   % 'constante' ou 'variable'
choix_mu1 = 'variable';


if strcmp(choix_mu1, 'constante')
  % Coefficient constant
  val = 1;
end

if strcmp(choix_mu1,'variable')
% %Interface penchee / mu1 variable et stratifie le long de l'interface
% Source dans Omega2
% % 
p=5;
mu2=1;
C= 2*mu2;
A= 4*C/5.;
alpha = (2*p+1)*pi/2;
a=1.;
val = C +  A*cos(alpha*(y-a*x));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%25
