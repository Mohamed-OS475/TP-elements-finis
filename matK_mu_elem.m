function [Kel] = matK_mu_elem(S1, S2, S3, Reftri)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matK_mu_elem :
% Calcule la matrice de raideur elementaire en P1 lagrange avec un
% coefficient variable.
%
% SYNOPSIS [Kel] = matK_mu_elem(S1, S2, S3, Reftri)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%         Reftri : reference du triangle (1 si appartient a Omega_1
%                     et 2 si dans Omega_2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exact (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par
%          les coordonnees barycentriques
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps)
  error('l aire d un triangle est nulle!!!');
end;


% calcul de la matrice de raideur
% -------------------------------
% Attention a prendre en compte le coefficient variable mu et le sous-domaine
% Omega
Kel = zeros(3,3);

%for i=1:3
  %for j=1:3
	% A COMPLETER
    %Kel(i,j) = 
  %end; % j
%end; % i

w0 = 1/6;
somme = (mu_1((4*x1 + x2 + x3) / 6, (4*y1 + y2 + y3) / 6) + mu_1((x1 + 4*x2 + x3) / 6, (y1 + 4*y2 + y3) / 6) + mu_1((x1 + x2 + 4*x3) / 6, (y1 + y2 + 4*y3) / 6)) * (Reftri == 1) +...
    (mu_2((4*x1 + x2 + x3) / 6, (4*y1 + y2 + y3) / 6) + mu_2((x1 + 4*x2 + x3) / 6, (y1 + 4*y2 + y3) / 6) + mu_2((x1 + x2 + 4*x3) / 6, (y1 + y2 + 4*y3) / 6)) * (Reftri == 2);
Kel = w0 * somme * (norm * norm') / abs(D);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%25
