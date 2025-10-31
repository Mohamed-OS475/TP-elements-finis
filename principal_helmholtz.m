clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINCIPAL_HELMHOLTZ;                                                        %
%                                                                             %
% Une routine pour la mise en oeuvre des EF P1 Lagrange pour resoudre         %
% l'equation de Helmholtz stationnaire suivante                               %
%                                                                             %
% -\rho \omega^2  P - div(\mu \grad P) = S,   dans \Omega=\Omega_1 U \Omega_2 %
%                                                                             %
% ou S est la source, omega la pulsation,                                     %
% \rho constante, et \mu variable :                                           %
%                       \mu = | \mu_1,  dans \Omega_1                         %
%                             | \mu_2,  dans \Omega_2.                        %
%                                                                             %
% A l'equation de Helmholtz s'ajoute l'une des conditions suivantes :         %
%     1) une condition de Dirichlet :                                         %
%                         P = P_\Gamma,   sur le bord,                        %
%        ou P_\Gamma represente la pression exterieure ;                      %
%                                                                             %
%     2) une condition de Fourier :                                           %
%                        dP/dn + beta P = 0 sur le bord,                      %
%        ou beta est le coefficient de Fourier ;                              %
%                                                                             %
%     3) des conditions mixtes de Dirichlet-Fourier.                          %
%                                                                             %
% La routine calcule et affiche la solution, qui peut ensuite ??tre compar??e   %
% ?? une solution de r??f??rence.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------- %
% Donnees du probleme %
% ------------------- %
nom_maillage = 'maillages/geom-2-couches.msh'; %'maillages/geomCarre.msh'; %'maillages/geomPacman.msh';
validation = 'non';           % Pour valider 'oui' ou 'non' les calculs

% Choix des conditions aux limites
CL = 'Fourier'; % 'Dirichlet', 'Neumann' ou 'Fourier'ou 'Mixtes'

% Donnees du probleme
rho = 1.0;        % Coefficient du probleme
omega = 5.0;      % Frequence
PP_Gamma = 1.0;   % Pression exterieure
alpha = -rho * omega^2;

if (strcmp(CL, 'Fourier') || strcmp(CL, 'fourier') ||...
    strcmp(CL, 'Mixtes')  || strcmp(CL, 'mixtes'))
    beta = -1i * omega;   % Condition de Fourier
    %beta = 0;           % Condition de Neumann
end

if strcmp(validation,'oui')
    mu_1 = 1;
    mu_2 = 1;
    PP_Gamma = 0;
end

% -------------------------------- %
% Lecture du maillage et affichage %
% -------------------------------- %
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, ...
 Reftri, Nbaretes, Numaretes, Refaretes] = lecture_msh(nom_maillage);

% ---------------------- %
% Calcul des matrices EF %
% ---------------------- %
% Declarations %
% ------------ %
MM = sparse(Nbpt, Nbpt); % matrice de masse
KK = sparse(Nbpt, Nbpt); % matrice de rigidite
LL = zeros(Nbpt, 1);     % vecteur second membre

% Boucle sur les triangles %
% ------------------------ %
for l = 1:Nbtri
   % Coordonnees des sommets du triangles
  % A COMPLETER
  S1=Coorneu(Numtri(l, 1),:);
  S2=Coorneu(Numtri(l, 2),:);
  S3=Coorneu(Numtri(l, 3),:);
  % calcul des matrices elementaires du triangle l 
  
       % Matrice de masse elementaire    
   Mel=matM_elem(S1, S2, S3);
   
    % Matrice de rigidite elementaire
    % LA ROUTINE matK_mu_elem.m DOIT ETRE MODIFIEE
   Kel=matK_mu_elem(S1, S2, S3,Reftri(l));
  
 

    % On fait l'assemblage de la matrice globale
    % A COMPLETER
    for i=1:3
      I = Numtri(l, i);
      for j=1:3
          J = Numtri(l, j);
          MM(I,J) = MM(I,J) + Mel(i,j);
          KK(I,J) = KK(I,J) + Kel(i,j);
      end % for j
  end % for i
    
end % for l

% Matrice de masse surfacique pour les conditions de Fourier %
% ---------------------------------------------------------- %
if (strcmp(CL, 'Fourier') || strcmp(CL, 'fourier') ||...
    strcmp(CL, 'Mixtes')  || strcmp(CL, 'mixtes'))
    % On trouve les aretes sur lesquelles la condition de Fourier est imposee
    idFourier = find(Refaretes == 2);   % Indices des aretes pour Fourier
    NbFourier = length(idFourier);      % Nombre d'aretes sur lesquelles Fourier est imposee
    SS = sparse(Nbpt, Nbpt);            % Matrice de masse surfacique

    % Assemblage sur les aretes de reference 2
    for l = 1:NbFourier
        % Matrice de masse surfacique elementaire
        % LA ROUTINE matS_elem.m DOIT ETRE COMPLETEE
        [Sel] = matS_elem(Coorneu(Numaretes(idFourier(l),1),:),...
                          Coorneu(Numaretes(idFourier(l),2), :), Refaretes(idFourier(l)));

        % On fait l'assemblage de la matrice de masse surfacique globale
        % A COMPLETER
        for i = 1:2
          I = Numaretes(l,i);
          for j = 1:2
              J = Numaretes(l,j);
              SS(I,J) = SS(I,J)+ Sel(i,j);
          end 
        end 
        
    end % for l
end % for if

% ---------- %
% Matrice EF %
% ---------- %
AA = alpha*MM + KK;

if (strcmp(CL, 'Fourier') || strcmp(CL, 'fourier') ||...
    strcmp(CL, 'Mixtes')  || strcmp(CL, 'mixtes'))
  % On rajoute la contribution de la matrice de masse surfacique
  % dans le cas des condtions de 'Fourier' ou 'Mixtes'
  % COMPLETER
  AA = AA + beta*SS;
end

% ------------------------- %
% Calcul du second membre F %
% ------------------------- %
% A COMPLETER EN UTILISANT LA ROUTINE f.m
%  /!\ ATTENTION : f prend aussi omega en argument /!\

FF = f(Coorneu(:,1), Coorneu(:,2), omega);
LL = MM*FF;

% ------------------------------- %
% Pseudo-elimination et inversion %
% ------------------------------- %
if (strcmp(CL, 'Dirichlet') || strcmp(CL, 'dirichlet') ||...
    strcmp(CL, 'Mixtes')    || strcmp(CL, 'mixtes'))
  % Conditions de Dirichlet ou conditions mixtes
  % tilde_AA ET tilde_LL SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
  % APRES PSEUDO_ELIMINATION
  % ECRIRE LA ROUTINE elimine.m ET INSERER L'APPEL A CETTE ROUTINE
  % A UN ENDROIT APPROPRIE
  [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu);
  UU = tilde_AA \ tilde_LL;
  PP = UU + PP_Gamma;
else
  % Conditions de Neumann ou de Fourier
  PP = AA \ LL;
end

if strcmp(validation, 'oui')
  % ------------------------------------------------------------------ %
  % Pour tracer la norme de l'inverse de A(omega) en fonction de omega %
  % ------------------------------------------------------------------ %
  % A COMPLETER
end

% ---------- %
% Validation %
% ---------- %
if strcmp(validation, 'oui')
  % Solution de reference

  %u(x, y) = sin(3πx) sin(4πy)
  PP_exact = PP_Gamma + sin(3*pi*Coorneu(:,1)) .* sin(4*pi*Coorneu(:,2));

  %f(x, y) = Asin(ωr)



  % Calcul de l'erreur L2
    diff = PP - PP_exact;

    % Calcul de l erreur L2
    norm_L2 = sqrt(diff' * MM * diff);
  % Calcul de l erreur H1
  semi_norm_H1 = sqrt(diff' * KK * diff);

  %Affichage erreurs
fprintf('erruer L2 = %e\n', norm_L2);
fprintf('erreur semi_H1 =  %e\n', semi_norm_H1);
  % attention de bien changer le terme source (dans FF)
end

% % ------------- %
% % Visualisation %
% % ------------- %
if (strcmp(CL, 'Dirichlet') || strcmp(CL, 'dirichlet'))
  % La solution est reelle : on l'affiche donc dans une seule figure
  figure;
  affiche(real(PP), Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));
  set(gca, 'DataAspectRatio',[1 1 1]);
end

if (strcmp(CL, 'Fourier') || strcmp(CL, 'fourier') ||...
    strcmp(CL, 'Mixtes')  || strcmp(CL, 'mixtes'))
  % La solution est COMPLEXE : on affiche donc ses parties reelle et
  % imaginaire dans des figures
  figure;
  affiche(real(PP), Numtri, Coorneu, sprintf('%s - Real(P) - %s', CL, nom_maillage));
  set(gca, 'DataAspectRatio',[1 1 1]);

  figure;
  affiche(imag(PP), Numtri, Coorneu, sprintf('%s - Imag(P) - %s', CL, nom_maillage));
  set(gca, 'DataAspectRatio',[1 1 1]);
end

% Affichage de la solution exacte
if (strcmp(validation, 'oui'))
    figure;
    affiche(real(PP_exact), Numtri, Coorneu, sprintf('Solution exacte - %s - %s', CL, nom_maillage));
    set(gca, 'DataAspectRatio',[1 1 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%25
