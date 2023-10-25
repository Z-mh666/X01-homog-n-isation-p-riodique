% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions periodiques
% sur le maillage nom_maillage.msh
%
% | -div(A grad u ) u + u= f,   dans \Omega
% |         u periodique,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------

nom_maillage = 'geomCarre_per2.msh';  %(0,1)^2
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre
eta = 1e-4;
e = 1;       %valeur d'epsilon

% boucle sur les triangles
% ------------------------
for i=1:Nbtri
  % Coordonnees des sommets du triangles
  % A COMPLETER
  k = Numtri(i,1);
  l = Numtri(i,2);
  m = Numtri(i,3);
  
  S1 = Coorneu(k,:);
  S2 = Coorneu(l,:);
  S3 = Coorneu(m,:);
  % calcul des matrices elementaires du triangle l 
  
  Kel = matK2_elem(S1,S2,S3,e);
           
  Mel = matM_elem(S1,S2,S3);
  %Mel = matM2_elem(S1,S2,S3,eta); %penalisation
    
  % On fait l'assemmblage de la matrice globale et du second membre
  % A COMPLETER
  MM(k,k) += Mel(1,1);
  MM(l,l) += Mel(2,2);
  MM(m,m) += Mel(3,3);
  MM(k,l) += Mel(1,2);
  MM(k,m) += Mel(1,3);
  MM(l,m) += Mel(2,3);
  MM(l,k) += Mel(2,1);
  MM(m,k) += Mel(3,1);
  MM(m,l) += Mel(3,2);
  
  KK(k,k) += Kel(1,1);
  KK(l,l) += Kel(2,2);
  KK(m,m) += Kel(3,3);
  KK(k,l) += Kel(1,2);
  KK(k,m) += Kel(1,3);
  KK(l,m) += Kel(2,3);
  KK(l,k) += Kel(2,1);
  KK(m,k) += Kel(3,1);
  KK(m,l) += Kel(3,2);

end % for l

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
 
% pour la methode de penalisation
%
FF = Coorneu(:,1);
FF2 = Coorneu(:,2);
LL = -KK*FF;
LL2 = -KK*FF2;
%


% Projection sur l espace V_p
% ——————————————————�?
% matrice de projection 
PP = zeros(Nbpt,Nbpt);
PP(1,1) = 1;
PP(2,1) = 1;
PP(3,1) = 1;
PP(4,1) = 1;
num = 2;
for i=5:(Nbaretes/2+2)
  PP(i,num) = 1;
  num += 1;
endfor
num = 2;
for i=(Nbaretes/2+Nbaretes/4+1):-1:(Nbaretes/2+3)
  PP(i,num) = 1;
  num += 1;
endfor
for i=Nbaretes:-1:(Nbaretes/2+Nbaretes/4+2)
  PP(i,num) = 1;
  num += 1;
endfor
num = Nbaretes+1;
for i=(Nbaretes+1):Nbpt
  PP(i,num) = 1;
  num += 1;
endfor

AA = eta*MM+KK;
AAp = PP'*AA*PP;
LLp = PP'*LL;

% inversion
% ----------
UUp = AAp\LLp;

% Expression de la solution dans toute la base
% ——————�?
UU = PP*UUp;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Periodique - %s', nom_maillage));

%calcul du tenseur homogeneise et visualisation
%
LLp2 = PP'*LL2;
UUp2 = AAp\LLp2;
UU2 = PP*UUp2;
Aeff = tenseur_homo(KK,UU,UU2,Coorneu(:,1),Coorneu(:,2));
affiche(UU2, Numtri, Coorneu, sprintf('Periodique - %s', nom_maillage));
%

validation = 'non';
% validation
% ----------
if strcmp(validation,'non')
UU_exact = Coorneu(:,1)*0;
% Calcul de l erreur L2
% A COMPLETER
errL2 = sqrt(UU'*MM*UU+UU_exact'*MM*UU_exact-2*UU'*MM*UU_exact)

errH1 = sqrt(UU'*KK*UU+UU_exact'*KK*UU_exact-2*UU'*KK*UU_exact)

% attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

