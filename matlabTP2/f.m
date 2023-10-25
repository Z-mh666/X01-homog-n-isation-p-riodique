function val = f(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER

%Pour la premi¨¨re partie

%val = (2*pi^2)*sin(pi*x).*sin(pi*y);  %pour A = idendite;

%val = (3*pi^2)*sin(pi*x).*sin(pi*y);  %pour A = [1 0;0 2]; 

%pour A = [2+sin(2*pi*x) 0; 0 4]
%val = (1+6*pi^2)*sin(pi*x).*sin(pi*y)-2*pi^2*cos(2*pi*x).*cos(pi*x).*sin(pi*y)...
%      +pi^2*sin(2*pi*x).*sin(pi*x).*sin(pi*y);

%pour A = [2+sin(2*pi*x) 0; 0 4+sin(2*pi*x)]
%val = (1+6*pi^2)*sin(pi*x).*sin(pi*y)-2*pi^2*cos(2*pi*x).*cos(pi*x).*sin(pi*y)...
%      +2*pi^2*sin(2*pi*x).*sin(pi*x).*sin(pi*y);

%Pour le probl¨¨me homog¨¦n¨¦is¨¦

%val = (3*pi^2)*sin(pi*x).*sin(pi*y);  %pour Aeff = [1 0;0 2]

%val = (sqrt(3)+4)*pi^2*sin(pi*x).*sin(pi*y);  %pour Aeff = [sqrt(3) 0;0 4]

%pour Aeff = [4*sqrt(3) 0;0 2*sqrt(15)]
val = (4*sqrt(3)+2*sqrt(15))*pi^2*sin(pi*x).*sin(pi*y);  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
