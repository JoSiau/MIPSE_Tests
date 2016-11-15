clc;
clear all;
close all;

% ~\src\formulations\g2elab\mipse\method\formulation\test\testhmatrix\Tests
% We have to go up to ~

T=load('.\..\..\..\..\..\..\..\..\..\Tps.out');

L = (T);
ord = polyfit(L(:,1),L(:,2),1);

figure(1);
loglog(T(:,1),T(:,2),'-bo',T(:,1),T(:,1).*log10(T(:,1)).*10^(-2),'-ko');
xlabel('Nbr de points');ylabel('Temps en secondes');title('loglog des temps de calculs pour assembler une Hmatrice');
legend('Java',['X.log_{10}(X)*1e' num2str(ord(2)) ]);
