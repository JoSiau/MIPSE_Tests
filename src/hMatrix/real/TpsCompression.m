clc;
clear all;
close all;

% ~\src\formulations\g2elab\mipse\method\formulation\test\testhmatrix\Tests
% We have to go up to ~

T=load('.\..\..\..\..\..\..\..\..\..\Tps.out');
tInd = T(:,1); tACA = T(:,2); tSVD = T(:,3); tTotal = tACA+tSVD;
L = log10([T(:,1),tTotal]);
ord = polyfit(L(:,1),L(:,2),1);

figure(1);
loglog(tInd,tACA,'-bo',tInd,tSVD,'-ro',tInd,tInd.*log10(tInd).*10^-5,':k',tInd,tTotal,'-mo');
xlabel('Nombre de points');ylabel('Temps en secondes');title('loglog des temps de calculs pour les compressions');
legend('ACA','Rk-SVD',['X.log_{10}(X)*1e' num2str(-5) ],'Compression complète');
%% 
% figure(2);hold on;
% % plot(T(:,1),T(:,2),'-bo');
% plot(T(:,1),T(:,1).*log10(T(:,1)),'-ko',T(:,1),ord(2)+ord(1)*T(:,1),'-rs');
% legend('X.log10(X)',[num2str(ord(2)) '+' num2str(ord(1)) '*X']);
% title(['Ordre de la methode: ' num2str(ord(1))]);