clc;
clear all;
close all;

% ~\src\formulations\g2elab\mipse\method\formulation\test\testhmatrix\Tests
% We have to go up to ~

T=load('.\..\..\..\..\..\..\..\..\..\Tps.out');

L = log10(T);

ord = polyfit(L(:,1),L(:,2),1);
figure(1);hold on;
% plot(T(:,1),T(:,2),'-bo');
plot(T(:,1),T(:,1).*log10(T(:,1)),'-ko',T(:,1),ord(2)+ord(1)*T(:,1),'-rs');
legend('X.log10(X)',[num2str(ord(2)) '+' num2str(ord(1)) '*X']);
title(['Ordre de la methode: ' num2str(ord(1))]);