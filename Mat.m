clear all
close all

i = 1;

%Format Needs to have number of nodes last

cos_10 = plotCSV('cos_10.csv');

cos_100 = plotCSV('cos_100.csv');

cos_1000 = plotCSV('cos_1000.csv');


[der_1_2_10,der_1_2_100,der_1_2_1000] = plotCSV3('der_1_2', 'der_1_2');

[der_1_2_err_10,der_1_2_err_100,der_1_2_err_1000] = plotCSV3('der_1_2_err', 'der_1_2_err');


[der_1_4_10,der_1_4_100,der_1_4_1000] = plotCSV3('der_1_4', 'der_1_4');

[der_1_4_err_10,der_1_4_err_100,der_1_4_err_1000] = plotCSV3('der_1_4_err', 'der_1_4_err');


[der_1_6_10,der_1_6_100,der_1_6_1000] = plotCSV3('der_1_6', 'der_1_6');

[der_1_6_err_10,der_1_6_err_100,der_1_6_err_1000] = plotCSV3('der_1_6_err', 'der_1_6_err');



sin_10 = plotCSV('sin_10.csv', 'sin_10');

sin_100 = plotCSV('sin_100.csv', 'sin_100');

sin_1000 = plotCSV('sin_1000.csv', 'sin_1000');


[der_2_2_10,der_2_2_100,der_2_2_1000] = plotCSV3('der_2_2', 'der_2_2');

[der_2_2_err_10,der_2_2_err_100,der_2_2_err_1000] = plotCSV3('der_2_2_err', 'der_2_2_err');


[der_2_4_10,der_2_4_100,der_2_4_1000] = plotCSV3('der_2_4', 'der_2_4');

[der_2_4_err_10,der_2_4_err_100,der_2_4_err_1000] = plotCSV3('der_2_4_err', 'der_2_4_err');


[der_2_6_10,der_2_6_100,der_2_6_1000] = plotCSV3('der_2_6', 'der_2_6');

[der_2_6_err_10,der_2_6_err_100,der_2_6_err_1000] = plotCSV3('der_2_6_err', 'der_2_6_err');

%[der_2_6_10,der_2_6_100,der_2_6_1000] = plotCSV3('der_2_6', 'der_2_6');

[der_2_2_cons_err_10,der_2_2_cons_err_100,der_2_2_cons_err_1000] = plotCSV3('der_2_2_cons_err', 'der_2_2_err');
[der_2_4_cons_err_10,der_2_4_cons_err_100,der_2_4_cons_err_1000] = plotCSV3('der_2_4_cons_err', 'der_2_4_err');
[der_2_6_cons_err_10,der_2_6_cons_err_100,der_2_6_cons_err_1000] = plotCSV3('der_2_6_cons_err', 'der_2_6_err');

[der_2_2_cons_10,der_2_2_cons_100,der_2_2_cons_1000] = plotCSV3('der_2_2_cons', 'der_2_2');
[der_2_4_cons_10,der_2_4_cons_100,der_2_4_cons_1000] = plotCSV3('der_2_4_cons', 'der_2_4');
[der_2_6_cons_10,der_2_6_cons_100,der_2_6_cons_1000] = plotCSV3('der_2_6_cons', 'der_2_6');


figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Cell Domain, First Order Derivative')
plot (cos_10(:,1),cos_10(:,2))
plot (der_1_2_10(:,1),der_1_2_10(:,2))
plot (der_1_4_10(:,1),der_1_4_10(:,2))
plot (der_1_6_10(:,1),der_1_6_10(:,2))

legend('Analytical','Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Cell Domain, First Order Derivative, Error')
plot (der_1_2_err_10(:,1),abs(der_1_2_err_10(:,2)))
plot (der_1_4_err_10(:,1),abs(der_1_4_err_10(:,2)))
plot (der_1_6_err_10(:,1),abs(der_1_6_err_10(:,2)))
legend('Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Cell Domain, First Order Derivative')
plot (cos_100(:,1),cos_100(:,2))
plot (der_1_2_100(:,1),der_1_2_100(:,2))
plot (der_1_4_100(:,1),der_1_4_100(:,2))
plot (der_1_6_100(:,1),der_1_6_100(:,2))

legend('Analytical','Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Cell Domain, First Order Derivative, Error')
plot (der_1_2_err_100(:,1),abs(der_1_2_err_100(:,2)))
plot (der_1_4_err_100(:,1),abs(der_1_4_err_100(:,2)))
plot (der_1_6_err_100(:,1),abs(der_1_6_err_100(:,2)))
legend('Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Cell Domain, First Order Derivative')
plot (cos_1000(:,1),cos_1000(:,2))
plot (der_1_2_1000(:,1),der_1_2_1000(:,2))
plot (der_1_4_1000(:,1),der_1_4_1000(:,2))
plot (der_1_6_1000(:,1),der_1_6_1000(:,2))

legend('Analytical','Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Cell Domain, First Order Derivative, Error')
plot (der_1_2_err_1000(:,1),abs(der_1_2_err_1000(:,2)))
plot (der_1_4_err_1000(:,1),abs(der_1_4_err_1000(:,2)))
plot (der_1_6_err_1000(:,1),abs(der_1_6_err_1000(:,2)))
legend('Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('Second Order Accuracy Error, First Order Derivative')
plot (der_1_2_err_10(:,1),abs(der_1_2_err_10(:,2)))
plot (der_1_2_err_100(:,1),abs(der_1_2_err_100(:,2)))
plot (der_1_2_err_1000(:,1),abs(der_1_2_err_1000(:,2)))
legend('10 Nodes','100 Nodes','1000 Nodes')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('Second Order Accuracy Error, Second Order Derivative')
plot (der_2_2_err_10(:,1),abs(der_2_2_err_10(:,2)))
plot (der_2_2_err_100(:,1),abs(der_2_2_err_100(:,2)))
plot (der_2_2_err_1000(:,1),abs(der_2_2_err_1000(:,2)))
legend('10 Nodes','100 Nodes','1000 Nodes')
xlabel('X Axis')
ylabel('Y Axis')
hold off

%% Second Derivative
figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Cell Domain, Second Order Derivative, Error')
plot (der_2_2_err_10(:,1),abs(der_2_2_err_10(:,2)))
plot (der_2_4_err_10(:,1),abs(der_2_4_err_10(:,2)))
plot (der_2_6_err_10(:,1),abs(der_2_6_err_10(:,2)))
legend('Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Cell Domain, Second Order Derivative')
plot (sin_100(:,1),-sin_100(:,2))
plot (der_2_2_100(:,1),der_2_2_100(:,2))
plot (der_2_4_100(:,1),der_2_4_100(:,2))
plot (der_2_6_100(:,1),der_2_6_100(:,2))

legend('Analytical','Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Cell Domain, Second Order Derivative, Error')
plot (der_2_2_err_100(:,1),abs(der_2_2_err_100(:,2)))
plot (der_2_4_err_100(:,1),abs(der_2_4_err_100(:,2)))
plot (der_2_6_err_100(:,1),abs(der_2_6_err_100(:,2)))
legend('Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Cell Domain, Second Order Derivative')
plot (sin_1000(:,1),-sin_1000(:,2))
plot (der_2_2_1000(:,1),der_2_2_1000(:,2))
plot (der_2_4_1000(:,1),der_2_4_1000(:,2))
plot (der_2_6_1000(:,1),der_2_6_1000(:,2))

legend('Analytical','Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Cell Domain, Second Order Derivative, Error')
plot (der_2_2_err_1000(:,1),abs(der_2_2_err_1000(:,2)))
plot (der_2_4_err_1000(:,1),abs(der_2_4_err_1000(:,2)))
plot (der_2_6_err_1000(:,1),abs(der_2_6_err_1000(:,2)))
legend('Second Order','Fourth Order','Sixth Order')
xlabel('X Axis')
ylabel('Y Axis')
hold off


figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('Second Order Accuracy Error, Second Order Derivative')
plot (der_2_2_err_10(:,1),abs(der_2_2_err_10(:,2)))
plot (der_2_2_err_100(:,1),abs(der_2_2_err_100(:,2)))
plot (der_2_2_err_1000(:,1),abs(der_2_2_err_1000(:,2)))
legend('10 Nodes','100 Nodes','1000 Nodes')
xlabel('X Axis')
ylabel('Y Axis')
hold off

%% Plot each error

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Nodes, Second Order Accuracy Error, Second Order Derivative')
loglog (der_2_2_err_10(:,1),abs(der_2_2_err_10(:,2)))
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Nodes, Second Order Accuracy Error, Second Order Derivative')
plot (der_2_2_err_100(:,1),abs(der_2_2_err_100(:,2)))
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Nodes, Second Order Accuracy Error, Second Order Derivative')
semilogy (der_2_2_err_1000(:,1),abs(der_2_2_err_1000(:,2)))
xlabel('X Axis')
ylabel('Y Axis')
hold off
%% Plot cons v non conservative error
figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Nodes, Second Order Error, Conservative Vs Non-Conservative')
loglog (der_2_2_cons_err_10(:,1),abs(der_2_2_cons_err_10(:,2)))
loglog (der_2_2_err_10(:,1),abs(der_2_2_err_10(:,2)))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Nodes, Second Order Error, Conservative Vs Non-Conservative')
loglog (der_2_2_cons_err_100(:,1),abs(der_2_2_cons_err_100(:,2)))
loglog (der_2_2_err_100(:,1),abs(der_2_2_err_100(:,2)))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Nodes, Second Order Error, Conservative Vs Non-Conservative')
loglog (der_2_2_cons_err_1000(:,1),abs(der_2_2_cons_err_1000(:,2)))
loglog (der_2_2_err_1000(:,1),abs(der_2_2_err_1000(:,2)))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Nodes, Fourth Order Error, Conservative Vs Non-Conservative')
loglog (der_2_4_cons_err_10(:,1),abs(der_2_4_cons_err_10(:,2)))
loglog (der_2_4_err_10(:,1),abs(der_2_4_err_10(:,2)))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Nodes, Fourth Order Error, Conservative Vs Non-Conservative')
loglog (der_2_4_cons_err_100(:,1),abs(der_2_4_cons_err_100(:,2)))
loglog (der_2_4_err_100(:,1),abs(der_2_4_err_100(:,2)))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Nodes, Fourth Order Error, Conservative Vs Non-Conservative')
loglog (der_2_4_cons_err_1000(:,1),abs(der_2_4_cons_err_1000(:,2)))
loglog (der_2_4_err_1000(:,1),abs(der_2_4_err_1000(:,2)))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Nodes, Sixth Order Error, Conservative Vs Non-Conservative')
loglog (der_2_6_cons_err_10(:,1),abs(der_2_6_cons_err_10(:,2)))
loglog (der_2_6_err_10(:,1),abs(der_2_6_err_10(:,2)))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Nodes, Sixth Order Error, Conservative Vs Non-Conservative')
loglog (der_2_6_cons_err_100(:,1),abs(der_2_6_cons_err_100(:,2)))
loglog (der_2_6_err_100(:,1),abs(der_2_6_err_100(:,2)))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
axes('XScale', 'linear', 'YScale', 'log')
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Nodes, Sixth Order Error, Conservative Vs Non-Conservative')
loglog (der_2_6_cons_err_1000(:,1),abs(der_2_6_cons_err_1000(:,2)))
loglog (der_2_6_err_1000(:,1),abs(der_2_6_err_1000(:,2)))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

%% Plot cons v non conservative 
figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Nodes, Second Order, Conservative Vs Non-Conservative')
loglog (der_2_2_cons_10(:,1),der_2_2_cons_10(:,2))
loglog (der_2_2_10(:,1),der_2_2_10(:,2))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Nodes, Second Order, Conservative Vs Non-Conservative')
loglog (der_2_2_cons_100(:,1),der_2_2_cons_100(:,2))
loglog (der_2_2_100(:,1),der_2_2_100(:,2))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Nodes, Second Order, Conservative Vs Non-Conservative')
loglog (der_2_2_cons_1000(:,1),der_2_2_cons_1000(:,2))
loglog (der_2_2_1000(:,1),der_2_2_1000(:,2))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Nodes, Fourth Order, Conservative Vs Non-Conservative')
loglog (der_2_4_cons_10(:,1),der_2_4_cons_10(:,2))
loglog (der_2_4_10(:,1),der_2_4_10(:,2))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Nodes, Fourth Order, Conservative Vs Non-Conservative')
loglog (der_2_4_cons_100(:,1),der_2_4_cons_100(:,2))
loglog (der_2_4_100(:,1),der_2_4_100(:,2))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Nodes, Fourth Order, Conservative Vs Non-Conservative')
loglog (der_2_4_cons_1000(:,1),der_2_4_cons_1000(:,2))
loglog (der_2_4_1000(:,1),der_2_4_1000(:,2))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('10 Nodes, Sixth Order, Conservative Vs Non-Conservative')
loglog (der_2_6_cons_10(:,1),der_2_6_cons_10(:,2))
loglog (der_2_6_10(:,1),der_2_6_10(:,2))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('100 Nodes, Sixth Order, Conservative Vs Non-Conservative')
loglog (der_2_6_cons_100(:,1),der_2_6_cons_100(:,2))
loglog (der_2_6_100(:,1),der_2_6_100(:,2))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off

figure
hold on
xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi','3\pi'});
title('1000 Nodes, Sixth Order, Conservative Vs Non-Conservative')
loglog (der_2_6_cons_1000(:,1),der_2_6_cons_1000(:,2))
loglog (der_2_6_1000(:,1),der_2_6_1000(:,2))
legend('Conservative','Non Conservative')
xlabel('X Axis')
ylabel('Y Axis')
hold off
function [A,B,C] = plotCSV3(file,name)

%     figure('Name', name);
% 
%     hold on

    A = readmatrix(append(file,'_10.csv'));
    B = readmatrix(append(file,'_100.csv'));
    C = readmatrix(append(file,'_1000.csv'));
%     n = size(A);
%     n2 = size(B);
%     n3 = size(C);
%     plot (A(3:n-3,1),A(3:n-3,2));
%     plot (B(3:n2-3,1),B(3:n2-3,2));
%     plot (C(3:n3-3,1),C(3:n3-3,2));

%    hold off


end


function A = plotCSV(file, name)

    %figure('Name', name);

    A = readmatrix(file);
    
    %plot (A(:,1),A(:,2));


end