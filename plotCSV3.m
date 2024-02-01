function plotCSV3(file, name)

    figure('Name', name);

    hold on

    A = readmatrix(append(file,'_10.csv'));
    B = readmatrix(append(file,'_100.csv'));
    C = readmatrix(append(file,'_1000.csv'));


    plot (A(:,1),A(:,2));
    plot (B(:,1),B(:,2));
    plot (C(:,1),C(:,2));

    hold off


end