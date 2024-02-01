function plotCSV(file, name)

    figure('Name', name);

    A = readmatrix(file);

    plot (A(:,1),A(:,2));


end