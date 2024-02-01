clear all
close all

i = 1;

%Format Needs to have number of nodes last

plotCSV('cos_10.csv', 'cos_10');

plotCSV('cos_100.csv', 'cos_100');

plotCSV('cos_1000.csv', 'cos_1000');


plotCSV3('der_1_2', 'der_1_2');

plotCSV3('der_1_2_err', 'der_1_2_err');


plotCSV3('der_1_4', 'der_1_4');

plotCSV3('der_1_4_err', 'der_1_4_err');


plotCSV3('der_1_6', 'der_1_6');

plotCSV3('der_1_6_err', 'der_1_6_err');



plotCSV('sin_10.csv', 'sin_10');

plotCSV('sin_100.csv', 'sin_100');

plotCSV('sin_1000.csv', 'sin_1000');


plotCSV3('der_2_2', 'der_2_2');

plotCSV3('der_2_2_err', 'der_2_2_err');


plotCSV3('der_2_4_10', 'der_2_4');

plotCSV3('der_2_4_err', 'der_2_4_err');


plotCSV3('der_2_6', 'der_2_6');

plotCSV3('der_2_6_err', 'der_2_6_err');

