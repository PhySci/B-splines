% Загружаем тестовые данные
newData1 = load('-mat', 'matlab.mat');
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
clear i newData1 vars

xData=date_interes;
yData=price_interes;
pweights=weight_interes;

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.003355272021954153; % Сплайн параметр
opts.Weights = pweights;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure(1);
h = plot( fitresult, xData, yData );
legend( h, 'price_interes vs. date_interes with weight_interes', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel date_interes
ylabel price_interes
grid on

sp2 = spaps(xData, yData, -opts.SmoothingParam, pweights, 2);

figure(2)
subplot(211);
    plot(xData, yData, 'rx', xData,fnval(xData,sp2),xData,fitresult(xData));

subplot(212);
    plot(xData, fnval(xData,sp2)-fitresult(xData));

%M = [xData yData pweights fitresult(xData)];
%dlmwrite('data.csv',M,'precision',7);