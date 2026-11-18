clear; clc; close all;

%% Постоянные
m = 2;
k = 17;

%% Палитра
maroon     = [0.4627, 0.3294, 0.2314]; % коричневый
red        = [0.7451, 0.5843, 0.4353]; % светло-коричневый
cream      = [0.8000, 0.8000, 0.8000]; % серая сетка
background = [1.0000, 1.0000, 1.0000]; % белый фон

%% Параметры построения
p_values = [-1, -28/323, 0, 1];
t = linspace(0, 2, 801);
styles = {'-', '--', '-.', ':'};
colors = [maroon; red; maroon; red];

% –––––

figure(1);
set(gcf, 'Color', background, 'Position', [80 100 1300 500]);

%% Левый график: растущие решения при разных p
ax_left = subplot(1, 2, 1);
hold on;
norm_min = Inf;
norm_max = 0;

for j = 1:length(p_values)
    p = p_values(j);
    A = [-3 k 0; m -2 p; m+k 0 -1];
    [V, D] = eig(A);
    values = diag(D);
    fprintf('\np = %.10f; собственные числа:\n', p);
    disp(values);

    [~, idx] = max(real(values));
    lambda = values(idx);
    v = V(:, idx);

    X = real(v * exp(lambda*t));
    X = X / norm(real(v));
    norms = sqrt(sum(X.^2, 1));

    semilogy(t, norms, styles{j}, ...
        'Color', colors(j,:), 'LineWidth', 2);

    norm_min = min(norm_min, min(norms));
    norm_max = max(norm_max, max(norms));
end
hold off;

% Такие же отступы по логарифмической оси
log_span = log(norm_max/norm_min);
set(gca, 'YScale', 'log');
xlim([0 2]);
ylim([norm_min*exp(-0.05*log_span), norm_max*exp(0.05*log_span)]);
grid on;
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12, 'Color', 'black');
ylabel('$\|X(t)\|$', 'Interpreter', 'latex', 'FontSize', 12, 'Color', 'black');
title('Растущие решения при разных p', 'FontSize', 13, ...
    'FontWeight', 'normal', 'Color', 'black');
legend({'$p=-1$', '$p=-28/323$', '$p=0$', '$p=1$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 10, ...
    'Color', background, 'EdgeColor', cream, 'TextColor', 'black');
set(gca, 'Color', background, 'GridColor', cream, 'GridAlpha', 1, ...
    'GridLineStyle', '-', 'LineWidth', 0.7, ...
    'XColor', 'black', 'YColor', 'black', ...
    'XTick', 0:0.25:2, 'YTick', 10.^(0:5), ...
    'YMinorTick', 'off', 'YMinorGrid', 'off', ...
    'FontSize', 10, 'Box', 'off', 'TickLength', [0 0]);

%% Правый график: разные начальные условия при p = 0
lambda_plus = (-5+sqrt(137))/2;
growth = exp(lambda_plus*t);
decay = exp(-t);
zero = zeros(size(t));

ax_right = subplot(1, 2, 2);
plot(t, growth, '-', 'Color', maroon, 'LineWidth', 2);
hold on;
plot(t, decay, '--', 'Color', red, 'LineWidth', 2);
plot(t, zero, ':', 'Color', maroon, 'LineWidth', 2);
hold off;

grid on;
xlim([0 2]);
ylim([-0.05*max(growth), 1.05*max(growth)]);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12, 'Color', 'black');
ylabel('$\|X(t)\|$', 'Interpreter', 'latex', 'FontSize', 12, 'Color', 'black');
title('Разные начальные условия при p = 0', 'FontSize', 13, ...
    'FontWeight', 'normal', 'Color', 'black');
legend('Растущее: e^{\lambda_+ t}', 'Затухающее: e^{-t}', ...
    'Нулевое решение', 'Location', 'southeast', 'FontSize', 9, ...
    'Color', background, 'EdgeColor', cream, 'TextColor', 'black');
set(gca, 'Color', background, 'GridColor', cream, 'GridAlpha', 1, ...
    'GridLineStyle', '-', 'LineWidth', 0.7, ...
    'XColor', 'black', 'YColor', 'black', ...
    'XTick', 0:0.25:2, 'YTick', 0:100:800, ...
    'FontSize', 10, 'Box', 'off', 'TickLength', [0 0]);

%% Расположение двух основных графиков — как в Jupyter
set(ax_left, 'Units', 'normalized', ...
    'Position', [0.06955 0.15011 0.40226 0.72467]);
set(ax_right, 'Units', 'normalized', ...
    'Position', [0.56460 0.15011 0.40226 0.72467]);

%% Увеличенный участок внутри правого графика
pos = get(ax_right, 'Position');
zoom_pos = [pos(1)+0.09*pos(3), pos(2)+0.49*pos(4), ...
            0.46*pos(3), 0.39*pos(4)];

ax_zoom = axes('Position', zoom_pos);
plot(t, decay, '--', 'Color', red, 'LineWidth', 2);
hold on;
plot(t, zero, ':', 'Color', maroon, 'LineWidth', 2);
hold off;

grid on;
xlim([0 2]);
ylim([-0.08 1.08]);
title('Увеличение около нуля', 'FontSize', 9, ...
    'FontWeight', 'normal', 'Color', 'black');
set(gca, 'Color', background, 'GridColor', cream, 'GridAlpha', 1, ...
    'GridLineStyle', '-', 'LineWidth', 0.6, ...
    'XColor', 'black', 'YColor', 'black', ...
    'XTick', [0 1 2], 'YTick', [0 0.5 1], ...
    'FontSize', 8, 'Box', 'on', 'TickLength', [0 0]);

% –––––

%% Запуск файла с Octave (раскомментировать)
%if exist('OCTAVE_VERSION', 'builtin')
%    disp('Графики построены.');
%    pause;
%end