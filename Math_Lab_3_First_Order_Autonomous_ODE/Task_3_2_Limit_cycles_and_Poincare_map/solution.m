clear; clc; close all;

%% Постоянные
m = 2;
k = 17;

%% Палитра
brown       = [0.4627, 0.3294, 0.2314];
light_brown = [0.7451, 0.5843, 0.4353];
beige       = [0.8549, 0.7765, 0.6667];
gray        = [0.8000, 0.8000, 0.8000];
background  = [1.0000, 1.0000, 1.0000];

% –––––

%% Нужный цикл и уравнение для отклонения
n = m + 1;
R = sqrt(n * pi);
T = 2 * pi / k;

fprime_R = 2 * m * n * pi * (-1)^n;
mu_exact = exp(fprime_R * T);

% Здесь delta = r - R — отклонение радиуса от выбранного цикла.
rhs_delta = @(t, delta) m * (-1)^n * (R + delta) ...
    .* sin(2 * R * delta + delta.^2);

opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-18, ...
    'MaxStep', T / 200);

fprintf('R = %.12f\n', R);
fprintf('T = %.12f\n', T);
fprintf('Точное значение mu = %.12e\n', mu_exact);

%% Численное отображение Пуанкаре
initial_delta = linspace(-0.15, 0.15, 121).';

% Компоненты этого вектора — независимые начальные отклонения.
[~, delta_all] = ode45(rhs_delta, [0 T], initial_delta, opts);
p_delta = delta_all(end, :).';

%% Производная отображения: уравнение в вариациях
[~, q] = ode45(@(t, q) fprime_R * q, [0 T], 1, opts);
mu_num = q(end);

fprintf('Численное значение mu = %.12e\n', mu_num);
assert(abs(mu_num / mu_exact - 1) < 1e-7);

% Проверяем производную центральной разностью.
for eps0 = [1e-2 1e-3 1e-4]
    [~, delta_pair] = ode45(rhs_delta, [0 T], [-eps0; eps0], opts);
    mu_difference = (delta_pair(end, 2) - delta_pair(end, 1)) ...
        / (2 * eps0);

    fprintf('Разность при eps = %g: %.12e\n', eps0, mu_difference);
end

% Повторяем расчёт с более строгими настройками точности.
opts_fine = odeset('RelTol', 1e-12, 'AbsTol', 1e-19, ...
    'MaxStep', T / 400);
[~, delta_fine] = ode45(rhs_delta, [0 T], initial_delta, opts_fine);
refinement_error = max(abs(delta_fine(end, :).' - p_delta));

fprintf('Изменение отображения при уточнении: %.3e\n', refinement_error);
assert(refinement_error < 1e-12);

%% Траектории для двух начальных условий
delta_start = [-0.15 0.15];
colors = [light_brown; brown];
styles = {'--', '-'};

t_two_periods = linspace(0, 2 * T, 3201).';
[~, delta_two_periods] = ode45(rhs_delta, t_two_periods, ...
    delta_start.', opts);

t_period = linspace(0, T, 1601).';
[~, delta_period] = ode45(rhs_delta, t_period, delta_start.', opts);

%% Рисунок: три графика
fig = figure(1);
set(fig, 'Color', background, 'Position', [80 100 1550 550], ...
    'Name', 'Предельные циклы и отображение Пуанкаре', ...
    'NumberTitle', 'off');

%% Первый график: циклы, близкие траектории и сечение
ax1 = subplot(1, 3, 1);
hold on;

phi = linspace(0, 2 * pi, 1001);

for j = 1:3
    nn = n + j - 2;
    radius = sqrt(nn * pi);

    if nn == n
        color = brown;
        style = '-';
        width = 2.3;
    elseif nn < n
        color = beige;
        style = '--';
        width = 1.5;
    else
        color = beige;
        style = ':';
        width = 1.5;
    end

    h_cycle(j) = plot(radius * cos(phi), radius * sin(phi), ...
        'Color', color, 'LineStyle', style, 'LineWidth', width);
end

for j = 1:2
    radius = R + delta_two_periods(:, j);

    h_trajectory(j) = plot(radius .* cos(k * t_two_periods), ...
        radius .* sin(k * t_two_periods), ...
        'Color', colors(j, :), 'LineStyle', styles{j}, 'LineWidth', 1.2);
end

h_section = plot([R - 0.15 R + 0.15], [0 0], ...
    'Color', brown, 'LineWidth', 4);

axis equal;
xlim([-4 4]);
ylim([-4 4]);
xlabel('x');
ylabel('y');
title('Третий цикл и сечение');
hold off;

leg1 = legend([h_cycle h_trajectory h_section], ...
    {sprintf('Цикл %d', n - 1), sprintf('Цикл %d', n), ...
     sprintf('Цикл %d', n + 1), '\delta_0 = -0.15', ...
     '\delta_0 = +0.15', 'Сечение'}, ...
    'Location', 'southoutside', 'NumColumns', 2);

% –––––

%% Второй график: отображение Пуанкаре
ax2 = subplot(1, 3, 2);

plot(initial_delta, 1e6 * p_delta, ...
    'Color', brown, 'LineWidth', 2.2);
hold on;

plot(initial_delta, 1e6 * mu_exact * initial_delta, '--', ...
    'Color', light_brown, 'LineWidth', 2);

xlim([-0.15 0.15]);
xticks([-0.15 -0.10 -0.05 0 0.05 0.10 0.15]);
xlabel('\delta_0 = r_0 - R');
ylabel('10^6 (P(R + \delta_0) - R)');
title('Отображение Пуанкаре');
hold off;

leg2 = legend('Численное отображение', '\mu\delta_0', ...
    'Location', 'southoutside');

% –––––

%% Третий график: сближение с циклом за один оборот
ax3 = subplot(1, 3, 3);

for j = 1:2
    semilogy(t_period, abs(delta_period(:, j)), ...
        'Color', colors(j, :), 'LineStyle', styles{j}, 'LineWidth', 2);
    hold on;
end

xlim([0 T]);
xlabel('t');
ylabel('|r(t) - R|');
title('Сближение за один оборот');
hold off;

leg3 = legend('\delta(0) = -0.15', '\delta(0) = +0.15', ...
    'Location', 'southoutside');

%% Общее оформление осей
axes_list = [ax1 ax2 ax3];

for j = 1:3
    ax = axes_list(j);

    set(ax, 'Color', background, 'Box', 'off', ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'GridColor', gray, 'GridAlpha', 1, 'GridLineStyle', '-', ...
        'Layer', 'bottom', 'LineWidth', 0.7, ...
        'FontName', 'Helvetica', 'FontSize', 10, ...
        'TickLength', [0 0], 'XMinorTick', 'off', 'YMinorTick', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');

    set(get(ax, 'XLabel'), 'FontSize', 12);
    set(get(ax, 'YLabel'), 'FontSize', 12);
    set(get(ax, 'Title'), 'FontSize', 13, 'FontWeight', 'normal');
end

%% Размещение графиков и легенд
set(ax1, 'Position', [0.055 0.31 0.260 0.61]);
set(ax2, 'Position', [0.385 0.31 0.252 0.61]);
set(ax3, 'Position', [0.725 0.31 0.250 0.61]);

legend_list = [leg1 leg2 leg3];

for j = 1:3
    set(legend_list(j), 'Color', background, 'EdgeColor', gray, ...
        'FontSize', 9, 'Interpreter', 'tex', 'Units', 'normalized');
end

set(leg1, 'Position', [0.055 0.045 0.260 0.15]);
set(leg2, 'Position', [0.385 0.085 0.252 0.10]);
set(leg3, 'Position', [0.725 0.085 0.250 0.10]);

drawnow;
fprintf('0 < mu < 1: третий цикл орбитально асимптотически устойчив.\n');

% –––––

%% Запуск файла с Octave (раскомментировать)
if exist('OCTAVE_VERSION', 'builtin')
    disp('Графики построены.');
    pause;
end