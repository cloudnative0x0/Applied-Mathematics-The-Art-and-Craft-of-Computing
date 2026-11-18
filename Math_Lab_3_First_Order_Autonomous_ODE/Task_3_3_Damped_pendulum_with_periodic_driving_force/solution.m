clear; clc; close all;

%% Постоянные
m = 2;
k = 17;

%% Палитра
brown       = [0.4627, 0.3294, 0.2314];
light_brown = [0.7451, 0.5843, 0.4353];
beige       = [0.8549, 0.7765, 0.6667];
cream       = [0.8000, 0.8000, 0.8000];
background  = [1.0000, 1.0000, 1.0000];

%% Шрифт для MATLAB и Octave
font_name = 'sans-serif';

% –––––

%% Физические параметры
l = 1;
m0 = 1;
g = 9.81;
b = k/50;
F = m/10;
c = b/m0;
omega0 = sqrt(g/l);
z0 = [pi/4;0];
alpha = (0:0.05:5).';
t_grid = (0:0.01:100).';
opts = odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.02);
opts_ref = odeset('RelTol',1e-11,'AbsTol',1e-13,'MaxStep',0.01);
amp = zeros(size(alpha));
amp70 = amp;
amp_ref = amp;
residual = amp;
condition = amp;
harmonic_ratio = amp;
period_error = nan(size(alpha));
angles = zeros(numel(t_grid),numel(alpha));

%% Восстановление полного периода
phase = linspace(0, 2 * pi, 4001).';

B = [
    sin(phase), cos(phase), ...
    sin(3 * phase), cos(3 * phase)
];

%% Численный расчёт для всех частот
frequency_count = numel(alpha);

% Сначала идут все углы, затем все угловые скорости.
% Каждая пара уравнений относится к своей частоте и не связана с другими.
initial_state = [
    z0(1) * ones(frequency_count, 1);
    z0(2) * ones(frequency_count, 1)
];

rhs = @(t, state) [
    state(frequency_count + 1:end);
    -c * state(frequency_count + 1:end) ...
    - (g / l) * sin(state(1:frequency_count)) ...
    + F * sin(alpha * t)
];

fprintf('Начат расчёт для %d частоты/частот...\n', frequency_count);
[t, states] = ode45(rhs, t_grid, initial_state, opts);
angles = states(:, 1:frequency_count);
fprintf('Основной расчёт завершён. Проверяем точность...\n');

[tr, states_ref] = ode45(rhs, t_grid, initial_state, opts_ref);
angles_ref = states_ref(:, 1:frequency_count);
fprintf('Уточнённый расчёт завершён. Определяем амплитуды...\n');

%% Частотный эксперимент
for j = 2:frequency_count
    a = alpha(j);

    % Три оценки: [60, 100], [70, 100] и уточнённый расчёт.
    for check = 1:3
        if check == 1
            mask = t >= 60;
            tt = t(mask);
            y = angles(mask, j);
        elseif check == 2
            mask = t >= 70;
            tt = t(mask);
            y = angles(mask, j);
        else
            mask = tr >= 60;
            tt = tr(mask);
            y = angles_ref(mask, j);
        end

        M = [
            sin(a * tt), cos(a * tt), ...
            sin(3 * a * tt), cos(3 * a * tt)
        ];

        q = M \ y;
        restored = B * q;
        amplitude = (max(restored) - min(restored)) / 2;

        if check == 1
            amp(j) = amplitude;
            residual(j) = max(abs(M * q - y));
            condition(j) = cond(M);
            harmonic_ratio(j) = hypot(q(3), q(4)) / hypot(q(1), q(2));
        elseif check == 2
            amp70(j) = amplitude;
        else
            amp_ref(j) = amplitude;
        end
    end

    period = 2 * pi / a;

    if period < 40
        tt = linspace(60, 100 - period, 1501).';
        before = interp1(t, angles(:, j), tt, 'linear');
        after = interp1(t, angles(:, j), tt + period, 'linear');
        period_error(j) = max(abs(after - before));
    end
end

%% Проверяем точность установившейся амплитуды
window_error = max(abs(amp-amp70));
solver_error = max(abs(amp-amp_ref));
finite_period_error = period_error(isfinite(period_error));
fprintf('Максимальная невязка периодического описания: %.3e рад\n',max(residual));
fprintf('Изменение окна [60,100] -> [70,100]: %.3e рад\n',window_error);
fprintf('Ужесточение допусков ode45: %.3e рад\n',solver_error);
fprintf('Повторяемость через период: %.3e рад\n',max(finite_period_error));
fprintf('Наибольшее число обусловленности: %.4f\n',max(condition));
fprintf('Наибольшее отношение третьей гармоники к первой: %.3e\n',max(harmonic_ratio));
assert(max(residual)<1e-4);
assert(window_error<1e-4);
assert(solver_error<1e-6);
assert(max(finite_period_error)<1e-4);
assert(max(condition)<10);
assert(max(harmonic_ratio)<1e-3);

%% Линейная оценка и максимум на заданной сетке
A_linear = F./sqrt((omega0^2-alpha.^2).^2+c^2*alpha.^2);
A_linear(1) = 0;
[max_amp,jmax] = max(amp);
alpha_res_linear = sqrt(omega0^2-c^2/2);
fprintf('\nМаксимум НА СЕТКЕ: alpha=%.2f рад/с; A=%.9f рад = %.6f град\n', ...
    alpha(jmax),max_amp,max_amp*180/pi);
fprintf('Линейная резонансная частота: %.9f рад/с\n',alpha_res_linear);
fprintf('\n alpha      A, рад       A_linear, рад\n');
for j = [1 2 21 41 61 62 63 64 65 81 101]
    fprintf('%6.2f    %.9f    %.9f\n',alpha(j),amp(j),A_linear(j));
end


%% Амплитудно-частотная характеристика
fig1 = figure(1);
set(fig1, 'Color', background, 'Position', [60 100 1400 520], ...
    'Name', 'Амплитуда и влияние трения', 'NumberTitle', 'off');

ax1 = subplot(1, 2, 1);

h_nonlinear = plot(alpha(2:end), amp(2:end), ...
    'Color', brown, 'LineWidth', 2.4);
hold on;

h_linear = plot(alpha(2:end), A_linear(2:end), '--', ...
    'Color', light_brown, 'LineWidth', 1.8);

h_zero = plot(0, 0, 'o', ...
    'Color', brown, 'MarkerFaceColor', brown, ...
    'MarkerSize', 5, 'Clipping', 'off');

h_maximum = plot(alpha(jmax), max_amp, 'o', ...
    'MarkerFaceColor', background, 'MarkerEdgeColor', brown, ...
    'MarkerSize', 7, 'LineWidth', 1.8);

plot([alpha(jmax) alpha(jmax)], [0 max_amp], ':', ...
    'Color', beige, 'LineWidth', 1.1, 'HandleVisibility', 'off');

text(3.40, 0.202, ...
    {sprintf('\\alpha = %.2f рад/с', alpha(jmax)), ...
     sprintf('A = %.5f рад', max_amp)}, ...
    'Color', brown, 'FontName', font_name, 'FontSize', 10, ...
    'VerticalAlignment', 'top', 'Interpreter', 'tex');

xlim([0 5]);
ylim([0 0.22]);
set(ax1, 'XTick', 0:1:5, 'YTick', 0:0.05:0.20);

xlabel('\alpha, рад/с');
ylabel('A, рад');
title('Амплитуда установившихся колебаний');
hold off;

leg1 = legend([h_nonlinear h_linear h_zero h_maximum], ...
    {'Нелинейная модель', 'Линейное приближение', ...
     '\alpha = 0', 'Максимум на сетке'}, ...
    'Location', 'southoutside', 'NumColumns', 2);


%% Влияние трения: линейное приближение
ax2 = subplot(1, 2, 2);
hold on;

damping_values = [c / 2, c, 2 * c];
damping_colors = [brown; light_brown; beige];
damping_styles = {'-', '--', '-.'};
damping_labels = cell(1, 3);

for j = 1:3
    damping = damping_values(j);

    response = F ./ sqrt( ...
        (omega0^2 - alpha(2:end).^2).^2 ...
        + damping^2 * alpha(2:end).^2);

    plot(alpha(2:end), response, ...
        'Color', damping_colors(j, :), ...
        'LineStyle', damping_styles{j}, 'LineWidth', 2.2);

    damping_labels{j} = sprintf('b = %.2f', damping * m0);
end

xlim([0 5]);
ylim([0 0.42]);
set(ax2, 'XTick', 0:1:5, 'YTick', 0:0.1:0.4);

xlabel('\alpha, рад/с');
ylabel('A, рад');
title('Влияние трения');

text(0.97, 0.96, 'Линейная модель', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontName', font_name, ...
    'FontSize', 10, 'Color', [0.47 0.42 0.38]);

hold off;
leg2 = legend(damping_labels, ...
    'Location', 'southoutside', 'NumColumns', 3);


%% Переходный процесс и участок оценки амплитуды
fig2 = figure(2);
set(fig2, 'Color', background, 'Position', [90 50 1200 730], ...
    'Name', 'Колебания во времени', 'NumberTitle', 'off');

targets = [1 3.1 5];
trajectory_colors = [brown; light_brown; brown];
trajectory_styles = {'-', '-', '--'};
trajectory_labels = cell(1, 3);

for j = 1:3
    [~, index] = min(abs(alpha - targets(j)));
    trajectory_labels{j} = sprintf('\\alpha = %.2f рад/с', alpha(index));
end

ax3 = subplot(2, 1, 1);
hold on;

% Светлая область показывает участок, используемый для оценки амплитуды.
window_color = 0.88 * background + 0.12 * beige;
patch([60 100 100 60], [-0.85 -0.85 0.85 0.85], window_color, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');

for j = 1:3
    [~, index] = min(abs(alpha - targets(j)));

    plot(t_grid, angles(:, index), ...
        'Color', trajectory_colors(j, :), ...
        'LineStyle', trajectory_styles{j}, 'LineWidth', 1);
end

plot([60 60], [-0.85 0.85], ':', ...
    'Color', light_brown, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');

text(61.3, 0.76, 't \geq 60 с: участок оценки амплитуды', ...
    'Color', brown, 'FontName', font_name, 'FontSize', 9, ...
    'VerticalAlignment', 'top', 'Interpreter', 'tex');

xlim([0 100]);
ylim([-0.85 0.85]);
set(ax3, 'XTick', 0:20:100, 'YTick', -0.8:0.4:0.8);

xlabel('t, с');
ylabel('\theta, рад');
title('Переходный процесс и вынужденные колебания');
hold off;


%% Установившиеся колебания: последние 20 секунд
ax4 = subplot(2, 1, 2);
hold on;
mask = t_grid >= 80;

for j = 1:3
    [~, index] = min(abs(alpha - targets(j)));

    plot(t_grid(mask), angles(mask, index), ...
        'Color', trajectory_colors(j, :), ...
        'LineStyle', trajectory_styles{j}, 'LineWidth', 1.5);
end

xlim([80 100]);
ylim([-0.22 0.22]);
set(ax4, 'XTick', 80:5:100, 'YTick', -0.2:0.1:0.2);

xlabel('t, с');
ylabel('\theta, рад');
title('Установившиеся колебания: последние 20 секунд');
hold off;

% Одна общая легенда для обоих временных графиков.
leg3 = legend(trajectory_labels, ...
    'Location', 'southoutside', 'NumColumns', 3);


%% Общее оформление осей
axes_list = [ax1 ax2 ax3 ax4];
text_color = [0.24 0.22 0.20];

for j = 1:numel(axes_list)
    ax = axes_list(j);

    set(ax, 'Color', background, 'Box', 'off', ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'GridColor', cream, 'GridAlpha', 0.6, 'GridLineStyle', '-', ...
        'Layer', 'bottom', 'LineWidth', 0.7, ...
        'XColor', text_color, 'YColor', text_color, ...
        'FontName', font_name, 'FontSize', 10, ...
        'TickLength', [0 0]);

    set(get(ax, 'XLabel'), ...
        'FontName', font_name, 'FontSize', 11, 'Interpreter', 'tex');
    set(get(ax, 'YLabel'), ...
        'FontName', font_name, 'FontSize', 11, 'Interpreter', 'tex');

    % Заголовок слева, с одинаковым отступом на всех графиках.
    set(get(ax, 'Title'), ...
        'FontName', font_name, 'FontSize', 13, 'FontWeight', 'bold', ...
        'Color', text_color, 'Interpreter', 'none', ...
        'Units', 'normalized', 'Position', [0 1.08 0], ...
        'HorizontalAlignment', 'left');
end


%% Компактные легенды и расположение графиков
legend_list = [leg1 leg2 leg3];

for j = 1:numel(legend_list)
    set(legend_list(j), ...
        'Box', 'off', 'Color', background, ...
        'FontName', font_name, 'FontSize', 10, ...
        'Interpreter', 'tex', 'Units', 'normalized');
end

% Первый рисунок: одинаковая высота графиков, легенды внизу.
set(ax1, 'Position', [0.07 0.25 0.405 0.60]);
set(ax2, 'Position', [0.56 0.25 0.405 0.60]);

set(leg1, 'Position', [0.07 0.035 0.405 0.105]);
set(leg2, 'Position', [0.56 0.055 0.405 0.065]);

% Второй рисунок: два широких графика и общая легенда.
set(ax3, 'Position', [0.08 0.61 0.895 0.28]);
set(ax4, 'Position', [0.08 0.18 0.895 0.29]);

set(leg3, 'Position', [0.22 0.025 0.62 0.06]);

%% Показываем графики
set(fig1, 'Visible', 'on', 'WindowStyle', 'normal');
set(fig2, 'Visible', 'on', 'WindowStyle', 'normal');

figure(fig1);
drawnow;

figure(fig2);
drawnow;

% –––––

%% Удерживаем окна в Octave
%if exist('OCTAVE_VERSION', 'builtin')
%    waitfor(fig1);
%    waitfor(fig2);
%end