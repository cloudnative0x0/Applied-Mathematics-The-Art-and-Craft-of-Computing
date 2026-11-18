clear; clc; close all;

%% Постоянные
m = 2;
k = 17;

%% Палитра
brown       = [0.4627, 0.3294, 0.2314]; % коричневый
light_brown = [0.7451, 0.5843, 0.4353]; % светло-коричневый
beige       = [0.8549, 0.7765, 0.6667]; % бежевый
gray        = [0.8000, 0.8000, 0.8000]; % серая сетка
background  = [1.0000, 1.0000, 1.0000]; % белый фон

%% Равновесия, начальные условия и точная формула
R = hypot(m, k);
phi = atan2(k, m);
x_left = -phi;
x_stable = pi-phi;
x_right = 2*pi-phi;
fractions = [0.08 0.25 0.50 0.75 0.92];
initial_x = x_left + 2*pi*fractions;
t = linspace(0, 0.7, 701).';

rhs = @(~, x) m*sin(x) + k*cos(x);

fprintf('Параметры и контроль точности\n');
fprintf('%s\n', repmat('-',1,38));
fprintf('R = %.12f,   phi = %.12f\n', R, phi);
fprintf('Неустойчивые равновесия: %.9f; %.9f\n', x_left, x_right);
fprintf('Устойчивое равновесие:   %.9f\n', x_stable);
fprintf('Характерное время 1/R:    %.6f\n', 1/R);
fprintf('Начальные условия:  %s\n', sprintf('%.6f  ', initial_x));

%% Численное решение и проверка точной формулы
n_sol = numel(initial_x);
x_num = zeros(numel(t),n_sol);
x_num(1,:) = initial_x;
for step = 1:numel(t)-1
    h = t(step+1) - t(step);
    a = rhs(t(step),       x_num(step,:));
    b = rhs(t(step) + h/2, x_num(step,:) + h*a/2);
    c = rhs(t(step) + h/2, x_num(step,:) + h*b/2);
    d = rhs(t(step) + h,   x_num(step,:) + h*c);
    x_num(step+1,:) = x_num(step,:) + h*(a + 2*b + 2*c + d)/6;
end

x_exact = zeros(numel(t),n_sol);
for j = 1:n_sol
    x0 = initial_x(j);
    x_exact(:,j) = x_stable + 2*atan( ...
        tan((x0-x_stable)/2)*exp(-R*t));
end

errors = abs(x_num-x_exact);
max_error = max(errors(:));
fprintf('Максимальная ошибка Рунге: %.3e\n', max_error);
assert(max_error < 1e-8);
assert(abs(rhs(0,x_stable)) < 1e-12);

%% Графики
fig = figure('Color', background, 'Position', [100 100 1560 576]);
ax1 = subplot(1, 2, 1, 'Parent', fig);
ax2 = subplot(1, 2, 2, 'Parent', fig);

colors = [brown; light_brown; brown; light_brown; brown];
line_styles = {'-','-',':','--','--'};
labels = arrayfun(@(x0) sprintf('x(0)=%.3f', x0), initial_x, ...
    'UniformOutput',false);

% Левый график: точные решения и равновесия
hold(ax1,'on');
h_sol = zeros(1,n_sol);
for j = 1:n_sol
    h_sol(j) = plot(ax1, t, x_exact(:,j), 'LineStyle', line_styles{j}, ...
        'Color', colors(j,:), 'LineWidth', 2.0);
end
plot(ax1,[t(1) t(end)],[x_left x_left],'--','Color',gray,'LineWidth',1.2, ...
    'HandleVisibility','off');
plot(ax1,[t(1) t(end)],[x_right x_right],'--','Color',gray,'LineWidth',1.2, ...
    'HandleVisibility','off');
h_stable = plot(ax1,[t(1) t(end)],[x_stable x_stable],':', ...
    'Color',beige,'LineWidth',2.0);
plot(ax1,[0.62 0.67],[x_stable+0.38 x_stable],'-', ...
    'Color',brown,'LineWidth',0.8,'HandleVisibility','off');
text(ax1,0.62,x_stable+0.38,'x_s','Color',brown,'FontSize',12, ...
    'HorizontalAlignment','center');
title(ax1,'Точные решения x(t)');
xlabel(ax1,'t'); ylabel(ax1,'x(t)');
xlim(ax1,[t(1) t(end)]);
lg1 = legend(ax1,[h_sol h_stable],[labels {'устойчивое равновесие'}], ...
    'Location','southeast');
set(lg1,'FontSize',8,'Box','off');

% Правый график: для каждой траектории своя ошибка
hold(ax2,'on');
h_error = zeros(1,n_sol);
for j = 1:n_sol
    h_error(j) = semilogy(ax2, t, max(errors(:,j), 1e-16), ...
        'LineStyle', line_styles{j}, 'Color', colors(j,:), 'LineWidth', 1.8);
end

title(ax2,'Ошибка относительно точной формулы');
xlabel(ax2,'t'); ylabel(ax2,'Абсолютная ошибка');
xlim(ax2,[t(1) t(end)]);
lg2 = legend(ax2,h_error,labels,'Location','northeast');
set(lg2,'FontSize',8,'Box','off');

% Общие параметры осей
for ax = [ax1 ax2]
    grid(ax,'on');
    set(ax,'FontSize',11,'Box','off','LineWidth',0.75);
    if isprop(ax,'GridColor')
        set(ax,'GridColor',gray);
    end
    if isprop(ax,'GridAlpha')
        set(ax,'GridAlpha',1);
    end
end

% –––––

%% Удерживаем окна в Octave
%if exist('OCTAVE_VERSION', 'builtin')
%    disp('Графики построены.');
%    pause;
%end