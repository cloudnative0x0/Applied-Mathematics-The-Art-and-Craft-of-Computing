%% Постоянные
m = 2;
k = 17;

%% Палитра
maroon     = [0.2078, 0.3608, 0.4902]; % вещественная часть
red        = [0.7529, 0.4235, 0.5176]; % мнимая часть
cream      = [0.8000, 0.8000, 0.8000]; % серая сетка и равновесия
background = [1.0000, 1.0000, 1.0000]; % фон

%% Временной отрезок и начальные условия
T = 12 / k;

initial_y = [ ...
    -0.85 * pi, -0.50 * pi, -0.15 * pi, ...
     0.15 * pi,  0.50 * pi,  0.85 * pi, ...
     1.15 * pi,  1.50 * pi,  1.85 * pi];

initial_x = initial_y / k;

%% x' = sin^2(kx)
rhs1 = @(t, x) sin(k * x).^2;

figure(1);
set(gcf, 'Color', background);

subplot(1, 2, 1);
set(gca, 'Color', background);
hold on;

plots = [];
labels = cell(1, length(initial_x) + 1);

for j = 1:length(initial_x)
    x0 = initial_x(j);

    [t_left, x_left] = ode45(rhs1, [0, -T], x0);
    [t_right, x_right] = ode45(rhs1, [0, T], x0);

    t = [flipud(t_left); t_right(2:end)];
    x = [flipud(x_left); x_right(2:end)];

    if mod(j - 1, 3) == 0
        plots(end + 1) = plot(t, x, '-', ...
            'Color', maroon, 'LineWidth', 1.5);
    elseif mod(j - 1, 3) == 1
        plots(end + 1) = plot(t, x, '--', ...
            'Color', red, 'LineWidth', 1.5);
    else
        plots(end + 1) = plot(t, x, ':', ...
            'Color', red, 'LineWidth', 1.5);
    end

    labels{j} = sprintf('kx(0)/\\pi = %.2f', k * x0 / pi);
end

stationary = plot(nan, nan, '--', 'Color', cream, 'LineWidth', 0.9);

for j = -1:2
    yline(j * pi / k, '--', 'Color', cream, ...
        'LineWidth', 0.9, 'HandleVisibility', 'off');
end

grid on;
set(gca, 'GridColor', cream);
xlabel('t');
ylabel('x(t)');
title('x'' = sin^2(kx)');
labels{end} = 'x = j\pi/k';
legend([plots, stationary], labels, 'Location', 'southoutside');
hold off;

%% x' = |sin(kx)|
rhs2 = @(t, x) abs(sin(k * x));

subplot(1, 2, 2);
set(gca, 'Color', background);
hold on;

plots = [];
labels = cell(1, length(initial_x) + 1);

for j = 1:length(initial_x)
    x0 = initial_x(j);

    [t_left, x_left] = ode45(rhs2, [0, -T], x0);
    [t_right, x_right] = ode45(rhs2, [0, T], x0);

    t = [flipud(t_left); t_right(2:end)];
    x = [flipud(x_left); x_right(2:end)];

    if mod(j - 1, 3) == 0
        plots(end + 1) = plot(t, x, '-', ...
            'Color', maroon, 'LineWidth', 1.5);
    elseif mod(j - 1, 3) == 1
        plots(end + 1) = plot(t, x, '--', ...
            'Color', red, 'LineWidth', 1.5);
    else
        plots(end + 1) = plot(t, x, ':', ...
            'Color', red, 'LineWidth', 1.5);
    end

    labels{j} = sprintf('kx(0)/\\pi = %.2f', k * x0 / pi);
end

stationary = plot(nan, nan, '--', 'Color', cream, 'LineWidth', 0.9);

for j = -1:2
    yline(j * pi / k, '--', 'Color', cream, ...
        'LineWidth', 0.9, 'HandleVisibility', 'off');
end

grid on;
set(gca, 'GridColor', cream);
xlabel('t');
ylabel('x(t)');
title('x'' = |sin(kx)|');
labels{end} = 'x = j\pi/k';
legend([plots, stationary], labels, 'Location', 'southoutside');
hold off;

%% x' = sqrt(|sin(kx)|)
rhs3 = @(t, x) sqrt(abs(sin(k * x)));

figure(2);
set(gcf, 'Color', background);

subplot(1, 2, 1);
set(gca, 'Color', background);
hold on;

sqrt_initial_y = [0, pi / 4, pi / 2, 3 * pi / 4, pi];
sqrt_initial_x = sqrt_initial_y / k;
plots = [];
labels = cell(1, length(sqrt_initial_x) + 1);

for j = 1:length(sqrt_initial_x)
    x0 = sqrt_initial_x(j);

    [t_left, x_left] = ode45(rhs3, [0, -T], x0);
    [t_right, x_right] = ode45(rhs3, [0, T], x0);

    t = [flipud(t_left); t_right(2:end)];
    x = [flipud(x_left); x_right(2:end)];

    if mod(j, 2) == 1
        color = maroon;
    else
        color = red;
    end

    plots(end + 1) = plot(t, x, 'Color', color, 'LineWidth', 1.5);
    labels{j} = sprintf('kx(0)/\\pi = %.2f', k * x0 / pi);
end

stationary = plot(nan, nan, '--', 'Color', cream, 'LineWidth', 0.9);

for j = -3:4
    yline(j * pi / k, '--', 'Color', cream, ...
        'LineWidth', 0.9, 'HandleVisibility', 'off');
end

grid on;
set(gca, 'GridColor', cream);
xlabel('t');
ylabel('x(t)');
title('x'' = sqrt(|sin(kx)|)');
labels{end} = 'x = j\pi/k';
legend([plots, stationary], labels, 'Location', 'southoutside');
hold off;

%% x' = sign(sin(kx)) sqrt(|sin(kx)|)
rhs4 = @(t, x) sign(sin(k * x)) .* sqrt(abs(sin(k * x)));

subplot(1, 2, 2);
set(gca, 'Color', background);
hold on;

plots = [];
labels = cell(1, length(initial_x) + 1);

for j = 1:length(initial_x)
    x0 = initial_x(j);

    [t_left, x_left] = ode45(rhs4, [0, -T], x0);
    [t_right, x_right] = ode45(rhs4, [0, T], x0);

    t = [flipud(t_left); t_right(2:end)];
    x = [flipud(x_left); x_right(2:end)];

    if mod(j - 1, 3) == 0
        plots(end + 1) = plot(t, x, '-', ...
            'Color', maroon, 'LineWidth', 1.5);
    elseif mod(j - 1, 3) == 1
        plots(end + 1) = plot(t, x, '--', ...
            'Color', red, 'LineWidth', 1.5);
    else
        plots(end + 1) = plot(t, x, ':', ...
            'Color', red, 'LineWidth', 1.5);
    end

    labels{j} = sprintf('kx(0)/\\pi = %.2f', k * x0 / pi);
end

stationary = plot(nan, nan, '--', 'Color', cream, 'LineWidth', 0.9);

for j = -1:2
    yline(j * pi / k, '--', 'Color', cream, ...
        'LineWidth', 0.9, 'HandleVisibility', 'off');
end

grid on;
set(gca, 'GridColor', cream);
xlabel('t');
ylabel('x(t)');
title('x'' = sign(sin(kx)) sqrt(|sin(kx)|)');
labels{end} = 'x = j\pi/k';
legend([plots, stationary], labels, 'Location', 'southoutside');
hold off;

%% Запуск файла с Octave (раскомментировать)
if exist('OCTAVE_VERSION', 'builtin')
    disp('Графики построены.');
    pause;
end
