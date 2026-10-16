clear;
clc;
close all;

%% Постоянные
m = 2;
k = 17;

%% Палитра
dark_purple = [0.2980, 0.1373, 0.3961]; % вещественная часть / основной график
red         = [0.7529, 0.4235, 0.5176]; % мнимая часть / вспомогательные линии
mauve       = [0.6196, 0.4941, 0.7098]; % сетка

%% Первые значения разностного уравнения
months = 80;

% –––––

% x(n + 1) хранит число пар, родившихся в месяце n
x = zeros(1, months + 1);
x(1:5) = [1, 0, 0, m, k];

% Вычисление последовательности по рекуррентной формуле
for n = 5:months
    x(n + 1) = x(n) + m * x(n - 2) + (k - m) * x(n - 3) + (1 - k) * x(n - 4);
end

% Вывод первых одиннадцати значений
fprintf('Месяц n    Число родившихся пар x(n)\n');
for n = 0:10
    fprintf('%7d %28.0f\n', n, x(n + 1));
end

% –––––

%% Характеристический многочлен и оценка Коши
P = @(s) s.^5 - s.^4 - m * s.^2 + (m - k) * s + k - 1;
dP = @(s) 5 * s.^4 - 4 * s.^3 - 2 * m * s + m - k;

% Коэффициенты записаны от старшей степени к свободному члену
coef = [1, -1, 0, -m, m - k, k - 1];
R = 1 + max(abs(coef(2:end) / coef(1)));

fprintf('\nОценка Коши: все корни лежат в круге |s| <= %g\n', R);

%% Графики характеристического многочлена
% График на отрезке, заданном оценкой Коши
figure(1);
s = linspace(-R, R, 6001);
plot(s, P(s), 'Color', red, 'LineWidth', 1.5);
hold on;
plot([-R, R], [0, 0], 'Color', dark_purple, 'LineWidth', 0.8);
grid on;
xlabel('s');
ylabel('P(s)');
title('Отрезок по оценке Коши');
legend('P(s)');
xlim([-R, R]);
hold off;

% Увеличенный участок, содержащий вещественные корни
figure(2);
s = linspace(-3, 3, 6001);
plot(s, P(s), 'Color', red, 'LineWidth', 1.5);
hold on;
plot([-3, 3], [0, 0], 'Color', dark_purple, 'LineWidth', 0.8);
grid on;
xlabel('s');
ylabel('P(s)');
title('Увеличенный участок [-3, 3]');
legend('P(s)');
xlim([-3, 3]);
hold off;

% Нулевые изолинии действительной и мнимой частей P(s)
real_part = linspace(-3, 3, 601);
imaginary_part = linspace(-3, 3, 601);
[X, Y] = meshgrid(real_part, imaginary_part);
Z = X + 1i * Y;

figure(3);
contour(X, Y, real(P(Z)), [0, 0], 'r-', 'LineWidth', 1.3);
hold on;
contour(X, Y, imag(P(Z)), [0, 0], 'b--', 'LineWidth', 1.3);
axis equal;
grid on;
xlabel('Re s');
ylabel('Im s');
title('Нулевые изолинии');
legend('Re P(s) = 0', 'Im P(s) = 0');
hold off;

% –––––

%% Метод Ньютона для пяти начальных приближений
starts = [-2, 0.8, 2, 2i, -2i];
tol_f = 1e-11;
tol_s = 1e-12;
max_iter = 100;

roots_newton = zeros(size(starts));
histories = cell(size(starts));

% Для каждой начальной точки сохраняется вся последовательность итераций.
fprintf('\nНачальная точка        Корень                              |s|          |P(s)|    Итерации\n');

for q = 1:numel(starts)
    s = starts(q);
    history = s;
    converged = false;

    for iteration = 1:max_iter
        derivative = dP(s);

        if abs(derivative) < 1e-14
            error('Производная слишком близка к нулю.');
        end

        s_new = s - P(s) / derivative;
        history(end + 1) = s_new;

        if abs(P(s_new)) < tol_f && abs(s_new - s) < tol_s
            converged = true;
            break;
        end

        s = s_new;
    end

    if ~converged
        error('Метод Ньютона не сошёлся за заданное число итераций.');
    end

    roots_newton(q) = s_new;
    histories{q} = history;

    fprintf('%9.4f%+9.4fi   % .10f%+.10fi    %.10f   %.2e     %2d\n', ...
        real(starts(q)), imag(starts(q)), real(s_new), imag(s_new), ...
        abs(s_new), abs(P(s_new)), numel(history) - 1);
end

% –––––

%% Проверка найденных корней
% Проверяем различность корней, невязку и ненулевую производную.
for q = 1:numel(roots_newton)
    for p = 1:q-1
        assert(abs(roots_newton(q) - roots_newton(p)) > 1e-6);
    end
end
assert(max(abs(P(roots_newton))) < 1e-11);
assert(min(abs(dP(roots_newton))) > 1e-6);

% –––––

% Восстанавливаем коэффициенты многочлена по найденным корням.
coefficients = [1, -1, 0, -m, m - k, k - 1];
coefficient_error = max(abs(poly(roots_newton) - coefficients));
fprintf('Ошибка восстановления коэффициентов: %.2e\n', coefficient_error);
assert(coefficient_error < 1e-10);

fprintf('\nИтерации для положительного растущего корня:\n');
for j = 1:numel(histories{3})
    s = histories{3}(j);
    fprintf('%2d: s = %.12f; |P(s)| = %.3e\n', j - 1, real(s), abs(P(s)));
end

% –––––

%% Общее решение и главный корень
% Матрица системы для коэффициентов A_1, ..., A_5.
V = zeros(5, 5);
for j = 1:5
    V(:, j) = roots_newton(j) .^ (0:4).';
end

initial = [1; 0; 0; m; k];
A = V \ initial;
[~, dominant_index] = max(abs(roots_newton));
rho = real(roots_newton(dominant_index));
C = real(A(dominant_index));

% Характеристики главного корня.
fprintf('\nНаибольший модуль: rho = %.12f\n', rho);
fprintf('Коэффициент: C = %.12f\n', C);
fprintf('Асимптотический прирост за месяц: %.4f%%\n', (rho - 1) * 100);
assert(C > 0);

% Сравнение рекуррентного решения с формулой через найденные корни.
n = 0:months;
V_all = zeros(months + 1, 5);
for j = 1:5
    V_all(:, j) = roots_newton(j) .^ n.';
end
from_roots = V_all * A;

relative_error = max(abs(from_roots.' - x) ./ max(1, x));
assert(relative_error < 1e-10);

% –––––

%% Графики роста решения
figure(4);

% Число рождений и его главная асимптотика.
subplot(1, 2, 1);
positive = x > 0;
semilogy(n(positive), x(positive), 'Color', red, 'LineWidth', 1.5);
hold on;
semilogy(n, C * rho.^n, '--', 'Color', dark_purple, 'LineWidth', 1.5);
grid on;
xlabel('Месяц n');
ylabel('x(n), логарифмическая шкала');
title('Рост числа рождений');
legend('Разностное уравнение', 'C\rho^n');
hold off;

% Отношение соседних значений последовательности.
subplot(1, 2, 2);
plot(n(5:end), x(5:end) ./ x(4:end-1), 'Color', red, 'LineWidth', 1.5);
hold on;
plot([n(5), n(end)], [rho, rho], '--', 'Color', dark_purple, 'LineWidth', 1.5);
grid on;
xlabel('Месяц n');
ylabel('Отношение соседних значений');
title('Приближение к постоянному множителю');
legend('x(n) / x(n-1)', sprintf('rho = %.6f', rho));
hold off;

% Численная проверка для нескольких месяцев.
for month = [20, 40, 60, 80]
    fprintf('n=%2d: x(n)/rho^n = %.9f; x(n)/x(n-1) = %.9f\n', ...
        month, x(month + 1) / rho^month, x(month + 1) / x(month));
end

%% Запуск файла с Octave (раскомментировать)
 if exist('OCTAVE_VERSION', 'builtin')
     disp('Графики построены.');
     pause;
 end
