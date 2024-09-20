clc;  % Borra el área de trabajo
clear;  % Borra las variables almacenadas
format long;  % Utiliza la máxima capacidad de la máquina

fprintf('                     FACTORIZACION LU CHOLESKY\n\n\n');

% Solicita la matriz A y el vector b
A = [10, 2, -1;
     -3, -6, 2;
     1, 1, 5];


b = [12; 18; -6];
[n, m] = size(A);
C = [A, b];  % Matriz aumentada [Ab]
disp(C);

% Inicialización de las matrices L y U
L = zeros(n);
U = zeros(n);

if n == m
    for k = 1:n
        suma1 = 0;
        for p = 1:k-1
            suma1 = suma1 + L(k, p) * U(p, k);
        end
        L(k, k) = sqrt(A(k, k) - suma1);
        U(k, k) = L(k, k);  % Inicio del método
        for i = k+1:n
            suma2 = 0;
            for q = 1:k-1
                suma2 = suma2 + L(i, q) * U(q, k);
            end
            L(i, k) = (A(i, k) - suma2) / L(k, k);  % Obtención de la matriz L
        end
        for j = k+1:n
            suma3 = 0;
            for r = 1:k-1
                suma3 = suma3 + L(k, r) * U(r, j);
            end
            U(k, j) = (A(k, j) - suma3) / L(k, k);  % Obtención de la matriz U
        end
    end

    producto = det(L) * det(U);  % Cálculo del determinante

    if producto ~= 0
        % Inicialización de los vectores z y x
        z = zeros(n, 1);
        x = zeros(n, 1);

        for i = 1:n
            suma = 0;
            for p = 1:i-1
                suma = suma + L(i, p) * z(p);
            end
            z(i) = (b(i) - suma) / L(i, i);  % Obtención del vector Z
        end

        for i = n:-1:1
            suma = 0;
            for p = i+1:n
                suma = suma + U(i, p) * x(p);
            end
            x(i) = (z(i) - suma) / U(i, i);  % Solución de las variables
        end
    else
        fprintf('\nEl determinante es igual a cero, por lo tanto el sistema tiene infinitas o ninguna solución\n');
    end
end

fprintf('\n Matriz Ab:\n');
disp(C);
fprintf('\n Matriz L:\n');
disp(L);
fprintf('\n Matriz U:\n');
disp(U);
fprintf('\n El vector Z:\n');
disp(z);

fprintf('\n\nLa solución de X1 hasta Xn es:\n');
for i = 1:n
    xi = x(i);
    fprintf('\nX%g=', i);
    disp(xi);
end

