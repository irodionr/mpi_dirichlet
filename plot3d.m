clear;
load('u.csv');

a = 2;
b = 1;
N = 100;
h = 1/N;

x = (0:a*h:a*(N-1)*h);
y = (0:b*h:b*(N-1)*h);
[xx, yy] = meshgrid(x, y);

figure;
surf(xx, yy, u);
xlabel('x');
ylabel('y');
zlabel('u(x, y)');
