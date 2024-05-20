% Learning Convolution

%% Symbolic approach
clear variables

syms x y h tau
%x = 0:0.5:10;
y = x;        % use -1 and -2 to ensure convolution exist 
h = 3-x;
%q = conv(y, h)
ytau = subs(y, x, tau);     
htau = subs(h, x, x-tau); 
q = int(htau * ytau, tau, -6, 6)   % this assumes x \in [0 \inf]
fplot(q, [0 10])

%% Discrete approach

dx = 0.1;
x = 0:dx:3; 
y = x;
h = 3-x;
q = conv(y, h);
figure; plot((0:(numel(q)-1))*dx,q*dx)
xlim([0 10])
