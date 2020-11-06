length = 200; 
width = 100;

num_steps = 100000;

global e_x e_y
e_x = zeros(9,1);
e_y = zeros(9,1);

% i really need to know how this is meant to be done cuz it cant be this
e_x(2) = 1; e_x(3) = -1;
e_y(4) = 1; e_y(5) = -1;

e_x(6) = 1; e_y(6) = 1;
e_x(7) = 1; e_y(7) = -1;
e_x(8) = -1; e_y(8) = 1;
e_x(9) = -1; e_y(9) = -1;

global omega
% 2-9 are directions, 1 = stationary
omega = zeros(9,1); 
omega(1) = 4/9;
omega(2:5) = 1/9;
omega(6:9) = 1/36;

u0 = 0.05;

[f, ux, uy, rho] = initialize(length, width, u0);

% stream_val determines how quickly collisions bring density to equilibrium
% also 1/relaxation constant?
stream_val = 1.78;

global barriers
barriers = zeros(length,width);
barriers(50, 40:60) = 1;

for step = 1:num_steps
    [f, ux, uy, rho] = collide_step(f, length, width, stream_val);
    f = stream_step(f, length, width);
    
    vel_squared = ux.*ux + uy.*uy;
    
    curl = zeros(length,width);
    [X,Y] = meshgrid(1:length, 1:width);
    for curlx = 2:length-1
        for curly = 2:width-1
            curl(curlx, curly) = uy(curlx+1, curly) - uy(curlx-1, curly) + ...
                ux(curlx, curly+1) - ux(curlx, curly-1);
        end
    end
    vel = sqrt(vel_squared);
    plot_val = surf(X,Y, transpose(curl));
    set(plot_val,'LineStyle','none');
    view(2);
    pause(0.002);
end

% get better function name
function f_ret = calculate(ux, uy, ex, ey, rho, prob_dist, l, w)
    f_ret = zeros(9, l, w);
    for i = 1:9
        edotu_i = ex(i) * ux + ey(i) * uy;
        usquared = ux .* ux + uy .* uy;
        % ignore c for now
        s_i = prob_dist(i) * (1 + 3*edotu_i + 4.5*(edotu_i .* edotu_i) - 1.5*usquared);
        f_ret(i,:,:) = s_i .* rho;
    end
end

function [f_init, ux_init, uy_init, rho_init] = initialize(l, w, u0)
    global e_x e_y omega
    rho_init = ones(l, w);
    ux_init = ones(l, w) * u0;
    uy_init = zeros(l, w);
    f_init = calculate(ux_init, uy_init, e_x, e_y, rho_init, omega, l, w);
end

% matlab is kinda weird with variable storage so i wanna be careful
% that's why the names are so convoluted here
function [f_ret, ux_ret, uy_ret, rho_ret] = collide_step(f_arg, l, w, stream_val)
    global e_x e_y omega
    % feels kinda bad declaring a new array each time but what can ya do
    rho_ret = zeros(l,w);
    for i = 1:9
        rho_ret = rho_ret + squeeze(f_arg(i,:,:));
    end
    ux_ret = zeros(l,w);
    uy_ret = zeros(l,w);
    for i = 1:9
        ux_ret = ux_ret + e_x(i)*squeeze(f_arg(i,:,:));
        uy_ret = uy_ret + e_y(i)*squeeze(f_arg(i,:,:));
    end
    ux_ret = ux_ret ./ rho_ret;
    uy_ret = uy_ret ./ rho_ret;
    
    f_ret = calculate(ux_ret, uy_ret, e_x, e_y, rho_ret, omega, l, w);
    f_ret = f_arg + stream_val * (f_ret - f_arg);
end

function f_new = stream_step(f, l, w)
    global barriers
    % going to do this manually, sorry matlab i have failed you
    f_new = zeros(9, l, w);
    % copy edges
    for i = 1:9
        for a = 1:l
            f_new(i,a,1) = f(i,a,1);
            f_new(i,a,w) = f(i,a,w);
        end
        for b = 1:w
            f_new(i,1,b) = f(i,1,b);
            f_new(i,l,b) = f(i,l,b);
        end
    end
    % move center
    for i = 2:l-1
        for j = 2:w-1
            f_new(1,i,j) = f(1,i,j); % stationary
            f_new(2,i,j) = f(2,i-1,j); % east
            f_new(3,i,j) = f(3,i+1,j); % west
            f_new(4,i,j) = f(4,i,j-1); % south
            f_new(5,i,j) = f(5,i,j+1); % north
            f_new(6,i,j) = f(6,i-1,j-1); % southeast?
            f_new(7,i,j) = f(7,i-1,j+1); % ive lost track
            f_new(8,i,j) = f(8,i+1,j-1);
            f_new(9,i,j) = f(9,i+1,j+1);
        end
    end
    
    % should handle barriers here
    for i = 2:l-1
        for j = 2:w-1
            if(barriers(i,j) == 1)
                f_new(3,i-1,j) = f_new(2,i,j);
                f_new(2,i+1,j) = f_new(3,i,j);
                f_new(5,i,j-1) = f_new(4,i,j);
                f_new(4,i,j+1) = f_new(5,i,j);
                f_new(9,i-1,j-1) = f_new(6,i,j);
                f_new(8,i-1,j+1) = f_new(7,i,j);
                f_new(7,i+1,j-1) = f_new(8,i,j);
                f_new(6,i+1,j+1) = f_new(9,i,j);
            end
        end
    end
end

