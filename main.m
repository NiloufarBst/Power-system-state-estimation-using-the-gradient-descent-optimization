%%%%%%%% Niloufar Besharatzad _ Nonlinear optimization project %%%%%%%
%% Power System State Estimation
clc; 
clear; 
close all;

x_true= [1.05; 0.98; 0.1];  %true states
x = [1; 1; 0]; %Initial states

%powers in true states 
P12 = (x_true(1)*x_true(2)*cos(x_true(3))-x_true(1)^2)+5*x_true(1)*x_true(2)*sin(x_true(3));
Q12 = 5*(x_true(1)*x_true(2)*cos(x_true(3))-x_true(1)^2)-x_true(1)*x_true(2)*sin(x_true(3));
P21 = (x_true(1)*x_true(2)*cos(x_true(3))-x_true(2)^2)-5*x_true(1)*x_true(2)*sin(x_true(3));
Q21=5*(x_true(1)*x_true(2)*cos(x_true(3))-x_true(2)^2)+x_true(1)*x_true(2)*sin(x_true(3));

%parameters
kmax=100000; 
alpha0=1;
eta= 1.1;
eps=0.4 ;

grad_norms=[];
err=[];

for k=1:kmax
    V1=x(1); 
    V2=x(2); 
    ang2=x(3);

    h1=(V1*V2*cos(ang2)-V1^2)+5*V1*V2*sin(ang2);
    h2=5*(V1*V2*cos(ang2)-V1^2)-V1*V2*sin(ang2);
    h3=(V1*V2*cos(ang2)-V2^2)-5*V1*V2*sin(ang2);
    h4=5*(V1*V2*cos(ang2)-V2^2)+V1*V2*sin(ang2);

    %Residuals
    r1=P12-h1;
    r2=Q12-h2;
    r3=P21-h3;
    r4=Q21-h4;

    %partial derivations
    dh1 =[V2*cos(ang2)-2*V1+5*V2*sin(ang2), V1*cos(ang2)+5*V1*sin(ang2), -V1*V2*sin(ang2)+5*V1*V2*cos(ang2)];
    dh2 =[5*V2*cos(ang2)-10*V1-V2*sin(ang2), 5*V1*cos(ang2)-V1*sin(ang2), -5*V1*V2*sin(ang2)-V1*V2*cos(ang2)];
    dh3 =[V2*cos(ang2)-5*V2*sin(ang2), V1*cos(ang2)-2*V2-5*V1*sin(ang2), -V1*V2*sin(ang2)-5*V1*V2*cos(ang2)];
    dh4 =[5*V2*cos(ang2)+V2*sin(ang2), 5*V1*cos(ang2)-10*V2+V1*sin(ang2),-5*V1*V2*sin(ang2)+V1*V2*cos(ang2)];

    gk = -2*(r1*dh1+r2*dh2+r3*dh3+ r4*dh4)';
    
    grad_norms(end+1) = norm(gk);
    err(end+1) = norm(x-x_true)^2;

    % if norm(gk)<=1e-9 || k==kmax
    %     break;
    % end

    %backtracking line search using armijo's rule
    alpha = alpha0;
    f = r1^2 + r2^2 + r3^2 + r4^2;
    while true
        x_new = x-alpha*gk;
        V1k = x_new(1);
        V2k = x_new(2); 
        angk = x_new(3);

        h1k = (V1k*V2k*cos(angk)-V1k^2) + 5*V1k*V2k*sin(angk);
        h2k = 5*(V1k*V2k*cos(angk)-V1k^2)-V1k*V2k*sin(angk);
        h3k = (V1k*V2k*cos(angk)-V2k^2)-5*V1k*V2k*sin(angk);
        h4k = 5*(V1k*V2k*cos(angk)-V2k^2) + V1k*V2k*sin(angk);

        r1k=P12-h1k;
        r2k=Q12-h2k;
        r3k=P21-h3k;
        r4k= Q21-h4k;

        f_k = r1k^2 + r2k^2 + r3k^2 + r4k^2;
        g_k2=norm(gk)^2 ;
        if f_k<=f-eps*alpha*g_k2
            break;
        end
        alpha =alpha/eta; %make alpha smaller to satify the condition
    end

    x= x-alpha*gk;
end

disp('X =')
disp((x))

disp('||xk - x*||^2 =')
disp(err(k))

disp('||g(xk)|| =')
disp(grad_norms(k))


%% Plot(i)  
figure;
semilogy(grad_norms);
xlabel('k'); 
ylabel('||g(xk)||');
title('gradient norm at each iteration');
grid on;

%% Plot(ii)  
figure;
semilogy(err);
xlabel('k');
ylabel('||xk - x*||^2');
title('error at each iteration'); 
grid on;



%% checking hessian at optimal point 
x_opt = x;

V1= x_opt(1);
V2=x_opt(2);
ang2=x_opt(3);

dh1 =[V2*cos(ang2)-2*V1+5*V2*sin(ang2), V1*cos(ang2)+5*V1*sin(ang2), -V1*V2*sin(ang2)+5*V1*V2*cos(ang2)];
dh2 =[5*V2*cos(ang2)-10*V1-V2*sin(ang2), 5*V1*cos(ang2)-V1*sin(ang2), -5*V1*V2*sin(ang2)-V1*V2*cos(ang2)];
dh3 =[V2*cos(ang2)-5*V2*sin(ang2), V1*cos(ang2)-2*V2-5*V1*sin(ang2), -V1*V2*sin(ang2)-5*V1*V2*cos(ang2)];
dh4 =[5*V2*cos(ang2)+V2*sin(ang2), 5*V1*cos(ang2)-10*V2+V1*sin(ang2),-5*V1*V2*sin(ang2)+V1*V2*cos(ang2)];

F_hess=2*(dh1'*dh1 + dh2'*dh2 + dh3'*dh3 + dh4'*dh4);
disp('F(x*)=')
disp(F_hess)

eig_F_hess=eig(F_hess) %eigenvalues of F(x) 

if all(eig_F_hess>0)
    disp('all eigenvalues of F(x)>0 -->  F(x) is positive definite --> x* relative local minimum');
end

