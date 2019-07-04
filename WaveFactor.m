%Code by Akhil Kandhari, axk751@case.edu
%This allows to maximize the wavefactor such that a given wave requires minimum power for motion

%This code uses fmincon to find the optimized values for w,m and b given constraint equations and
%boundary conditions

clear all
clc

global R L t n nu;
R = 23/2; %Radius of segment
t = 0.3175;
L = 20; %Length of segment
nu = 1.3;


%w; %Number of waves
%m; %Pair of moving segments
%b; %Number of bridged segments

%Set initial guess values
w_guess = 1;
m_guess = 1;
b_guess = 0;

x0 = [w_guess, m_guess, b_guess]; %Initial values
A =[]; %Linear equality 
B = [];
Aeq = []; 
Beq = [];
LB = [1.0,1.0,0.0]; %Lower bound for waveform variables
UB = [];
%UB = [inf, n, inf]; %Upper bound for waveform variables

%Call the solver to find optimized result
%help fmincon for understanding the syntax below

for num = 6:150
    n = num;
    UB = [floor(n/2),floor(n/2),floor(n/2)];
    xopt = fmincon(@objective, x0, A, B, Aeq, Beq, LB, UB, @constraint);
    %xopt = floor(xopt); %Rounds the optimized waveform found
    [min_COT,wopt,mopt,bopt] = CalcMinCOT(xopt);
    xopt = [wopt,mopt,bopt];
    num_waves(num) = xopt(1);
    num_pairs(num) = xopt(2);
    num_bridge(num) = xopt(3);
    COT(num) = CalcCot(xopt); %Calculated COT based on optimum waveform
    percentage(num) = AnchoringSegments(xopt);
end

figure
plot(6:150, COT(6:150));

figure
plot(6:150,percentage(6:150));

%[X,Y] = meshgrid(num_waves,num_pairs);
%mesh(X,Y,COT)


%Define function that needs to be minimzed
function MinFunc = objFunc(x)
    global R L t n nu;
    w = x(1);
    m = x(2);
    b = x(3);
    
    %Function based on the COT equation derived for thin cylinder
    MinFunc = (1/nu)*(((n*pi^3*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));
    %(1/nu)*(((n*pi^3*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/(384/5)))*(n/(w*m*(m+b)));
    %(1/nu)*(((n*pi^3*R^3)/((n-w*(2*m+b))*4*L*t^2)))*(n/(w*m*(m+b)));
    %(1/nu)*(((n*pi^3*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));
end

%Define constraint equations
function ConstEq = bridgeSegments(x)
    global n;
    w = x(1);
    m = x(2);
    b = x(3);
    
    ConstEq = (w)*(2*m+b)-n;
end

%Define objective function for minimzation
function obj = objective(x)
    obj = objFunc(x);
end

%Define constraint for optimization
function [c, ceq] = constraint(x)
    c = bridgeSegments(x);
    ceq = [];
end

function [min_COT,wopt,mopt,bopt] = CalcMinCOT(x)
    global n;
    w = [floor(x(1)),ceil(x(1))];
    m = [floor(x(2)),ceil(x(2))];
    b = [floor(x(3)),ceil(x(3))];
    min_COT = inf;
    wopt = w(1);
    mopt = m(1);
    bopt = b(1);
    for i=1:2
        for j = 1:2
            for k = 1:2
                c = bridgeSegments([w(i),m(j),b(k)]);
                linconst = 2*m(j)+b(k)-floor(n/2)+1;
                if(c<0 && linconst <0)
                    COT = CalcCot([w(i),m(j),b(k)]);
                    if (COT<min_COT)
                        min_COT = COT;
                        wopt = w(i);
                        mopt = m(j);
                        bopt = b(k);
                    end
                end
            end
        end
    end
end

function COT = CalcCot(x)
    global R L t n nu;
    w = x(1);
    m = x(2);
    b = x(3);
    COT = (1/nu)*(((n*pi^3*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b))); 
end

function perc_anchoring = AnchoringSegments(x)
    
    global R L t n nu;
    w = x(1);
    m = x(2);
    b = x(3);
    num_anchoring = n-w*(2*m+b);
    num_moving = w*(2*m+b);
    perc_anchoring = (num_anchoring/n)*100; 
end

