function MILPproblem = buildMILPproblem(model)
% Constructs Thermodynamically feasible MILP problem from
% pseudo-thermodynamic model.
% Definitions
m     = model.N.m;
n1    = model.N.n1;
n2    = model.N.n2;
i1    = model.N.i1;
i2    = model.N.i2;
St    = transpose(model.S);
rev   = model.rev;
u_std = model.u_std;
R     = model.R;
T     = model.T;
K     = 1e6;                                     % Must be at least 10*Fmax
k1    = 1e-3;                                    % Represents min flux
k2    = 1e-3;                                    % Represents min Dgr

% A matrix construction
A_up   = [model.A zeros(m+i1+i2,3*i1)];
I1     = eye(i1,n1);
I2     = I1(:,rev);
I3     = eye(i1);
Z1     = zeros(i1,n1);
Z2     = zeros(i1,n2);
Z3     = zeros(i1,m);
Z4     = zeros(i1,i1);
B3     = R*T*St(1:i1,:);
A_down = [Z1,Z2,Z3,I3,I3,I3;
    (1/K)*I1,-(1/K)*I2,Z3,-I3,Z4,Z4;
    -(1/K)*I1,(1/K)*I2,Z3,Z4,-I3,Z4;
    I1,-I2,Z3,-K*I3,k1*I3,-k1*I3;
    -I1,I2,Z3,k1*I3,-K*I3,-k1*I3;
    Z1,Z2,B3,K*I3,Z4,Z4;
    Z1,Z2,-B3,Z4,K*I3,Z4;
    Z1,Z2,-B3,Z4,Z4,K*I3;
    Z1,Z2,B3,Z4,Z4,K*I3];

MILPproblem.A = [A_up; A_down];

% Other fields
c3                  = K * ones(i1,1) + St(1:i1,:) * u_std;
c4                  = K * ones(i1,1) - St(1:i1,:) * u_std;
MILPproblem.b       = [model.b; ones(i1,1); zeros(4*i1,1);
    c4-k2; c3-k2; c3; c4];                                % k2 is added in order to simulate a strict inequality
MILPproblem.c       = zeros(n1+n2+m+3*i1,1);
MILPproblem.lb      = [model.lb; zeros(3*i1,1)];
MILPproblem.ub      = [model.ub; ones(3*i1,1)];
MILPproblem.osense  = 1;                      % Minimize
MILPproblem.csense  = strcat(repmat('E',1,m),repmat('L',1,i1+i2), ...
    repmat('E',1,i1),repmat('L',1,8*i1));
MILPproblem.vartype = strcat(repmat('C',1,n1+n2+m),repmat('B',1,3*i1));

save('MP.mat',"MILPproblem")

end