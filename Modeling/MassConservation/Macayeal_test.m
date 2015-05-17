% This script represents a finite-element model of an ice sheet %
% Warning: element areas have cancelled for this % particular exercise only.

% Load mesh
xy = load('points.txt');
index = load('elements.txt');
index=index(:,:)+1;


N=16;
nodes=length(xy);
nel=length(index);
row=zeros(4050,1);
col=zeros(4050,1);

% Model parameters
g=9.81;
rho=910;
phi=zeros(3,3);
Ao=1/31556926 * 1e-16;
a=0.3/31556926;
L=1500e3/2;
Z=( 5*a*L^ 4/( 2 * Ao * (rho*g)^ 3 ) )^ (1/8); r=linspace(0,1,N)';
sexact=( 4 * ( (1/2).^ (4/3) - (r/2).^ (4/3) ) ).^ (3/8); hold off; plot(L*r,Z*sexact,'g*'); hold on 
nsteps=1500;
dt=10*31556926*a/Z;

% Initialize at zero ice thickness
sn=zeros(N,N);
s=zeros(nodes,1);

% Create the interpolation functions.
alpha=zeros(nel,3);
beta=zeros(nel,3);
for n=1:nel
    [lowtri,uptri]=lu([[xy(index(n,1),1) xy(index(n,2),1) xy(index(n,3),1)]'... 
    [xy(index(n,1),2) xy(index(n,2),2) xy(index(n,3),2)]' ones(3,1)]); 
    phi(:,1)=uptri\(lowtri\ [1 0 0]');
    phi(:,2)=uptri\(lowtri\ [0 1 0]');
    phi(:,3)=uptri\(lowtri\ [0 0 1]');
    for     k=1:3
        alpha(n,k)=phi(1,k);
        beta(n,k)=phi(2,k);
    end
end

for time=1:nsteps
    time
    value=zeros(nodes,1);
    R=zeros(nodes,1);
    d=zeros(nel,1);
    count=0;
    for n=1:nel
    % Effective diffusivity
        for l=1:3
            for k=1:3
                d(n)=d(n)+s(index(n,l))*s(index(n,k))*(alpha(n,l)*alpha(n,k)+beta(n,l)*beta(n,k)); 
            end
        end
        d(n)=d(n)*((s(index(n,1))+s(index(n,2))+s(index(n,3)))/3)^ 5;
        %
        % Load matrix and right-hand-side vector:
        %
        R(index(n,1))=R(index(n,1))+1/3;
        R(index(n,2))=R(index(n,2))+1/3;
        R(index(n,3))=R(index(n,3))+1/3;
        for l=1:3
            for k=1:3
                count=count+1;
                row(count)=index(n,k);
                col(count)=index(n,l);
                if l == k
                    R(index(n,k))=R(index(n,k))+s(index(n,l))/(6*dt);
                    value(count)=1/(6*dt);
                else
                    R(index(n,k))=R(index(n,k))+s(index(n,l))/(12*dt);
                    value(count)=1/(12*dt);
                end
                value(count)=value(count)+d(n)*(alpha(n,l)*alpha(n,k)+beta(n,l)*beta(n,k)); end
            end
        % End loop over elements
    end
    A=sparse(row,col,value);
    % Boundary condition at terminus
    for i=1:31
        A(Bound(i),Bound(i))=1.e12; % is a trick to use large number
        R(Bound(i))=0;
    end
    % Solve the system for new thickness
    s=A\R;
    if rem(time,20) == 1
        for i=1:16
            for j=1:16
                sn(i,j)=s(gamma(i,j));
            end
        end
        plot(L*r,Z*sn(:,1));
    end
% End time-stepping loop
end