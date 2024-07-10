/*
ToDo (in no particular order):
1)  Implement a command line parser for the variables (which will be renamed) a,b,e, L_0, Tend, Delt
2)  Design the structure of the Gas phase problem
    a)  Should have a Rhs/Jac function (the equivalent of a virtual function)
    b)  Should have a solver (will start with RK4 and implement more)
    c)  Should have visualization (see the Crate "plotters")
3)  Refactor to allow for more complex chemical mechanisms (maybe a yaml plugin?)
*/

struct KapProb{
    active : bool,
    state  : Vec<f32>
}
impl KapProb{
    fn rhs(&self){
        println!("This is the rhs function");
    }
    fn jac(&self){
        println!("This is the jac function");
    }
    fn solve(&self){
        println!("This is the attempt to solve the problem");
    }
    fn printState(&self){
        println!("The vector is {:?}", &self.state);
    }
    fn step(&self){
        println!("This impl will perform a step");
    }
}


fn main() {
    println!("Hello, world!");
    let test_problem = KapProb{ active : true, state : vec![1.0,0.0,0.0]};
    test_problem.rhs();
    test_problem.jac();
    test_problem.solve();
    test_problem.printState();


}
/* 
function [y,TIME]=KapilaSimple(a,b,e,l_0,Tend,Delt)
%Typical run kapila(1.0,0.5,1e-2,0.5,10,1e-2);
%a is the first species amount
%b controls z, the second species.  Modulated by epsilon
%l_0 thermal conductivity
%e is epsilon, the inverse reaction parameter. As it goes to zero more
%explosive
L=l_0*exp(1/e);
z=e*b; %for stiffest use z=5e-3;
y0=[1,a,z]; % theta(0),y(0),z(0)
y=y0;

F=@(t,y)[rhs(t,y,e,L)]';
Fepic = @(y) F(y(end),y')';

NumSteps=round(Tend/Delt);
options=odeset('RelTol', 1e-10,'AbsTol',1e-10);
EpicOption.N_step=1;
EpicOption.RHS_function=F;
h=Delt;
start=tic;
for i=1:NumSteps
%t=y(i,end);
yIn=y(i,:);
%[tt,yUp]=ode15s(F, [yIn(end),yIn(end)+h],yIn,options);
epic(Exp_Int_Leja,Fepic,[y(end),y(end)+h],y',EpicOption);
yUp=yUp(end,:);
for j=1:3
    if(yUp(j)<0)
        yUp(j)=0;
    end
end
y(i+1,:)=yUp;
  fprintf(['Completed time step:', num2str( (i)*h), '\n'])
end
TIME=toc(start);
y=y';
time=0:h:Tend;
figure(1)
plot(time,y,"LineWidth",1.5,"MarkerSize",15)
xlabel("Time")
ylabel("Non-Dimensional quantity")
title('states over time')
legend({'Temp','y','z'});

for i=1:length(y)
   J=jac(y(end,i),y(1:3,i),e,L);
   eig_(i,:)=sort(eig(J));
end

figure(2)
semilogy(time,abs(real(eig_)),"LineWidth",1.5,"MarkerSize",15);
title('real part of eigs over time')
xlabel('time')
ylabel('real part')
legend({'E1','E2','E3'});

%figure(3)
%plot(time,imag(eig_));
%xlabel('time')
%ylabel('complex part')
%title('complex part of eigs over time')
%legend({'E1','E2','E3'});
%%

function yd=rhs(t,y,e,L)

yd=zeros(size(y)-1);

omega=L*y(2).*y(3).*exp(-1./(e*y(1)));
yd(1)=y(3);
yd(2)=-omega;
yd(3)=-yd(2)-yd(1);


% y = [theta,y,z]
function J=jac(~,y,e,L)

n=length(y);
J=zeros(n,n);

%omega=l*y(2).*y(3).*exp(-1./(e*y(1)));

J(1,3)=1;
J(2,1)=-(L*y(2)*y(3)*exp(-1/(e*y(1))))/(e*y(1)^2);
J(2,2)=-L*y(3)*exp(-1/(e*y(1)));
J(2,3)=-L*y(2)*exp(-1/(e*y(1)));
J(3,1)=-J(2,1);
J(3,2)=-J(2,2);
J(3,3)=-J(2,3)-1;

*/