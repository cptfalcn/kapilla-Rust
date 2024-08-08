/*
ToDo (in no particular order):
1)  Implement a command line parser for the variables (which will be renamed) a,b,e, L_0, Tend, Delt
2)  Design the structure of the Gas phase problem
    a)  Should have a Rhs/Jac function (the equivalent of a virtual function)
    b)  Should have a solver (will start with RK4 and implement more)
    c)  Should have visualization (see the Crate "plotters")
3)  Refactor to allow for more complex chemical mechanisms (maybe a yaml plugin?)
*/


//Everyone should have an some time integration capabilities
//All the pertinent functions should be supplied by the impl of the structure
//rhs:  The right-hand-side forcing term function f(vec)=vec
//jac:  the Jacobian of the rhs function
//

/*
These public functions will be used to help perform the basic vector operations needed by the stepper
*/
pub fn elementwise_addition(vec_a: Vec<f32>, vec_b: Vec<f32>) -> Vec<f32> {
    vec_a.into_iter().zip(vec_b).map(|(a, b)| a + b).collect()
}

pub fn elementwise_scale(mut vec_a: Vec<f32>, scale: f32) -> Vec<f32> {
    for i in vec_a.iter_mut(){
        *i *= scale;
    }
    return vec_a;
}


/*==============================================
A time stepper needs to be implemented.
It should return the solution, in case a specific outside "class" needs it.
=================================================*/
pub trait IntegratorStep{
    fn step(&self) -> Vec<f32>;
}


        
/*==============================================
A Right hand side function should be implemented
It should return the result as a vector.
=================================================*/
pub trait RightHandSide{
    fn rhs(& self, time : f32, state: Vec<f32>)->Vec<f32>;
}

/*==============================================
A jacobian calculation should be implemented
=================================================*/
pub trait JacobianMatrix{
    fn jac(&self) -> Vec<f32>;
}





//==================================================================
//The Kapila Reaction kinetics problem structure and implementations
//==================================================================
struct KapProb{
    state       : Vec<f32>,
    end_time    : f32,  //Final integration time.
    start_time  : f32,  //The time we are starting at.
    curr_time   : f32,  //This is the "t" variable 
    delta_t     : f32   //The time step size
}

//An implementation to print the state.
impl KapProb{

    fn print_state(&self){
        println!("The vector is {:?}", &self.state);
    }
}


impl RightHandSide for KapProb{
    fn rhs(&self, time: f32, state : Vec<f32>)-> Vec<f32>{
        return elementwise_scale(self.state.clone(),2.0);
    }
}

pub trait Integrate{ 
    fn time_integration(& mut self)->Vec<f32>;
}

//Make it so various stepping methods can be implemented
pub trait RK2Step{
    fn rk2_step(& self)->Vec<f32>;
}

pub trait RK4Step{
    fn rk4_step(& self)->Vec<f32>;
}


/*==============================================
This implements the time integration for the problem.
The following are implemented for the problem:
rhs
start_time
end_time
delta_t
step
=================================================*/
impl Integrate for KapProb{

    fn time_integration(& mut self)->Vec<f32>{
    //Set a loop
    let step_max = ((self.end_time- self.start_time)/self.delta_t) as i32;
    println!("Performing {:?} time interation steps", &step_max);
    println!("The starting time is: {:?}", self.start_time);
    println!("The ending time is: {:?}", self.end_time);
    println!("Starting the time integration procees, the state is {:?}", &self.state);
    let mut step_ct = 0;
    while step_ct < step_max {
        self.state = self.rk4_step();
        step_ct += 1;
        self.curr_time+=self.delta_t;
    }
    println!("At the end of the simulation we are at time {:?}", self.curr_time);
    return self.state.clone();
    }
}

//Here are various step options, which can be called by the IntegratorStep
impl RK2Step for KapProb{
    fn rk2_step(& self) -> Vec<f32>{
        return elementwise_addition(self.state.clone(), 
        elementwise_scale(  self.rhs(self.curr_time.clone(), self.state.clone()) , self.delta_t) );
    }
}

//RK4

impl RK4Step for KapProb{
    fn rk4_step(& self) -> Vec<f32>{
        let k1 =    self.rhs(self.curr_time.clone(), self.state.clone());
        
        let y2 =    elementwise_addition(self.state.clone(), elementwise_scale(k1.clone(), self.delta_t/2.0));
        let k2 =    self.rhs(self.curr_time.clone()+ self.delta_t/2.0, y2);
        
        let y3 =    elementwise_addition(self.state.clone(), elementwise_scale(k2.clone(), self.delta_t/2.0));
        let k3 =    self.rhs(self.curr_time.clone() + self.delta_t/2.0, y3);

        let y4 =    elementwise_addition(self.state.clone(), elementwise_scale(k3.clone(), self.delta_t));
        let k4 =    self.rhs(self.curr_time.clone() + self.delta_t, y4);

        let k1k2 =  elementwise_addition(k1.clone(), elementwise_scale(k2.clone(),2.0) );
        let k3k4 =  elementwise_addition(k4.clone(), elementwise_scale(k3.clone(),2.0) );
        let ks   =  elementwise_addition(k1k2.clone(), k3k4.clone()); 
        let ans  =  elementwise_addition(self.state.clone(), elementwise_scale(ks, self.delta_t/6.0));
        return ans;
    }
    
}



/*=================
Main
=================*/

fn main() {
    println!("Hello, Testing ground!");

    let mut test_problem = KapProb{ state : vec![1.0,2.0,3.0], end_time : 1.0, 
                                    start_time : 0.0, curr_time : 0.0, delta_t : 0.1 };
    test_problem.time_integration();
    test_problem.print_state();
    

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