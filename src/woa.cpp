#include <iostream>
#include <woa_config.hpp>


// set the precision of the solution
#define CONSTRAINS 10
#define OUT_PRECISION 3
#define FUNCTION  levi
#define MAX_MIN_STATUS MIN_FIT



//schaffer function，constrains=-10~10，minimun f（0,0）=0, max -f=0
double schaffer(double px,double py)
{
    double numer=pow(sin(sqrt(px*px+py*py)),2)-0.5;
    double denom=pow(1+0.001*(px*px+py*py),2);
    double fitness=0.5+numer*(1.0/denom);
    return fitness;
}


//eggholder function，constrains=512，minimun f(512, 404.2319)=-959.6407
double eggholder(double px,double py)
{
    double py1=py+47.0;
    double t1=(-1.0)*py1*sin(sqrt(fabs(px/2.0+py1)));
    double t2=(-1.0)*px*sin(sqrt(fabs(px-py1)));
    double fitness=t1+t2;
    return fitness;
}


//booth function,constrains=10，minimun f(1, 3)=0
double booth(double px,double py)
{
    double t1=pow(px+2.0*py-7,2);
    double t2=pow(2.0*px+py-5,2);
    double fitness=t1+t2;
    return fitness;
}


//matyas function,constrains=10，minimum f(0, 0)=0
double matyas(double px,double py)
{
    double t1=0.26*(px*px+py*py);
    double t2=(-0.48*px*py);
    double fitness=t1+t2;
    return fitness;
}

//cross_in_tray function,constrains=10，minimum f(1.34941, -1.34941)=-2.06261,f(1.34941, 1.34941)=-2.06261,f(-1.34941, 1.34941)=-2.06261,f(-1.34941, -1.34941)=-2.06261
double cross_in_tray(double px,double py)
{
    double a=exp(fabs(100.0-(sqrt(px*px+py*py)/M_PI)));
    double b=fabs(sin(px)*sin(py)*a)+1;
    double fitness=-0.0001*pow(b,0.1);
    return fitness;
}


//levi function,constrains=10，minimum f(1,1)=0.0
double levi(double px,double py)
{
    double t1=pow(sin(3.0*M_PI*px),2);
    double t2=(pow(px-1,2)*pow(1+sin(3.0*M_PI*py),2));
    double t3=(pow(py-1,2)*pow(1+sin(2.0*M_PI*py),2));
    double fitness=t1+t2+t3;
    return fitness;
}


double fitness_cal(double px,double py, double (*pf)(double px,double py))
{
	return  (*pf)(px,py);
}
 


int main(int argc, char** argv)
{

    std::vector<double> constrain(4);
    constrain[0]=-1*CONSTRAINS;
    constrain[1]=CONSTRAINS;
    constrain[2]=-1*CONSTRAINS;
    constrain[3]=CONSTRAINS;
    woa::woa woa_a(2.0,0.5,50,20,constrain,OUT_PRECISION,MAX_MIN_STATUS);
    std::cout<<woa_a.nsols;

   for(size_t k=0;k<woa_a.ngens;k++)
    {
        std::cout<<"cur_genration:"<<k<<"\t"<<"a:"<<woa_a.a<<std::endl;

        for(size_t p=0;p<woa_a.nsols;p++)
        {
            std::cout<<"p"<<p<<std::endl;
            woa_a.solution[p].fitness=fitness_cal(woa_a.solution[p].px,woa_a.solution[p].py,FUNCTION);
            std::cout<<woa_a.solution[p].px<<"\t"<<woa_a.solution[p].py<<"\t"<<woa_a.solution[p].fitness<<"\t"<<std::endl;
        }
        woa_a.optimize();
    }

   for(size_t p=0;p<woa_a.nsols;p++)
   {
      woa_a.solution[p].fitness=fitness_cal(woa_a.solution[p].px,woa_a.solution[p].py,FUNCTION);
      std::cout<<woa_a.solution[p].px<<"\t"<<woa_a.solution[p].py<<"\t"<<woa_a.solution[p].fitness<<"\t"<<std::endl;;
   }

    for(size_t m=0;m<woa_a.best_solution.size();m++)
    {
        std::cout<<woa_a.best_solution[m].px<<"\t"<<woa_a.best_solution[m].py<<"\t"<<woa_a.best_solution[m].fitness<<std::endl;
    
    }
   return 0;

}






