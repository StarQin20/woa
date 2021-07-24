#ifndef WOA_CONFIG_HPP
#define WOA_CONFIG_HPP

#include <stdlib.h>
#include <math.h>
#include <vector>
#include <time.h>


#define MAX_FIT 0
#define MIN_FIT 1


namespace woa
{

 //Define the struct of the agent
 struct sol_type
{
    double px;
    double py;
    double fitness;
};

//Vector operation class under eigen library is not included
class dvector
{
    public:
    double v1;
    double v2;
    dvector()
    {}
    dvector(double x,double y)
    {
        v1=x;
        v2=y;
    }
    ~dvector()
    {
    }
    
    //Numerical approximation
    void round(int precision)
    {
        
        v1=((int)(v1*pow10(precision)+0.5))/pow10(precision);
        v2=((int)(v2*pow10(precision)+0.5))/pow10(precision);
    }

    //Vector numbers are multiplied separately
    void vmulti(dvector r)
    {
        v1=v1*r.v1;
        v2=v2*r.v2;
    }

    void nmulti(double n)
    {
        v1=v1*n;
        v2=v2*n;
    }

    void vminus(dvector r)
    {
        v1=v1-r.v1;
        v2=v2-r.v2;
    }

    void nminus(double n)
    {
        v1=v1-n;
        v2=v2-n;
    }

    //Find the modulus of a vector
    double normlize()
    {
        return sqrt(v1*v1+v2*v2);
    }

};


class woa
{
private:
    
public:
    //a is the search range, b is the helix parameter, nsols is the population number, ngens is the genration number, a_step is the search interval
    double a;
    double b;
    int nsols;
    int ngens;
    double a_step;
    //store the best id of current population
    int cur_best_id;
    int precision;
    bool maxmin;

    //constarin is a range constraint and solution is a population
    std::vector<double> constrain;
    std::vector<sol_type> solution;

    //best_ solution stores the optimal values of each generation
    std::vector<sol_type> best_solution;
    woa(double an,double bn,int sols,int gens,std::vector<double> con,int pre,bool maxmin_status);
    ~woa();

    void initSols(void);
    void getMaxMinSol(void);
    void optimize(void);

    dvector generateA(void);
    dvector generateC(void);
};



woa::woa(double an,double bn,int sols,int gens,std::vector<double> con,int pre,bool maxmin_status)
{
    a=an;
    b=bn;
    nsols=sols;
    ngens=gens;
    a_step=an/ngens;
    cur_best_id=-1;
    constrain.clear();
    for(size_t k=0;k<=con.size();k++)
    {
        constrain.push_back(con[k]);
    }
    best_solution.clear();
    precision=pre;
    maxmin=maxmin_status;
    initSols();
}

void woa::initSols(void)
{
    solution.clear();
    sol_type tmp_sol;
    srand( (unsigned)time( NULL ) );
    for(size_t k=0;k<nsols;k++)
    {
        tmp_sol.px=((rand()%(int((constrain[1]-constrain[0])/pow(0.1,precision))))*pow(0.1,precision)+constrain[0]);
        tmp_sol.py=((rand()%(int((constrain[3]-constrain[2])/pow(0.1,precision))))*pow(0.1,precision)+constrain[2]);
        //set the fitness to 0
        tmp_sol.fitness=0;
        solution.push_back(tmp_sol);
    }
}

void woa::getMaxMinSol(void)
{
    double tmp_fitness=solution[0].fitness;
    if(maxmin=MAX_FIT)
    {
        for(size_t k=0;k<nsols;k++)
        {
            if(solution[k].fitness>tmp_fitness)
            {
                cur_best_id=k;
                tmp_fitness=solution[k].fitness;
            }
        }
    }
    else if(maxmin=MIN_FIT)
    {
        for(size_t k=0;k<nsols;k++)
        {
            if(solution[k].fitness<tmp_fitness)
            {
                cur_best_id=k;
                tmp_fitness=solution[k].fitness;
            }
        }
    }

    if(cur_best_id!=-1)
    {
        best_solution.push_back(solution[cur_best_id]);
        std::cout<<"cur_best"<<"\t"<<solution[cur_best_id].px<<"\t"<<solution[cur_best_id].py<<"\t"<<solution[cur_best_id].fitness<<std::endl;
    }
    else
    {
        std::cout<<"wrong id number!";
    }

}


void woa::optimize(void)
{
    //find the best first
    getMaxMinSol();
    sol_type cur_best_sol=best_solution[best_solution.size()-1];
    dvector best_sol(cur_best_sol.px,cur_best_sol.py);
    dvector cur_sol(0,0);
    dvector result_sol(0,0);
    dvector random_sol(0,0);
    double poss;
    double tmp_rand1=0;
    double tmp_rand2=0;

    //generate A number
    dvector A;
    double A_norm;
    dvector C;
    double D;

    //For each individual
    for(size_t k=0;k<nsols;k++)
    {
        //Skip the best value
        if(k==cur_best_id)
            continue;
        //for the non-optimal value
        else
        {
            //Fetch individual value
            cur_sol.v1=solution[k].px;
            cur_sol.v2=solution[k].py;
            
            poss=rand()*1.0f/(RAND_MAX);
            //if the probability is greater than 0.5, it is the encirclement stage
            if(poss>0.5)
            {
                A=generateA();
                A_norm=A.normlize();
                //if the modulus of A is less than 1, it is considered exploitation
                if(A_norm<1.0)
                {
                    C=generateC();
                    C.vmulti(best_sol);
                    C.vminus(cur_sol);
                    D=C.normlize();
                    A.nmulti(D);
                    result_sol.v1=best_sol.v1-A.v1;
                    result_sol.v2=best_sol.v2-A.v2;
                }
                //if the module length of a is greater than 1, it is considered exploration
                else
                {
                    tmp_rand1=(rand()%nsols);
                    random_sol.v1=solution[tmp_rand1].px;
                    random_sol.v2=solution[tmp_rand1].py;
                    C=generateC();
                    C.vmulti(random_sol);
                    C.vminus(cur_sol);
                    D=C.normlize();
                    A.nmulti(D);
                    result_sol.v1=random_sol.v1-A.v1;
                    result_sol.v2=random_sol.v2-A.v2;
                }
            }
            //if the probability is less than 0.5, it is the attack stage
            else
            {
                    D=sqrt((best_sol.v1-cur_sol.v1)*(best_sol.v1-cur_sol.v1)+(best_sol.v2-cur_sol.v2)*(best_sol.v2-cur_sol.v2));
                    tmp_rand1=2.0*(rand()*1.0f/(RAND_MAX))-1.0;
                    tmp_rand2=2.0*(rand()*1.0f/(RAND_MAX))-1.0;
                    result_sol.v1=D*exp(b*tmp_rand1)*cos(2.0*M_PI*tmp_rand1)+best_sol.v1;
                    result_sol.v2=D*exp(b*tmp_rand2)*cos(2.0*M_PI*tmp_rand2)+best_sol.v2;
            }
            //Population renewal and constraints
            result_sol.v1=result_sol.v1<constrain[0]?constrain[0]:result_sol.v1;
            result_sol.v1=result_sol.v1>constrain[1]?constrain[1]:result_sol.v1;
            result_sol.v2=result_sol.v2<constrain[2]?constrain[2]:result_sol.v2;
            result_sol.v2=result_sol.v2>constrain[3]?constrain[3]:result_sol.v2;
            result_sol.round(precision);

            solution[k].px=result_sol.v1;
            solution[k].py=result_sol.v2;
        }
        
    }
       a=a-a_step;
}



 dvector woa::generateA(void)
 {
     double rand1=rand()*1.0f/(RAND_MAX);
     double rand2=rand()*1.0f/(RAND_MAX);

     dvector r(rand1,rand2);
     //A=2*a*r-a
     r.nmulti(a);
     r.nmulti(2.0);
     r.nminus(a);
     return r;
 }

dvector woa::generateC(void)
{
    double rand1=rand()*1.0f/(RAND_MAX);
    double rand2=rand()*1.0f/(RAND_MAX);

    dvector r(rand1,rand2);
    r.nmulti(2.0);
    return r;
}


woa::~woa()
{
}


}
#endif
