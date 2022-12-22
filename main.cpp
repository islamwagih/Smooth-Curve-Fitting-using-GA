#include <iostream>
#include <vector>
#include <random>
#include <time.h>
#include <algorithm>
using namespace std;

//random double between min and max
double randomBetween(double minD = -10, double maxD = 10)
{
    double random = (double)rand() / RAND_MAX;
    return minD + random * (maxD - minD);
}

//point representation
struct Point
{
    double x, y;
    Point(double x=0.0, double y=0.0):x(x), y(y){}
};

//chromosome representation
struct Chromosome
{
    vector<double> genes;
    double fitness;
    Chromosome(int degreeOfPolynomial)
    {
        genes = vector<double>(degreeOfPolynomial + 1);
        for (int i = 0; i < degreeOfPolynomial + 1; i++)
        {
            //return a random double from -10 to 10
            genes[i] = randomBetween();
        }
        fitness = 0;
    }
};



double difference(Chromosome* chrom, Point* point)
{
    double result = 0.0;
    for (int i = 0; i < chrom->genes.size(); i++)
    {
        double power = 1;
        for (int j = 0; j < i; j++)
        {
            power *= point->x;
        }
        result += (chrom->genes[i] * power);
    }
    result -= point->y;
    return result * result;
}

//goal is to minimize the value
double fitness(Chromosome* chrom, vector<Point>& points)
{
    double fit = 0.0;
    for (Point point : points)
    {
        //difference between y of the point and y of the curve
        fit += difference(chrom, &point);
    }
    return fit;
}

vector<Chromosome*> crossover(Chromosome* p1, Chromosome* p2)
{
    int size = p1->genes.size();
    int edge1 = (rand() % size) + 1; //1 -> n
    if (edge1 == size) { edge1--; } //decrease if it's n
    int edge2 = (rand() % size) + 1; //1 -> n
    if (edge2 == size) { edge2--; } //decrease if it's n
    //handle if edge1 == edge2
    while (edge1 == edge2)
    {
        edge1 = (rand() % size) + 1; //1 -> n
        if (edge1 == size) { edge1--; } //decrease if it's n
    }
    //handle if edge1 is bigger than edge2
    if (edge1 > edge2) { swap(edge1, edge2); }
    //new children
    Chromosome* firstSpring = new Chromosome(size - 1);
    Chromosome* secSpring = new Chromosome(size - 1);
    //apply 2 points crossover
    for (int i = 0; i < edge1; i++)
    {
        firstSpring->genes[i] = p1->genes[i];
        secSpring->genes[i] = p2->genes[i];
    }
    for (int i = edge1; i <= edge2; i++)
    {
        firstSpring->genes[i] = p2->genes[i];
        secSpring->genes[i] = p1->genes[i];
    }
    for (int i = edge2 + 1; i < size; i++)
    {
        firstSpring->genes[i] = p1->genes[i];
        secSpring->genes[i] = p2->genes[i];
    }

    vector<Chromosome*> childs(2);
    childs[0] = firstSpring;
    childs[1] = secSpring;

    return childs;
}

/*
Generate random number ri1 ϵ [0, 1]

y = ∆L if ri1 ≤ 0.5

y = ∆U if ri1 > 0.5

Let ∆(t,y)

= value of mutation at generation t

= y(1-r(1-t/T)^b)

where:

r = random number ϵ [0, 1]

t = current generation

T = maximum number of
generations

b = dependency factor ≈ 1…5
*/

//to do
Chromosome* mutation(Chromosome* chrom, double mutationChance, int currentGeneration, int maximumGenerations)
{
    double random = randomBetween(0.0, 1.0);
    if (random <= mutationChance)
    {
        //do a mutation
        for (int i = 0; i < chrom->genes.size(); i++)
        {
            double r1 = randomBetween(0.0, 1.0);
            int y;
            if (r1 <= 0.5) { y = -10; }
            else { y = 10; }

            int dependencyFactor = (rand() % 5) + 1;
            double delta = 1.0 - (currentGeneration / maximumGenerations);
            delta = pow(delta, dependencyFactor);
            double r2 = randomBetween(0.0, 1.0);
            delta = y*(1.0 - pow(r2, delta));

            if (r1 <= 0.5)
            {
                chrom->genes[i] -= delta;

            }
            else
            {
                chrom->genes[i] += delta;

            }
        }

    }
    return chrom;
}

//selection strategy
int tournamentSelection(vector<Chromosome*>& chroms, vector<Point>& points)
{
    //two random indicies
    int first = rand() % chroms.size();
    int sec = rand() % chroms.size();
    while (first == sec)
    {
        first = rand() % chroms.size();
    }
    //compare the fitness of the selected two
    double firstFit = chroms[first]->fitness;
    double secFit = chroms[sec]->fitness;

    //return the index of the best
    return (firstFit > secFit) ? sec : first;
}

// smalles in index of 0
bool comparator(Chromosome* a, Chromosome* b)
{
    return a->fitness < b->fitness;
}


int main()
{
    //initialize a random seed
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    srand(time(NULL));
    int testCases, maxNumOfGenerations = 1e3;
    double mutationChance = 0.2;
    cin >> testCases;
    int t = 1;
    while(testCases--)
    {

        int numberOfPoints, degree, currGeneration = 0;
        cin >> numberOfPoints >> degree;
        vector<Point> points;
        while(numberOfPoints--)
        {
            double x, y;
            cin >> x >> y;
            points.push_back(Point(x, y));
        }

        int populationSize = 500;
        vector<Chromosome*> population(populationSize);
        for(int i=0;i<populationSize;i++)
        {
            population[i] = new Chromosome(degree);
            population[i]->fitness = fitness(population[i], points);
        }

        while(currGeneration <= maxNumOfGenerations)
        {
            vector<Chromosome*> childs;
            for(int i=0;i<populationSize;i+=2)
            {
                //selection
                int fatherInd = tournamentSelection(population, points);
                int motherInd = tournamentSelection(population, points);
                while(fatherInd == motherInd)
                {
                    fatherInd = tournamentSelection(population, points);
                }

                //crossover & mutaion
                vector<Chromosome*> currChilds = crossover(population[fatherInd], population[motherInd]);
                for(Chromosome* child:currChilds)
                {
                    childs.push_back(mutation(child, mutationChance, currGeneration, maxNumOfGenerations));
                    child->fitness = fitness(child, points);
                }
            }


            for(Chromosome* child:childs)
            {
                population.push_back(child);
            }

            sort(population.begin(), population.end(), comparator);
            int newSize = population.size();
            for(int i=newSize-1;i>=populationSize;i--)
            {
                delete population[i];
                population.pop_back();
            }
            currGeneration++;
        }

        //the best chromosome at index 0
        cout<<"test case "<<t++<<"# best chromosome: ";
        for(int i=0;i<population[0]->genes.size();i++)
        {
            cout<<population[0]->genes[i]<<' ';
        }
        cout<<"fitness "<<population[0]->fitness<<endl;

        //free the memory
        for(int i=0;i<populationSize;i++)
        {
            delete population[i];
        }
    }


    return 0;
}
