#include "barnesTree.cpp"

int main(int argc, char **argv)
{
        if (argc < 3)
        {
                cout << "\n wrong usage";
                cout<<"\n Correct usage ./Barneshut <no_of_particles> <no_of_steps>";
                exit(0);
        }

		int test;
        #pragma omp parallel
        {
                #pragma omp master
                {
                        #ifdef _OPENMP
                        //printf("\nnum of threads = %d\n",omp_get_num_threads());
                        test = omp_get_num_threads();
                        #endif
                }
        }
		cout<<"\n num of threads = "<<test;
        long int noOfParticle = atoi(argv[1]);
        int noOfSteps = atoi(argv[2]);
		double startTime, endTime;

        barnesTree *tree = new barnesTree(noOfParticle, noOfSteps, THETA);

		startTime = omp_get_wtime();
	    //Build the tree now;
        bool retFlag = tree->run();
		endTime = omp_get_wtime();

        if (retFlag == false)
        {
                cout << "\n failed to run the barnes hut tree";
                exit(1);
        }
		else
		{
			cout<<"\nOutput->";
			cout<<"\n no of particles = "<< noOfParticle<<"\n no of steps = "<<noOfSteps;
			cout.precision(15);
			cout<<"\n Total Time taken = "<<endTime - startTime;
		}
        return 0;
}
