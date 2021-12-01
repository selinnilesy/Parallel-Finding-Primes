#include <iostream>
#include <vector>
#include <cmath>
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[]) {
    int rank, i, n = atoi(argv[1]), size;
    int numbers[n];
    bool isPrime[n];
    double t1, t2;
    //Since 1 and 2 are outliers, we defined them at first.
    isPrime[0]=0;
    isPrime[1]=1;
    numbers[0]=1;
    numbers[1]=2;
    //We eliminated even numbers
    for( i=2; i<n; i++){
        if((i%2)==1){
            isPrime[i]=0;
        }else{
            isPrime[i]=1;
        }
        numbers[i]=(i+1);
    }

    MPI_Status status;
    MPI_Init(&argc, &argv); 
	MPI_Comm_rank( MPI_COMM_WORLD, &rank); // rank = process id 
    MPI_Comm_size(MPI_COMM_WORLD, &size); // size = total number of processes

    i = ceil(sqrt(n)); // Until k*k > n.

    int k=1;
        t1 = MPI_Wtime();      
        while(1){
            if(rank==0){    //master
            int first_unm;
            // find first unmarked k, broadcast_all. k=2
            for( first_unm=k+1; first_unm<i ;first_unm++ ){
                if(isPrime[first_unm-1]) break;
            }
            k= first_unm;
            for(int second_unm=2*k-1; second_unm<i; second_unm=second_unm+k){
                isPrime[second_unm]=false;
            }
            MPI_Bcast(&k, 1,  MPI_INT, 0,  MPI_COMM_WORLD) ;
            if(first_unm==i) break;
            }
            else{
                MPI_Bcast(&k, 1, MPI_INT, rank, MPI_COMM_WORLD);
                bool flag = false;
                for (int j = i+((n-i+1)* (rank-1) / (size-1));     j < i+((n-i+1)* rank / (size-1)); j=j+k){
                    if(flag){
                        isPrime[j]  = false;
                    }
                    else{
                    int mod = numbers[j]%k;
                        if(mod){
                            int fwd_index = j + (k - mod);                                   // new index (forward) to find the first false-prime
                            if(fwd_index < i-1+((n-i+1)* rank / (size-1)) && fwd_index < n ) j = fwd_index;            // beware of the range before changing i.
                            else break;
                            if(isPrime[j]==false){
                                continue;
                            }
                            if(!(numbers[j] % k) ) {
                                flag=true;
                                isPrime[j]  = false;
                            }
                        }
                        else{
                          isPrime[j] = false;
                          flag = true;
                        }
                    }
                }
                if(k==i) break;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        t2 = MPI_Wtime();
        if(rank==0){
            for (int t1 = 1; t1 < size; t1++)
            {
                
                bool copyPrime[n];
                
                int finish, start;
                MPI_Recv(&start, 1, MPI_INT, t1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&finish, 1, MPI_INT, t1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&copyPrime, n, MPI_C_BOOL, t1, 3, MPI_COMM_WORLD, &status);
                
                
                for(int t2 = start; t2 < finish; t2++){
                    isPrime[t2]=copyPrime[t2];
                }
            }

            for(int i=0; i<n; i++){
                cout << "Number " << i+1 << " prime bool value is: " << isPrime[i] << "rank:" << rank <<endl;
            }

            cout << t2-t1 <<endl;
            
        }else{
            int finish = i+((n-i+1)* rank / (size-1));
            int start = i+((n-i+1)* (rank-1) / (size-1));
            MPI_Send(&start, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&finish, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		    MPI_Send(&isPrime, n, MPI_C_BOOL, 0, 3, MPI_COMM_WORLD);
        }
    MPI_Finalize();
    
    return 0;
}
