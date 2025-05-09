#include <omp.h>
#include <stdio.h>



int main(){

    #pragma omp parallel
    {

        for(int i = 0 ; i < 2 ; i++){
            printf("i is %d\n", i);
        }
    }

}
