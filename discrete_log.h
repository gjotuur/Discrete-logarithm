#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

//We work with multiplicative groups here - so we create a structure to save them (generator, module - prime number, order = p-1)
//Since the size of such a struct is comparably small, we can use just stack, without direct memory allocation
typedef struct 
{
    uint64_t module;
    uint64_t generator;
    uint64_t order;

}AlgebraicGroup;


void init_group(AlgebraicGroup* group, uint64_t module, uint64_t generator);


int main(){
    AlgebraicGroup z_5;
    init_group(&z_5, 5, 2);

    printf("\nGroup is: Z%llu*\nChosen generator is: %llu\nGroup order is: %llu", z_5.module, z_5.generator, z_5.order);
    return 0;
}

//Simple _init through pointer
void init_group(AlgebraicGroup* group, uint64_t module, uint64_t generator){
    group->generator = generator;
    group->module = module;
    group->order = module-1;
}