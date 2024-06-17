#include "funcs.c"

int main(argc, argv)

    int argc;
char *argv[];

{
    PreprocData *pd = walloc(1,sizeof(PreprocData));
    
    pd->P = 100; 
    pd->total = 1; 
    pd->LEVEL_NUM = 10; 

    loadGraphData(pd); 

    initPreprocData(pd); 

    preProc(pd); 

    calcuRandomPairs(200,pd); 
  
    exit(0);

}
