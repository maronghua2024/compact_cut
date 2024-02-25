#include "funcs.c"

int main(argc, argv)

    int argc;
char *argv[];

{
    PreprocData *pd = walloc(1,sizeof(PreprocData));
    //300是mode 3, 100 以内是混合mode 1 和 2
    pd->P = 100; // P% percent of total passes in mode 1, the remaining in mode 2
    pd->total = 1; //number of total passes
    pd->LEVEL_NUM = 10; 

    loadGraphData(pd); //load graph data from standard input

    initPreprocData(pd); //init data structure

    preProc(pd); // preproc by traversing the graph for pd->total times

    calcuRandomPairs(200,pd); // randomly choose 100 node pairs and calcu their min-cut and output
  
    exit(0);

}
