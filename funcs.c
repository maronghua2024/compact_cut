
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <values.h>
#include <math.h>
#include <time.h>

#include "types.h"  /* type definitions */
#include "parser.c" /* parser */
#include "timer.c"        /* timing routine */

#define MAX_LONG LONG_MAX

#define min(s, t) (s) < (t) ? (s) : (t)
#define max(s, t) (s) > (t) ? (s) : (t)













void *walloc(unsigned int num, unsigned int size)
{
  void *ptr = calloc(num, size);
  assert(ptr != NULL);
  return ptr;
}


cType mrand(RandomData* rd)
{
  return rd->randNums[rd->randNumIdx++];
}

RandomData* initrand(cType len)
{
  RandomData *rd = walloc(1, sizeof(RandomData));

  srand((int)(timer() * 1000));

  rd->randNums = (cType *)walloc(len, sizeof(cType));
  rd->randNumIdx = 0;

  for (int i = 0; i < len; i++)
  {
    rd->randNums[i] = (cType)rand();
  }

  return rd;
}


void HeapAdjustDown(sType *idx, edgeP * edges ,int start,int end)  
{  
    sType tempIdx = idx[start];  

    int i = 2*start+1;      
    
    
    
    while(i<=end)  
    {  
        if(i+1<=end && edges[idx[i+1]].tmp > edges[idx[i]].tmp )    
            i++;  

        if(edges[idx[i]].tmp <= edges[tempIdx].tmp )   
            break;  

        idx[start] = idx[i];

        start = i;  
        i = 2*start+1;  
    }  

    idx[start] = tempIdx;  
}  
  
void HeapSort(sType *idx, edgeP * edges, int len)  
{  

    int i;  
    for(i=(len-1)/2;i>=0;i--){  
        HeapAdjustDown(idx,edges,i,len-1);  
    }

    for(i=len-1;i>0;i--)
    {  
        
        sType temp = idx[i];  
        idx[i] = idx[0];  
        idx[0] = temp;  

        HeapAdjustDown(idx,edges,0,i-1);  
    }  

}  
  

void deOrderEdgeByRandomCap(nodeP *np,PreprocData *pd)
{

  cType cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    
    pedges[i].tmp = 1000-pedges[i].cap+1; 
  }

    assert(cnt<4 || idxs[2]!=idxs[3]);
    HeapSort(idxs,pedges,cnt);
    assert(cnt<4 || idxs[2]!=idxs[3]);  
}


void aOrderEdgeByAvgCV(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    long acv = pp->avgCV;
    if(acv == 0 ){
      pp->tmp = MAX_LONG;
    }
    else{
      pedges[i].tmp = mrand(pd->rd) % acv; 
    }
  }

    HeapSort(idxs,pedges,cnt);
}


void aOrderEdgeByDegree(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;
  nodeP* nodes = pd->gd->nodes;
  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    cType zn = pp->endNode;
    pedges[i].tmp =  (nodes+zn)->nIdx;
  }

    HeapSort(idxs,pedges,cnt);
}

























/*
  the function to add more data to a traversal tree to accelerate the searching in the tree
  The idea is to precalcuate minimal cv value of a span of nodes in a traversal tree, e.g., when SPAN_LEN = 100, and a node has depth of 200, then the algorithm will pre-calculate the minimal cv value of the nodes between the node (dep=200) and an acestor(dep=101)
  upid is the id of the last SPAN node, mcv is the min cv among all previous nodes in the recent SPAN
  lastDepMCV is the depth of the node depth that has the minimal cv in the span
  lastJointNodeId is the last ancestor node id that has more than one child nodes
  lastJointMCV is the cv of lastJoineNodeId

  改版后预处理算法要有变化：
      对每个节点的值，都要改进下考虑当前出发节点z的考虑下游分支的更小的cv'=cv2+-cof

  (1)处理时：
      在非段头节点中还要考虑父亲cv'的值，如果更小，则更新指向父亲
      跨段的cv2+-cof直接放到下段mcv中，因为solve中确认跨段了才使用当前段的mcv
  (2)求解时：
      看solve的注释


*/

long gRoot = 0;
void buildAcc(PreprocData *pd, cType curN, cType upid, cType upid_LC, long mcv, long* applyAdj, cType* depMap)
{
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);
  int cnt = (nodes + curN)->nIdx;
  edgeP *pedges = (nodes + curN)->edges;

  long curCV = pd->gpcv[curN];
  cType curDep = pd->gpdep[curN];

  assert(curDep == 0 || curCV >0);

  depMap[curDep] = curN;
  pd->gpaccup[curN] = upid; 

  int SPAN_LEN = pd->SPAN_LEN;
  int LN = pd->LEVEL_NUM;
  /******************************************/
  
  pd->gpaccup[curN] = upid; 
  pd->gpaccposmcv[curN] = upid_LC; 
  
  if(curDep % SPAN_LEN == 0){
    upid = curN;
    upid_LC = 0; 
  }
  else{
    

  }


  /****************************************/
  cType ances; 
  for (int i = 1; i <= pd->LEVEL_NUM; i++)
  {
    if (curDep  < (unsigned int)i)
    {
      break;
    }

    ances = depMap[curDep - i];
    assert(ances != 0);
    if (applyAdj[ances] == 0) 
    {
      if (pd->gcutadjsign[curN][i] == 1 ) 
      {
        
        applyAdj[ances] = pd->gcutadj[curN][i];
        pd->gcutadj[curN][i] = -1 * pd->gcutadj[curN][i];
      }
    }

  }

  
  if(curDep >= LN ){
      if(curDep % SPAN_LEN == LN){
        upid_LC = curN; 
        if(LN > 0){
          
          pd->gpaccposmcv[curN] = upid_LC; 
        }
      }
      else if(curDep % SPAN_LEN == (LN+1)){
        mcv = MAX_LONG;
      }

      assert(applyAdj[ances] <= 0);
      assert(pd->gcutadj[ances][0] <= 0);
      assert(pd->gpcv[ances]+pd->gcutadj[ances][0] > 0);
      
      ances = depMap[curDep-LN];
      
      
      
      pd->gpoh[curN] = pd->gpcv[ances] + pd->gcutadj[ances][0] - applyAdj[ances]; 
      mcv = min(mcv, pd->gpoh[curN]);
      pd->gpaccmcv[curN] = mcv; 

  }
  else{
    

  }


  while (cnt > 0)
  {
    cType zn = pedges->endNode;
    if (pd->gpfa[zn] == curN)
    {
      buildAcc(pd,zn, upid, upid_LC, mcv, applyAdj, depMap);
    }
    pedges++;
    cnt--;
  }

  
  for (int i = 1; i <= pd->LEVEL_NUM; i++)
  {
    if (curDep  < (unsigned int)i)
    {
      break;
    }
    cType ances = depMap[curDep - i];
    
    
    assert(ances != 0);
    
    if (pd->gcutadjsign[curN][i] == 1 && pd->gcutadj[curN][i] > 0)
    {
      pd->gcutadj[curN][i] = -1 * pd->gcutadj[curN][i];
      assert(applyAdj[ances] == pd->gcutadj[curN][i]);
      applyAdj[ances] = 0;
    }
  }

}

void calcuTotalCap(PreprocData *pd)
{
  nodeP *nodes = pd->gd->nodes;

  for (cType curN = 1; curN <= pd->gd->N; curN++)
  {
    nodeP *np = nodes + curN;
    edgeP *pedges = np->edges;
    int cnt = np->nIdx;
    np->totalCap = 0;
    for (int ni = 0; ni < cnt; ni++)
    {
      edgeP *eh = pedges + ni;
      np->totalCap += eh->cap;
    }
  }
}

/*

算法整体步骤：相对之前treem算法的改进
(1)预处理
	
	sum_f = 0
	对f的每个为访问的child n:
		在遍历n前，在f上设置oh_f并置0
		在遍历计算中，每次访问祖先，除了计算mc，还更新oh_f (加上就可以)
		n返回后
			此时知道mc_n, oh_f(oh_f就是n这一支连到f的边的值，要包含f_n父子的边)
			这时计算n这一支去掉对f的mc的影响cof_n = oh_f - (mc_n-oh_f) = 2* oh_f - mc_n
				
			如果cof_n是负值 
				就加到sum_f上,即sum_f 加上 负值，sum_f指的是f的所有减少割值的子n去掉，总共减少的割值
				如果不减少，这个n就不去掉
	
	得到mc2_f = mc_f + sum_f 
				
(2)计算时: solve 和 build的时候
	
	每次向上回溯，n回溯到f时，如果cof_n是负值，说明如果要保留n这一支，目前f最优割就包含n，即此时f用于计算的mc应该取(mc2_f + -cof_n)，即加上n这一支减掉的值
  现在问题来了：
    f的mc2和访问哪一支有关系，buildAcc预处理咋做？
    这样就意味着，不同的底层上来，每个节点的mc还不一样，导致 节点段 中最小值还不一样
    buildAcc记录的是向上的，所以可以记录的呀

*/


void markCut(cType curN, PreprocData *pd)
{
  
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);

  
  assert((nodes + curN)->nIdx > 0);

  short *curS = pd->gps + curN;

  assert(*curS == 0);

  long *curCV = pd->gpcv + curN;
  cType *curDep = pd->gpdep + curN;

  *curS = 1;
  *curCV = 0;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;

  if (pedges == NULL)
  {
    *curS = 2;
    return;
  }

  if (np->orderedEdges == NULL)
  {
    np->orderedEdges = (sType *)walloc(cnt + 1, sizeof(sType));
    for (int i = 0; i < cnt; i++)
    {
      np->orderedEdges[i] = i;
    }
  }

  long cap;
  sType *idxs = np->orderedEdges;

  if (pd->mode == 1)
  {
    deOrderEdgeByRandomCap(np,pd); 
  }
  else if (pd->mode == 2)
  {
    aOrderEdgeByAvgCV(np,pd);
  }
  else if (pd->mode == 3){
    
    aOrderEdgeByDegree(np,pd);
  }

  

  
  
  cType fa = curN;
  for(int i=1; i<=pd->LEVEL_NUM; i++ ){
    if(pd->gpdep[fa] == 0){
      break;
    }
    
    fa = pd->gpfa[fa];
    
    
    pd->gcutadj[curN][i] = pd->gpoh[fa]; 
  }


  for (int ni = 0; ni < cnt; ni++)
  {
    

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;

    assert(zn != 0);
    assert(zn != curN);

    short zs = pd->gps[zn];

    assert(!(pd->gpfa[zn] == curN && zs == 2));

    
    if (zs == 1) 
    {
      
      cap = eh->cap;
      *curCV += cap;
      pd->gpcv[zn] -= cap;
      pd->gpoh[zn] += cap;
    }
    else if (zs == 0) 
    {
      
      

      pd->gpfa[zn] = curN;
      pd->gpdep[zn] = *curDep + 1;

      
      markCut(zn,pd);

      
      
      assert(pd->gpdep[zn] == pd->gpdep[curN] + 1);
      
      *curCV += pd->gpcv[zn];




		 
         
      
      
      
      
      
      
      

    }
    else
    {
      
      
      assert(pd->gpdep[curN] < pd->gpdep[zn]);
      
    }

  }

  
  
  
  
  
  
  
  fa = curN;
  cType faBelow = 0; 
  int actual_ln = -1;
  for(int i=1; i<=pd->LEVEL_NUM; i++ ){
    
  
    if(pd->gpdep[fa] == 0){
      actual_ln = i-1;
      break;
    }

    fa = pd->gpfa[fa];
    
    pd->gcutadj[curN][i] = pd->gpoh[fa] - pd->gcutadj[curN][i];
    
    faBelow += pd->gcutadj[curN][i];
    
    
    pd->gcutadj[curN][i] = faBelow - (*curCV - faBelow);    
  
  }

  
  if(actual_ln < 0){
    actual_ln = pd->LEVEL_NUM;
  }

  
  

  
  
  
  

  


  
  
  
  long tempCumSum = 0;
  pd->gcutadj[curN][0] = 0; 
  for (int i = (actual_ln < pd->LEVEL_NUM ? actual_ln : actual_ln - 1); i >= 0; i--) 
  {
    
    tempCumSum = 0;

    for (int ni = 0; ni < cnt; ni++)
    {
      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      assert(zn != 0);
      assert(zn != curN);

      if (pd->gpfa[zn] != curN)
      {
        continue;
      }

      if (pd->gcutadj[zn][i + 1] < 0)
      {
        tempCumSum += pd->gcutadj[zn][i + 1];
      }
    }

      

    if (pd->gcutadj[curN][i] < tempCumSum) 
    {
      
      
      
      pd->gcutadjsign[curN][i] = 1; 
      
    }
    else
    {


      
      
      
      pd->gcutadj[curN][i] = tempCumSum;
      pd->gcutadjsign[curN][i] = 0; 
      
    }

  }


  
  
  
  
  if(*curCV == 0 && curN != gRoot){
    printf("cv=0 error: (cnt is %d) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n",cnt,curN, gRoot,*curCV,pd->gpdep[curN]);
    assert(1==2);
  }

  

  if(pd->mode == 1){

    for (int ni = 0; ni < cnt; ni++)
    {
      

      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      

      assert(zn != 0);
      assert(zn != curN);
      short zs = pd->gps[zn];
    

      
      if (zs == 1 && pd->gpdep[zn] != *curDep - 1)
      {
          cType weight = eh-> w;
          if(eh->avgCV == 0){
            eh->avgCV = MAX_LONG;
          }
          eh->avgCV = min(eh->avgCV, *curCV);
          eh->w = weight+1;
          
          edgeP *reh = eh->rev;

          if(reh->avgCV == 0){
            reh->avgCV = MAX_LONG;
          }

          weight = reh-> w;
          reh->avgCV = min(reh->avgCV, *curCV);
          reh->w = weight+1;

      }
      
    }
  }

    
  assert(pd->gpdep[curN] == 0 || *curCV >0);

  
  *curS = 2;

}


/*
 (2)求解时：
      深度大的出发节点直接cv2，
      另外一个不是这个的祖先时，另外一个也可以用cv2
      [这个先不优化有点复杂]如果t是s祖先，
        如果cv2正好割掉下面，直接可以用cv2
        如果不是，也可以计算割掉对应树后这个t的割

      PS: t是s祖先，可以用cv2,t不是s的祖先也可以cv2，两个cv2都可以用
      问题是，如果t是s祖先，t的cv2对应的割并不能切断s这条支线，这样就不对了
        如果t是s祖先，其实我们算的是去掉这个支后的t的最小割
      那就简化：
        只要t不是s祖先，就可以用t的cv2

      后面的逐个逼近，需要用cv'
      
*/

void traceUp2LN(NodePropArr* pnp, cType startDep, cType curN, cType curDep, long *adj, cType LN)
{

  for (int i = 1; i <= LN; i++)
  {
    if ( curDep < (unsigned int)i || (curDep+LN) <= (startDep+i) )
    {
      break; 
    }

    if (pnp->pacc_cut_adjust_sign[curN][i] == 1)
    {
      adj[startDep+i-(curDep)] = pnp->pacc_cut_adjust[curN][i]; 
    }
  }
}




long isOptimalExclude(NodePropArr* pnp, cType v, cType v_top, cType LN_v, cType *pDep, cType *pFa, cType LN){
  
  cType minDep = 0;
  
  while(v != LN_v){
    
    
    if(pDep[v] - pDep[LN_v] <= LN && pnp->pacc_cut_adjust_sign[v][pDep[v] - pDep[LN_v]]==1){
      
      minDep = pDep[v];    
    }
    v = pFa[v];
  }

  

  
  
  
  
  
  
  /*(2) LN_s(含)及之上的祖先，则需要看不含一方的最高点，如果不含的最高点在LN_s(不含)以下，可以认为不含，因为去掉t不影响s。
        如果在之上就麻烦了，去掉t也去掉s了，这个就不好办了

  */
  
  
  if(minDep > 0){
    if(minDep <= pDep[v_top]){
      return 2; 
    }
    else{
      return 1;
    }
  }

  return 0;
}



long solveMaxFlowAccVER4(cType root, long minCandi, NodePropArr* pnp, cType s, cType t, int SPAN_LEN, aType LN)
{
  cType *pDep = pnp->pdep;
  long *pCV = pnp->pcv;
  cType *pFa = pnp->pfa;
  long *poh = pnp->poh;
  cType *paccup = pnp->pacc_upid;
  cType *paccup_LC = pnp->pacc_pos_upmincv;
  long *paccmcv = pnp->pacc_upmincv;
  cType cycCnt = 0;
  
  assert(s != t);
  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  printf("\nstart: s %ld dep %ld, t %ld dep %ld\n",s,pDep[s],t,pDep[t]);

  assert(pDep[s] >= pDep[t]);

  long *adj = (long *)walloc(LN+1, sizeof(long));
  long *adj2 = (long *)walloc(LN+1, sizeof(long));
  memset(adj, 0, (LN+1) *  sizeof(long));
  memset(adj2, 0, (LN+1) *  sizeof(long));

  long mcv = MAX_LONG;
  cType LN_s = s , LN_t = t ;
  int alsoT = 0;

  if(pDep[LN_s] == pDep[LN_t]){
    alsoT = 1;
  }

  assert(pDep[LN_t] <= pDep[LN_s]);

  
    while (pDep[s] - pDep[LN_s] < LN) 
    {
		cycCnt ++;
      
      assert(pDep[LN_t] <= pDep[LN_s]);

      if (alsoT > 0)
      {
        assert(pDep[LN_t] == pDep[LN_s]);
        
        mcv = min(mcv, pCV[LN_t] + pnp->pacc_cut_adjust[LN_t][0] - adj2[pDep[t] - pDep[LN_t]]); 
        traceUp2LN(pnp, pDep[t], LN_t, pDep[LN_t], adj2, LN);
        LN_t = pFa[LN_t];
      }

      mcv = min(mcv, pCV[LN_s] + pnp->pacc_cut_adjust[LN_s][0] - adj[pDep[s] - pDep[LN_s]]);
      traceUp2LN(pnp, pDep[s], LN_s, pDep[LN_s], adj, LN);
      LN_s = pFa[LN_s];

      if (pDep[LN_s] == pDep[LN_t])
      {
        if (LN_s == LN_t)
        {
          
          goto end;
        }

        
        alsoT = 1;
      }

      
      
      
      
      
    }
    
    

    assert(pDep[LN_s] > 0);
    assert(pDep[s] - pDep[LN_s] == LN);
    assert(pDep[LN_t] > 0 || LN_t == root);
    assert(pDep[LN_t] <= pDep[LN_s]);

  
  
  

  
    while (1 == 1)
    {
      cycCnt ++;
      if (LN_s == LN_t)
      {
        goto end;
      }

      
      assert(pDep[LN_s] > 0 || LN_s == root);
      assert(LN_s > 0);
      assert(pDep[LN_t] <= pDep[LN_s]);

      if (pDep[paccup[LN_s]] >= pDep[LN_t])
      {
        
        mcv = min(mcv, paccmcv[s]);
      }
      else
      {
        break;
      }

      
      
      if (paccup_LC[LN_s] == 0)
      {
        assert(pDep[LN_s] % SPAN_LEN < LN);
        for (int i = 0; i < pDep[LN_s] - pDep[paccup[LN_s]]; i++)
        {
          s = pFa[s];
        }
      }
      else
      {
        s = paccup_LC[LN_s];
        
        assert(pDep[s] % SPAN_LEN == LN);
      }

      LN_s = paccup[LN_s];

  
  assert(pDep[s] - pDep[LN_s] == LN);

    }

  
  assert(pDep[LN_t] <= pDep[LN_s]);
  assert(pDep[paccup[LN_s]] < pDep[LN_t]);

  
  
  assert(pDep[LN_s] - pDep[LN_t] < SPAN_LEN);



  

  while (pDep[LN_s] > pDep[LN_t])
  {
	  cycCnt ++;
    mcv = min(mcv, poh[s]);
    s = pFa[s];
    LN_s = pFa[LN_s];
  }

  assert(pDep[LN_s] == pDep[LN_t]);
  if(LN_s == LN_t){
    goto end;
  }
  

  
  
  
  while(pDep[t] - pDep[LN_t] < LN)
  {
	  cycCnt ++;
      assert(pDep[LN_t] == pDep[LN_s]);
      if (LN_s == LN_t)
      {
        goto end;
      }
      mcv = min(mcv, pCV[LN_t] + pnp->pacc_cut_adjust[LN_t][0] - adj2[pDep[t] - pDep[LN_t]]);
      traceUp2LN(pnp, pDep[t], LN_t, pDep[LN_t], adj2, LN);
      LN_t = pFa[LN_t];      

      mcv = min(mcv, poh[s]);
      s = pFa[s];
      LN_s = pFa[LN_s];
  }

  if (LN_s == LN_t)
  {
    goto end;
  }

  assert(pDep[LN_t] == pDep[LN_s]);
  assert(pDep[t] - pDep[LN_t] == LN);
  
  assert(pDep[s] - pDep[LN_s] == LN);

  
  
  while(paccup[LN_s] != paccup[LN_t]){
    
	cycCnt ++;
    mcv  = min(mcv, paccmcv[s]);

    if(paccup_LC[LN_s] == 0){
      assert(pDep[LN_s] %SPAN_LEN < LN);
      for(int i=0; i<pDep[LN_s] - pDep[paccup[LN_s]]; i++){
        s = pFa[s];
      }
    }
    else{
      s = paccup_LC[LN_s];
      assert(pDep[s] %SPAN_LEN >= LN);
    }


    LN_s = paccup[LN_s];  

    mcv  = min(mcv, paccmcv[t]);
    if(paccup_LC[LN_t] == 0){
      assert(pDep[LN_t] %SPAN_LEN < LN);
      for(int i=0; i<pDep[LN_t] - pDep[paccup[LN_t]]; i++){
        t = pFa[t];
      }
    }
    else{
      t = paccup_LC[LN_t];
      assert(pDep[t] %SPAN_LEN >= LN);
    }

    LN_t = paccup[LN_t];  

  }

  assert(paccup[LN_s] == paccup[LN_t]);

  
  
  while(LN_s != LN_t){
	  cycCnt ++;
    mcv = min(mcv, poh[s]);
    s = pFa[s];
    LN_s = pFa[LN_s];

    mcv = min(mcv, poh[t]);
    t = pFa[t];
    LN_t = pFa[LN_t];        
  }

end:
  
  assert(LN_s == LN_t);
  printf("mcv : %ld, final root: (LN_s/LN_t is %ld, dep is %ld ) \n",mcv,LN_s,pDep[LN_s]);

  if(t == LN_t || s == LN_s){
    
    goto end2;
  }

  cType V = LN_s;
  
  
  

  
  for (; pDep[LN_s] - pDep[V] < LN; V = pFa[V])
  {
	  cycCnt ++;
    int inc_s = isOptimalExclude(pnp, s, LN_s,V, pDep, pFa,LN);
    int inc_t = isOptimalExclude(pnp, t, LN_t,V, pDep, pFa,LN);

    if(inc_s == 2 || inc_t == 2){
      
      continue;
    }

    
    if (inc_s + inc_t == 2)
    {
      
      continue;
    }

    if (inc_s + inc_t == 0)
    {
      
      continue;
    }

    
    if (inc_s == 1)
    {
      mcv = min(mcv, pCV[V] + pnp->pacc_cut_adjust[V][0]);
    }

    if (inc_t == 1)
    {
      mcv = min(mcv, pCV[V] + pnp->pacc_cut_adjust[V][0]);
    }

    if(pDep[V] == 0){
      break;
    }

    
  }

end2:
  free(adj);
  free(adj2);
  return mcv;
  

}














































































































































































void loadGraphData(PreprocData *pd){
  pd->gd = walloc(1,sizeof(GraphData));  

  parse(&(pd->gd->N), &(pd->gd->M), &(pd->gd->nodes));

  printf("c nodes:       %10ld\nc arcs:        %10ld\nc\n", pd->gd->N, pd->gd->M);
}


void initPreprocData(PreprocData *pd){
  pd->rd = NULL;
  pd->SPAN_LEN = (int)(sqrt(pd->gd->N));

  pd->roots = (cType *)walloc(pd->total + 2, sizeof(cType));
  pd->allResults = walloc(pd->total+2, sizeof(NodePropArr));

  NodePropArr * allResults = pd->allResults;
  cType len = pd->gd->N + 2;

  calcuTotalCap(pd);

  int LN = pd->LEVEL_NUM+1;
  for (int i = 0; i < pd->total; i++)
  {
    allResults[i].pfa = (cType *)walloc(len, sizeof(cType));
    allResults[i].pdep = (cType *)walloc(len, sizeof(cType));
    allResults[i].pcv = (long *)walloc(len, sizeof(long));
    allResults[i].poh = (long *)walloc(len, sizeof(long));
    allResults[i].pcof = (long *)walloc(len, sizeof(long));
    allResults[i].ps = (short *)walloc(len, sizeof(short));
    allResults[i].pacc_upid = (cType *)walloc(len, sizeof(cType));
    allResults[i].pacc_upmincv = (long *)walloc(len, sizeof(long));
    allResults[i].pacc_pos_upmincv = (cType *)walloc(len, sizeof(cType));

    long* ptr = (long *)walloc(len*LN, sizeof(long));
    memset(ptr, 0, len * LN * sizeof(long));
    allResults[i].pacc_cut_adjust = (long **)walloc(len, sizeof(long*));
    for(int j = 0; j<len; j++){
      allResults[i].pacc_cut_adjust[j] = ptr+j*LN;
    }


    aType* ptr2 = (aType *)walloc(len*LN, sizeof(aType));
    memset(ptr2, 0, len * LN * sizeof(aType));
    allResults[i].pacc_cut_adjust_sign = (aType **)walloc(len, sizeof(aType*));
    for(int j = 0; j<len; j++){
      allResults[i].pacc_cut_adjust_sign[j] = ptr2+j*LN;
    }    

    memset(allResults[i].pfa, 0, len * sizeof(cType));
    memset(allResults[i].pdep, 0, len * sizeof(cType));
    memset(allResults[i].pcv, 0, len * sizeof(long));
    memset(allResults[i].poh, 0, len * sizeof(long));
    memset(allResults[i].pcof, 0, len * sizeof(long));
    memset(allResults[i].ps, 0, len * sizeof(short));
    memset(allResults[i].pacc_upid, 0, len * sizeof(cType));
    memset(allResults[i].pacc_upmincv, 0, len * sizeof(long));
    
    

  }  
}



void preProc(PreprocData *pd){
  double tm;
  double totalProcTime = 0;
  NodePropArr *allResults = pd->allResults;
  
  cType root;

  cType len = pd->gd->N + 2;
  long *apply_adj = (long *)walloc(len, sizeof(long));
  cType *depth_map = (cType *)walloc(len, sizeof(cType));

  for (int ipass = 0; ipass < pd->total; ipass++)
  {
    if(pd->rd != NULL){
      free(pd->rd);
      pd->rd = NULL;
    }
    pd->rd = initrand(pd->gd->M*2);    
    
    pd->gpfa = allResults[ipass].pfa;
    pd->gpdep = allResults[ipass].pdep;
    pd->gpcv = allResults[ipass].pcv;
    pd->gpoh = allResults[ipass].poh;
    pd->gpcof = allResults[ipass].pcof;
    pd->gps = allResults[ipass].ps;
    pd->gpaccup = allResults[ipass].pacc_upid;
    pd->gpaccmcv = allResults[ipass].pacc_upmincv;
    pd->gpaccposmcv = allResults[ipass].pacc_pos_upmincv;
    pd->gcutadj = allResults[ipass].pacc_cut_adjust;
    pd->gcutadjsign = allResults[ipass].pacc_cut_adjust_sign;

    if (pd->P == 300)
    {
      pd->mode = 3;
    }
    else
    {
      pd->mode = ipass < pd->P * pd->total / 100 ? 1 : 2;
    }
    
    
    
    
    
    
    root = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % pd->gd->N);
    
    pd->roots[ipass] = root;
    pd->gpdep[root] = 0;
    printf("pass %d before markCut: root is %ld\n",ipass,root);
    fflush(stdout);
    
    
    tm = timer();
    gRoot = root;
    
    markCut(root,pd);
    pd->gpcv[root] = MAX_LONG;
    pd->gpoh[root] = MAX_LONG;


    
    printf("c before buildAcc\n");
	  fflush(stdout);

    memset(apply_adj, 0, len *  sizeof(long));
    memset(depth_map, 0, len *  sizeof(cType));
    buildAcc(pd, root, root, 0,MAX_LONG, apply_adj,depth_map);
    printf("c after buildAcc\n");
    fflush(stdout);


    totalProcTime += timer() - tm;
    printf("c proctime for onepass: %10.06f\n", timer() - tm);
    if (ipass % 10 == 0)
    {
      printf("c the %d passes\n", ipass);
    }
  }

  free(apply_adj);
  free(depth_map);

  printf("c preprocess times %10.6f\n", totalProcTime);

}



void calcuRandomPairs(int numOfPairs, PreprocData *pd){
  double totalTime = 0;
  long mv = MAX_LONG;

  double curTime = 0;
  cType ns, nt;

  if(pd->rd != NULL){
    free(pd->rd);
    pd->rd = NULL;
  }
  pd->rd = initrand(pd->gd->M*2);  
	struct timespec time_start={0,0},time_end={0,0};
  for (int ipair = 0; ipair < numOfPairs;)
  {

    
    ns = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % (pd->gd->N));
    nt = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % (pd->gd->N));
    if (ns != nt)
    {
      
      mv = min((pd->gd->nodes+ns)->totalCap, (pd->gd->nodes+nt)->totalCap);
      
      clock_gettime(CLOCK_REALTIME,&time_start);
      for (int j=0; j < pd->total; j++)
      {
        cType root = pd->roots[j];
        
        long tmp = solveMaxFlowAccVER4(root,mv, pd->allResults+j, ns, nt,pd->SPAN_LEN,pd->LEVEL_NUM);
        
        if (mv > tmp)
        {
          mv = tmp;
        }
        
      }
	  clock_gettime(CLOCK_REALTIME,&time_end);
	  curTime = 10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec;
	  curTime = curTime/10e9;
	  totalTime += curTime;
	  ipair++;
      printf("c hi_treem_res(n,s,mflow,tm) %lu %lu %12.01f %f \n", ns, nt, 1.0 * mv, curTime);
    }
  }

  printf("c run ok! average time %10.6f\n", totalTime / numOfPairs);


}


