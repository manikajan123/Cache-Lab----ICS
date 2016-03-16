#include "cachelab.h"
#include <stdio.h>
#include <getopt.h> 
#include <stdlib.h> 
#include <unistd.h>
#include<string.h>

int verbose;
int s,b;
int line;
int hit=0,miss=0,eviction=0;
int counter;
struct cache
{
    long long tag;
    int valid;
    int t;
};

struct cache* ele;

void visit(long long ad)
{
    int f=(1<<b)-1;
    f=(1<<(s+b))-1-f;
    int ss=(ad&f)>>b;
    long long tag=ad>>(s+b);
    int i;
    int minn=100000000,mini=0;
    for(i=ss*line;i<(ss+1)*line;++i)
    {
        if(ele[i].valid==1&&ele[i].tag==tag)
        {
            hit++;
            ele[i].t = counter;
            counter++;
            return;
        }
        if(ele[i].t<minn)
        {
            minn=ele[i].t;
            mini=i;
        }
        if(ele[i].valid==0)
            break;
    }
    miss++;
    //有空
    if(ele[i].valid==0&&i<(ss+1)*line)
    {
        ele[i].tag=tag,ele[i].valid=1,ele[i].t=counter;
    }
    //没空替换
    else
    {
        ele[mini].tag=tag;
        ele[mini].t=counter;
        eviction++;
    }
    counter++;
}


int main(int argc, char* argv[])
{
    int opt;
    while((opt=getopt(argc,argv,"vs:E:b:t:"))!=-1)
    {
    	switch(opt)
    	{
    		case'v':
    			verbose=1;
    			break;
    		case 's':
    			s=atoi(optarg);
    			break;
    		case'E':
    			line=atoi(optarg);
    			break;
    		case 'b':
    			b=atoi(optarg);
    			break;
    		case 't':
    			freopen(optarg, "r", stdin);
    			break;
    	}
    }
    char op;
    long long ad;
    int n;
    ele = (struct cache*)malloc((1<<s)*line*sizeof(struct cache));
    memset(ele, 0, (1<<s)*line*sizeof(struct cache));
    while(scanf("%c",&op)!=EOF)
    {
        if(op!='L'&&op!='S'&&op!='M'&&op!='I')
            continue;
        scanf("%llx,%d",&ad, &n);
        switch (op)
        {
            case'L':
            case'S':
            {
                visit(ad);
                break;
            }
            case'M':
            {
                visit(ad);
                visit(ad);
                break;
            }
            case'I':
            break;
        }
    }
   printSummary(hit, miss, eviction);
    free(ele);
    return 0;
}
