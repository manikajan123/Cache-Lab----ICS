
/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"
#include "contracts.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. The REQUIRES and ENSURES from 15-122 are included
 *     for your convenience. They can be removed if you like.
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, a[8];
    int row,col;
    REQUIRES(M > 0);
    REQUIRES(N > 0);
    //32*32时
    if(M==32&&N==32){
	for(j=0;j<4;++j) {
	    for(i=0;i<4;++i)
	    {
		if(i==j)
		{
		    for(row=i*8;row<8*i+8;++row)
			for(col=j*8;col<8*j+8;++col)
			{
			    if(i!=3)
				B[col][row+8]=A[row][col];
			}
		    if(i!=3)
		    {
			for(row=i*8;row<8*i+8;++row)
			    for(col=j*8;col<8*j+8;++col)		
				B[col][row]=B[col][row+8];			
		    }
		}
		else {
		    for(row=i*8;row<8*i+8;++row)
			for(col=j*8;col<8*j+8;++col)
			{
			    B[col][row]=A[row][col];
			}
		}
	    }
	}
	for (row=3*8; row<4*8; ++row)
	    for (col=3*8; col<4*8; ++col)
		B[col][row]=A[row][col];
    }
    //M=64,N=64
    else if(M==64&&N==64){
	//用[0][0],[63][63]作为跳板
	//把[0][7]和[7][0]先借助[1][1]方块先写
	//A[0][7]
	//这里的两个for循环不可以合并，会产生额外的miss
	for(row=0;row<4;++row)
	    for(col=0;col<4;col++)
		//A[0][7]的左上角放到B[7][0]的左上角
		B[56+row][col]=A[col][56+row];
	for(row=0;row<4;++row)
	    for(col=0;col<4;col++)
		//A[0][7]的右上角放到B[1][1]的左下角
		B[8+4+row][8+col]=A[col][56+4+row];
	for(row=0;row<4;++row)
	    for(col=0;col<4;col++)
		//A[0][7]的左下角放到B[7][0]的右上角
		B[56+row][4+col]=A[4+col][row+56];
	for(row=0;row<4;++row)
	    for(col=0;col<4;col++)
		//A[0][7]的右下角放到B[1][1]的右下角
		B[8+4+row][8+4+col]=A[4+col][56+4+row];
	for(row=4;row<8;++row)
	    for(col=0;col<8;++col)
		B[56+row][col]=B[8+row][8+col];
	//A[7][0] B[0][7]
	for(row=0;row<4;++row)
	    for(col=0;col<4;col++)
		//A[7][0]的左上角放到B[0][7]的左上角
		B[row][56+col]=A[56+col][row];
	for(row=0;row<4;++row)
	    for(col=0;col<4;col++)
		//A[7][0]的右上角放到B[1][1]的左下角
		B[8+4+row][8+col]=A[56+col][4+row];
	for(row=0;row<4;++row)
	    for(col=0;col<4;col++)
		//A[7][0]的左下角放到B[0][7]的右上角
		B[row][56+4+col]=A[56+4+col][row];
	for(row=0;row<4;++row)
	    for(col=0;col<4;col++)
		//A[7][0]的右下角放到B[1][1]的右下角
		B[8+4+row][8+4+col]=A[56+4+col][4+row];
	for(row=4;row<8;++row)
	    for(col=0;col<8;++col)
		B[row][56+col]=B[8+row][8+col];
	//计算[1][1]到[6][6]
	for(i=1;i<=6;++i)
	    for(j=1;j<=6;++j)
	    {
		//左上角A还是变到B[0][0]左上角
		for(row=0;row<4;++row)
		    for(col=0;col<4;++col)
			B[row][col]=A[i*8+col][8*j+row];
		//右上角A变成B[7][7]左下角
		for(row=0;row<4;++row)
		    for(col=0;col<4;++col)
			B[60+row][56+col]=A[i*8+col][j*8+4+row];
		//左下角A变成B[0][0]右上角
		for(row=0;row<4;++row)
		    for(col=0;col<4;++col)
			B[row][4+col]=A[i*8+4+col][j*8+row];
		//右下角A变成B[7][7]右下角
		for(row=0;row<4;++row)
		    for(col=0;col<4;++col)
			B[60+row][60+col]=A[i*8+4+col][j*8+4+row];
		//然后将B[0][0],B[7][7]赋值给相应的B,
		for(row=0;row<4;++row)
		    for(col=0;col<8;++col)
			B[j*8+row][col+i*8]=B[row][col];
		for(row=4;row<8;++row)
		    for(col=0;col<8;++col)
			B[j*8+row][col+i*8]=B[row+56][col+56];
	    }


	//利用[0][0]来算最后一列和最后一行
	for(i=1;i<=6;++i)
	{
	    //先算A的最后一行[7][i]
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //左上角A的变成左上角B[i][7]
		    B[8*i+row][56+col]=A[56+col][i*8+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //右上角A变成左下角的B[0][0]
		    B[4+row][col]=A[56+col][i*8+4+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //左下角A变成右上角的B[i][7]
		    B[i*8+row][60+col]=A[60+col][i*8+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //右下角A变成右下角的B[0][0]
		    B[4+row][4+col]=A[60+col][i*8+4+row];
	    //将B[0][0]的值赋给B[i][7]
	    for(row=4;row<=7;++row)
		for(col=0;col<8;++col)
		    B[i*8+row][56+col]=B[row][col];
	    //先算A的最后一列[i][7]
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //左上角A的变成左上角B[7][i]
		    B[56+row][8*i+col]=A[8*i+col][56+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //右上角A变成左下角的B[0][0]
		    B[4+row][col]=A[8*i+col][56+4+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //左下角A变成右上角的B[7][i]
		    B[56+row][8*i+4+col]=A[8*i+4+col][56+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //右下角A变成右下角的B[0][0]
		    B[4+row][4+col]=A[8*i+4+col][56+4+row];
	    //将B[0][0]的值赋给B[7][i]
	    for(row=4;row<=7;++row)
		for(col=0;col<8;++col)
		    B[56+row][8*i+col]=B[row][col];
	}

	//再利用B[7][7]去计算第0行和第0列
	for(i=1;i<=6;++i)
	{
	    //先算A的第0行[0][i]
	    //B[i][0]
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //左上角A的变成左上角B[i][0]
		    B[8*i+row][col]=A[col][i*8+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //右上角A变成左下角的B[7][7]
		    B[56+4+row][56+col]=A[col][i*8+4+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //左下角A变成右上角的B[i][0]
		    B[i*8+row][4+col]=A[4+col][i*8+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //右下角A变成右下角的B[7][7]
		    B[56+4+row][56+4+col]=A[4+col][i*8+4+row];
	    //将B[7][7]的值赋给B[i][0]
	    for(row=4;row<=7;++row)
		for(col=0;col<8;++col)
		    B[i*8+row][col]=B[56+row][56+col];
	    //先算A的第0列[i][0]
	    //B[0][i]
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //左上角A的变成左上角B[0][i]
		    B[row][8*i+col]=A[8*i+col][row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //右上角A变成左下角的B[7][7]
		    B[56+4+row][56+col]=A[8*i+col][4+row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //左下角A变成右上角的B[0][i]
		    B[row][8*i+4+col]=A[8*i+4+col][row];
	    for(row=0;row<4;++row)
		for(col=0;col<4;++col)
		    //右下角A变成右下角的B[7][7]
		    B[56+4+row][56+4+col]=A[8*i+4+col][4+row];
	    //将B[7][7]的值赋给B[0][i]
	    for(row=4;row<=7;++row)
		for(col=0;col<8;++col)
		    B[row][8*i+col]=B[56+row][56+col];
	}

	//再计算4个角落A[0][0],A[7][7],A[0][7],A[7][0]
	for(row=0;row<8;++row)
	    for(col=0;col<8;++col)
	    {
		B[row][col]=A[col][row];
		B[56+row][56+col]=A[56+col][56+row];
	    }
    }
    else if(M==61&&N==67){
	//切成8*8的格子暴力
	//A[i][j]到B[j][i]
	for(i=0;i<=7;++i)//行数不超过7组，余5
	    for(j=0;j<=6;++j)//列数不超过8组，余3
	    {
		for(col=0;col<8;++col)
		{
		    for(row=0;row<8;++row)
		    {
			a[row]=A[i*8+row][j*8+col];
		    }
		    for(row=0;row<8;++row)
		    {
			B[j*8+col][i*8+row]=a[row];
		    }
		}
	    }
	//处理56 57 58 59 60行
	for(j=0;j<=6;++j)//列数分成这么多块
	{
	    for(row=64;row<67;++row)
	    {
		for(col=0;col<8;++col)
		{
		    a[col]=A[row][j*8+col];
		}
		for(col=0;col<8;++col)
		{
		    B[j*8+col][row]=a[col];
		}
	    }
	}	
	//处理64 65 66列
	for(i=0;i<=7;++i)
	{
	    for(col=56;col<=61;++col)
	    {
		for(row=0;row<8;++row)
		{
		    a[row]=A[row+i*8][col];
		}
		for(row=0;row<8;++row)
		{
		    B[col][row+8*i]=a[row];
		}
	    }
	}
	//还有5*3的格子直接倒？
	for(row=64;row<67;row++)
	    for(col=56;col<61;++col)
	    {
		B[col][row]=A[row][col];
	    }
    }
    ENSURES(is_transpose(M, N, A, B));
}

/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{

    ENSURES(is_transpose(M, N, A, B));
}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
	for (j = 0; j < M; ++j) {
	    if (A[i][j] != B[j][i]) {
		return 0;
	    }
	}
    }
    return 1;
}

