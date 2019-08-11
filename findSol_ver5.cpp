// input
// grids: 20 * 10, 0 / 1;
// cube : 4  * 4 , 0 / 1;
// cx   : x coodination of the cube
// return
// opt  : [int, steps of move(negative left or positive right); int, numbers of rotation, clockwise]

#include <bits/stdc++.h>
using namespace std;

#define pii pair<int, int>
#define mp make_pair
#define pb push_back
#define fi first
#define se second

#define a1 -4.500158825082766
#define a2 3.4181268101392694
#define a3 -3.2178882868487753
#define a4 -9.348695305445199
#define a5 -7.899265427351652
#define a6 -3.3855972247263626

void matrix_sort(int matrix_s[][10], int matrix[][10], int r_size, int *numL);
int count_RT(int matrix[][10], int r_size);
int count_CT(int matrix[][10], int r_size);
int count_NH(int matrix[][10], int r_size);
int count_WS(int matrix[][10], int r_size);
void rotate(int cube[][4], int n);

double fit(int grids[][10], int grids_sort[20][10], int cube[][4], int col, double &ret) // col 开始的那一列
{
    // 首先确定是否越界
    vector<pii> v;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (cube[i][j] > 0) {
                v.pb(mp(i, j + col));
                if (j + col < 0 || j + col >= 10) return -1e9;
            }

    if (v.size() != 4) {
        return -1e10;
    }

    // 下落
    for (int row = 0; row < 20; ++row) {
        int flag = 1;
        for (int i = 0; i < 4; ++i) {
            int x = v[i].fi + 1, y = v[i].se;
            if (x >= 20 || grids[x][y]) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            for (int i = 0; i < 4; ++i)
                ++v[i].fi;
        }
        else {
            break;
        }
    }

    int rmax = -1, rmin = 21;
    static int new_grids[20][10];
    static int raw_size_drop, LH, RE, RT, CT, NH, WS;

    // 拷贝，并倒过来
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 10; ++j)
            new_grids[20 - i - 1][j] = grids[i][j];
    }

    for (int i = 0; i < 4; ++i) {
        int x = 20 - v[i].fi - 1, y = v[i].se;
        new_grids[x][y] = 1;
        rmax = max(rmax, x);
        rmin = min(rmin, x);
    }

    LH = (rmax + rmin) / 2;
    matrix_sort(grids_sort, new_grids, rmax + 1, &RE);//用于计算对应落点下能消掉的行数以及消掉后的矩阵
    RT = count_RT(grids_sort, rmax - RE);
    CT = count_CT(grids_sort, rmax - RE);
    NH = count_NH(grids_sort, rmax - RE);
    WS = count_WS(grids_sort, rmax - RE);
    ret = LH * a1 + a2 * RE;
    cout << LH << ' ' << RE << ' ' << RT << ' ' << CT << ' ' << NH << ' ' << WS << endl;
    return a1 * LH + a2 * RE + a3 * RT + a4 * CT + a5 * NH + a6 * WS;
}

double findBestSolution2(int grids[][10], int cube[][4], int operation[2], int cx = 3)
{
    double bestScore = -1000, score;
    static int delta[4];
	static int grids_sort[20][10];
	static int cube_cp[4][4];
	
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
			cube_cp[i][j] = cube[i][j];
	}
	
    delta[0] = delta[1] = delta[2] = delta[3] = 0;

    int row = -1, cnt;
    for (int i = 0; i < 20; ++i) {
        cnt = 0;
        for (int j = 0; j < 10; ++j) {
            cnt += grids[i][j];
        }
        if (cnt == 0) {
            row = i;
        }
    }

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < 10; ++j) {
            grids[i][j] = 0;
        }
    }

    int rcnt = 4; // 最多需要旋转几次
    int nsize = 4; // 旋转块的大小
    int type = 0; // 特殊类型吗
	int midy;
    if (cube[3][0] && cube[3][1] && cube[2][0] && cube[2][1]) {
        // 2 * 2
        rcnt = 1;
    }
    else if (cube[2][1] && !(cube[3][0] == 0 && (cube[3][1] == 0 || cube[2][0] == 0))) {
        // S型或者T型
        nsize = 3;
        
        if (cube[3][0] == 1 && cube[3][1] == 1 && cube[3][2] == 1) {
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < 4; ++j) {
					if (i == 3) cube[i][j] = 0;
					else cube[i][j] = cube[i + 1][j];
				}
			}
		} else if (cube[3][0] == 1 && cube[2][0] == 1 && cube[1][0] == 1) {
			for (int i = 0; i < 4; ++i) {
				for (int j = 3; j >= 0; --j) {
					if (j == 0) cube[i][j] = 0;
					else cube[i][j] = cube[i][j - 1];
				}
			}
			cx -= 1;
		} else if (cube[3][1] == 1 && cube[2][1] == 1 && cube[1][1] == 1) {
		} else if (cube[2][0] == 1 && cube[2][1] == 1 && cube[2][2] == 1) {
		} else {
			// 是S型了！
			type = 2;
			
			int sumy = 0;
			for (int i = 1; i < 4; ++i) {
				for (int j = 0; j < 3; ++j) {
					if (cube[i][j]) {
						sumy += j;
					}
				}
			}
			midy = sumy / 4;
		}
        
    }
    else if (cube[3][3] || cube[0][0]) {
        // 1 * 4
        rcnt = 2;
        if (cube[3][3]) delta[1] = 1;
        else delta[1] = -1;
    }
    else {
        // L型
        type = 3;
        nsize = 3;
    }

    for (int rs = 0; rs < rcnt; ++rs) {

        // L型 或者 S型 特判
        if (type == 3 || type == 2) {
            if ((cube[3][0] | cube[2][0] | cube[1][0]) == 0) {
                delta[rs] = -1;
            }
        }

        for (int dx = -6; dx <= 6; ++dx) {
        	cout << dx << ' ' << rs << endl;
        	static double tmp;
            score = fit(grids, grids_sort, cube, cx + dx + delta[rs], tmp);
            //cout << score << ' ';
            if (score > bestScore) {
	            bestScore = score;
	            operation[0] = dx;
	        	operation[1] = rs;
            }
        }

        // 旋转90°
        rotate(cube, nsize);
    }
	
	
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
			cube[i][j] = cube_cp[i][j];
	}
	
	return bestScore;
}

void findBestSolution(int grids[][10], int cube[][4], int cx, int cube_next[][4], int operation[])
{
    double bestScore = -1000, score;
    static int delta[4];
	static int grids_sort[20][10], new_grids[20][10];
    delta[0] = delta[1] = delta[2] = delta[3] = 0;
    operation[0] = -1; operation[1] = -1;

    int row = -1, cnt;
    for (int i = 0; i < 20; ++i) {
        cnt = 0;
        for (int j = 0; j < 10; ++j) {
            cnt += grids[i][j];
        }
        if (cnt == 0) {
            row = i;
        }
    }

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < 10; ++j) {
            grids[i][j] = 0;
        }
    }

    int rcnt = 4; // 最多需要旋转几次
    int nsize = 4; // 旋转块的大小
    int type = 0; // 特殊类型吗
	int midy;
    if (cube[3][0] && cube[3][1] && cube[2][0] && cube[2][1]) {
        // 2 * 2
        rcnt = 1;
    }
    else if (cube[2][1] && !(cube[3][0] == 0 && (cube[3][1] == 0 || cube[2][0] == 0))) {
        // S型或者T型
        nsize = 3;
        
        if (cube[3][0] == 1 && cube[3][1] == 1 && cube[3][2] == 1) {
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < 4; ++j) {
					if (i == 3) cube[i][j] = 0;
					else cube[i][j] = cube[i + 1][j];
				}
			}
		} else if (cube[3][0] == 1 && cube[2][0] == 1 && cube[1][0] == 1) {
			for (int i = 0; i < 4; ++i) {
				for (int j = 3; j >= 0; --j) {
					if (j == 0) cube[i][j] = 0;
					else cube[i][j] = cube[i][j - 1];
				}
			}
			cx -= 1;
		} else if (cube[3][1] == 1 && cube[2][1] == 1 && cube[1][1] == 1) {
		} else if (cube[2][0] == 1 && cube[2][1] == 1 && cube[2][2] == 1) {
		} else {
			// 是S型了！
			type = 2;
			
			int sumy = 0;
			for (int i = 1; i < 4; ++i) {
				for (int j = 0; j < 3; ++j) {
					if (cube[i][j]) {
						sumy += j;
					}
				}
			}
			midy = sumy / 4;
		}
        
    }
    else if (cube[3][3] || cube[0][0]) {
        // 1 * 4
        rcnt = 2;
        if (cube[3][3]) delta[1] = 1;
        else delta[1] = -1;
    }
    else {
        // L型
        type = 3;
        nsize = 3;
    }

    for (int rs = 0; rs < rcnt; ++rs) {

        // L型 或者 S型 特判
        if (type == 3 || type == 2) {
            if ((cube[3][0] | cube[2][0] | cube[1][0]) == 0) {
                delta[rs] = -1;
            }
        }

        for (int dx = -6; dx <= 6; ++dx) {
        	
        	if (!(dx == -3  && rs == 0 || dx == -2  && rs == 0)) 
        		continue;
        	
        	static double tmp;
            score = fit(grids, grids_sort, cube, cx + dx + delta[rs], tmp);
            if (score <= -900) continue;
            //cout << score << endl;
            
            for (int i = 0; i < 20; ++i) {
		        for (int j = 0; j < 10; ++j)
		            new_grids[20 - i - 1][j] = grids_sort[i][j];
		    }
            
            static int opt[2];
        	opt[0] = opt[1] = -1;
			score = tmp + findBestSolution2(new_grids, cube_next, opt);
            cout << dx << ' ' << rs << ' ' << opt[0] << ' ' << opt[1] << ' ' << score << endl;
            cout << tmp << endl;
            
            if (score > bestScore) {
                bestScore = score;
                operation[0] = dx;
                operation[1] = rs;
            }
        }

        // 旋转90°
        rotate(cube, nsize);
    }
}

void rotate(int cube[][4], int n)
{
    static int tmp[4][4];
    int dx = 4 - n;

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            tmp[dx + j][n - 1 - i] = cube[dx + i][j];

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            cube[dx + i][j] = tmp[dx + i][j];
}


void matrix_sort(int matrix_s[][10],int matrix[][10],int r_size,int *numL)   //1126黄颖浩测试完毕，单独可用
{
    *numL=0;
    int rr;
    for(int ii=0;ii<20;ii++)
    {
        for(int kk=0;kk<10;kk++){matrix_s[ii][kk]=matrix[ii][kk];}
    }
    rr=r_size;
    for(int i=0;i<rr;i++)
    {
        bool flag=1;
        for(int f=0;f<10;f++){if(matrix_s[i][f]==0)flag=0;}
        if(flag)
        {
            *numL=*numL+1;
            for(int j=i;j<rr;j++)
            {
                for(int ii=0;ii<10;ii++){matrix_s[j][ii]=matrix_s[j+1][ii];}
            }
            rr=rr-1;
            i=i-1;
        }
    }
}

int count_RT(int matrix[][10],int r_size)     //1126黄颖浩测试完毕，单独可用
{
    int RT_num=0;
    for(int i=0;i<r_size;i++)
    {
        for(int j=0;j<9;j++)
        {
            if(matrix[i][j]-matrix[i][j+1]!=0){RT_num+=1;}
        }
        if(matrix[i][0]==0){RT_num+=1;}
        if(matrix[i][9]==0){RT_num+=1;}
    }
    return RT_num;
}

int count_CT(int matrix[][10],int r_size)     //1126黄颖浩测试完毕，单独可用
{
    int CT_num=0;
    for(int i=0;i<10;i++)
    {
        for(int j=0;j<r_size;j++)
        {
            if(matrix[j][i]-matrix[j+1][i]!=0){CT_num+=1;}
        }
        if(matrix[0][i]==0){CT_num+=1;}
    }
    return CT_num;
}

int count_NH(int matrix[][10],int r_size)     //1126黄颖浩测试完毕，单独可用
{
    int NH_num=0;
    for(int i=0;i<10;i++)
    {
        int roof=0;
        for(int j=r_size-1;j>-1;j--)
        {
            if(matrix[j][i]==1){roof=j;break;}
        }
        if(roof!=0)
        {
            for(int j=0;j<roof;j++)
                {
                    if(matrix[j][i]==0)NH_num+=1;
                }
        }
    }
    return NH_num;
}

int count_WS(int matrix[][10],int r_size)     //1126黄颖浩测试完毕，单独可用
{
	int cnt = 0;
	for (int i = 0; i < r_size; ++i) {
		for (int j = 0; j < 10; ++j) {
			if (!matrix[i][j] && (j == 0 || matrix[i][j - 1]) && (j == 9 || matrix[i][j + 1])) {
				cnt += 1;
			}
		}
	}
	return cnt; 
	/*
    int Well_sum=0,Well_depth=0;
    for(int i=1;i<9;i++)
    {
        for(int j=r_size-1;j>-1;j--)
        {
            if(matrix[j][i]==0 && matrix[j][i+1]==1 && matrix[j][i-1]==1){Well_depth+=1;}
            else if(matrix[j][i]==1){Well_sum=Well_sum + 0.5*Well_depth*(Well_depth+1);  Well_depth=0;}
            if(j==0){Well_sum=Well_sum + 0.5*Well_depth*(Well_depth+1);  Well_depth=0;}
        }
        //cout<<i<<' '<<Well_sum<<endl;
    }
    for(int j=r_size-1;j>-1;j--)
    {
        if(matrix[j][0]==0 && matrix[j][1]==1){Well_depth+=1;}
        if(matrix[j][0]==1 || j==0){Well_sum=Well_sum + 0.5*Well_depth*(Well_depth+1);  Well_depth=0;}
    }
    //cout<<"0 "<<Well_sum<<endl;
    for(int j=r_size-1;j>-1;j--)
    {
        if(matrix[j][9]==0 && matrix[j][8]==1){Well_depth+=1;}
        if(matrix[j][9]==1 || j==0){Well_sum=Well_sum + 0.5*Well_depth*(Well_depth+1);  Well_depth=0;}
    }
    //cout<<"9 "<<Well_sum<<endl;
    return Well_sum;
    */
}

// test locally

int main() {
    int  operation[2]={1,2};
    
    int grids_origin[20][10]={
        {0,0,0,0,0,0,0,0,0,0},//0
        {0,0,0,0,0,0,0,0,0,0},//1
        {0,0,0,0,0,0,0,0,0,0},//2
        {0,0,0,0,0,0,0,0,0,0},//3
        {0,0,0,0,0,0,0,0,0,0},//4
        {0,0,0,0,0,0,0,0,0,0},//5
        {0,0,0,0,0,0,0,0,0,0},//6
        {0,0,0,0,0,0,0,0,0,0},//7
        {0,0,0,0,0,0,0,0,0,0},//8
        {0,0,0,0,0,0,0,0,0,0},//9
        {0,0,0,0,0,0,0,0,0,0},//10
        {0,0,0,0,0,0,0,0,0,0},//11
        {0,0,0,0,0,0,0,0,0,0},//12
        {0,0,0,0,0,0,0,0,0,0},//13
        {0,0,0,0,0,0,0,0,0,0},//14
        {0,0,0,0,0,0,0,0,0,0},//15
        {0,0,0,0,0,0,0,0,0,0},//13
        {0,0,0,0,0,0,0,0,0,0},//14
        {0,0,0,1,1,0,0,0,0,0},//15
        {1,1,1,1,1,1,0,0,0,0},//19
        };
        
    int cube[4][4] = {
    	{0,0,0,0},
    	{0,0,0,0},
    	{0,1,1,0},
    	{1,1,0,0}
	};
	
	int cube_next[4][4] = {
    	{0,0,0,0},
    	{0,0,0,0},
    	{1,1,0,0},
    	{1,1,0,0}
	};
    
    findBestSolution(grids_origin, cube, 3, cube_next, operation);
    printf("%d %d\n", operation[0], operation[1]);
    return 0;
}
