#include <iostream>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <time.h>
#include <windows.h>
#include<fstream>

#define n 4		//一个网络的规模，即为4x4mesh排列
#define N n*n	//网络中的节点总数，即染色体中的基因位数
#define M 1000		//种群中的环组数目，即种群规模
#define R n*(n-1)/2	//每个环组中包含的环的数目
#define adjust_index 0.5	//均值和方差之间权重的调节因子
#define TM N*(N-1)/2
#define S 4*(n-1)

using namespace std;

double random(double start, double end)
{
	return start + (end - start)*rand() / (RAND_MAX + 1.0);
}

struct Ring{
	int a[R][N];	//每行表示一个环
	int length[R];	//length表示环的长度
}ring;

struct Group{
	int a[R][N];
	int length[R];
	double fit;
	double rfit;//相对的fit值，即所占的百分比  
	double cfit;//积累概率
	double occupation_index;
	double mean;
	int mean_rank;
	double var;
	int var_rank;
	double index_perlink;
	int perlink;
}group[M];

double inject[N][N] = { { 0 } };		//二维数组inject中存放了所有节点的注入概率，通过改变每个节点的注入概率，可以表示不同的通信模式
										//inject为全1的情况为均匀流量

double occupation[N][N] = { { 0 } };

void ring_creating();
void group_output();
int judge_rings(Ring ring);
void group_creating();
void occupation_index(Group group[M]);
void initial_inject();
int location(int start, int t, int r);
void initial_probability(double probability[R][N][N]);
void initial_num(int num[R][N][N]);
void initial_data(double data[R][N][N]);
void Sort_mean();
void Sort_var();
void Sort_perlink();

int main(int argc, char *argv[]) {

	int i, j;
	ofstream fout;
	clock_t start, finish;
	double totaltime;
	start = clock();

	group_creating();
	initial_inject();
	occupation_index(group);
	//group_output();

	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n此程序的运行时间为" << totaltime << "秒！" << endl;
	
	fout.open("output.txt");

	for (int t = 0; t < M; t++)
	{
		if (group[t].occupation_index < 7)
		{
			for (i = 0; i<R; i++)
			{
				for (j = 0; j<N; j++)
				{
					fout << group[t].a[i][j];
				}
				fout << "	length:" << group[t].length[i] << endl;
			}
			fout << "fit：" << group[t].fit << endl;
			fout << "mean: " << group[t].mean << "	rank: " << group[t].mean_rank << endl;
			fout << "var: " << group[t].var << "	rank: " << group[t].var_rank << endl;
			fout << "occupation index: " << group[t].occupation_index << endl;
			fout << "index_perlink: " << group[t].index_perlink << "	rank: " << group[t].perlink << endl << endl;
		}
	}

	fout << flush;
	fout.close();

	system("pause");
	return 0;
}

void ring_creating()
{
	int i;
	int j;
	int len = 0;

	double random(double, double);
	srand(unsigned(time(0) * random(0,10000)));

	//产生环
	ring.a[R][M] = { 0 };
	for (i = 0; i < R; i++)
	{
		for (j = 0; j < N; j++)
		{
			ring.a[i][j] = int(random(0, 2));	//1表示环包含节点（i,j）,0表示环不包含该节点
		}
	}

	//计算环长度
	for (i = 0; i<R; i++)
	{
		for (j = 0; j<N; j++)
		{
			if (ring.a[i][j] == 1)
				len += 1;
		}
		ring.length[i] = len;
		len = 0;
	}

	for (i = 0; i < R; i++)
	{
		if (ring.length[i] > S)
			ring_creating();
	}
}

void group_output()
{
	int i;
	int j;
	int r;

	//输出生成环的二进制编码，环的长度
	for (r = 0; r < M; r++)
	{
		for (i = 0; i<R; i++)
		{
			for (j = 0; j<N; j++)
			{
				cout << group[r].a[i][j];
			}
			cout << "	length:" << group[r].length[i] << endl;
		}
		cout << "fit：" << group[r].fit << endl;
		cout << "mean: " << group[r].mean << "	rank: " << group[r].mean_rank << endl;
		cout << "var: " << group[r].var << "	rank: " << group[r].var_rank << endl;
		cout << "occupation index: " << group[r].occupation_index << endl;
		cout << "index_perlink: " << group[r].index_perlink <<"	rank: "<< group[r].perlink << endl << endl;
	}
}

int judge_rings(Ring ring)		//判断环组是否能满足不跨环全连通的约束条件
{
	int i;
	int j;
	int k;
	int r = 0;
	int flag_s[N][N] = { {0} };				//用一个数组的每一位表示一对节点之间的是否可以实现不跨环通信，1表示可以，0表示不行

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			for (k = 0; k < R; k++)
			{
				if ((ring.a[k][i] == 1) && (ring.a[k][j] == 1))
				{
					flag_s[i][j] = 1;
				}
			}
		}
	}

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			if ((flag_s[i][j] == 0) && (i != j))
				r++;
			//cout<<flag_s[i][j];
		}
		//cout << endl;
	}
	if (r == 0)
		return 1;
	else
		return 0;
}

void occupation_index(Group group[M])
{
	int i;
	int j;
	int r;
	double s = 0;
	int t;
	double probability[R][N][N];
	int num[R][N][N];
	int a[R] = { 0 };
	int p;
	int q;
	double data[R][N][N];
	int tn;
	int tm, ts = 0;
	double mean[M] = { 0 };
	double var[M] = { 0 };

	for (t = 0; t < M; t++)
	{
		initial_probability(probability);
		initial_num(num);
		initial_data(data);
		tn = 0;
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				s = 0;
				for (r = 0; r < R; r++)
				{
					if ((group[t].a[r][i] == 1) && (group[t].a[r][j] == 1) && (i != j))
					{
						a[r] = 1;
						s++;
					}
				}

				if (s != 0)
				{
					for (r = 0; r < R; r++)
					{
						if (a[r] == 1)
						{
							for (p = 0; p < N; p++)
							{
								if (group[t].a[r][p] == 1)
								{
									q = location(p, t, r);
									probability[r][p][q] += inject[i][j] / (2 * s);
									num[r][p][q]++;
								}
							}
						}
					}
				}
			}
		}

		for (i = 0; i < R; i++)
		{
			for (j = 0; j < N; j++)
			{
				for (int k = 0; k< N; k++)
				{
					if (num[i][j][k] != 0)
					{
						data[i][j][k] += probability[i][j][k] / num[i][j][k];
					}
				}
			}
		}

		tm = 0;
		for (i = 0; i < R; i++)
		{
			for (j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
						mean[t] += data[i][j][k];
						tm++;
				}
			}
		}
		group[t].mean = mean[t] / tm;

		ts = 0;
		for (i = 0; i < R; i++)
		{
			for (j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					var[t] += (data[i][j][k] - mean[t])*(data[i][j][k] - mean[t]);
					ts++;
				}
			}
		}
		group[t].var = var[t] / ts;

		//group[t].occupation_index = 1000*mean[t] + 1000*var[t];
		//cout << group[t].occupation_index << endl;
	}

	Sort_mean();
	Sort_var();

	for (i = 0; i < M; i++)
	{
		group[i].occupation_index = group[i].mean_rank + group[i].var_rank;
		group[i].index_perlink = group[i].occupation_index / group[i].fit;
	}
	Sort_perlink();
}

void group_creating()		//创建种群
{
	//int s = -1;
	int i;
	int j;
	int r;
	double sum = 0;

	for (r = 0; r < M; r++)
	{
		while (1)
		{
			ring_creating();
			if (judge_rings(ring) == 1)
				break;
		}

		for (i = 0; i < R; i++)
		{
			for (j = 0; j < N; j++)
			{
				group[r].a[i][j] = ring.a[i][j];
			}
			group[r].length[i] = ring.length[i];
			group[r].fit += ring.length[i];
		}
		sum += group[r].fit;
		//Sleep(random(2000, 3000));
	}

	for (r = 0; r < M; r++)
	{
		group[r].rfit = group[r].fit / sum;
		group[r].cfit = 0;//将其初始化为0  
	}
}

void initial_inject()
{
	int i;
	int j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			inject[i][j] = 1.0;
		}
	}
}

int location(int start, int t, int r)
{
	int i;

	for (i = start+1; i < N; i++)
	{
		if (group[t].a[r][i] == 1)
			return i;
	}

	if (i == N)
		location(-1, t, r);
}

void initial_probability(double probability[R][N][N])
{
	int i;
	int j;
	int s;

	for (s = 0; s < R; s++)
	{
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				probability[s][i][j] = 0;
			}
		}
	}
}

void initial_num(int num[R][N][N])
{
	int i;
	int j;
	int s;

	for (s = 0; s < R; s++)
	{
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				num[s][i][j] = 0;
			}
		}
	}
}

void initial_data(double data[R][N][N])
{
	int i;
	int j;
	int s;

	for (j = 0; j < R; j++)
	{
		for (s = 0; s < N; s++)
		{
			for (i = 0; i < N; i++)
			{
				data[j][s][i] = 0;
			}
		}
	}
}

/*均值排序*/
void Sort_mean()
{
	int j;
	int i;
	double min;
	int min_index = 0;
	int a[M];

	for (i = 0; i < M; i++)
	{
		a[i] = 1;
	}

	for (i = 0; i< M; i++)
	{
		min = 1000;
		for (j = 0; j < M; j++)
		{
			if ((group[j].mean < min)&&(a[j] != 0))
			{
				min_index = j;
				min = group[j].mean;
			}
		}
		group[min_index].mean_rank = i + 1;
		a[min_index] = 0;
	}
}

/*方差排序*/
void Sort_var()
{
	int j;
	int i;
	double min;
	int min_index = 0;
	int a[M];

	for (i = 0; i < M; i++)
	{
		a[i] = 1;
	}

	for (i = 0; i< M; i++)
	{
		min = 1000;
		for (j = 0; j < M; j++)
		{
			if ((group[j].var < min) && (a[j] != 0))
			{
				min_index = j;
				min = group[j].var;
			}
		}
		group[min_index].var_rank = i + 1;
		a[min_index] = 0;
	}
}

/*方差排序*/
void Sort_perlink()
{
	int j;
	int i;
	double min;
	int min_index = 0;
	int a[M];

	for (i = 0; i < M; i++)
	{
		a[i] = 1;
	}

	for (i = 0; i< M; i++)
	{
		min = 1000;
		for (j = 0; j < M; j++)
		{
			if ((group[j].index_perlink < min) && (a[j] != 0))
			{
				min_index = j;
				min = group[j].index_perlink;
			}
		}
		group[min_index].perlink = i + 1;
		a[min_index] = 0;
	}
}
