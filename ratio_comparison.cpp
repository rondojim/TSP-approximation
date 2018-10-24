#include <stdio.h>
#include <algorithm>
#include <vector>
#include <string.h>
#include <stack>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <time.h>

using namespace std;

ifstream in("input.txt");
ofstream out("results(k=10)_new.txt");
	
typedef pair<int, int> pii;
#define x first
#define y second

struct edge{
	int u, v;
	double w;
	bool operator<(const edge &rhs) const{
		return w < rhs.w;
	}
};

const int MAXN = 22;
pii A[MAXN];
int T, N, S, p[MAXN], vis[MAXN], WW, best[MAXN][(1<<MAXN)], O;
pair<int, int> bM[(1<<MAXN)];
double dp[MAXN][(1<<MAXN)], w[MAXN][MAXN], M[(1<<MAXN)], F;
vector<int> adj[MAXN], adj2[MAXN], odds, path, path2;
edge E[MAXN * MAXN];

int find(int u){
	if(p[u] == u) return u;
	return p[u] = find(p[u]);
}

int sq(int a){ return a * a; }
double dist(pii a, pii b){ return sqrt(1.0 * sq(a.x - b.x) + 1.0 * sq(a.y - b.y)); }

double solve_DP(int i, int mask){
	if(__builtin_popcount(mask) == 1) return dist(A[1], A[i]);
	double &res = dp[i][mask];
	if(res > 0) return res;
	res = 2000000000;
	for(int j=0; j<N; ++j)
		if((mask & (1<<j)) && (j + 1) != i){
			double ww = dist(A[i], A[j + 1]) + solve_DP(j + 1, mask ^ (1<<(i - 1)));
			if(res > ww){
				res = ww;
				best[i][mask] = j;
			}
		}
	return res;
}

void print(int i, int mask){
	if(__builtin_popcount(mask) == 1){
		out<<A[1].x<<" "<<A[1].y<<" "<<A[i].x<<" "<<A[i].y<<endl;
		return;
	}
	print(best[i][mask] + 1, mask ^ (1<<((i - 1))));
	out<<A[i].x<<" "<<A[i].y<<endl;
}

void init(){
	for(int i=1; i<=N; ++i){
		p[i] = i;
		for(int j=1; j<=N; ++j){
			E[++S] = (edge){i, j, dist(A[i], A[j])};
			w[i][j] = dist(A[i], A[j]);
		}
	}
}

void mst(){
	double res = 0;
	for(int i=1; i<=S; ++i){
		int u = find(E[i].u), v = find(E[i].v);
		if(u != v){
			p[u] = v;
			res += E[i].w;
			adj[E[i].u].push_back(E[i].v);
			adj[E[i].v].push_back(E[i].u);
		}
	}
}

void find_odds(){
	for(int i=1; i<=N; ++i)
		if(adj[i].size() & 1) odds.push_back(i);
	O = odds.size();
}

double solve_matching(int mask){
	if(mask == 0) return 0;
	double &res = M[mask];
	if(res > 0) return res;
	res = 2000000000;
	for(int i=0; i<O; ++i)
		for(int j=i+1; j<O; ++j)
			if((mask & (1<<i)) && (mask & (1<<j))){
				double ww = dist(A[odds[i]], A[odds[j]]) + solve_matching((mask ^ (1<<i)) ^ (1<<j));
				if(res > ww){
					res = ww;
					bM[mask] = {i, j};
				}
			}
	return res;
}

void find_edge(int mask){
	if(mask == 0) return;
	adj[odds[bM[mask].x]].push_back(odds[bM[mask].y]);
	adj[odds[bM[mask].y]].push_back(odds[bM[mask].x]);
	find_edge((mask ^ (1<<bM[mask].x)) ^ (1<<bM[mask].y));
}

void best_matching(){
	solve_matching((1<<O) - 1);
	find_edge((1<<O) - 1);
}

void euler_toor(){
	for(int i=1; i<=N; ++i) adj2[i] = adj[i];
	stack<int> S;
	int pos = 1;
	while(!S.empty() || adj2[pos].size() > 0){
		if(adj2[pos].empty()){
			path.push_back(pos);
			pos = S.top();
			S.pop();
		}
		else{
			S.push(pos);
			int v = adj2[pos].back();
			adj2[pos].pop_back();
			for(int i=0; i<adj2[v].size(); ++i){
				if(adj2[v][i] == pos){
					adj2[v].erase(adj2[v].begin() + i);
					break;
				}
			}
			pos = v;
		}
	}
	path.push_back(pos);
}

double hamiltonian_path(){
	double res = 0;
	int l = 0, r = 1;
	vis[path[0]] = 1;
	path2.push_back(path[0]);
	while(r < path.size()){
		if(!vis[path[r]]){
			path2.push_back(path[r]);
			res += w[path[l]][path[r]];
			l = r;
			vis[path[l]] = 1;
			r = l + 1;
		}
		else r++;
	}
	res += w[path[l]][1];
	path2.push_back(1);
	return res;
}

void clear(){
	memset(dp, -1, sizeof(dp));
	memset(M, -1, sizeof(M));
	memset(vis, 0, sizeof(vis));
	for(int i=1; i<MAXN; ++i){ adj[i].clear(); adj2[i].clear(); }
	odds.clear();
	path.clear();
	path2.clear();
}

double solve_Christofides(){
	init();
	sort(E + 1, E + S + 1);
	double mn = 2000000000;
	// change 
	for(int k=1; k<=10; ++k){
		clear();
		for(int o=1; o<=N; ++o) p[o] = o;
		mst();
		find_odds();
		best_matching();
		euler_toor();
		double x = hamiltonian_path();
		if(k == 1) F = x;
		mn = min(mn, x);
		int i = 1;
		while(i <= S){
			int j = i + 1;
			while(E[i].w == E[j].w && j <= S) j++;
			j--;
			random_shuffle(E + i, E + j);
			i = j + 1;
		}
	}
	return mn;
}

int main(){
	clock_t beg = clock(), ed;
	in>>T;
	int X = T, same = 0;
	double mx = -1, maxratio = -1;
	vector<pii> vecmx;
	int cnt = 0;
	double sum = 0, sum2 = 0;
	while(T--){
		clear();
		in>>N;
		S = 0;
		for(int i=1; i<=N; ++i) in>>A[i].x>>A[i].y;
		init();
		double result_Ch = solve_Christofides();
		clear();
		double result_DP = solve_DP(1, (1<<N) - 1);
		double ratio = result_Ch / result_DP, ratioF = F / result_DP;
		if(fabs(ratio - 1) < 0.00001) same++;
		//printf("%lf\n", ratio);
		if(maxratio < ratio) maxratio = ratio;
		if(mx < (ratioF - ratio)){
			mx = ratioF - ratio;
			vecmx.clear();
			for(int i=1; i<=N; ++i) vecmx.push_back(A[i]);
		}
		sum += ratio;
		sum2 += ratioF;
		out.precision(6);
		out<<fixed<<ratioF<<" "<<fixed<<ratio<<endl;
		printf("%d --> {%lf} %lf %lf\n", ++cnt, mx, ratioF, ratio);
		double seconds = double(ed - beg) / CLOCKS_PER_SEC;
	}
	clear();
	out<<"Best Average ratio: "<<sum/X<<", Other Average: "<<sum2/X<<endl;
	out<<"Maximum difference: "<<mx<<endl;
	out<<"Maximum ratio: "<<maxratio<<endl;
	out<<vecmx.size()<<endl;
	int sz = vecmx.size();
	N = sz;
	for(int i=0; i<N; ++i) out<<vecmx[i].x<<" "<<vecmx[i].y<<endl;
	ed = clock();
	double seconds = double(ed - beg) / CLOCKS_PER_SEC;
	out<<"Seconds: "<<seconds<<endl;
	out<<"Number of experiments with ratio 1: "<<same<<endl;
	return 0;
}
