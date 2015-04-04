#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <queue>
#include <vector>
#define X first
#define Y second
#define MAXN 100
#define MAXP 10000
#define INF 1 << 30
using namespace std;

struct Graph
{
	int fst[MAXP];
	int u[MAXP], v[MAXP], cap[MAXP], nxt[MAXP];
	int rev[MAXP];
	int s, t;
	int dis[MAXP];
	int _Index;
    
	int find(int cnt, int Min)
	{
		int del;
		if (cnt == t) return Min;
		for (int i = fst[cnt]; i != -1; i = nxt[i]) {
			if (cap[i] && dis[v[i]] == dis[cnt] + 1) {
				del = find(v[i], min(Min, cap[i]));
				if(del) {
					cap[i] -= del;
					cap[rev[i]] += del;
					return del;
				}
			}
		}
		return 0;
	}
    
	bool BFS()
	{
		int cnt ;
		static queue<int> q;
		memset(dis, -1, sizeof(dis));
		q.push(s);
		dis[s] = 0;
        
		while(!q.empty()) {
			cnt = q.front();
			q.pop();
			for (int i = fst[cnt]; i != -1; i = nxt[i]) {
				if (dis[v[i]] == -1 && cap[i]) {
					dis[v[i]] = dis[cnt] + 1;
					q.push(v[i]);
				}
			}
		}
		
		return dis[t] != -1;
	}
    
	int MaxFlow()
	{
		int ret = 0, del;
		while (BFS()) {
			while (del = find(s, INF)) ret += del;
		}
		return ret;
	}
    
	void addedge(int a, int b, int c) {
		_Index ++;
		u[_Index] = a, v[_Index] = b, cap[_Index] = c;
		nxt[_Index] = fst[a];
		fst[a] = _Index;
		rev[_Index] = _Index + 1;
		
		_Index ++;
		u[_Index] = b, v[_Index] = a, cap[_Index] = 0;
		nxt[_Index] = fst[b];
		fst[b] = _Index;
		rev[_Index] = _Index - 1;
	}
    
	void Init(int _s, int _t)
	{
		_Index = 0;
		s = _s, t = _t;
		memset(fst, -1, sizeof(fst));
	}
};

Graph FlowGraph;
int h, w, c, m, nw, nc, nm;
int id_h[MAXN];
int id_w[MAXN], id_c[MAXN], id_m[MAXN];
int id_w2[MAXN], id_c2[MAXN], id_m2[MAXN];
vector<int> gw[MAXN], gc[MAXN], gm[MAXN];

bool input()
{
	int x, a;
	scanf("%d%d%d%d%d%d%d", &h, &w, &c, &m, &nw, &nc, &nm);
	if (h == -1) return false;
    
	for (int i = 1; i <= w; i++) {
		scanf("%d", &x);
		gw[i].clear();
		for (int j = 1; j <= x; j++) {
			scanf("%d", &a);
			gw[i].push_back(a);
		}
	}
	for (int i = 1; i <= c; i++) {
		scanf("%d", &x);
		gc[i].clear();
		for (int j = 1; j <= x; j++) {
			scanf("%d", &a);
			gc[i].push_back(a);
		}
	}
	for (int i = 1; i <= m; i++) {
		scanf("%d", &x);
		gm[i].clear();
		for (int j = 1; j <= x; j++) {
			scanf("%d", &a);
			gm[i].push_back(a);
		}
	}
	return true;
}

void init()
{
	int _Tot = 2;
	FlowGraph.Init(1, 2);
	for (int i = 1; i <= h; i++) {
		id_h[i] = ++_Tot;
		FlowGraph.addedge(1, id_h[i], 1);
	}
	for (int i = 1; i <= w; i++) {
		id_w[i] = ++_Tot;
		id_w2[i] = ++_Tot;
		FlowGraph.addedge(id_w[i], id_w2[i], 1);
	}
	for (int i = 1; i <= c; i++) {
		id_c[i] = ++_Tot;
		id_c2[i] = ++_Tot;
		FlowGraph.addedge(id_c[i], id_c2[i], 1);
	}
	for (int i = 1; i <= m; i++) {
		id_m[i] = ++_Tot;
		id_m2[i] = ++_Tot;
		FlowGraph.addedge(id_m[i], id_m2[i], 1);
	}
    
	id_w[0] = ++_Tot, id_w2[0] = ++_Tot;
	id_c[0] = ++_Tot, id_c2[0] = ++_Tot;
	id_m[0] = ++_Tot, id_m2[0] = ++_Tot;
	FlowGraph.addedge(id_w[0], id_w2[0], nw);
	FlowGraph.addedge(id_c[0], id_c2[0], nc);
	FlowGraph.addedge(id_m[0], id_m2[0], nm);
	
	//BuildGraph
    
	int cnt;
    
	for (int i = 1; i <= h; i++) FlowGraph.addedge(id_h[i], id_w[0], 1);
	for (int i = 1; i <= w; i++) {
		for (int j = 0; j < gw[i].size(); j++) {
			cnt = gw[i][j];
			FlowGraph.addedge(id_h[cnt], id_w[i], 1);
		}
		FlowGraph.addedge(id_w2[i], id_c[0], 1);
	}
    
	for (int i = 1; i <= c; i++) {
		for (int j = 0; j < gc[i].size(); j++) {
			cnt = gc[i][j];
			FlowGraph.addedge(id_w2[cnt], id_c[i], 1);
		}
		FlowGraph.addedge(id_w2[0], id_c[i], 1);
		FlowGraph.addedge(id_c2[i], id_m[0], 1);
	}
    
	for (int i = 1; i <= m; i++) {
		for (int j = 0; j < gm[i].size(); j++) {
			cnt = gm[i][j];
			FlowGraph.addedge(id_c2[cnt], id_m[i], 1);
		}
		FlowGraph.addedge(id_m2[i], 2, 1);
		FlowGraph.addedge(id_c2[0], id_m[i], 1);
	}
	FlowGraph.addedge(id_m2[0], 2, INF);

	#ifdef DEBUG
	printf("[%d]\n", FlowGraph._Index);
	#endif
}

void solve()
{
	printf("%d\n", FlowGraph.MaxFlow());
}

int main()
{
	freopen("attackonworldtree.in", "r", stdin);
	freopen("attackonworldtree.out", "w", stdout);
    
	while (input()) {
		init();
		solve();
	}
    
    #ifdef DEBUG
	printf("\nRuntime : %d ms\n", clock());
	#endif
    
	fclose(stdin);
	fclose(stdout);
	return 0;
}
