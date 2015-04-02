#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <queue>
//#define DEBUG
#define X first
#define Y second
#define err(...) fprintf(stderr, __VA_ARGS__)
#define INF 2e9
#define MAXN 100010
#define MAXM 400010
using namespace std;

typedef long long LL;

#define Input
int n, m;
int pCity = 0;
bool p[MAXN];

#define Output
int q;
int Hash[MAXN]; //映射
int _Hash = 0;

struct SegmentTree
{
	struct node{
		int l, r;
		LL maxv;
	};

	int N;
	node nodes[MAXN << 2];

	LL query(int cnt, int L, int R)
	{
		int l = nodes[cnt].l, r = nodes[cnt].r;
		if (L == l && R == r) {
			return nodes[cnt].maxv;
		}
		int mid = (l + r) >> 1;
		if (R <= mid) return query(cnt << 1, L, R);
		else if (L > mid) return query(cnt << 1 | 1, L, R);
		else return max(query(cnt << 1, L, mid), query(cnt << 1 | 1, mid + 1, R));
	}

	void maintain(int cnt)
	{
		if(nodes[cnt].l != nodes[cnt].r) {
			nodes[cnt].maxv = max(nodes[cnt << 1].maxv, nodes[cnt << 1 | 1].maxv);
		}
	}

	void BuildTree(int cnt, int l, int r, LL val[])
	{
		nodes[cnt].l = l, nodes[cnt].r = r;
		if(l == r) {
			nodes[cnt].maxv = val[l];
			return ;
		}
		int mid = (l + r) >> 1;
		BuildTree(cnt << 1, l, mid, val);
		BuildTree(cnt << 1 | 1, mid + 1, r, val);
		maintain(cnt);
	}

	void Init(int _n, LL val[])
	{
		N = _n;
		BuildTree(1, 1, N, val);
	}

};

struct UnionSet
{
	int p[MAXN],rank[MAXN];

	void Init(int _n)
	{
		for(int i = 1; i <= _n; i++) {
			p[i] = i;
			rank[i] = 0;
		}
	}

	int find(int u)
	{
		return (p[u] == u) ? (u) : (p[u] = find(p[u]));
	}

	bool Union(int x, int y)
	{
		int t1 = find(x), t2 = find(y);
		if (t1 == t2) return false;
		if (rank[t1] > rank[t2]) p[t2] = t1;
		else if (rank[t1] == rank[t2]) p[t2] = t1, rank[t1]++;
		else p[t1] = t2;
		return true;
	}
};

struct Edge
{
	int u, v;
	LL w;
	bool operator < (const Edge& rhs) const
	{
		return w < rhs.w;
	}
};

struct Graph
{
	int N;
	int fst[MAXN], pre[MAXN];
	int u[MAXM], v[MAXM], nxt[MAXM];
	LL w[MAXM];
	LL dis[MAXN];
	int _;

	struct Node
	{
		LL dist;
		int u;
		bool operator < (const Node& rhs) const
		{
			return dist > rhs.dist;
		}
	};

	void Init(int _n)
	{
		N = _n;
		memset(fst, -1, sizeof(fst));
		_ = 0;
	}

	void Dijkstra()
	{
		Node cnt,tmp;
		static priority_queue< Node > q;
		static bool vis[MAXN];
		memset(vis, 0, sizeof(vis));
		memset(pre, -1, sizeof(pre));
		for(int i = 1; i <= n; i++) {
			if(p[i]) {
				pre[i] = i;
				dis[i] = 0;
				tmp.u = i;
				tmp.dist = 0;
				q.push(tmp);
			}
			else {
				dis[i] = INF;
			}
		}

		while (!q.empty()) {
			cnt = q.top();
			q.pop();
			if (!vis[cnt.u]) {
				vis[cnt.u] = true;
				for (int i = fst[cnt.u]; i != -1; i = nxt[i]) {
					if (dis[v[i]] > dis[cnt.u] + w[i]) {
						pre[v[i]] = pre[cnt.u];
						dis[v[i]] = dis[cnt.u] + w[i];
						tmp.u = v[i];
						tmp.dist = dis[v[i]];
						q.push(tmp); 
					}
				}	
			}
			
		}
	}

	void MinwTree(Graph& to)
	{
		static Edge DisData[MAXM];
		static UnionSet us;
		us.Init(N);
		int x, y;
		for(int i = 1; i <= _; i++) {
			DisData[i].u = u[i];
			DisData[i].v = v[i];
			DisData[i].w = w[i];
		}
		sort(DisData + 1, DisData + _ + 1);
		for(int i = 1; i <= _; i ++) {
			if(us.Union(DisData[i].u, DisData[i].v)) {
				x = Hash[DisData[i].u], y = Hash[DisData[i].v];
				if(x == -1) x = Hash[DisData[i].u] = ++_Hash;
				if(y == -1) y = Hash[DisData[i].v] = ++_Hash;
				to.addedge(x, y, DisData[i].w);
				to.addedge(y, x, DisData[i].w);
			}
		}
	}

	void addedge(int uu, int vv, LL ww) 
	{	
		_ ++;
		u[_] = uu,v[_] = vv,w[_] = ww;
		nxt[_] = fst[uu];
		fst[uu] = _;
	}
	
	void prt()
	{
		for(int i = 1; i <= N; i++) {
			printf("%d", i);
			for(int j = fst[i]; j != -1; j = nxt[j]) {
				printf(" - %d", v[j]);
			}	
			printf("\n");
		}
	}
};

struct Tree_Cut
{
	int N;
	LL _link;
	int father[MAXN], son[MAXN], top[MAXN], w[MAXN], size[MAXN], dep[MAXN];
	LL W[MAXN];
	LL Max[MAXN];
	SegmentTree sg;

	void dfs1(Graph& res, int cnt, int pre, int depth)
	{
		LL Max = -1;
		int MaxInd = cnt;
		size[cnt] = 1;
		father[cnt] = pre;
		dep[cnt] = depth;

		for(int i = res.fst[cnt]; i != -1; i = res.nxt[i]) {
			if(res.v[i] != pre) {
				dfs1(res, res.v[i], cnt, depth + 1);
				size[cnt] += size[res.v[i]];
				if(size[res.v[i]] > Max) {
					Max = size[res.v[i]];
					MaxInd = res.v[i];
				}
			}
		}

		son[cnt] = MaxInd;
	}

	void dfs2(Graph& res, int cnt, int pre)
	{
		if(son[cnt] != cnt) {
			for(int i = res.fst[cnt]; i != -1; i = res.nxt[i]) {
				if(res.v[i] != pre && res.v[i] == son[cnt]) {
					top[res.v[i]] = top[cnt];
					w[res.v[i]] = ++_link;
					W[_link] = res.w[i];
					Max[res.v[i]] = max(Max[cnt], res.w[i]);
					dfs2(res, res.v[i], cnt);
				}
			}
		}

		for(int i = res.fst[cnt]; i != -1; i = res.nxt[i]) {
			if(res.v[i] != pre && res.v[i] != son[cnt]) {
				w[res.v[i]] = ++_link;
				W[_link] = res.w[i];
				Max[res.v[i]] = res.w[i];
				dfs2(res, res.v[i], cnt);
			}
		}
	}

	LL query(int x, int y)
	{
		int t1 = top[x], t2 = top[y];
		LL ret = -INF;
		while (t1 != t2) {
			if(dep[t1] < dep[t2]) {
				swap(x, y);
				swap(t1, t2);
			}
			ret = max(ret, Max[x]);
			x = father[t1];
			t1 = top[x];
		}

		if(dep[x] > dep[y]) swap(x, y);
		if(x != y) {
			ret = max(ret, sg.query(1, w[son[x]], w[y]));
		}
		return ret;
	}

	void Init(Graph& res, int Root)
	{
		N = res.N;
		_link = 0;
		dfs1(res, Root, Root, 1);
		for(int i = 1; i <= N; i++) {
			top[i] = i;
		}
		dfs2(res, Root, Root);
		sg.Init(_link, W);
	}
};

#define Vars
static Graph InitGraph,TmpGraph,MinTree;
static Tree_Cut Tc;

void readInt(int &x)
{
	x = 0;
	char t = getchar();
	while(t > '9' || t < '0') t = getchar();
	while(t <= '9' && t >= '0') {
		x = (x << 1) + (x << 3) + t - '0';
		t = getchar();
	}
}

void readLL(LL &x)
{
	x = 0;
	char t = getchar();
	while(t > '9' || t < '0') t = getchar();
	while(t <= '9' && t >= '0') {
		x = (x << 1) + (x << 3) + t - '0';
		t = getchar();
	}
}

void input()
{
	int a,b;
	LL c;
	char t;
	scanf("%d%d", &n, &m);
	getchar();
	
	InitGraph.Init(n);

	for(int i = 1; i <= n; i++) {
		t = getchar();
		while(t != '0' && t != '1') t = getchar();
		p[i] = (t == '1');
		if (p[i]) {
			pCity ++;
		}
	}

	for(int i = 1; i <= m; i++) {
		readInt(a);
		readInt(b);
		readLL(c);
		InitGraph.addedge(a, b, c);
		InitGraph.addedge(b, a, c);
	}
	
	memset(Hash, -1, sizeof(Hash));
}

void BuildNewGraph()
{
	int x, y;
	for(int i = 1; i <= InitGraph._; i++)
	{
		x = InitGraph.u[i], y = InitGraph.v[i];
		InitGraph.w[i] += InitGraph.dis[x] + InitGraph.dis[y];
		InitGraph.u[i] = InitGraph.pre[x], InitGraph.v[i] = InitGraph.pre[y];
	}
}

void output()
{
	int q;
	int x, y;
	readInt(q);
	for(int i = 1; i <= q; i++)
	{
		readInt(x);
		readInt(y);
		printf("%I64d\n", Tc.query(Hash[x], Hash[y]));
	}
}

int main()
{
	freopen("travel.in","r",stdin);
	freopen("travel.out","w",stdout);

	input();
	InitGraph.Dijkstra();
	BuildNewGraph();
	MinTree.Init(pCity);
	InitGraph.MinwTree(MinTree);
	Tc.Init(MinTree, 1);
	output();

	#ifdef DEBUG
	printf("\nRuntime : %d ms\n", clock());
	#endif

	fclose(stdin);
	fclose(stdout);
	return 0;
}
