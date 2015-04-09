#include <cstdio>
#include <cstring>
#include <algorithm>
#define err(...) fprintf(stderr, __VA_ARGS__)
using namespace std;
typedef long long LL;
const int maxn = 400010;
const int maxm = 1200300;
const long long inf = 21390621430000;

int stack[30], top;

void readInt(int& x)
{
	char t = getchar();
	while (t < '0' || t > '9') t = getchar();
	x = 0;
	while (t >= '0' && t <= '9') {
		x = (x << 1) + (x << 3) + t - '0';
		t = getchar();
	}
}


void readLL(LL& x)
{
	char t = getchar();
	while (t < '0' || t > '9') t = getchar();
	x = 0;
	while (t >= '0' && t <= '9') {
		x = (x << 1) + (x << 3) + t - '0';
		t = getchar();
	}
}


void printLL(LL x)
{
	top = -1;
	while (x) {
		stack[++top] = x % 10;
		x /= 10;
	}
	for (; top >= 0; top--) {
		putchar(stack[top] + '0');
	}
	putchar('\n');	
}


struct Graph
{
	LL w[maxn];
	LL f[maxn], d[maxn], Sum[maxn];
	int fst[maxn];
	int u[maxm], v[maxm], nxt[maxm];
	int _index;


	Graph() {
		memset(fst, -1 ,sizeof(fst));
		_index = 0;
	}

	void addedge(int x, int y) {
		_index ++;
		u[_index] = x, v[_index] = y;
		nxt[_index] = fst[x];
		fst[x] = _index;
	}
	
	void InitDfs(int cnt, int pre)
	{
		LL sum = inf;
		f[cnt] = w[cnt];
		Sum[cnt] = inf;
		for (register int i = fst[cnt]; i != -1; i = nxt[i]) {
			if (v[i] != pre) {
				InitDfs(v[i], cnt);
				if (Sum[cnt] == inf) Sum[cnt] = 0;
				Sum[cnt] += f[v[i]];
			}
		}
		if (Sum[cnt] < f[cnt]) {
			f[cnt] = Sum[cnt];
		}
		d[cnt] = w[cnt] - Sum[cnt];
	}
};

struct SegmentTree
{
	// val -> Dx = Wx - Fx
	struct Node{
		int l, r;
		LL add;
		LL minv, maxv;
	};

	Node nodes[maxn << 2];

	void maintain(int cnt)
	{
		int l = cnt << 1, r = cnt << 1 | 1;
		if (l != r) {
			nodes[cnt].maxv = max(nodes[l].maxv, nodes[r].maxv);
			nodes[cnt].minv = min(nodes[l].minv, nodes[r].minv);
		}
	}

	void BuildTree(int cnt, int l, int r, LL val[])
	{
		int mid = (l + r) >> 1;
		nodes[cnt].l = l, nodes[cnt].r = r;
		nodes[cnt].add = 0;
		if (l == r) {
			nodes[cnt].maxv = nodes[cnt].minv = val[l];
			return ;
		}
		BuildTree(cnt << 1, l, mid, val);
		BuildTree(cnt << 1 | 1, mid + 1, r, val);
		maintain(cnt);
	}

	void pushdown(int cnt)
	{
		int l = cnt << 1, r = cnt << 1 | 1;
		if (nodes[cnt].add != 0) {
			nodes[l].add += nodes[cnt].add;
			nodes[r].add += nodes[cnt].add;
			nodes[l].maxv += nodes[cnt].add, nodes[l].minv += nodes[cnt].add;
			nodes[r].maxv += nodes[cnt].add, nodes[r].minv += nodes[cnt].add;
			nodes[cnt].add = 0;
		}
	}

	LL queryMin(int cnt, int L, int R)
	{
		LL ret;
		int l = nodes[cnt].l, r = nodes[cnt].r;
		int mid = (l + r) >> 1;
		if (L == l && R == r) return nodes[cnt].minv;
		pushdown(cnt);
		if (R <= mid) ret = queryMin(cnt << 1, L, R);
		else if (L > mid) ret = queryMin(cnt << 1 | 1, L, R);
		else ret = min(queryMin(cnt << 1, L, mid), queryMin(cnt << 1 | 1, mid + 1, R));
		maintain(cnt);
		return ret;
	}

	LL query(int cnt, int x)
	{
		LL ret;
		int l = nodes[cnt].l, r = nodes[cnt].r;
		if (l == x && r == x) return nodes[cnt].maxv;
		pushdown(cnt);
		int mid = (l + r) >> 1;
		if (x <= mid) ret = query(cnt << 1, x);
		else ret = query(cnt << 1 | 1, x);
		maintain(cnt);
		return ret;
	}

	void update(int cnt, int L, int R, LL v)
	{
		int l = nodes[cnt].l, r = nodes[cnt].r;
		if (L == l && R == r) {
			nodes[cnt].add += v;
			nodes[cnt].maxv += v;
			nodes[cnt].minv += v;
			return ;
		}
		pushdown(cnt);
		int mid = (l + r) >> 1;
		if (R <= mid) update(cnt << 1, L, R, v);
		else if (L > mid) update(cnt << 1 | 1, L, R, v);
		else {
			update(cnt << 1, L, mid, v);
			update(cnt << 1 | 1, mid + 1, R, v);
		}
		maintain(cnt);
	}
	
	LL Div(int cnt, int L, int R, LL v) //find the last element that lower than v 
	{
		LL ret = -1;
		int l = nodes[cnt].l, r = nodes[cnt].r;
		int mid = (l + r) >> 1;
		if (nodes[cnt].minv >= v) return ret;
		if (l == r) {
			if (nodes[cnt].minv < v) return l;
			else return -1;
		}
		pushdown(cnt);
		if (R <= mid) ret = Div(cnt << 1, L, R, v);
		else if (L > mid) ret = Div(cnt << 1 | 1, L, R ,v);
		else {
			ret = Div(cnt << 1 | 1, mid + 1, R, v);
			if (ret == -1) ret = Div(cnt << 1, L, mid, v);
		} 
		maintain(cnt);
		return ret;		
	}

	void init(int _n, LL val[])
	{
		BuildTree(1, 1, _n, val);
	}
};

Graph InitGraph;

struct HLD
{
	int n, _link;
	int father[maxn], top[maxn], size[maxn], dep[maxn], w[maxn], son[maxn], ind[maxn];
	LL W[maxn];
	SegmentTree sg;

	void dfs1(int cnt, int pre, int depth, Graph& res)
	{
		int Maxind = cnt, Max = -1;
		size[cnt] = 1;
		father[cnt] = pre;
		dep[cnt] = depth;

		for (register int i = res.fst[cnt]; i != -1; i = res.nxt[i]) {
			if (res.v[i] != pre) {
				dfs1(res.v[i], cnt, depth + 1, res);
				size[cnt] += size[res.v[i]];
				if (size[res.v[i]] > Max) {
					Max = size[res.v[i]];
					Maxind = res.v[i];
				}
			}
		}	

		son[cnt] = Maxind;
	}

	void dfs2(int cnt, Graph& res)
	{
		if (son[cnt] != cnt) {
			for (register int i = res.fst[cnt]; i != -1; i = res.nxt[i]) {
				if (res.v[i] == son[cnt]) {
					top[res.v[i]] = top[cnt];
					w[res.v[i]] = ++_link;
					ind[_link] = res.v[i];
					W[_link] = res.d[res.v[i]];
					dfs2(res.v[i], res);		
				}
			}
		}

		for (register int i = res.fst[cnt]; i != -1; i = res.nxt[i]) {
			if (res.v[i] != son[cnt] && res.v[i] != father[cnt]) {
				w[res.v[i]] = ++_link;
				ind[_link] = res.v[i];
				W[_link] = res.d[res.v[i]];
				dfs2(res.v[i], res);
			}
		}
	}

	LL query(int x) 
	{
		return sg.query(1, w[x]);
	}

	void init(int _n, Graph& res)
	{
		_link = 0;
		res.addedge(1, 0);
		res.addedge(0, 1);
		dfs1(0, 0, 1, res);

		for (int i = 0; i <= _n; i++) {
			top[i] = i;
		}

		w[1] = ++_link;
		W[_link] = InitGraph.d[1];
		dfs2(1, res);
		sg.init(_n, W);
		ind[1] = 1;
		father[1] = 0;
	}

	void Modify(int x, LL v)
	{
		int cnt, aim, nxt;
		LL delta, tmp;
		bool goon = false;
		LL Dx = query(x);
		
		InitGraph.w[x] += v;
		sg.update(1, w[x], w[x], v);
		if (Dx >= 0) return;
		if (Dx + v >= 0) delta = -Dx;
		else delta = v;

		cnt = father[x];
		if(cnt == 0) return;
		while (1) {
			if (cnt == 0) break;
			goon = true;
			aim = sg.Div(1, w[top[cnt]], w[cnt], delta);
			if (aim == -1) aim = w[top[cnt]];
			
			tmp = sg.query(1, aim);
			sg.update(1, aim, w[cnt], -delta);
			if(tmp < 0) goon = false;
			else if(tmp < delta) delta = tmp;

			if (!goon) break;
			else cnt = father[ind[aim]];
		}
	}
};

int n, q;
HLD Treecut;

void input()
{
	int a, b;
	readInt(n);
	for (register int i = 1; i <= n; i++) {
		readLL(InitGraph.w[i]);
	}
	for (register int i = 1; i <= n - 1; i++) {
		readInt(a);
		readInt(b);
		InitGraph.addedge(a, b);
		InitGraph.addedge(b, a);
	}
	readInt(q);
}

void solve()
{
	char command;
	int x;
	LL y;
	LL tmp, ans;

	InitGraph.InitDfs(1, 1);
	Treecut.init(n, InitGraph);
	while (q--) {
		command = getchar();
		while (command != 'Q' && command != 'C') command = getchar();
		if (command == 'Q') {
			readInt(x);
			tmp = Treecut.query(x);
			ans = min(InitGraph.w[x], InitGraph.w[x] - tmp);
			printLL(ans);
		} else {
			readInt(x);
			readLL(y);
			Treecut.Modify(x, y);
		}
	}
}

int main()
{
	freopen("flood.in", "r", stdin);
	freopen("flood.out", "w", stdout);

	input();
	solve();

	fclose(stdin);
	fclose(stdout);
	return 0;
}
