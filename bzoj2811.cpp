#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace std;

const int maxn = 100010;
const int inf = 0x3f3f3f3f;

int n, k, m;

void read(int &x)
{
	int f = 1;
	char t = getchar();
	while (t < '0' || t > '9') 
	{
		if (t == '-') f = -1;
		t = getchar();
	}
	x = 0;
	while (t >= '0' && t <= '9') {
		x = (x << 3) + (x << 1) + t - '0';
		t = getchar();
	}
}

struct SegmentTree
{
	struct Node
	{
		int l, r;
		int set;
	};

	Node nodes[maxn << 2];		

	void BuildTree(int cnt, int l, int r, int v)
	{
		int mid = (l + r) >> 1;
		nodes[cnt].l = l, nodes[cnt].r = r;
		if (cnt == 1) {
			nodes[cnt].set = v;
		} else {
			nodes[cnt].set = -1;
		}
		if (l == r) return;
		if (l <= mid) BuildTree(cnt << 1, l, mid, v);
		if (mid + 1 <= r) BuildTree(cnt << 1 | 1, mid + 1, r, v);
	}

	void pushdown(int cnt)
	{
		if (nodes[cnt].set != -1)
		{
			nodes[cnt << 1].set = nodes[cnt << 1 | 1].set = nodes[cnt].set;
			nodes[cnt].set = -1;
		}
	}

	void Modify(int cnt, int L, int R, int v)
	{
		int l = nodes[cnt].l, r = nodes[cnt].r;
		int mid = (l + r) >> 1;
		if (l == L && r == R) {
			nodes[cnt].set = v;
			return;
		}
		pushdown(cnt);
		if (R <= mid) Modify(cnt << 1, L, R, v);
		else if (L > mid) Modify(cnt << 1 | 1, L, R, v);
		else {
			Modify(cnt << 1, L, mid, v);
			Modify(cnt << 1 | 1, mid + 1, R, v);
		}
	}

	void Query(int cnt, bool Can[], int &tot)
	{
		int l = nodes[cnt].l, r = nodes[cnt].r;
		if (l == r) {
			Can[l] = nodes[cnt].set;
			if (Can[l]) {
				tot ++;
			}
			return;
		}
		pushdown(cnt);
		Query(cnt << 1, Can, tot);
		Query(cnt << 1 | 1, Can, tot);
	}
};

SegmentTree status;

int f[maxn], sum[maxn], h[maxn], pos[maxn];
bool Can[maxn];
int tot_Can = 0;

vector<int> ans;
vector< pair<int, int> > srcquerys, tmpquerys, querys;

void input()
{
	int x, y, v;
	read(n);
	read(k);
	read(m);
	status.BuildTree(1, 1, n, 1);
	for (register int i = 1; i <= m; i++) {
		read(x);
		read(y);
		read(v);
		if (!v) {
			status.Modify(1, x, y, v);
		} else {
			srcquerys.push_back(make_pair(x, y));
		}
	}
	status.Query(1, Can, tot_Can);
}

void solve()
{
	if (tot_Can == k) {
		for (register int i = 1; i <= n; i++) {
			if (Can[i]) {
				printf("%d\n", i);
			}
		}
	} else {
		bool flag;
		int l, r, cnt, x, y, _index = 0;
		sum[0] = h[0] = 0;
		for (register int i = 1; i <= n; i++) {
			sum[i] = sum[i - 1];
			h[i] = h[i - 1];
			if (Can[i]) {
				sum[i] ++;
				h[i] = ++_index;
				pos[_index] = i;
			}
		}
		//去除包含区间
		sort(srcquerys.begin(), srcquerys.end());
		for (register int i = srcquerys.size() - 1, r = n + 1; i >= 0; i--) {
			if (srcquerys[i].second < r) {
				r = srcquerys[i].second;
				tmpquerys.push_back(srcquerys[i]);
			}
		}
		//如果一个区间只有一个忍者，那必须是答案
		for (register int i = 0; i < tmpquerys.size(); i++) {
			if (sum[tmpquerys[i].second] - sum[tmpquerys[i].first - 1] == 1) {
				cnt = pos[h[tmpquerys[i].second]];
				if (Can[cnt]) {
					ans.push_back(cnt);
					Can[cnt] = false;
				}
			}
		}
		sum[0] = pos[0] = 0;
		r = 0;
		for (register int i = 1; i <= n; i++) {
			sum[i] = sum[i - 1];
			pos[i] = pos[i - 1];
			if (Can[i]) {
				sum[i] ++;
				pos[i] = i;
			}
		}
		//包含这些答案点的区间不用考虑
		sort(ans.begin(), ans.end());
		for (register int i = 0; i < tmpquerys.size(); i++) {
			cnt = lower_bound(ans.begin(), ans.end(), tmpquerys[i].first) - ans.begin();
			if (cnt >= ans.size() || ans[cnt] > tmpquerys[i].second) {
				if (sum[tmpquerys[i].second] - sum[tmpquerys[i].first - 1]) {
					querys.push_back(tmpquerys[i]);
				}
			}
		}
		sort(querys.begin(), querys.end());
		//f[i], 代表从i到最后一个区间，至少要取多少个点才能满足每个区间都至少有一个点
		//因为区间已经排序，左端点、右端点已经有序，因此对于一个i区间的最右点，只会被一段连续的区间包含
		//那么假设i区间的最右点被i到j - 1区间包含，f[i] = f[j] + 1
		//二分查找j，即第一个left[j] > right[i]的区间
		for (register int i = querys.size() - 1; i >= 0; i--) {
			x = pos[querys[i].second];
			cnt = lower_bound(querys.begin(), querys.end(), make_pair(x + 1, 0)) - querys.begin();
			f[i] = f[cnt] + 1;
		}
		k -= ans.size();
		cnt = 0, r = -1;
		for (register int i = 0; i < querys.size(); i++) {
			if (querys[i].first > r) {
				r = pos[querys[i].second];
				x = pos[r - 1];
				y = lower_bound(querys.begin(), querys.end(), make_pair(x + 1, 0)) - querys.begin();
				if (cnt + 1 + f[y] > k) ans.push_back(r);
				cnt ++;
			}
		}
		sort(ans.begin(), ans.end());
		if (ans.size() == 0) {
			printf("-1\n");
		} else {
			for (register int i = 0; i < ans.size(); i++) {
				printf("%d\n", ans[i]);
			}
		}
	}
}

int main()
{
	#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	#endif

	input();	
	solve();

	#ifndef ONLINE_JUDGE
	fclose(stdin);
	fclose(stdout);
	#endif
	return 0;
}
