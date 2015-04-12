#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cmath>
using namespace std;

const int maxn = 200100;
const double pi = acos(-1);
const int inf = 1 << 30;

struct point
{
	int x, y;

	bool operator < (const point& rhs) const
	{
		if (x == rhs.x) return y < rhs.y;
		else return x < rhs.x;
	}

	double tan (const point& rhs)
	{
		if (x == rhs.x) return 0;
		else return (double) (rhs.y - y) / (rhs.x - x);
	}
} mts[maxn];


point R[maxn], L[maxn];
int n, q;
double ans[maxn];
bool vis[maxn];
pair<int, int> query[maxn];

void input()
{
	int Max = -1;
	int x, y;
	scanf("%d%d", &n, &q);
	mts[0].x = 0, mts[0].y = -1;
	mts[n + 1].y = -1;
	for (int i = 1; i <= n; i++) {
		scanf("%d%d", &mts[i].x, &mts[i].y);
		Max = max(mts[i].x, Max);
	}
	mts[n + 1].x = Max;
	for (int i = 1; i <= q; i++) {
		scanf("%d%d", &x, &y);
		query[i] = make_pair(x, y);
		R[i] = mts[n + 1], L[i] = mts[0];
		if (!x) vis[y] = true;
	}
	sort(mts + 1, mts + n + 1);
	for (int i = 1; i <= n; i++) {
		if (!vis[i]) {
			R[i] = mts[n + 1], L[i] = mts[0];
			query[++q] = make_pair(0, i);
		}
	}
}

struct Q
{
	int idx, pos;
	bool operator < (const Q& rhs) const {
		return mts[idx] < mts[rhs.idx];
	}
} Qs[maxn];

point Ds[maxn], stack[maxn];
int top;

int cross(point& a, point& b, point& c)
{
	int x1 = a.x - c.x, y1 = a.y - c.y, x2 = b.x - c.x, y2 = b.y - c.y;
	return x1 * y2 - y1 * x2;
}

void pushin(point& p, int d)
{                            
	// d = 1, 从左向右 , 反之 -1
	if (d == 1) {
		while (top > 0 && cross(p, stack[top - 1], stack[top]) >= 0)
			top --;
		stack[++top] = p;
	} else {
		while (top > 0 && cross(stack[top - 1], p, stack[top]) >= 0)
			top --;
		stack[++top] = p;
	}
}

int find(int ind, int d)
{
	//寻找凸壳最近的点 d = 1, 从左到右, 反之 -1
	int l = 0, r = top, mid;
	int ret = top;
	point& p = mts[ind];
	if (d == 1) {

		while (r >= l) {
			mid = (l + r) >> 1;
			if (cross(p, stack[mid - 1], stack[mid]) < 0) {
				l = mid + 1;
				ret = mid;
			} else {
				r = mid - 1;
			}
		}

	} else {

		while (r >= l) {
			mid = (l + r) >> 1;
			if (cross(stack[mid + 1], p, stack[mid]) < 0) {
				r = mid - 1;
				ret = mid;
			} else {
				l = mid + 1;
			}
		}
	}
	return ret;
}

int pos[maxn];

int cmp(const int& a, const int& b)
{
	return mts[a] < mts[b];
}

void solve(int l, int r)
{
	//前半部分的询问，后半部分的加点。
	int mid, tot_Q = 0, tot_D = 0, cnt;
	top = -1;
	mid = (l + r) >> 1;
	for (int i = l; i <= mid; i++) {
		if (query[i].first) {
			cnt = query[i].second;
			Qs[++tot_Q].idx = cnt;
			Qs[tot_Q].pos = i;
		}
	}
	for (int i = mid + 1; i <= r; i++) {
		if (!query[i].first) {
			cnt = query[i].second;
			Ds[++tot_D] = mts[cnt];
		}
	}

	sort(Qs + 1, Qs + tot_Q + 1);
	sort(Ds + 1, Ds + tot_D + 1);

	int k;
	for (int i = 1, k = 1; i <= tot_Q; i++) {
		while (k <= tot_D && Ds[k].x < mts[Qs[i].idx].x) {
			pushin(Ds[k], 1);
			k++;
		}
		if (top >= 0) {
			cnt = find(Qs[i].idx, 1);
			//是否比以前的答案要优
			if (cross(mts[Qs[i].idx], L[Qs[i].pos], stack[cnt]) <= 0) {
				L[Qs[i].pos] = stack[cnt];
			}
		}
	}
	top = -1;
	for (int i = tot_Q, k = tot_D; i >= 1; i--) {
		while (k >= 1 && Ds[k].x > mts[Qs[i].idx].x) {
			pushin(Ds[k], -1);
			k--;
		}
		if (top >= 0) {
			cnt = find(Qs[i].idx, -1);
			if (cross(R[Qs[i].pos], mts[Qs[i].idx], stack[cnt]) <= 0) {
				R[Qs[i].pos] = stack[cnt];
			}
		}
	}
	if (l < mid) solve(l, mid);
	if (r > mid + 1) solve(mid + 1, r);
}

void getanwser()
{
	double t1, t2;
	int ind;
	for (int i = 1; i <= q; i++) {
		if (query[i].first) {
			ind = query[i].second;
			t1 = pi / 2 - atan(mts[ind].tan(R[i]));
			t2 = pi / 2 + atan(mts[ind].tan(L[i]));
			printf("%.6lf\n", t1 + t2);
		}
	}
}

int main()
{
	freopen("lhxsb.in", "r", stdin);
	freopen("lhxsb.out", "w", stdout);

	input();
	solve(1, q);
	getanwser();

	fclose(stdin);
	fclose(stdout);
	return 0;
}
