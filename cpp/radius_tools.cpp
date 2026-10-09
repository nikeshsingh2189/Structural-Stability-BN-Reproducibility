// radius_tools.cpp -- exhaustive verification, chain scaling and clique stress tests.
//   g++ -O3 -std=c++17 -o radius_tools radius_tools.cpp
//   ./radius_tools verify            exhaustive check for all DAGs with n<=6 (cardinality + weighted costs 0..15)
//   ./radius_tools chain             brute-force enumeration vs shortest path on X->C1->...->Cs->Y
//   ./radius_tools clique            sparse auxiliary graph vs explicit moral clique on X->Ci->Y, i=1..q
#include <bits/stdc++.h>
#include <chrono>
using namespace std;
using u32 = uint32_t; using u64 = uint64_t;
static const int MAXN = 32;
struct Dag { int n; u32 par[MAXN]; };

static u32 ancestors(const Dag& g, u32 S) {
  u32 A = S, fr = S;
  while (fr) { u32 nf = 0; for (int v = 0; v < g.n; v++) if (fr >> v & 1) nf |= g.par[v]; nf &= ~A; A |= nf; fr = nf; }
  return A;
}
static void moral(const Dag& g, u32 A, u32* adj) {
  for (int v = 0; v < g.n; v++) adj[v] = 0;
  for (int v = 0; v < g.n; v++) if (A >> v & 1) {
    u32 ps = g.par[v] & A;
    for (int p = 0; p < g.n; p++) if (ps >> p & 1) { adj[v] |= 1u << p; adj[p] |= 1u << v; adj[p] |= ps & ~(1u << p); }
  }
}
static bool dsep(const Dag& g, int X, int Y, u32 Z) {
  u32 adj[MAXN]; moral(g, ancestors(g, (1u << X) | (1u << Y) | Z), adj);
  u32 seen = 1u << X, fr = 1u << X;
  while (fr) { u32 nf = 0; for (int v = 0; v < g.n; v++) if (fr >> v & 1) nf |= adj[v]; nf &= ~seen & ~Z; if (nf >> Y & 1) return false; seen |= nf; fr = nf; }
  return true;
}
// min conditioned-vertex cost of an X-Y path in H_C (vertex costs; cost[v] for v in C, else 0); -1 = no path
static long long formula(const Dag& g, int X, int Y, u32 C, const int* cost) {
  u32 adj[MAXN]; moral(g, ancestors(g, (1u << X) | (1u << Y) | C), adj);
  const long long INF = LLONG_MAX / 4; long long d[MAXN]; for (int i = 0; i < g.n; i++) d[i] = INF; d[X] = 0;
  bool ch = true;
  while (ch) { ch = false; for (int v = 0; v < g.n; v++) if (d[v] < INF) for (int w = 0; w < g.n; w++) if (adj[v] >> w & 1) {
      long long nd = d[v] + ((C >> w & 1) ? cost[w] : 0); if (nd < d[w]) { d[w] = nd; ch = true; } } }
  return d[Y] >= INF ? -1 : d[Y];
}
static u64 splitmix(u64 x) { x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL; x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31); }

static void verify() {
  long long tDag = 0, tQ = 0, tMis = 0, tMisW = 0, tFin = 0, tInf = 0, tNM = 0, tR2 = 0;
  printf("%2s %8s %10s %9s %9s %9s %9s %9s %8s\n", "n", "DAGs", "valid", "mismatch", "mismatchW", "finite", "inf", "nonmono", "rho>=2");
  for (int n = 2; n <= 6; n++) {
    vector<pair<int,int>> pr; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pr.push_back({i, j});
    int np = pr.size(); long long dag = 0, q = 0, mis = 0, misW = 0, fin = 0, inf = 0, nm = 0, r2 = 0; long long rc[8] = {0};
    for (long long m = 0; m < (1LL << np); m++) {
      Dag g; g.n = n; for (int v = 0; v < n; v++) g.par[v] = 0;
      for (int e = 0; e < np; e++) if (m >> e & 1) g.par[pr[e].second] |= 1u << pr[e].first;
      dag++;
      for (int X = 0; X < n; X++) for (int Y = X + 1; Y < n; Y++) {
        u32 others = 0; for (int v = 0; v < n; v++) if (v != X && v != Y) others |= 1u << v;
        for (u32 C = others;; C = (C - 1) & others) {
          if (dsep(g, X, Y, C)) {
            q++;
            int cost[MAXN]; for (int v = 0; v < n; v++) cost[v] = (int)(splitmix(((u64)m << 24) ^ (X << 16) ^ (Y << 12) ^ (u64)C * 131 ^ v) % 16);
            int unit[MAXN]; for (int v = 0; v < n; v++) unit[v] = 1;
            vector<u32> subs; for (u32 R = C;; R = (R - 1) & C) { subs.push_back(R); if (!R) break; }
            int best = INT_MAX; long long bestW = LLONG_MAX; unordered_map<u32, char> conn;
            for (u32 R : subs) {
              bool c = !dsep(g, X, Y, C & ~R); conn[R] = c;
              if (c) { best = min(best, __builtin_popcount(R)); long long w = 0; for (int v = 0; v < n; v++) if (R >> v & 1) w += cost[v]; bestW = min(bestW, w); }
            }
            bool nmono = false;
            for (u32 R1 : subs) if (conn[R1]) { for (u32 R2 : subs) if ((R2 & R1) == R1 && R2 != R1 && !conn[R2]) { nmono = true; break; } if (nmono) break; }
            if (nmono) nm++;
            long long f = formula(g, X, Y, C, unit), fW = formula(g, X, Y, C, cost);
            long long b = best == INT_MAX ? -1 : best, bW = bestW == LLONG_MAX ? -1 : bestW;
            if (f != b) mis++;
            if (fW != bW) misW++;
            if (b < 0) inf++; else { fin++; if (b >= 2) r2++; if (b < 8) rc[b]++; }
          }
          if (!C) break;
        }
      }
    }
    printf("%2d %8lld %10lld %9lld %9lld %9lld %9lld %9lld %8lld", n, dag, q, mis, misW, fin, inf, nm, r2);
    if (n == 6) printf("   rho=1..4: %lld %lld %lld %lld", rc[1], rc[2], rc[3], rc[4]);
    printf("\n");
    tDag += dag; tQ += q; tMis += mis; tMisW += misW; tFin += fin; tInf += inf; tNM += nm; tR2 += r2;
  }
  printf("Total: DAGs=%lld valid=%lld mismatch=%lld weighted_mismatch=%lld finite=%lld inf=%lld nonmono=%lld rho>=2=%lld\n", tDag, tQ, tMis, tMisW, tFin, tInf, tNM, tR2);
}

static void chain() {
  printf("%4s %10s %14s %16s %4s\n", "s", "subsets", "brute_ms", "shortest_path_us", "rho");
  for (int s : {4, 6, 8, 10, 12, 14, 16, 18, 20}) {
    Dag g; g.n = s + 2; for (int v = 0; v < g.n; v++) g.par[v] = 0;   // X=0, C_i=i, Y=s+1
    for (int v = 1; v < g.n; v++) g.par[v] = 1u << (v - 1);
    int X = 0, Y = s + 1; u32 C = 0; for (int i = 1; i <= s; i++) C |= 1u << i;
    auto t0 = chrono::steady_clock::now(); int best = INT_MAX;
    for (u32 R = C;; R = (R - 1) & C) { if (!dsep(g, X, Y, C & ~R)) best = min(best, __builtin_popcount(R)); if (!R) break; }
    double brute = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
    int unit[MAXN]; for (int v = 0; v < g.n; v++) unit[v] = 1;
    const int reps = 10000; long long acc = 0; auto t1 = chrono::steady_clock::now();
    for (int r = 0; r < reps; r++) acc += formula(g, X, Y, C, unit);
    double sp = chrono::duration<double, micro>(chrono::steady_clock::now() - t1).count() / reps;
    printf("%4d %10llu %14.3f %16.3f %4d %s\n", s, 1ULL << s, brute, sp, best, (acc / reps == best) ? "" : "MISMATCH");
  }
}

static void clique() {
  printf("%9s %22s %10s %10s %10s\n", "q", "explicit_moral_edges", "aux_verts", "aux_edges", "BFS_ms");
  for (long long q : {1000LL, 10000LL, 100000LL, 1000000LL}) {
    // vertices: X=0, C_i=1..q, Y=q+1, aux=q+2 ; edges X-Ci, Ci-Y (2q skeleton) and aux-Ci (q)
    long long V = q + 3; vector<int> deg(V, 0);
    auto add = [&](int a, int b, vector<pair<int,int>>& E) { E.push_back({a, b}); };
    vector<pair<int,int>> E; E.reserve(3 * q);
    for (int i = 1; i <= q; i++) { add(0, i, E); add(i, q + 1, E); add(q + 2, i, E); }
    vector<int> off(V + 1, 0); for (auto& e : E) { off[e.first + 1]++; off[e.second + 1]++; }
    for (long long i = 0; i < V; i++) off[i + 1] += off[i];
    vector<int> adj(2 * E.size()), pos(off.begin(), off.end() - 1);
    for (auto& e : E) { adj[pos[e.first]++] = e.second; adj[pos[e.second]++] = e.first; }
    vector<char> isC(V, 0); for (int i = 1; i <= q; i++) isC[i] = 1;
    auto t0 = chrono::steady_clock::now();
    vector<int> d(V, INT_MAX); deque<int> dq; d[0] = 0; dq.push_back(0);      // 0-1 BFS, entering a conditioned vertex costs 1
    while (!dq.empty()) { int v = dq.front(); dq.pop_front();
      for (int k = off[v]; k < off[v + 1]; k++) { int w = adj[k]; int c = isC[w] ? 1 : 0; if (d[v] + c < d[w]) { d[w] = d[v] + c; if (c) dq.push_back(w); else dq.push_front(w); } } }
    double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
    if (d[q + 1] != 1) { printf("UNEXPECTED radius %d\n", d[q + 1]); }
    printf("%9lld %22lld %10lld %10lld %10.3f   (radius=%d)\n", q, q * (q + 3) / 2, q + 3, 3 * q, ms, d[q + 1]);
  }
}

int main(int argc, char** argv) {
  string mode = argc > 1 ? argv[1] : "verify";
  if (mode == "verify") verify(); else if (mode == "chain") chain(); else if (mode == "clique") clique();
  else { fprintf(stderr, "usage: %s verify|chain|clique\n", argv[0]); return 1; }
}
