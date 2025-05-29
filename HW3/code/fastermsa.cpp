/*  fastermsa.cpp  
 *
 *  Compile:  g++ -O3 fastermsa.cpp -lpsapi -o fastermsa.exe  
 *  Run:      ./fastermsa.exe 
 */

#include <array>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <ctime>

/* ----- platform-independent peak-RSS helper --------------------------- */
#if defined(_WIN32)
    #ifndef NOMINMAX      // avoid redefinition warning
    #define NOMINMAX
    #endif
    #include <windows.h>
    #include <psapi.h>
    std::size_t getPeakRSS()
    {
        PROCESS_MEMORY_COUNTERS info;
        GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
        return static_cast<std::size_t>(info.PeakWorkingSetSize); // bytes
    }
#else
    #include <sys/resource.h>
    std::size_t getPeakRSS()
    {
        struct rusage ru{};
        getrusage(RUSAGE_SELF, &ru);
    #if defined(__APPLE__) && defined(__MACH__)
        return static_cast<std::size_t>(ru.ru_maxrss);            // already bytes
    #else
        return static_cast<std::size_t>(ru.ru_maxrss) * 1024;     // kB → bytes
    #endif
    }
#endif


/* -------------------------------------------------------------------------- */
/*  Scoring matrix and helpers                                                */
/* -------------------------------------------------------------------------- */
using namespace std;

constexpr array<array<int,5>,5> S = {{
    /* A   C   G   T   –  */
    {{  5,-4,-4,-4,-8 }},   // A
    {{ -4, 5,-4,-4,-8 }},   // C
    {{ -4,-4, 5,-4,-8 }},   // G
    {{ -4,-4,-4, 5,-8 }},   // T
    {{ -8,-8,-8,-8, 0 }}    // –
}};
inline int DELTA(int a,int b,int c) noexcept {
    // single branch-free look-up, inlined by the compiler
    return S[a][b] + S[b][c] + S[c][a];
}
/* translate A,C,G,T,- into 0…4 once and keep it integer thereafter         */
unordered_map<char,int> base_map = {{'A',0},{'C',1},{'G',2},{'T',3},{'-',4}};

/* -------------------------------------------------------------------------- */
/*  Result container                                                          */
/* -------------------------------------------------------------------------- */
struct Result {
    int    aln_score = 0;
    string aln[3];            // three aligned sequences
};

/* -------------------------------------------------------------------------- */
/*  Flat-buffer 2-layer DP for a 3-D slice                                    */
/* -------------------------------------------------------------------------- */
static vector<int> score_half(int m,int n,int p,
                              const vector<int>& X,const vector<int>& Y,const vector<int>& Z,
                              int xs,int ys,int zs,      // slice start indices
                              bool reverse)              // false = prefix, true = suffix
{
    const int pitch = p + 1;                      // row stride
    vector<int> buf(2 * (n+1) * (p+1));          // rolling buffers
    int* cur = buf.data();
    int* nxt = cur + (n+1)*(p+1);

    auto valX = [&](int i){ return X[ reverse ? xs + m - i     : xs + i - 1 ]; };
    auto valY = [&](int j){ return Y[ reverse ? ys + n - j     : ys + j - 1 ]; };
    auto valZ = [&](int k){ return Z[ reverse ? zs + p - k     : zs + k - 1 ]; };

    /* ---------- layer 0 (i == 0) -------------------------------------- */
    cur[0] = 0;
    for(int j=1;j<=n;++j) cur[j*pitch]      = cur[(j-1)*pitch] + DELTA(4,valY(j),4);
    for(int k=1;k<=p;++k) cur[k]            = cur[k-1]         + DELTA(4,4,valZ(k));

    for(int j=1;j<=n;++j)
        for(int k=1;k<=p;++k){
            int a = cur[(j-1)*pitch + k-1] + DELTA(4,valY(j),valZ(k));
            int b = cur[(j-1)*pitch + k  ] + DELTA(4,valY(j),4);
            int c = cur[(j  )*pitch + k-1] + DELTA(4,4,valZ(k));
            cur[j*pitch+k] = std::max({a,b,c});
        }

    /* ---------- layers 1 … m ------------------------------------------ */
    for(int i=1;i<=m;++i){
        nxt[0] = cur[0] + DELTA(valX(i),4,4);

        /* first column (vary j, k=0) */
        for(int j=1;j<=n;++j){
            int a = cur[(j-1)*pitch+0] + DELTA(valX(i),valY(j),4);
            int b = cur[(j  )*pitch+0] + DELTA(valX(i),4,4);
            int c = nxt[(j-1)*pitch+0] + DELTA(4,valY(j),4);
            nxt[j*pitch+0] = std::max({a,b,c});
        }
        /* first row (vary k, j=0) */
        for(int k=1;k<=p;++k){
            int a = cur[0*pitch+k-1] + DELTA(valX(i),4,valZ(k));
            int b = cur[0*pitch+k  ] + DELTA(valX(i),4,4);
            int c = nxt[0*pitch+k-1] + DELTA(4,4,valZ(k));
            nxt[0*pitch+k] = std::max({a,b,c});
        }
        /* inner j-k grid */
        for(int j=1;j<=n;++j)
            for(int k=1;k<=p;++k){
                int cd  = cur[(j-1)*pitch + k-1] + DELTA(valX(i),valY(j),valZ(k));
                int fij = cur[(j-1)*pitch + k  ] + DELTA(valX(i),valY(j),4);
                int fjk = nxt[(j-1)*pitch + k-1] + DELTA(4,valY(j),valZ(k));
                int fki = cur[(j  )*pitch + k-1] + DELTA(valX(i),4,valZ(k));
                int ei  = cur[(j  )*pitch + k  ] + DELTA(valX(i),4,4);
                int ej  = nxt[(j-1)*pitch + k  ] + DELTA(4,valY(j),4);
                int ek  = nxt[(j  )*pitch + k-1] + DELTA(4,4,valZ(k));
                nxt[j*pitch+k] = std::max({cd,fij,fjk,fki,ei,ej,ek});
            }
        std::swap(cur,nxt);                   // roll the buffers
    }
    /* return the final layer (cur) ------------------------------------- */
    return vector<int>(cur, cur + (n+1)*(p+1));
}

/* -------------------------------------------------------------------------- */
/*  Simple 3-D DP for small blocks (≤ 1 in x-dimension)                        */
/* -------------------------------------------------------------------------- */
static Result* basic_3Dalignment(const string& X,const string& Y,const string& Z,
                                 int sx,int ex,int sy,int ey,int sz,int ez)
{
    int m = ex-sx, n = ey-sy, p = ez-sz;
    auto idx3 = [&](int i,int j,int k){ return (i*(n+1)*(p+1)) + (j*(p+1)) + k; };

    vector<int> score((m+1)*(n+1)*(p+1));
    vector<uint8_t> trace((m+1)*(n+1)*(p+1));

    /* --- initialise edges ------------------------------------------------ */
    score[idx3(0,0,0)] = 0;
    for(int i=1;i<=m;++i){
        score[idx3(i,0,0)] = score[idx3(i-1,0,0)] +
                             DELTA(base_map[X[sx+i-1]],4,4);
        trace[idx3(i,0,0)] = 5;
    }
    for(int j=1;j<=n;++j){
        score[idx3(0,j,0)] = score[idx3(0,j-1,0)] +
                             DELTA(4,base_map[Y[sy+j-1]],4);
        trace[idx3(0,j,0)] = 6;
    }
    for(int k=1;k<=p;++k){
        score[idx3(0,0,k)] = score[idx3(0,0,k-1)] +
                             DELTA(4,4,base_map[Z[sz+k-1]]);
        trace[idx3(0,0,k)] = 7;
    }
    /* --- initialise faces ------------------------------------------------ */
    for(int i=1;i<=m;++i)
        for(int j=1;j<=n;++j){
            int a = score[idx3(i-1,j-1,0)] + DELTA(base_map[X[sx+i-1]],
                                                   base_map[Y[sy+j-1]],4);
            int b = score[idx3(i-1,j,0)]   + DELTA(base_map[X[sx+i-1]],4,4);
            int c = score[idx3(i  ,j-1,0)] + DELTA(4,base_map[Y[sy+j-1]],4);
            int best = std::max({a,b,c});
            score[idx3(i,j,0)] = best;
            trace[idx3(i,j,0)] = (best==a)?2:(best==b)?5:6;
        }
    for(int j=1;j<=n;++j)
        for(int k=1;k<=p;++k){
            int a = score[idx3(0,j-1,k-1)] + DELTA(4,base_map[Y[sy+j-1]],base_map[Z[sz+k-1]]);
            int b = score[idx3(0,j-1,k  )] + DELTA(4,base_map[Y[sy+j-1]],4);
            int c = score[idx3(0,j  ,k-1)] + DELTA(4,4,base_map[Z[sz+k-1]]);
            int best = std::max({a,b,c});
            score[idx3(0,j,k)] = best;
            trace[idx3(0,j,k)] = (best==a)?3:(best==b)?6:7;
        }
    for(int k=1;k<=p;++k)
        for(int i=1;i<=m;++i){
            int a = score[idx3(i-1,0,k-1)] + DELTA(base_map[X[sx+i-1]],4,base_map[Z[sz+k-1]]);
            int b = score[idx3(i-1,0,k  )] + DELTA(base_map[X[sx+i-1]],4,4);
            int c = score[idx3(i  ,0,k-1)] + DELTA(4,4,base_map[Z[sz+k-1]]);
            int best = std::max({a,b,c});
            score[idx3(i,0,k)] = best;
            trace[idx3(i,0,k)] = (best==a)?4:(best==b)?5:7;
        }

    /* --- full cube ------------------------------------------------------- */
    for(int i=1;i<=m;++i)
        for(int j=1;j<=n;++j)
            for(int k=1;k<=p;++k){
                int cd  = score[idx3(i-1,j-1,k-1)] + DELTA(base_map[X[sx+i-1]],
                                                           base_map[Y[sy+j-1]],
                                                           base_map[Z[sz+k-1]]);
                int fij = score[idx3(i-1,j-1,k  )] + DELTA(base_map[X[sx+i-1]],
                                                           base_map[Y[sy+j-1]],4);
                int fjk = score[idx3(i  ,j-1,k-1)] + DELTA(4,base_map[Y[sy+j-1]],
                                                           base_map[Z[sz+k-1]]);
                int fki = score[idx3(i-1,j  ,k-1)] + DELTA(base_map[X[sx+i-1]],4,
                                                           base_map[Z[sz+k-1]]);
                int ei  = score[idx3(i-1,j  ,k  )] + DELTA(base_map[X[sx+i-1]],4,4);
                int ej  = score[idx3(i  ,j-1,k  )] + DELTA(4,base_map[Y[sy+j-1]],4);
                int ek  = score[idx3(i  ,j  ,k-1)] + DELTA(4,4,base_map[Z[sz+k-1]]);
                int best = std::max({cd,fij,fjk,fki,ei,ej,ek});
                uint8_t tr = 1;
                if(best==fij) tr=2; else if(best==fjk) tr=3; else if(best==fki) tr=4;
                else if(best==ei) tr=5; else if(best==ej) tr=6; else if(best==ek) tr=7;
                score[idx3(i,j,k)] = best;
                trace[idx3(i,j,k)] = tr;
            }

    /* --- traceback ------------------------------------------------------- */
    Result* res = new Result;
    res->aln_score = score[idx3(m,n,p)];
    string& A = res->aln[0];
    string& B = res->aln[1];
    string& C = res->aln[2];
    int i=m,j=n,k=p;
    while(i||j||k){
        uint8_t t = trace[idx3(i,j,k)];
        switch(t){
            case 1: A += X[sx+--i]; B += Y[sy+--j]; C += Z[sz+--k]; break;
            case 2: A += X[sx+--i]; B += Y[sy+--j]; C += '-';      break;
            case 3: A += '-';       B += Y[sy+--j]; C += Z[sz+--k]; break;
            case 4: A += X[sx+--i]; B += '-';       C += Z[sz+--k]; break;
            case 5: A += X[sx+--i]; B += '-';       C += '-';       break;
            case 6: A += '-';       B += Y[sy+--j]; C += '-';       break;
            case 7: A += '-';       B += '-';       C += Z[sz+--k]; break;
        }
    }
    std::reverse(A.begin(),A.end());
    std::reverse(B.begin(),B.end());
    std::reverse(C.begin(),C.end());
    return res;
}

/* -------------------------------------------------------------------------- */
/*  Divide-and-conquer alignment (Hirschberg-like in 3-D)                     */
/* -------------------------------------------------------------------------- */
static Result* optimalAlignment3D(const string& X,const string& Y,const string& Z,
                                  const vector<int>& Xi,const vector<int>& Yi,const vector<int>& Zi,
                                  int sx,int ex,int sy,int ey,int sz,int ez)
{
    int m = ex - sx;
    if(m == 1)                       // base case => plain DP
        return basic_3Dalignment(X,Y,Z,sx,ex,sy,ey,sz,ez);

    int mid = m / 2;

    /* ---------- score prefix (sx … sx+mid) -------------------------- */
    vector<int> pre = score_half(mid,           // length of left half
                                 ey-sy, ez-sz,
                                 Xi,Yi,Zi,
                                 sx,sy,sz,false);

    /* ---------- score suffix (sx+mid … ex) -------------------------- */
    vector<int> suf = score_half(m - mid,       // length of right half
                                 ey-sy, ez-sz,
                                 Xi,Yi,Zi,
                                 sx+mid,sy,sz,true);

    int n = ey - sy, p = ez - sz, pitch = p + 1;
    int best = std::numeric_limits<int>::min();
    int bestJ = 0, bestK = 0;

    for(int j=0;j<=n;++j)
        for(int k=0;k<=p;++k){
            int val = pre[j*pitch + k] + suf[(n-j)*pitch + (p-k)];
            if(val > best){
                best = val; bestJ = j; bestK = k;
            }
        }

    /* ---------- recurse on the two halves --------------------------- */
    Result* left  = optimalAlignment3D(X,Y,Z,Xi,Yi,Zi,
                                       sx, sx+mid,
                                       sy, sy+bestJ,
                                       sz, sz+bestK);

    Result* right = optimalAlignment3D(X,Y,Z,Xi,Yi,Zi,
                                       sx+mid, ex,
                                       sy+bestJ, ey,
                                       sz+bestK, ez);

    Result* res = new Result;
    res->aln_score    = left->aln_score + right->aln_score;
    res->aln[0] = left->aln[0] + right->aln[0];
    res->aln[1] = left->aln[1] + right->aln[1];
    res->aln[2] = left->aln[2] + right->aln[2];
    delete left; delete right;
    return res;
}

/* -------------------------------------------------------------------------- */
/*  Main driver                                                               */
/* -------------------------------------------------------------------------- */
int main()
{
    /* ----- I/O file names ------------------------------------------------ */
    std::string num = "5"; // change this to 1, 2, 3, etc. for different test files
    std::string input_file  = "cpp_test" + num + ".txt";
    std::string output_file = "fast_cpp_result" + num + ".txt";

    /* ----- read sequences ------------------------------------------------ */
    std::string X, Y, Z;
    std::ifstream in(input_file);
    if (!in) { std::cerr << "Cannot open " << input_file << '\n'; return 1; }
    std::getline(in, X); std::getline(in, Y); std::getline(in, Z); in.close();
    std::cout << "Strings: " << X.size() << ", " << Y.size() << ", " << Z.size() << '\n';

    /* ----- encode once into 0…4 ----------------------------------------- */
    std::vector<int> Xi(X.size()), Yi(Y.size()), Zi(Z.size());
    for (size_t i = 0; i < X.size(); ++i) Xi[i] = base_map[X[i]];
    for (size_t i = 0; i < Y.size(); ++i) Yi[i] = base_map[Y[i]];
    for (size_t i = 0; i < Z.size(); ++i) Zi[i] = base_map[Z[i]];

    /* ----- align & time -------------------------------------------------- */
    std::clock_t t0 = std::clock();
    std::size_t   mem_before = getPeakRSS();

    Result* res = optimalAlignment3D(X, Y, Z,
                                     Xi, Yi, Zi,
                                     0, X.size(),
                                     0, Y.size(),
                                     0, Z.size());

    std::clock_t t1 = std::clock();
    std::size_t  mem_after = getPeakRSS();
    double secs = double(t1 - t0) / CLOCKS_PER_SEC;
    double peakMiB = mem_after / (1024.0 * 1024.0);

    /* ----- analyse alignment -------------------------------------------- */
    int perfect = 0;
    for (size_t i = 0; i < res->aln[0].size(); ++i)
        if (res->aln[0][i] == res->aln[1][i] && res->aln[0][i] == res->aln[2][i])
            ++perfect;

    /* ----- console report ----------------------------------------------- */
    std::cout << "Optimal alignment score: " << res->aln_score << '\n';
    for (int h = 0; h < 3; ++h)
        std::cout << "Length " << res->aln[h].size()
                  << " | preview: " << res->aln[h].substr(0, 10) << "...\n";
    std::cout << "Perfectly matched nucleotides: " << perfect << '\n';
    std::cout << "Runtime (seconds): "          << secs     << '\n';
    std::cout << "Peak memory: " << peakMiB << " MiB\n";

    /* ----- write output file -------------------------------------------- */
    std::ofstream out(output_file);
    out << "Optimal alignment score: "       << res->aln_score      << '\n'
        << "Length of the alignment: "       << res->aln[0].size()  << '\n'
        << "Perfectly matched nucleotides: " << perfect            << '\n'
        << "Runtime (seconds): "             << secs               << '\n'
        << "Peak memory (MiB): "             << peakMiB            << '\n'
        << "Actual aligned sequences:\n"
        << res->aln[0] << '\n' << res->aln[1] << '\n' << res->aln[2] << '\n';
    out.close();

    delete res;
    return 0;
}
