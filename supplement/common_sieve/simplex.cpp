// Directed simplex upper sum for the two-step positive Buchstab leaves.
// All exponent comparisons and all rounding are integer.
#ifdef NDEBUG
#error "This arithmetic certificate requires assertions."
#endif
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using Int=unsigned __int128;
constexpr int denominator=10000, max_mesh=1024, max_rank=20;
constexpr Int scale=1000000000000000ULL;

Int ceil_div(Int n,Int d) { assert(d); return n/d+(n%d!=0); }
Int checked_product(Int a,Int b) { assert(!b || a<=~Int(0)/b); return a*b; }
Int checked_sum(Int a,Int b) { assert(a<=~Int(0)-b); return a+b; }
std::string integer_string(Int x) {
    if (!x) return "0";
    std::string s;
    while (x) { s.push_back(char('0'+x%10)); x/=10; }
    std::reverse(s.begin(),s.end());
    return s;
}

class Simplex {
    int mesh,rank,left,right,reflection,least;
    using Bits=std::bitset<max_mesh+1>;
    std::array<int,max_rank> lower{};
    std::array<std::array<Bits,max_rank+1>,max_rank+1> subsets{};
    std::array<Bits,max_rank+1> automatically_good{};
    std::array<Int,max_rank+1> factorial{};
    std::array<std::array<Int,max_rank+1>,max_rank+1> eulerian{};
    Int total_upper=0;
    std::uint64_t active=0;

    bool inside_window(int lo,int hi) const {
        return denominator*lo>=left*mesh && denominator*hi<=right*mesh;
    }

    bool good_cell(int h) const {
        for(int count=1;count<rank;count++) {
            int lo_add=std::max(0,h-rank+count),hi_add=std::min(h,count);
            int lo=(left*mesh+denominator-1)/denominator-lo_add;
            int hi=right*mesh/denominator-hi_add;
            for(int s=std::max(0,lo);s<=std::min(mesh,hi);s++)
                if(subsets[rank][count][s]) return true;
        }
        return false;
    }

    unsigned int leaves(unsigned int mask,int raw_sum,int count,int smallest,int h) const {
        if(count>=rank) return 0;
        int lo=raw_sum+lower[smallest]+std::max(0,h-rank+count);
        int hi=raw_sum+lower[smallest]+std::min(h,count)+1;
        unsigned int stop=0,expand=0;
        if(denominator*hi>=reflection*mesh && smallest==0) stop=1;
        if(denominator*lo<reflection*mesh && count+2<rank) {
            for(int i=1;i<smallest;i++) if(!(mask&(1u<<i)))
                for(int j=0;j<i;j++) if(!(mask&(1u<<j)))
                    expand+=leaves(mask|(1u<<i)|(1u<<j),raw_sum+lower[i]+lower[j],count+2,j,h);
        }
        return std::max(stop,expand);
    }

    void finish(int sum) {
        int h=mesh-sum;
        if(h<=0 || h>=rank || good_cell(h)) return;
        unsigned int count=0;
        for(int i=1;i<rank;i++) if(2*lower[i]<mesh)
            for(int j=0;j<i;j++)
                count+=leaves((1u<<i)|(1u<<j),lower[i]+lower[j],2,j,h);
        if(!count) return;
        Int weight=ceil_div(checked_product(checked_product(scale,mesh),eulerian[rank-1][h-1]),
                            factorial[rank-1]);
        for(int i=0;i<rank;i++) weight=ceil_div(weight,lower[i]);
        for(int i=0;i<rank;) {
            int last=i+1;
            while(last<rank && lower[last]==lower[i]) last++;
            weight=ceil_div(weight,factorial[last-i]);
            i=last;
        }
        total_upper=checked_sum(total_upper,checked_product(weight,count));
        active++;
    }

    void enumerate(int depth,int total,int minimum) {
        int remain=rank-depth;
        if(remain==0) { finish(total); return; }
        int maximum=(mesh-1-total)/remain;
        if(remain==1) minimum=std::max(minimum,mesh-rank+1-total);
        for(int v=minimum;v<=maximum;v++) {
            lower[depth]=v;
            subsets[depth+1][0]=subsets[depth][0];
            bool reject=false;
            for(int n=1;n<=depth+1;n++) {
                auto shifted=subsets[depth][n-1]<<v;
                if((shifted&automatically_good[n]).any()) { reject=true; break; }
                subsets[depth+1][n]=subsets[depth][n]|shifted;
            }
            if(reject) continue;
            for(int n=depth+2;n<=rank;n++) subsets[depth+1][n].reset();
            enumerate(depth+1,total+v,v);
        }
    }

public:
    Simplex(int n,int r,int a,int b):mesh(n),rank(r),left(a),right(b),
        reflection(denominator-a),least((b-a)*n/denominator) {
        assert(0<least && mesh<=max_mesh && rank<=max_rank);
        factorial[0]=1;
        eulerian[0][0]=1;
        for(int i=1;i<=rank;i++) {
            factorial[i]=checked_product(i,factorial[i-1]);
            for(int j=0;j<i;j++)
                eulerian[i][j]=(j?checked_product(i-j,eulerian[i-1][j-1]):0)
                    +checked_product(j+1,eulerian[i-1][j]);
        }
        for(int i=1;i<rank;i++)
            for(int s=0;s<=mesh;s++)
                if(inside_window(s,s+i)) automatically_good[i].set(s);
        subsets[0][0].set(0);
    }

    Int run() {
        enumerate(0,0,least);
        std::cout<<"rank="<<rank<<" upper="<<integer_string(total_upper)
                 <<" active="<<active<<"\n"<<std::flush;
        return total_upper;
    }
};

int main(int argc, char** argv) {
    if (argc != 6) return 2;
    int a = std::stoi(argv[1]), b = std::stoi(argv[2]), n = std::stoi(argv[3]);
    int first = std::stoi(argv[4]), last = std::stoi(argv[5]);
    assert(0 < a && a < b && b < 10000 && 3 <= first && first <= last && last <= max_rank);
    Int total = 0;
    for (int r=first; r<=last; ++r) total=checked_sum(total,Simplex(n,r,a,b).run());
    std::cout << "total=" << integer_string(total) << "\n";
}
