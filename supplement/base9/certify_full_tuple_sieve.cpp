// Rigorous outer-cell upper sum for the original Maynard positive leaves,
// after deleting every complete prime tuple admitting a Type II subset.
// Exact integer arithmetic only.  Parameter denominator 10000.
#ifdef NDEBUG
#error "Compile this mathematical certificate with assertions enabled."
#endif
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
using U=unsigned __int128;using I=long long;
constexpr int MAXN=2048,MAXK=24;
constexpr I END=10000,AN=3155,BN=3691,EN=BN-AN,AC=END-AN,BC=END-BN;
constexpr unsigned long long FP=1000000000000000000ULL;
std::string str(U u){if(!u)return "0";std::string r;while(u){r.push_back('0'+u%10);u/=10;}std::reverse(r.begin(),r.end());return r;}
U cd(U x,U y){assert(y);return x/y+(x%y!=0);}
U mul(U x,U y){assert(!y||x<=(~U(0))/y);return x*y;}
U add(U x,U y){assert(x<=(~U(0))-y);return x+y;}
void inc(unsigned long long &x){assert(x<std::numeric_limits<unsigned long long>::max());++x;}
std::string fixed(U u){auto s=str(u);if(s.size()<=18)s=std::string(19-s.size(),'0')+s;return s.substr(0,s.size()-18)+"."+s.substr(s.size()-18);}
struct Range{I lo,hi;};
struct Result{int k;U upper=0;unsigned long long nodes=0,cells=0,active=0,good_early=0,good_final=0;};
class Certificate{
 int N,k,mn;std::array<int,MAXK> j{};
 std::array<std::array<std::bitset<MAXN>,MAXK>,MAXK> dp{};
 std::array<std::bitset<MAXN>,MAXK> good{};
 std::array<U,MAXK> fact{};
 std::array<std::array<U,MAXK>,MAXK> euler{};
 Result out;
 bool gt(Range r,I rhs)const{return END*r.hi>rhs*N;}
 bool lt(Range r,I rhs)const{return END*r.lo<rhs*N;}
 Range sum(int h,std::initializer_list<int> inds,bool doubledlast=false)const{
  int n=inds.size();I b=0;int last=-1;for(int i:inds){b+=j[i];last=i;}if(doubledlast)b+=j[last];
  return {b+std::max(0,h-k+n),b+std::min(h,n)+(doubledlast?1:0)};
 }
 unsigned count(int h){
  unsigned ans=0;
  for(int x=1;x<k;x++){
   auto uv=sum(h,{x,0});auto u2v=sum(h,{x,0},true);
   bool low=END*j[x]<AN*N;
   bool d2=low&&gt(uv,BN)&&lt(uv,BC)&&gt(u2v,AC)&&lt(u2v,END);
   if(k==3&&d2)ans++; // I2: exactly three primes.
   if(k>=4&&d2){int w=(x==1?2:1);if(lt(sum(h,{x,0,w},true),END))ans++;} // I4; D3 empty.
   bool cen=END*(j[x]+1)>BN*N&&2*j[x]<N&&lt(u2v,END);
   if(cen&&gt(uv,AC))ans++; // I7
   if(cen&&gt(uv,BN)&&lt(uv,BC)&&gt(u2v,AC))ans++; // I8
  }
  if(k>=5)for(int x=1;x<k;x++)for(int y=x+1;y<k;y++)for(int z=y+1;z<k;z++){
   // marked t0=t, t_x=w, t_y=v, t_z=u, all remaining primes exceed t.
   auto uv=sum(h,{z,y});
   bool large=gt(uv,BN)&&lt(uv,BC)&&lt(sum(h,{z,y},true),AC)
     &&lt(sum(h,{z,y,x},true),END)&&lt(sum(h,{z,y,x,0},true),END);
   if(large&&END*j[z]<AN*N)ans++; // I5
   if(END*j[z]<AN*N&&lt(uv,AN))ans++; // I6
   if(large&&END*(j[z]+1)>BN*N&&2*j[z]<N)ans++; // I9
  }
  return ans;
 }
 bool final_good(int h){
  for(int sz=1;sz<k;sz++){
   int extra_lo=std::max(0,h-k+sz),extra_hi=std::min(h,sz);
   for(auto interval: {std::pair<I,I>{AN,BN},std::pair<I,I>{BC,AC}}){
    I lo=(interval.first*N+END-1)/END-extra_lo;
    I hi=interval.second*N/END-extra_hi;
    lo=std::max<I>(0,lo);hi=std::min<I>(N,hi);
    for(I s=lo;s<=hi;s++)if(dp[k][sz][s])return true;
   }
  }return false;
 }
 void finish(int total){
  int h=N-total;if(h<=0||h>=k)return;
  inc(out.cells);
  if(final_good(h)){inc(out.good_final);return;}
  unsigned c=count(h);if(!c)return;
  U weight=cd(mul(mul(FP,U(N)),euler[k-1][h-1]),fact[k-1]);
  for(int i=0;i<k;i++)weight=cd(weight,U(j[i]));
  for(int i=0;i<k;){int end=i+1;while(end<k&&j[end]==j[i])end++;weight=cd(weight,fact[end-i]);i=end;}
  out.upper=add(out.upper,mul(weight,c));inc(out.active);
 }
 void rec(int depth,int total,int start){
  inc(out.nodes);int rem=k-depth;
  if(rem==0){finish(total);return;}
  int high=(N-1-total)/rem;
  if(rem==1)start=std::max(start,N-k+1-total);
  for(int v=start;v<=high;v++){
   j[depth]=v;
   bool kill=false;
   dp[depth+1][0]=dp[depth][0];
   for(int sz=1;sz<=depth+1;sz++){
    auto news=dp[depth][sz-1]<<v;
    if((news&good[sz]).any()){kill=true;break;}
    dp[depth+1][sz]=dp[depth][sz]|news;
   }
   if(kill){inc(out.good_early);continue;}
   for(int sz=depth+2;sz<=k;sz++)dp[depth+1][sz].reset();
   rec(depth+1,total+v,v);
  }
 }
public:
 Certificate(int n,int K):N(n),k(K),mn(EN*n/END),out{K}{
  assert(mn>0&&N<MAXN&&k<MAXK);fact[0]=1;
  for(int i=1;i<MAXK;i++)fact[i]=mul(fact[i-1],i);
  euler[0][0]=1;
  for(int nn=1;nn<MAXK;nn++)for(int m=0;m<nn;m++){
   euler[nn][m]=add((m?mul(nn-m,euler[nn-1][m-1]):0),mul(m+1,euler[nn-1][m]));
  }
  for(int sz=1;sz<k;sz++)for(int s=0;s<N;s++){
   if((END*s>=AN*N&&END*(s+sz)<=BN*N)||(END*s>=BC*N&&END*(s+sz)<=AC*N))good[sz].set(s);
  }
  dp[0][0].set(0);
 }
 Result run(){rec(0,0,mn);return out;}
};
int main(int argc,char**argv){
 int N=128,klo=3,khi=END/EN;std::string path;
 for(int a=1;a<argc;a++){std::string z=argv[a];if(z=="--n")N=std::stoi(argv[++a]);else if(z=="--kmin")klo=std::stoi(argv[++a]);else if(z=="--kmax")khi=std::stoi(argv[++a]);else if(z=="--out")path=argv[++a];else throw std::runtime_error("bad option");}
 U sum=0;std::vector<Result>res;
 for(int k=klo;k<=khi;k++){
  auto r=Certificate(N,k).run();sum=add(sum,r.upper);res.push_back(r);
  std::cout<<"k="<<k<<" upper="<<fixed(r.upper)<<" cells="<<r.cells<<" active="<<r.active<<" nodes="<<r.nodes<<" cumulative="<<fixed(sum)<<"\n"<<std::flush;
 }
 std::ostringstream o;o<<"{\"scope\":\"EXACT OUTER SUM; requires accompanying domain proof\",\"a\":\"631/2000\",\"b\":\"3691/10000\",\"eta\":\"67/1250\",\"N\":"<<N<<",\"denominator\":\""<<FP<<"\",\"total_upper_numerator\":\""<<str(sum)<<"\",\"total_upper_decimal\":\""<<fixed(sum)<<"\",\"results\":[";
 for(unsigned z=0;z<res.size();z++){auto r=res[z];if(z)o<<",";o<<"{\"k\":"<<r.k<<",\"upper_numerator\":\""<<str(r.upper)<<"\",\"upper_decimal\":\""<<fixed(r.upper)<<"\",\"nodes\":"<<r.nodes<<",\"cells\":"<<r.cells<<",\"active\":"<<r.active<<",\"good_early\":"<<r.good_early<<",\"good_final\":"<<r.good_final<<"}";}
 o<<"]}\n";if(path.empty())std::cout<<o.str();else{std::ofstream f(path);f<<o.str();}
}
