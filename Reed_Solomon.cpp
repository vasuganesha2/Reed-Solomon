#include<bits/stdc++.h>
// #include"Int.h"
#pragma GCC optimize("03")

using namespace std;

void __print(int x) {cerr << x;}
void __print(long x) {cerr << x;}
void __print(long long x) {cerr << x;}
void __print(unsigned x) {cerr << x;}
void __print(unsigned long x) {cerr << x;}
void __print(unsigned long long x) {cerr << x;}
void __print(float x) {cerr << x;}
void __print(double x) {cerr << x;}
void __print(long double x) {cerr << x;}
void __print(char x) {cerr << '\'' << x << '\'';}
void __print(const char *x) {cerr << '\"' << x << '\"';}
void __print(const string &x) {cerr << '\"' << x << '\"';}
void __print(bool x) {cerr << (x ? "true" : "false");}

template<typename T, typename V>
void __print(const pair<T, V> &x) {cerr << '{'; __print(x.first); cerr << ','; __print(x.second); cerr << '}';}
template<typename T>
void __print(const T &x) {int f = 0; cerr << '{'; for (auto &i: x) cerr << (f++ ? "," : ""), __print(i); cerr << "}";}
void _print() {cerr << "]\n";}
template <typename T, typename... V>
void _print(T t, V... v) {__print(t); if (sizeof...(v)) cerr << ", "; _print(v...);}
#ifndef ONLINE_JUDGE
#define debug(x...) {cerr << "[" << #x << "] = ["; _print(x);}
#define reach cerr << "reached" << endl
#else
#define debug(x...)
#define reach 
#endif

mt19937_64 RNG(chrono::steady_clock::now().time_since_epoch().count());

const int base = 1000000000;
const int base_digits = 9;
 
struct Int {
   
    vector<int> a;
    int sign;
   
    /*<arpa>*/
   
    int size(){
        if(a.empty()) return 0;
        int ans=(a.size()-1)*base_digits;
        int ca=a.back();
        while(ca)
            ans++,ca/=10;
        return ans;
    }
   
    Int operator ^(const Int &v){
        Int ans=1,a=*this,b=v;
        while(!b.isZero()){
            if(b%2)
            ans*=a;
            a*=a,b/=2;
        }
        return ans;
    }
   
    string to_string(){
        stringstream ss;
        //ss << *this;
        string s;
        ss >> s;
        return s;
    }
   
    int sumof(){
        string s = to_string();
        int ans = 0;
        for(auto c : s)  ans += c - '0';
        return ans;
    }
   
    /*</arpa>*/
   
    Int() :
    sign(1) {
    }
 
    Int(long long v) {
        *this = v;
    }
 
    Int(const string &s) {
        read(s);
    }
 
    void operator=(const Int &v) {
        sign = v.sign;
        a = v.a;
    }
 
    void operator=(long long v) {
        sign = 1;
        a.clear();
        if (v < 0)
            sign = -1, v = -v;
        for (; v > 0; v = v / base)
            a.push_back(v % base);
    }
 
    Int operator+(const Int &v) const {
        if (sign == v.sign) {
            Int res = v;
            for (int i = 0, carry = 0; i < (int) max(a.size(), v.a.size()) || carry; ++i) {
                if (i == (int) res.a.size())
                    res.a.push_back(0);
                res.a[i] += carry + (i < (int) a.size() ? a[i] : 0);
                carry = res.a[i] >= base;
                if (carry)
                    res.a[i] -= base;
            }
            return res;
        }
        return *this - (-v);
    }
 
    Int operator-(const Int &v) const {
        if (sign == v.sign) {
            if (abs() >= v.abs()) {
                Int res = *this;
                for (int i = 0, carry = 0; i < (int) v.a.size() || carry; ++i) {
                    res.a[i] -= carry + (i < (int) v.a.size() ? v.a[i] : 0);
                    carry = res.a[i] < 0;
                    if (carry)
                    res.a[i] += base;
                }
                res.trim();
                return res;
            }
            return -(v - *this);
        }
        return *this + (-v);
    }
 
    void operator*=(int v) {
        if (v < 0)
            sign = -sign, v = -v;
        for (int i = 0, carry = 0; i < (int) a.size() || carry; ++i) {
            if (i == (int) a.size())
                a.push_back(0);
            long long cur = a[i] * (long long) v + carry;
            carry = (int) (cur / base);
            a[i] = (int) (cur % base);
            //asm("divl %%ecx" : "=a"(carry), "=d"(a[i]) : "A"(cur), "c"(base));
        }
        trim();
    }
 
    Int operator*(int v) const {
        Int res = *this;
        res *= v;
        return res;
    }
 
    void operator*=(long long v) {
        if (v < 0)
            sign = -sign, v = -v;
        for (int i = 0, carry = 0; i < (int) a.size() || carry; ++i) {
            if (i == (int) a.size())
                a.push_back(0);
            long long cur = a[i] * (long long) v + carry;
            carry = (int) (cur / base);
            a[i] = (int) (cur % base);
            //asm("divl %%ecx" : "=a"(carry), "=d"(a[i]) : "A"(cur), "c"(base));
        }
        trim();
    }
 
    Int operator*(long long v) const {
        Int res = *this;
        res *= v;
        return res;
    }
 
    friend pair<Int, Int> divmod(const Int &a1, const Int &b1) {
        int norm = base / (b1.a.back() + 1);
        Int a = a1.abs() * norm;
        Int b = b1.abs() * norm;
        Int q, r;
        q.a.resize(a.a.size());
 
        for (int i = a.a.size() - 1; i >= 0; i--) {
            r *= base;
            r += a.a[i];
            int s1 = r.a.size() <= b.a.size() ? 0 : r.a[b.a.size()];
            int s2 = r.a.size() <= b.a.size() - 1 ? 0 : r.a[b.a.size() - 1];
            int d = ((long long) base * s1 + s2) / b.a.back();
            r -= b * d;
            while (r < 0)
                r += b, --d;
            q.a[i] = d;
        }
 
        q.sign = a1.sign * b1.sign;
        r.sign = a1.sign;
        q.trim();
        r.trim();
        return make_pair(q, r / norm);
    }
 
    Int operator/(const Int &v) const {
        return divmod(*this, v).first;
    }
 
    Int operator%(const Int &v) const {
        return divmod(*this, v).second;
    }
 
    void operator/=(int v) {
        if (v < 0)
            sign = -sign, v = -v;
        for (int i = (int) a.size() - 1, rem = 0; i >= 0; --i) {
            long long cur = a[i] + rem * (long long) base;
            a[i] = (int) (cur / v);
            rem = (int) (cur % v);
        }
        trim();
    }
 
    Int operator/(int v) const {
        Int res = *this;
        res /= v;
        return res;
    }
 
    int operator%(int v) const {
        if (v < 0)
            v = -v;
        int m = 0;
        for (int i = a.size() - 1; i >= 0; --i)
            m = (a[i] + m * (long long) base) % v;
        return m * sign;
    }
 
    void operator+=(const Int &v) {
        *this = *this + v;
    }
   
    void operator-=(const Int &v) {
        *this = *this - v;
    }
   
    void operator*=(const Int &v) {
        *this = *this * v;
    }
   
    void operator/=(const Int &v) {
        *this = *this / v;
    }
   
    Int operator ++(){
        *this += 1;
        return *this;
    }
   
    Int operator ++(int){
        Int temp = *this;
        *this += 1;
        return temp;
    }
   
    Int operator --(){
        *this -= 1;
        return *this;
    }
   
    Int operator --(int){
        Int temp = *this;
        *this -= 1;
        return temp;
    }
 
    bool operator<(const Int &v) const {
        if (sign != v.sign)
            return sign < v.sign;
        if (a.size() != v.a.size())
           return a.size() * sign < v.a.size() * v.sign;
        for (int i = a.size() - 1; i >= 0; i--)
            if (a[i] != v.a[i])
                return a[i] * sign < v.a[i] * sign;
        return false;
    }
 
    bool operator>(const Int &v) const {
        return v < *this;
    }
   
    bool operator<=(const Int &v) const {
        return !(v < *this);
    }
   
    bool operator>=(const Int &v) const {
        return !(*this < v);
    }
   
    bool operator==(const Int &v) const {
        return !(*this < v) && !(v < *this);
    }
   
    bool operator!=(const Int &v) const {
        return *this < v || v < *this;
    }
 
    void trim() {
        while (!a.empty() && !a.back())
            a.pop_back();
        if (a.empty())
            sign = 1;
    }
 
    bool isZero() const {
        return a.empty() || (a.size() == 1 && !a[0]);
    }
 
    Int operator-() const {
        Int res = *this;
        res.sign = -sign;
        return res;
    }
 
    Int abs() const {
        Int res = *this;
        res.sign *= res.sign;
        return res;
    }
 
    long long longValue() const {
        long long res = 0;
        for (int i = a.size() - 1; i >= 0; i--)
            res = res * base + a[i];
        return res * sign;
    }
 
    friend Int gcd(const Int &a, const Int &b) {
        return b.isZero() ? a : gcd(b, a % b);
    }
   
    friend Int lcm(const Int &a, const Int &b) {
        return a / gcd(a, b) * b;
    }
 
    void read(const string &s) {
        sign = 1;
        a.clear();
        int pos = 0;
        while (pos < (int) s.size() && (s[pos] == '-' || s[pos] == '+')) {
            if (s[pos] == '-')
                sign = -sign;
            ++pos;
        }
       
        for (int i = s.size() - 1; i >= pos; i -= base_digits) {
            int x = 0;
            for (int j = max(pos, i - base_digits + 1); j <= i; j++)
                x = x * 10 + s[j] - '0';
            a.push_back(x);
        }
        trim();
    }
 
    friend istream& operator>>(istream &stream, Int &v) {
        string s;
        stream >> s;
        v.read(s);
        return stream;
    }
 
    friend ostream& operator<<(ostream &stream, const Int &v) {
        if (v.sign == -1)
           stream << '-';
        stream << (v.a.empty() ? 0 : v.a.back());
        for (int i = (int) v.a.size() - 2; i >= 0; --i)
            stream << setw(base_digits) << setfill('0') << v.a[i];
        return stream;
    }
 
    static vector<int> convert_base(const vector<int> &a, int old_digits, int new_digits) {
        vector<long long> p(max(old_digits, new_digits) + 1);
        p[0] = 1;
        for (int i = 1; i < (int) p.size(); i++)
            p[i] = p[i - 1] * 10;
        vector<int> res;
        long long cur = 0;
        int cur_digits = 0;
        for (int i = 0; i < (int) a.size(); i++) {
            cur += a[i] * p[cur_digits];
            cur_digits += old_digits;
            while (cur_digits >= new_digits) {
                res.push_back(int(cur % p[new_digits]));
                cur /= p[new_digits];
                cur_digits -= new_digits;
            }
        }
        res.push_back((int) cur);
        while (!res.empty() && !res.back())
            res.pop_back();
        return res;
    }
 
    typedef vector<long long> vll;
 
    static vll karatsubaMultiply(const vll &a, const vll &b) {
        int n = a.size();
        vll res(n + n);
        if (n <= 32) {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    res[i + j] += a[i] * b[j];
            return res;
        }
 
        int k = n >> 1;
        vll a1(a.begin(), a.begin() + k);
        vll a2(a.begin() + k, a.end());
        vll b1(b.begin(), b.begin() + k);
        vll b2(b.begin() + k, b.end());
 
        vll a1b1 = karatsubaMultiply(a1, b1);
        vll a2b2 = karatsubaMultiply(a2, b2);
 
        for (int i = 0; i < k; i++)
            a2[i] += a1[i];
        for (int i = 0; i < k; i++)
            b2[i] += b1[i];
 
        vll r = karatsubaMultiply(a2, b2);
       
        for (int i = 0; i < (int) a1b1.size(); i++)
            r[i] -= a1b1[i];
        for (int i = 0; i < (int) a2b2.size(); i++)
            r[i] -= a2b2[i];
 
        for (int i = 0; i < (int) r.size(); i++)
            res[i + k] += r[i];
        for (int i = 0; i < (int) a1b1.size(); i++)
            res[i] += a1b1[i];
        for (int i = 0; i < (int) a2b2.size(); i++)
            res[i + n] += a2b2[i];
        return res;
    }
 
    Int operator*(const Int &v) const {
        vector<int> a6 = convert_base(this->a, base_digits, 6);
        vector<int> b6 = convert_base(v.a, base_digits, 6);
       
        vll a(a6.begin(), a6.end());
        vll b(b6.begin(), b6.end());
       
        while (a.size() < b.size())
            a.push_back(0);
        while (b.size() < a.size())
            b.push_back(0);
        while (a.size() & (a.size() - 1))
            a.push_back(0), b.push_back(0);
           
        vll c = karatsubaMultiply(a, b);
        Int res;
        res.sign = sign * v.sign;
        for (int i = 0, carry = 0; i < (int) c.size(); i++) {
            long long cur = c[i] + carry;
            res.a.push_back((int) (cur % 1000000));
            carry = (int) (cur / 1000000);
        }
        res.a = convert_base(res.a, 6, base_digits);
        res.trim();
        return res;
    }
   
    //Added part.
   
    friend Int max(const Int &a,const Int &b){
        if(a>b){
            return a;
        }
        return b;
    }
   
    friend Int max(const Int &a,const int32_t &B){
        Int b = B;
        return max(a,b);
    }
   
    friend Int max(const Int &a,const int64_t &B){
        Int b = B;
        return max(a,b);
    }
   
    friend Int min(const Int &a,const Int &b){
        if(a<b){
            return a;
        }
        return b;
    }
   
    friend Int min(const Int &a,const int32_t &B){
        Int b = B;
        return min(a,b);
    }
   
    friend Int min(const Int &a,const int64_t &B){
        Int b = B;
        return min(a,b);
    }
};

Int rndm()
{
    Int a(1);
    int d = 10;

    while(d --)
    {
        a *= Int(10);
        a += Int((int)(RNG() % 10));
    }

    return a;
}



Int get_rndm(Int l, Int r)
{
    if(l == r) return l;
    return l + (rndm() % (r - l));
}


Int binpower(Int _base, Int e, Int mod) 
{
    Int result = 1;
    _base = _base % mod;
    while (e > Int(0)) {
        if (e % Int(2) == Int(1))
            result = (Int)result * _base % mod;
        _base = (Int)_base * _base % mod;
        e /= Int(2);
    }
    return result;
}



bool check_composite(Int n, Int a, Int d, int s) 
{
    Int x = binpower(a, d, n);
    if (x == Int(1) || x == n - Int(1))
        return false;
    for (int r = 1; r < s; r++) {
        x = x * x % n;
        if (x == n - Int(1))
            return false;
    }
    return true;
}

bool MillerRabin(Int n, int iter=5) 
{
    if (n < Int(4))
        return n == Int(2) || n == Int(3);

    int s = 0;
    Int d = n - Int(1);
    while ((d % Int(2)) == Int(0)) {
        d /= Int(2);
        s++;
    }

    for (int i = 0; i < iter; i++) {
        Int a = Int(2) + rndm() % (n - Int(3));
        if (check_composite(n, a, d, s))
            return false;
    }
    return true;
}





vector<Int> getK(int k)
{
    vector<Int> primes;
    set<Int> s;
    long long a = 11;
    
    while (primes.size() < k) 
    {
        if (MillerRabin(a) && s.find(a) == s.end())
        {
            primes.push_back(a);
            s.insert(a);
        }
        a++;
    }
    return primes;
}


    
Int M = Int("1" + string("0", 1000));
double mu=0.3;
int k=0;
vector<Int> primes;



Int extendedGCD(Int a, Int b, Int &x, Int &y) 
{
    if (a == Int(0)) 
    {
        x = Int(0);
        y = Int(1);
        return b;
    }
    Int x1, y1;
    Int gcd = extendedGCD(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return gcd;
}



Int modInverse(Int a, Int m) 
{
    Int x, y;
   
    Int gcd = extendedGCD(a, m, x, y);

    if (gcd != Int(1)) 
    {
        throw invalid_argument("Modular inverse does not exist.");
    }
    return (x % m + m) % m;
}



Int chineseRemainderTheorem(vector<Int>& a, vector<Int>& m) 
{
    Int M(1);
    for (int i = 0; i < m.size(); ++i) {
        M *= m[i];
    }
    Int x(0);
    for (int i = 0; i < m.size(); ++i) 
    {
        Int Mi = M / m[i];

        Int Ni = modInverse(Mi, m[i]);

        x += a[i] * Mi * Ni;
    }
    return x % M;
}



void GlobalSetup(Int mu, Int M) 
{
    k = 1000;
    primes = getK(k);
    sort(primes.rbegin(), primes.rend()); 
}

vector<Int> Transmit(vector<Int>& a) 
{
    // Choose l at random from the range [0, µ.k]
    int l = rand() % static_cast<int>(mu * k);

    unordered_set<long long> I;
    while (I.size() < l) 
    {
        long long index = rand() % k;
        I.insert(index);
    }
    
    // Construct bi
    vector<Int> b(k);

    for (long long i = 0; i < k; ++i) 
    {
        if (I.find(i) != I.end()) 
        {
            Int pi_minus_1 = primes[i] - 1;
            Int temp1 = a[i];
            while (temp1 == a[i]) 
            {
                temp1 = rndm() % pi_minus_1;
            }
            b[i] = temp1;
        }
        else 
        {
            b[i] = a[i];
        }
    }
    return b;
}

vector<Int> ReedSolomonSend(Int a) 
{
    vector<Int> alist;
    for (Int i : primes) 
    {
        alist.push_back(a % i);
    }
    return Transmit(alist); 
}


vector<vector<Int>> EEA(Int a, Int b) 
{
    Int r = a;
    Int rd = b;
    Int s = 1;
    Int sd = 0;
    Int t = 0;
    Int td = 1;
    vector<Int> rlist, slist, tlist;
    rlist.push_back(r);
    slist.push_back(s);
    tlist.push_back(t);

    while (rd != 0) {
        Int q = r / rd;
        Int rdd = r % rd;
        r = rd;
        rd = rdd;
        Int sm = s - q * sd;
        Int tm = t - q * td;
        s = sd;
        sd = sm;
        t = td;
        td = tm;
        rlist.push_back(r);
        slist.push_back(s);
        tlist.push_back(t);
    }
    return {rlist, slist, tlist};
}



Int ReedSolomonReceive(vector<Int>& b) 
{
    Int P = 1;
    Int n = 1;
    int l = mu * k;


    for (int i = 0; i < l; ++i) 
    {
        P *= primes[i];
    }
    
    for (int i = 0; i < k; ++i)
    {
        n *= primes[i];
    }



    Int bval = chineseRemainderTheorem(b, primes);


    vector<vector<Int>> a = EEA(n, bval);

    Int rst = M * P;
    Int tst = P;

    int j = 0;

    for (int i = 0; i < a[0].size(); i++) 
    {
        if (a[0][i] <= rst) 
        {
            j = i;
            break;
        }
    }

    Int p = a[0][j];
    Int q = a[2][j];

    if(p % q == 0) return p / q;
    return (Int)-1;
}



int32_t main()
{
    auto start = std::chrono::high_resolution_clock::now();
    GlobalSetup(mu,M);
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Output the duration
    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;

    // Do some computation or task here...

    // Get the current time point again
    Int a = rndm();
  
    vector<Int> b = ReedSolomonSend(a);
    cout << ReedSolomonReceive(b) << endl;
}