import random
import gmpy2
import math

def ab(a):
    return gmpy2.mpz(a)

def power(a,b):
    c=ab(1)
    n = a
    while(b > 0):
        if b & 1:
            c = c * a
        a = a * a
        b = b >> 1
    return c

def powm(a,b,mod1):
    c = ab(1)
    n = a
    while(b > 0):
        if b & 1:
            c = mod(c * a,mod1)
        a = mod(a * a, mod1)
        b = b >> 1
    return c


def begcd(a,b):
    r = a
    rd = b
    e = (0)
    while (r & 1 == 0 and rd & 1 == 0):
        r = r >> 1
        rd = rd >> 1
        e = e + 1
    ad = r
    bd = rd
    s = ab(1)
    t = ab(0)
    sd = ab(0)
    td = ab(1)
    while rd:
        while r & 1 == 0:
            r = r >> 1
            if s & 1 == 0 and t & 1 == 0:
                s = s >> 1
                t = t >> 1
            else:
                s = (s + bd) >> 1
                t = (t - ad) >> 1
        while rd & 1 == 0:
            rd = rd >> 1
            if sd & 1 == 0 and td & 1 == 0:
                sd = sd >> 1
                td = td >> 1
            else:
                sd = (sd + bd) >> 1
                td = (td - ad) >> 1
        if r > rd:
            r,rd = rd,r
            s,sd = sd,s
            t,td = td,t
        rd = rd - r
        sd = sd - s
        td = td - t
    
    ans = ab(1)
    ans = ans << e      
    ans = ans * r
    return ans,s,t

def modinv(a,mod):
    d,s,t = begcd(a,mod)
    assert(d == 1)
    return s % mod

def mod(a,b):
    ans=a - b * (a // b)
    if ans < 0:
        ans += b
    return ans

def crt(primes_list,a_list):
    n = ab(1)
    for i in primes_list:
        n = n * i
    k = len(primes_list)
    nstar = [ab(1) for i in range(k)]
    e =[ab(1) for i in range(k)]
    for i in range(k):
        nstar[i] = n // primes_list[i]
        b = mod(nstar[i],primes_list[i])
        t = modinv(b,primes_list[i])
        e[i] = nstar[i] * t
    ans = ab(0)
    for i in range(k):
        ans = mod(ans + a_list[i] * e[i],n)
    return ans

def rand(a,b):
    x = min(a,b)
    y = max(a,b)
    a = x
    b = y
    return ab(a + random.random() * (b - a))



def check_composite(a, h, t, n):
    x = pow(a, t, n)
    if x == 1 or x == n - 1:
        return False
    for _ in range(1, h):
        x = (x * x) % n
        if x == n - 1:
            return False
    return True

def miller_rabin(n, threshold):
    if (n % 2 == 0 and n > 2) or n == 1:
        return False
    if n == 2:
        return True
    nd = n - 1
    h = 0
    while nd & 1 == 0:
        nd >>= 1
        h += 1
    t = (n - 1) // (1 << h)
    for _ in range(threshold):
        a = random.randint(2, n - 1)
        if math.gcd(a, n) != 1:
            return False
        if check_composite(a, h, t, n):
            return False
    return True



def getK(k):
    primes=[]
    s=set()
    a=11
    while len(primes)<k:
        if miller_rabin(a, 10) and a not in s:
            primes.append(a)
            s.add(a)
        a+=1
    return primes
 
 
 

    
    
M=ab('1'+'0'*1000)
mu=0.3
k=1000
primes=[]





def GlobalSetup(mu, M):
    global k, primes
    k = 1000
    # primes=getkprimes(k,ab(gmpy2.isqrt(M)),M-1)
    primes=getK( k )
    primes.sort(reverse = 1)



def ReedSolomonSend(p):
    global primes
    List = []
    for i in primes:
        List.append(mod(p , i))
    return Transmit(List)


def Transmit(p):
    global mu, k, primes
    l = random.randint(0, int(mu * k))
    I = set()
    while len(I) < l:
        temp = random.randint(0, k - 1)
        if temp not in I:
            I.add(temp)
    b = [0] * k
    for i in range(len(p)):
        if i in I:
            temp1 = p[i]
            while temp1 == p[i]:
                temp1 = random.randint(0, primes[i] - 1)
            b[i] = temp1
        else:
            b[i] = p[i]
    return b



def EEA(a, b):
    r = a
    rd = b
    s = 1
    sd = 0
    t = 0
    td = 1
    r_values = [r]
    s_values = [s]
    t_values = [t]
    while rd:
        quotient, remainder = divmod(r, rd)
        r, s, t, rd, sd, td = rd, sd, td, remainder, s - quotient * sd, t - quotient * td
        r_values.append(r)
        s_values.append(s)
        t_values.append(t)
    return r_values, s_values, t_values


def reed_solomon_receive(b):
    global M, primes, mu, k
    P = 1
    n = 1
    l = int(mu * k)
    primes.sort(reverse=True) 
    for i in range(l):
        P *= primes[i]
    for i in range(k):
        n *= primes[i]
    # assert(n > 2 * M * P * P)
    b_val = crt(primes, b)
    r_list, s_list, t_list = EEA(n, b_val)
    rst = M * P
    tst = P
    j = 0
    for i in range(len(r_list)):
        if r_list[i] <= rst:
            j = i
            break
    rd = r_list[j]
    sd = s_list[j]
    td = t_list[j]
    
    try:
        if (rd % td) == 0:
            return rd // td
        else:
            return -1
    except:
        return -1


    
    
def main():
    GlobalSetup(mu,M)
   
    a = int(input("Enter a number: "))
    #mu = float(input("Enter a mu: "))
    b=ReedSolomonSend(a)
    print(a)
    #print(mu)
    print(reed_solomon_receive(b))
    
if __name__ == "__main__":
    main()
