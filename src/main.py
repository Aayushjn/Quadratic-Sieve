from math import fabs, ceil, sqrt, exp, log
from shanks import STonelli
from itertools import chain
from random import randint
import helpers

def is_probable_prime(a):
    """ Perform Rabin-Miller primality test to determine whether given number
        is prime. Return True if number is very likely to be a prime, and False
        if it is definitely composite
    """
    if a == 2:
        return True

    if a == 1 or a % 2 == 0:
        return False

    return rabin_miller_primality_test(a, 50)

def rabin_miller_primality_test(a, iterations):
    """ Rabin Miller primality test
    """
    r, s = 0, a - 1

    while s % 2 == 0:
        r += 1
        s //= 2

    for _ in range(iterations):
        n = randint(2, a - 1)
        x = pow(n, s, a)
        if x == 1 or x == a - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, a)
            if x == a - 1:
                break
        else:
            return False
    return True

def check_perfect_power(n):
    """ Check if the given integer is a perfect power. If yes, return (r, b)
        such that r^b = n, otherwise return None.
        Assume that global small_primes has already been initialised and that n
        does not have any prime factors from small_primes.
    """
    prime = small_primes[-1]
    for p in small_primes:
        pth_root = helpers.kth_iroot(n, p)
        if pth_root < prime:
            break
        if pth_root ** p == n:
            return (pth_root, p)
    return None

def mprint(M): #prints a matrix in readable form
    for row in M:
        print(row)
        
def prime_gen(n): # sieve of Eratosthenes, generates primes up to a bound n
    if n < 2:
        return []
    
    nums = []
    isPrime = []
    
    for i in range(0, n+1):#Creates list of numbers from 0 to n
        nums.append(i)
        isPrime.append(True)
        
    isPrime[0]=False
    isPrime[1]=False
    
    for j in range(2,int(n/2)):#tries all size gaps that make sense
        if isPrime[j] == True:
            for i in range(2*j,n+1,j):#starts from j+j, jumps by gap size j and crosses out that number
                isPrime[i] = False
                
    primes = []
    for i in range(0, n+1):#Adds leftovers
        if isPrime[i] == True:
            primes.append(nums[i])
            
    return primes


# def quad_residue(a, p): #quad_residue symbol of (a/p)
#     return pow(a, (p - 1) // 2, p)

def quad_residue(a,n):
    #checks if a is quad residue of n
    l=1
    q=(n-1)//2
    x = q**l
    if x==0:
        return 1
        
    a =a%n
    z=1
    while x!= 0:
        if x%2==0:
            a=(a **2) % n
            x//= 2
        else:
            x-=1
            z=(z*a) % n

    return z
        
def size_bound(N): # finds optimal factor base size and interval

    F = pow(exp(sqrt(log(N)*log(log(N)))),sqrt(2)/4)
    I = F**3
    #print(F,I)
    return int(F),int(I)


def find_base(N,B):
# generates a B-smooth factor base

    factor_base = []
    primes = prime_gen(B)
    #print(primes)
    
    for p in primes: # such that N is a quadratic residue mod p
        if quad_residue(N,p) == 1:
            factor_base.append(p)
    return factor_base

def find_smooth(factor_base,N,I):
# tries to find B-smooth numbers in sieve_seq, using sieving

    def sieve_prep(N,sieve_int):
    # generates a sequence from Y(x) = x^2 - N, starting at x = root 
        sieve_seq = [x**2 - N for x in range(root,root+sieve_int)]
        #sieve_seq_neg = [x**2 - N for x in range(root,root-sieve_int,-1)]
        return sieve_seq

    sieve_seq = sieve_prep(N,I)
    sieve_list = sieve_seq.copy() # keep a copy of sieve_seq for later
    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i,len(sieve_list),2): # found the 1st even term, now every other term will also be even
            while sieve_list[j] % 2 == 0: #account for powers of 2
                sieve_list[j] //= 2
    #print("")
    for p in factor_base[1:]: #not including 2
        residues = STonelli(N,p) #finds x such that x^2 = n (mod p). There are two start solutions
        
        for r in residues:
            for i in range((r-root) % p, len(sieve_list), p): # Now every pth term will also be divisible
                while sieve_list[i] % p == 0: #account for prime powers
                    sieve_list[i] //= p
    xlist = [] #original x terms
    smooth_nums = []
    indices = [] # index of discovery
    
    for i in range(len(sieve_list)):
        if len(smooth_nums) >= len(factor_base)+T: #probability of no solutions is 2^-T
            break
        if sieve_list[i] == 1: # found B-smooth number
            smooth_nums.append(sieve_seq[i])
            xlist.append(i+root)
            indices.append(i)

    return(smooth_nums,xlist,indices)

def build_matrix(smooth_nums,factor_base):
# generates exponent vectors mod 2 from previously obtained smooth numbers, then builds matrix

    def factor(n,factor_base):#trial division from factor base
        factors = []
        if n < 0:
            factors.append(-1)
        for p in factor_base:
            if p == -1:
                pass
            else:
                while n % p == 0:
                    factors.append(p)
                    n //= p
        return factors

    M = []
    factor_base.insert(0,-1)
    for n in smooth_nums:
        exp_vector = [0]*(len(factor_base))
        n_factors = factor(n,factor_base)
        #print(n,n_factors)
        for i in range(len(factor_base)):
            if factor_base[i] in n_factors:
                exp_vector[i] = (exp_vector[i] + n_factors.count(factor_base[i])) % 2

        #print(n_factors, exp_vector)
        if 1 not in exp_vector: #search for squares
            return True, n
        else:
            pass
        
        M.append(exp_vector)  
    #print("Matrix built:")
    #mprint(M)
    return(False, transpose(M))

    
def transpose(matrix):
#transpose matrix so columns become rows, makes list comp easier to work with
    new_matrix = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        new_matrix.append(new_row)
    return(new_matrix)

'''def optimize(M):
    for row in M: #matrix optimization; delete factors that only occur once
        if row.count(1) == 1:
            for r in M:
                del r[row.index(1)]
            del row

    return(M)'''
        
def gauss_elim(M):
#reduced form of gaussian elimination, finds rref and reads off the nullspace
#https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
    
    #M = optimize(M)
    marks = [False]*len(M[0])
    
    for i in range(len(M)): #do for all rows
        row = M[i]
        #print(row)
        
        for num in row: #search for pivot
            if num == 1:
                #print("found pivot at column " + str(row.index(num)+1))
                j = row.index(num) # column index
                marks[j] = True
                
                for k in chain(range(0,i),range(i+1,len(M))): #search for other 1s in the same column
                    if M[k][j] == 1:
                        for i in range(len(M[k])):
                            M[k][i] = (M[k][i] + row[i])%2
                break
    
    M = transpose(M)
        
    sol_rows = []
    for i in range(len(marks)): #find free columns (which have now become rows)
        if marks[i]== False:
            free_row = [M[i],i]
            sol_rows.append(free_row)
    
    if not sol_rows:
        return("No solution found. Need more smooth numbers.")
    print("Found {} potential solutions".format(len(sol_rows)))
    return sol_rows,marks,M

def solve_row(sol_rows,M,marks,K=0):
    solution_vec, indices = [],[]
    free_row = sol_rows[K][0] # may be multiple K
    for i in range(len(free_row)):
        if free_row[i] == 1: 
            indices.append(i)
    for r in range(len(M)): #rows with 1 in the same column will be dependent
        for i in indices:
            if M[r][i] == 1 and marks[r]:
                solution_vec.append(r)
                break
            
    solution_vec.append(sol_rows[K][1])       
    return(solution_vec)
    
def solve(solution_vec,smooth_nums,xlist,N):
    
    solution_nums = [smooth_nums[i] for i in solution_vec]
    x_nums = [xlist[i] for i in solution_vec]
    
    Asquare = 1
    for n in solution_nums:
        Asquare *= n
        
    b = 1
    for n in x_nums:
        b *= n

    a = helpers.isqrt(Asquare)
    
    factor = helpers.gcd(b-a,N)
    return factor


def QS(n,B,I):
#single polynomial version of quadratic sieve, given smoothness bound B and sieve interval I
    
    global N
    global root
    global T #tolerance factor
    N,root,K,T = n,int(sqrt(n)),0,1

    if is_probable_prime(N):
        return "prime"
    
    if isinstance(sqrt(N),int):
        return helpers.isqrt(N)
    
    #print(root)
    print("Attempting to factor {}...".format(N))
    #F,I = size_bound(N)
    
    print("Generating {}-smooth factor base...".format(B))
    factor_base = find_base(N,B) #generates a B-smooth factor base
    #print(factor_base)

    global F
    F = len(factor_base)
    
    print("Looking for {} {}-smooth relations...".format(F+T,B))
    smooth_nums,xlist,indices = find_smooth(factor_base, N,I)
    #finds B-smooth relations, using sieving and Tonelli-Shanks
    
    print("Found {} B-smooth numbers.".format(len(smooth_nums)))
   
    print(smooth_nums)
    
    if len(smooth_nums) < len(factor_base):
        return("Not enough smooth numbers. Increase the sieve interval or size of the factor base.")
    
    print("Building exponent matrix...")
    is_square, t_matrix = build_matrix(smooth_nums,factor_base)
    #builds exponent matrix mod 2 from relations
    
    if is_square == True:
        x = smooth_nums.index(t_matrix)
        factor = helpers.gcd(xlist[x]+sqrt(t_matrix),N)
        print("Found a square!")
        return factor, N/factor
    
    print("Performing Gaussian Elimination...")
    sol_rows,marks,M = gauss_elim(t_matrix) #solves the matrix for the null space, finds perfect square
    solution_vec = solve_row(sol_rows,M,marks,0)
    
    '''vec = [0]*len(smooth_nums) # prints solution vector
    for i in solution_vec:
        vec[i] = 1
    print("Solution vector found: " + str(vec))'''
    
    print("Solving congruence of squares...")
    #print(solution_vec)
    factor = solve(solution_vec,smooth_nums,xlist,N) #solves the congruence of squares to obtain factors

    for K in range(1,len(sol_rows)):
        if (factor == 1 or factor == N):
            print("Didn't work. Trying different solution vector...")
            solution_vec = solve_row(sol_rows,M,marks,K)
            factor = solve(solution_vec,smooth_nums,xlist,N)
        else:
            print("Found factors!")
            return factor, int(N/factor)
            
            
    return("Didn't find any nontrivial factors!")

def check_factor(n, i, factors):
    """ Checks whether 'i' is a factor of 'n' and adds 'i' to 'factors' if true
        by trial division
    """
    while n % i == 0:
        n //= i
        factors.append(i)
        if is_probable_prime(n):
            factors.append(n)
            n = 1
    return n

def find_small_primes(n, upper_bound):
    """ Perform trial division on the given number using all the primes up
        to the upper bound. Initialize the global variable 'small_primes' with
        a list of all the primes <= upper_bound.
    """
    print("Trial division and initializing small primes...")
    global small_primes
    is_prime = [True] * (upper_bound + 1)
    is_prime[0:2] = [False] * 2
    factors = []
    small_primes = []
    max_i = helpers.isqrt(upper_bound)
    rem = n
    for i in range(2, max_i + 1):
        if is_prime[i]:
            small_primes.append(i)
            rem = check_factor(rem, i, factors)
            if rem == 1:
                return factors, 1

            for j in range(i ** 2, upper_bound + 1, i):
                is_prime[j] = False

    for i in range(max_i + 1, upper_bound + 1):
        if is_prime[i]:
            small_primes.append(i)
            rem = check_factor(rem, i, factors)
            if rem == 1:
                return factors, 1

    print("Primes initialised!")
    return factors, rem

def find_prime_factors(n):
    """ Return one or more prime factors of the given number n.
        Assume that n is not a prime and does not have very small factors, and
        that the global small_primes has already been initialised. Do not
        return duplicate factors.
    """
    print("Checking whether {} is a perfect power ...".format(n))
    perfect_power = check_perfect_power(n)
    if perfect_power:
        print("{} is {}^{}".format(n, perfect_power[0], perfect_power[1]))
        factors = perfect_power[0]
    else:
        print("Not a perfect power")
        digits = len(str(n))
        print("Using Brent's variant of Pollard's rho factorization " + \
                "algorithm to factorise {} ({} digits)".format(n, digits))
        factors = [brent_factorise(n)]

    prime_factors = []
    for f in set(factors):
        for pf in find_all_prime_factors(f):
            prime_factors.append(pf)

    return prime_factors

def find_all_prime_factors(n):
    """ Return all prime factors of the given number n.
        Assume that n does not have very small factors and that the global
        small_primes has already been initialised.
    """
    rem = n
    factors = []

    while rem > 1:
        if is_probable_prime(rem):
            factors.append(rem)
            break

        for f in find_prime_factors(rem):
            print("Prime factor found: {}".format(f))
            assert is_probable_prime(f)
            assert rem % f == 0
            while rem % f == 0:
                rem //= f
                factors.append(f)

    return factors

def _pollard_brent_func(c, n, x):
    """ Return f(x) = (x^2 + c) % n
        Assume c < n
    """
    y = (x ** 2) % n + c
    if y >= n:
        y -= n

    assert y >= 0 and y < n
    return y

def brent_factorise(n, iterations=None):
    """ Perform Brent's variant of Pollard's rho factorization algorithm to
        attempt to find a non-trivial factor of the given number number, n.
        If iterations > 0, return None if no factors are found within its range
    """
    y, c, m = (randint(1, n - 1) for _ in range(3))
    r, q, g = 1, 1, 1
    i = 0
    while g == 1:
        x = y
        for _ in range(r):
            y = _pollard_brent_func(c, n, y)
        k = 0
        while k < r and g == 1:
            ys = y
            for _ in range(min(m, r - k)):
                y = _pollard_brent_func(c, n, y)
                q = (q * abs(x - y)) %  n
            g = helpers.gcd(q, n)
            k += m
        r *= 2
        if iterations:
            i += 1
            if i == iterations:
                return None

    if g == n:
        while True:
            ys = _pollard_brent_func(c, n, ys)
            g = helpers.gcd(abs(x - ys), n)
            if g > 1:
                break
    return g

def pollard_brent_iterator(n, factors):
    """ Iterator function for Brent's variant of Pollard's rho factorization
        algorithm to find all small prime factors. Restart every time a factor
        is found.
        Return 1 if all prime factors are found, or otherwise the remaining
        factor
    """
    rem = n
    while True:
        if is_probable_prime(n):
            factors.append(n)
            rem = 1
            break

        digits = len(str(n))
        if digits < 30:
            iterations = 20
        else:
            iterations = 25

            f = brent_factorise(rem, iterations)
            if f and f < rem:
                if is_probable_prime(f):
                    print("Brent's (Pollard's rho): Prime factor found: " + \
                        "{}".format(f))
                    factors.append(f)
                    rem //= f
                else:
                    print("Brent's (Pollard's rho): Composite factor " + \
                          "found: {}".format(f))
                    rem_f = pollard_brent_iterator(f, factors)
                    rem = (rem // f) * rem_f
            else:
                print("No more small factors found")
                break
    return rem

def factorise(n):
    if type(n) != int or n < 1:
        raise ValueError("Number must be a POSITIVE INTEGER")

    print("Factorizing {} ({} digits)...".format(n, len(str(n))))

    if n == 1:
        return []

    if is_probable_prime(n):
        return [n]

    factors, rem = find_small_primes(n, 10000)

    if factors:
        print("Prime factors found so far:")
        factors_temp = []
        for _ in factors:
            if _ not in factors_temp:
                factors_temp.append(_)
        print(*factors_temp, sep=', ')
    else:
        print("No small factors found!")

    if rem != 1:
        digits = len(str(rem))
        if digits > 30:
            print("Attempting Quick Pollard's rho (Brent's variation) to " + \
                  "find slightly larger factors...")
            rem = pollard_brent_iterator(rem, factors)
        if rem > 1:
            for f in find_all_prime_factors(rem):
                factors.append(f)

    factors.sort()
    assert helpers.product(factors) == n
    for p in factors:
        assert is_probable_prime(p)
    return factors
                   
def main():
    n = int(input("Enter the number to be factorized: "))

    if len(str(n)) > 50:
        print("Using Quadratic Sieve to factorise {} ({} digits)".format(n, len(str(n))))
        #B, interval = int(input("Enter B-value and interval: "))
        #print(QS(n, B, i))
        print(QS(1811706971,1000,1000))
    else:
        result = factorise(n)
        new_result = []

        for _ in result:
            if _ not in new_result:
                new_result.append(_)

        print("\nPrime factors: {}".format(new_result))

if __name__ == '__main__':
    main()
    
