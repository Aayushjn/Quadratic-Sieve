import math
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

	return z==1

'''def calculateB(n):
	factor_base=[]
	x=math.ceil(math.sqrt(n))
	return x
'''

def convertx2e(x,e):
	e = 0
	while x % 2 == 0:
		x /= 2
		e+=1
	return x

def order(p, b):
	if (gcd(p, b) != 1):
		printf("p and b are not co-prime.\n")
		return -1

	# Initializing k with first odd prime number
	k = 3;
	while (1):
		if (pow(b, k, p) == 1):
			return k
		k+=1
def gcd(x,y):
    if (y == 0):
        return x
    else:
        return gcd(y, x % y)

def sqrt_mod_prime(n, p):
	"""Return the square root of n modulo the prime p. Behaviour is
	undefined if n is not a quadratic residue mod p."""
	if n==0:
		return 0
	if p==2:
		return n
	if p%2==0:
		return None
	#Shanks algorithm
	#a and p should be coprime for finding the modular square root
	if gcd(n,p) != 1:
		return 0
 
	#If below expression return (p - 1)  then modular square root is not possible
	if pow(n,(p-1)//2,p)==(p - 1):
		return -1;

	# expressing p - 1, in terms of s * 2^e,  where s is an odd number
	e=0
	s = convertx2e(p - 1, e);
 
	#finding smallest q such that q ^ ((p - 1) / 2)
	#(mod p) = p - 1
	 #q - 1 is in place of  (-1 % p)
	q=2
	while(True):
		if pow(q, (p - 1) // 2, p) == (p - 1):
			break
		q+=1
	x = pow(n,(s + 1) // 2, p);
	b = pow(n, s, p);
	g = pow(q, s, p);
 
	r = e;
 
	# keep looping until b become 1 or m becomes 0
	while (True):
		for m in range(0,r):
			if (order(p, b) == -1):
				return -1;
 
				#finding m such that b^ (2^m) = 1
			if (order(p, b) == pow(2, m)):
				break;
		if m == 0:
			return x;
 
		#updating value of x, g and b according to algorithm
		x = (x * pow(g, pow(2, r - m - 1), p)) % p;
		g = pow(g, pow(2, r - m), p);
		b = (b * g) % p;
 
		if b == 1:
			return x;
		r = m
   

def siqs_factor_base_primes(n, nf):
	"""Compute and return nf factor base primes suitable for a Quadratic
	Sieve on the number n.
	"""
	small_primes=[] #List of small primes upto B.
	factor_base = []
	for p in small_primes:
		if quad_residue(n, p):
			t = sqrt_mod_prime(n%p,p)
			lp = round(log2(p))
			factor_base.append(FactorBasePrime(p, t, lp))
			if len(factor_base)>= nf:
				break
	return factor_base

def siqs_sieve(factor_base, m):
	"""Perform the sieving step of the SIQS. Return the sieve array."""
	sieve_array = [0]*(2*m+1)
	for fb in factor_base:
		if fb.soln1 is None:
			continue
		p = fb.p
		i_start_1 = -((m + fb.soln1) // p)
		a_start_1 = fb.soln1 + i_start_1 * p
		lp = fb.lp
		if p > 20:
			for a in range(a_start_1 + m, 2 * m + 1, p):
				sieve_array[a] += lp

			i_start_2 = -((m + fb.soln2) // p)
			a_start_2 = fb.soln2 + i_start_2 * p
			for a in range(a_start_2 + m, 2 * m + 1, p):
				sieve_array[a] += lp
	return sieve_array
def main():
	n=87463909928837
	#x=sieve(n)
	#print(x)
	#find quadratic residues from 1 to sqrt(n) 
	#B=sqrt(n)
	#for i in range(1,x):		
	#if quad_residue(i,x)==True:
		#print(i)
	print(quad_residue(2,113))
	print(sqrt_mod_prime(2,113))
if __name__=='__main__':
	main()