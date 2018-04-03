import math
def quad_residue(a,n):
    l=1
    q=(n-1)//2
    x = q ** l
    if x == 0:
        return 1
        
    a =a% n
    z=1
    while x != 0:
        if x % 2 == 0:
            a = (a ** 2) % n
            x //= 2
        else:
            x -= 1
            z = (z * a) % n

    return z==1

def sieve(n):
    factor_base=[]
    x=math.ceil(math.sqrt(n))
    return x

def main():
    n=87463
    x=sieve(n)
    for i in range(1,x):
        if quad_residue(i,x)==True:
            print(i)

if __name__=='__main__':
    main()