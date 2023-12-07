from zokrates_pycrypto.babyjubjub import Point as point
from zokrates_pycrypto.field import FQ

q = 2736030358979909402780800718157159386076813972158567259200215660948447373041
alpha = 21356

# n = 128
# k = 64

# n = 64
# k = 32

# n = 32
# k = 16

n = 16
k = 8


def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m
    
def mod(x):
    return modinv(x, q)

def div(x,y):
    return (x-y) % q

def mul(x,y):
    return (x*y) % q

def addmod(x,y):
    return (x+y)%q

#TSS reconstruct function
def reconstruct(shares):
    G = point.generator()
    k = len(shares)
    total_yj_delta = point.infinity()
    for j in range(k):
        xj, yj = shares[j]
        delta = 1
        for m in range(k):
            if(m != j):
                xm, _ = shares[m]
                delta = mul(delta, (mul(xm,mod(div(xm,xj)))))
        total_yj_delta = total_yj_delta.add(yj.mult(delta))
    return total_yj_delta

#generate [h]_t-1
def genPP():
    p = [point.infinity() for x in range(k)]
    for i in range(k):
        p[i] = point.generator().mult((alpha+i)%q)
    return p

#ComX^TSS
def genZ(p, x):
    buf = point.infinity()
    for i in range(k):
        buf = buf.add(p[i].mult(x**(i+1)%q))
    return buf

#GenY^TSS
def genY(ghx, r, zeta):
    return ghx.mult(r%q).add(zeta)

if __name__ == "__main__":
    G = point.generator()
    r = 21888242871839275222246405745257275088548364400416034343698204
    # r = 11111994
    secret = G.mult(r+1234)
    print("secret: ", secret)

    p = genPP()
    ghx = [point.infinity() for x in range(n)]
    shares = [point.infinity() for x in range(n)]
    gfx = [point.infinity() for x in range(n)]

    for i in range(n):
        ghx[i] = genZ(p, 37821381+i)
        gfx[i] = genY(ghx[i], r, secret)
        shares[i] = [37821381+i, gfx[i]]

    print(reconstruct(shares))
    print("abc")
