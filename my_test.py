import my_ml_kem
import subprocess, sys
from my_ml_kem import k, Rq, du, q
import numpy as np

tv_d = bytes.fromhex('8085c873120ca7c2325f22418210332b9de6664656a6fa736864134165b489b6')
tv_z = bytes.fromhex('c91460e5aebd987585c9e7175155ecaa02008c8d9ae00aedfa11205ecc8b42c5')
tv_m = bytes.fromhex('530e0f36e99c8035e47b17656dc8e6c24e2510ab587b93c6a9a9547816026612')

def xtimes(poly: Rq, p: int):
    for i in range(p):
        coeff = poly.coeff.pop(0)
        coeff = coeff * -1
        poly.coeff.append(coeff if coeff > 0 else coeff + q)

def get_pk_masked_ct(A_, c: bytearray, sk: bytearray, row: int, scalar: int, rot: int):
    """
    Mult-then-comp
    """
    ### Decode
    pk = sk[384*k:768*k+32]
    A_star = [Rq.intt(A_[row][i]) * scalar for i in range(k)]   # Mult
    xtimes(A_star[0], rot)  # Rotate
    xtimes(A_star[1], rot)  # Rotate

    t_ = [Rq.decode(pk[12*32*i:]) for i in range(k)]
    t_star = Rq.intt(t_[row]) * scalar  # Mult
    xtimes(t_star, rot) # Rotate

    ### pk mask
    u = np.array(Rq.polyvecDecodeDecomp(c))
    v = Rq.polyDecodeDecomp(c[32*du*k:])
    # v.coeff[0] = v.coeff[0] + 832
    # u[0].coeff[0] = u[0].coeff[0] + 100

    # A_star = Rq.polyvecCompEncode(A_star) 
    # t_star = t_star.polyCompEncode() 
    # c = A_star + t_star
    # A_star = Rq.polyvecDecodeDecomp(c)
    # A_star[0] *= scalar
    # A_star[1] *= scalar
    # xtimes(A_star[0], rot)
    # xtimes(A_star[1], rot)
    # t_star = Rq.polyDecodeDecomp(c[32*du*k:]) * scalar
    # xtimes(t_star, rot)
    u = A_star
    v = t_star

    ### Encode
    c1 = Rq.polyvecCompEncode(u)    ### 640 bytes
    c2 = v.polyCompEncode()         ### 128 bytes

    return c1+c2

def main():
    inst = my_ml_kem.my_ML_KEM()
    pk, sk = inst.cca_keygen(tv_z, tv_d)
    ct, K = inst.cca_enc(pk, tv_m)
    # inst.cca_dec(ct, sk)
    A_ = inst.genA(pk[-32:])
    pk_masked_ct = get_pk_masked_ct(A_, ct, sk, 0, 1, 0)

    cp = subprocess.run(['../bin/my_test_ml-kem-512_clean'], input = pk_masked_ct + pk,\
                        stdout=subprocess.PIPE)
    print(cp)
    # out = cp.stdout.decode('utf-8').split('/n')
    # print(out[0])

if __name__ == "__main__":
    main()

