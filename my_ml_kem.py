from __future__ import annotations
import hashlib
import itertools
import numpy as np
import copy
#import sha3

### 512
k = 2
eta1 = 3

### 768
# k = 3
# eta1 = 2

### 1024
# k = 4
# eta1 = 2

### General
n = 256
q = 3329
eta2 = 2
du = 10
dv = 4
tree = [
  0, 64, 32, 96, 16, 80, 48, 112, 8, 72, 40, 104, 24, 88, 56, 120,
  4, 68, 36, 100, 20, 84, 52, 116, 12, 76, 44, 108, 28, 92, 60, 124,
  2, 66, 34, 98, 18, 82, 50, 114, 10, 74, 42, 106, 26, 90, 58, 122,
  6, 70, 38, 102, 22, 86, 54, 118, 14, 78, 46, 110, 30, 94, 62, 126,
  1, 65, 33, 97, 17, 81, 49, 113, 9, 73, 41, 105, 25, 89, 57, 121,
  5, 69, 37, 101, 21, 85, 53, 117, 13, 77, 45, 109, 29, 93, 61, 125,
  3, 67, 35, 99, 19, 83, 51, 115, 11, 75, 43, 107, 27, 91, 59, 123,
  7, 71, 39, 103, 23, 87, 55, 119, 15, 79, 47, 111, 31, 95, 63, 127
]
zetas = [17**i % q for i in range(n)]
n_inv = pow(n, q-2, q)

class Rq:
    """
    Name:        Rq
    12 bits * 256 = 384 bytes
    """

    def __init__(self):
        self.coeff = [0 for x in range(n)]

    def __repr__(self):
        return str(list(map(hex, self.coeff)))

    def __getitem__(self, index):
        return self.coeff[index]

    def __eq__(self, other):
        return self.coeff == other.coeff

    def __add__(self, other):
        tmp = self.__class__()
        for i in range(n):
            tmp.coeff[i] = (self.coeff[i] + other.coeff[i]) % q ### reduction
            #上の%qで値が正になるので、ref実装に合わせるため[-q/2,q/2]の範囲に戻す。
            if tmp.coeff[i] > 1664:
                tmp.coeff[i] -= q 
        return tmp

    def __sub__(self, other):
        tmp = self.__class__()
        for i in range(n):
            tmp.coeff[i] = (self.coeff[i] - other.coeff[i]) % q ### reduction
            #上の%qで値が正になるので、ref実装に合わせるため[-q/2,q/2]の範囲に戻す。
            if tmp.coeff[i] > 1664:
                tmp.coeff[i] -= q 
        return tmp

    def __matmul__(self, other: Rq) -> Rq:
        ### Polynomial mult
        tmp = self.__class__()
        for i in range(0, n, 2):
            tmp.coeff[i] = self.coeff[i+1] * other.coeff[i+1]
            tmp.coeff[i] = tmp.coeff[i] * 17**(2*tree[i//2] + 1)
            tmp.coeff[i] += self.coeff[i] * other.coeff[i]
            tmp.coeff[i] = tmp.coeff[i] % q
            tmp.coeff[i+1] = self.coeff[i] * other.coeff[i+1]
            tmp.coeff[i+1] += self.coeff[i+1] * other.coeff[i]

        return tmp

    def __mul__(self, scalar: int) -> Rq:
        ### Scalar mult
        tmp = self.__class__()
        for i in range(n):
            tmp.coeff[i] = self.coeff[i] * scalar % q
        return tmp
    
    @classmethod
    def ntt(cls, poly_in):
        # Straight forward version
        poly_out = cls()
        for i in range(128):
            for j in range(128):
                zeta = zetas[(2*tree[i]+1)*j % n]
                poly_out.coeff[2*i] += poly_in.coeff[2*j] * zeta
                poly_out.coeff[2*i+1] += poly_in.coeff[2*j+1] * zeta
        
        # # FFT version
        # poly_out = copy.deepcopy(poly_in)
        # kk = 1
        # for len in [128, 64, 32, 16, 8, 4, 2]:
        #     for start in range(0, 256, 2*len):
        #         zeta = 17**tree[kk]
        #         kk += 1
        #         for j in range(start, start + len):
        #             t = zeta * poly_out.coeff[j + len]
        #             poly_out.coeff[j + len] = poly_out.coeff[j] - t
        #             poly_out.coeff[j] = poly_out.coeff[j] + t

        # Reduction
        for i in range(n):
            poly_out.coeff[i] %= q
            #上の%qで値が正になるので、ref実装に合わせるため[-q/2,q/2]の範囲に戻す。
            if poly_out.coeff[i] > 1664:
                poly_out.coeff[i] -= q 

            #bit_reverse = int(f'{i:08b}'[::-1], 2)

        return poly_out

    @classmethod
    def intt(cls, poly_in):
        # FFT version
        poly_out = copy.deepcopy(poly_in)
        kk = 127
        for len in [2, 4, 8, 16, 32, 64, 128]:
            for start in range(0, 256, 2*len):
                zeta = 17**tree[kk]
                kk -= 1
                for j in range(start, start + len):
                    t = poly_out.coeff[j]
                    poly_out.coeff[j] = t + poly_out.coeff[j+len]
                    poly_out.coeff[j+len] = zeta * (poly_out.coeff[j+len] - t)

        # Reduction
        for i in range(n):
            poly_out.coeff[i] *= 3303 ### Multiply every entry by 128^-1
            poly_out.coeff[i] %= q
            #上の%qで値が正になるので、ref実装に合わせるため[-q/2,q/2]の範囲に戻す。
            if poly_out.coeff[i] > 1664:
                poly_out.coeff[i] -= q 

            #bit_reverse = int(f'{i:08b}'[::-1], 2)

        return poly_out

    @classmethod
    def sample_uniform(cls, rho, j, i) -> Rq:
        poly = cls()
        s = hashlib.shake_128()
        s.update(rho + bytes([j]) + bytes([i]))
        xof_out = bytes.fromhex(s.hexdigest(672))

        sample_cnt = 0
        for i in range(0, 672, 3):
            val = xof_out[i:i+3]
            d1 = val[0] + 256 * (val[1] & 0xf)
            d2 = (val[1] >> 4) + 16 * val[2]
            if(d1 < q and sample_cnt < n):
                poly.coeff[sample_cnt] = d1
                sample_cnt += 1
            if(d2 < q and sample_cnt < n):
                poly.coeff[sample_cnt] = d2
                sample_cnt += 1
        
        return poly

    @classmethod
    def sample_cbd(cls, sigma, nonce, eta):

        ### PRF part
        poly = cls()
        s = hashlib.shake_256()
        s.update(sigma + bytes([nonce]))
        prf_out = s.digest(192)

        ### CBD part
        if eta == 3:
            for i in range(256//4):
                t = int.from_bytes(prf_out[i*3:i*3+3], 'little')    ### By 24bits
                d = t & 0x00249249;
                d += (t>>1) & 0x00249249
                d += (t>>2) & 0x00249249

                for j in range(4):
                    a = (d >> (6*j)) & 0x07
                    b = (d >> (6*j + 3)) & 0x07
                    poly.coeff[4*i + j] = a - b
                    #負数を利用

        elif eta == 2:
            for i in range(256//8):
                t = int.from_bytes(prf_out[i*4:i*4+4], 'little')    ### By 32bits
                d = t & 0x55555555;
                d += (t>>1) & 0x55555555

                for j in range(8):
                    a = (d >> (4*j)) & 0x03
                    b = (d >> (4*j + 2)) & 0x03
                    poly.coeff[8*i + j] = a - b
                    #負数を利用
        else:
            exit(-1)

        return poly

    def encode(self) -> bytearray:
        val = 0
        for i in range(0, n, 2):
            elem0 = self.coeff[i]
            elem1 = self.coeff[i+1]

            if elem0 < 0:
                elem0 += q
            if elem1 < 0:
                elem1 += q

            val <<= 8
            val += elem0 & 0xff
            val <<= 8
            val += ((elem1 & 0xf) << 4) | (elem0 >> 8)
            val <<= 8
            val += elem1 >> 4

        return val.to_bytes(384, 'big')

    @classmethod
    def decode(cls, barray) -> Rq:
        ### d = 12
        poly = cls()
        for i in range(256//2):
            b0 = barray[3*i]
            b1 = barray[3*i+1]
            b2 = barray[3*i+2]

            uh = b1 >> 4
            lh = b1 & 0x0f

            poly.coeff[2*i] = (lh << 8) | b0
            poly.coeff[2*i+1] = (b2 << 4) | uh
        
        return poly
        
    @classmethod
    def msgdecode(cls, msg: bytearray):
        ### d = 1
        poly = cls()
        for i in range(32):
            t = msg[i]
            for j in range(8):
                bit = (t >> j) & 0x01
                if bit:
                    poly.coeff[8*i+j] = 1665
                else:
                    poly.coeff[8*i+j] = 0
        
        return poly
    
    @classmethod
    def msgencode(cls, poly: Rq) -> int:
        ### d = 1
        msg = 0
        for i in range(32):
            t = 0
            for j in range(8):
                coeff = poly.coeff[8*i+j]
                coeff <<= 1
                coeff += 1665
                coeff *= 80635
                coeff >>= 28
                coeff &= 1
                t |= coeff << j
            msg <<= 8
            msg |= t
        return msg

    @classmethod
    def polyvecCompEncode(cls, polyvec):
        ### d = 10
        barray = bytearray()
        ### compress and encode
        for i in range(k):
            for j in range(256//4):
                list_t = [0] * 4
                for kk in range(4):
                    list_t[kk] = polyvec[i][j*4+kk]
                    if list_t[kk] < 0:
                        list_t[kk] += q
                    list_t[kk] = ((((list_t[kk] << 10) + 1665)*1290167) >> 32) & 0x3ff    ###  & 0x3ffいらない？ -> いる
                barray.append(list_t[0] & 0xff)
                barray.append(list_t[0] >> 8 | (list_t[1] << 2) & 0xff)
                barray.append(list_t[1] >> 6 | (list_t[2] << 4) & 0xff)
                barray.append(list_t[2] >> 4 | (list_t[3] << 6) & 0xff)
                barray.append(list_t[3] >> 2)

        return barray

    def polyCompEncode(self):
        ### d = 4
        barray = bytearray()
        ### compress and encode
        for i in range(256//8):
            list_t = [0] * 8
            for j in range(8):
                list_t[j] = self.coeff[i*8+j]
                if list_t[j] < 0:
                    list_t[j] += q
                list_t[j] = ((((list_t[j] << 4) + 1665)*80635) >> 28) & 0xf
            
            barray.append(list_t[0] | list_t[1] << 4)
            barray.append(list_t[2] | list_t[3] << 4)
            barray.append(list_t[4] | list_t[5] << 4)
            barray.append(list_t[6] | list_t[7] << 4)

        return barray

    @classmethod
    def polyvecDecodeDecomp(cls, barray: bytearray) -> polyvec:
        ### d = 10
        polyvec = [Rq() for x in range(k)]
        list_t = [0] * 4

        ### decode and decompress
        for i in range(k):
            for j in range(256//4):
                list_int= [x for x in barray[i*320+j*5:i*320+j*5+5]]
                list_t[0] = (list_int[0] >> 0) | (list_int[1] << 8)
                list_t[1] = (list_int[1] >> 2) | (list_int[2] << 6)
                list_t[2] = (list_int[2] >> 4) | (list_int[3] << 4)
                list_t[3] = (list_int[3] >> 6) | (list_int[4] << 2)

                for kk in range(4):
                    polyvec[i].coeff[4*j+kk] = ((list_t[kk] & 0x3ff) * q + 512) >> 10

        return polyvec

    @classmethod
    def polyDecodeDecomp(cls, barray: bytearray) -> Rq:
        ### d = 4
        poly = Rq()
        list_t = [0] * 4

        ### decode and decompress
        for i in range(256//2):
            integer = barray[i]
            poly.coeff[2*i + 0] = ((integer & 0xf) * q + 8) >> 4
            poly.coeff[2*i + 1] = ((integer >>  4) * q + 8) >> 4

        return poly


def hash_H(din: bytearray) -> str:
    s = hashlib.sha3_256()
    s.update(din)
    return s.digest()

def hash_G(din: bytearray) -> str:
    s = hashlib.sha3_512()
    s.update(din)
    return s.digest()

def hash_J(din: bytearray) -> str:
    s = hashlib.shake_256()
    s.update(din)
    K = s.digest(32)

    return K

class my_ML_KEM:
    """
    Name:        my_ML_KEM
    """

    def __init__(self):
        pass
        # if security == 768:
        #     k = 3
        #     eta1 = 2
        #     eta2 = 2
        # elif security == 512:
        #     self.k = 2
        #     self.eta1 = 3
        #     self.eta2 = 2
        # else:
        #     print("Invalid security level")
        #     exit(-1)

    def genA(self, rho, transpose = False) -> list[Rq][Rq]:
        A = []
        for i in range(k):
            list_temp = []
            for j in range(k):
                if transpose:
                    list_temp.append(Rq.sample_uniform(rho, i, j))
                else:
                    list_temp.append(Rq.sample_uniform(rho, j, i))
            A.append(list_temp)
        return A

    def unpack(self, arr) -> int:
        val = 0
        for poly in np.array(arr).flatten():
            for elem in poly:
                if elem < 0:
                    #ref実装のテストベクタに合わせるため16ビット符号付数に変換
                    elem += q
                val <<= 16
                # Little endian
                val += elem >> 8
                val += (elem & 0xff) << 8
        return val

    def cpa_keygen(self, d: bytearray) -> bytearray:
        Gout = hash_G((d + bytes([k])))
        rho = Gout[:32]
        sigma = Gout[32:]
        A = self.genA(rho)
        N = 0

        s = []
        for i in range(k):
            s.append(Rq.sample_cbd(sigma, N, eta1))
            N += 1

        e = []
        for i in range(k):
            e.append(Rq.sample_cbd(sigma, N, eta1))
            N += 1

        ntt_s = [Rq.ntt(s[i]) for i in range(k)]

        ntt_e = [Rq.ntt(e[i]) for i in range(k)]

        # Matrix-vector multiplication
        As = [Rq() for x in range(k)]
        for i in range(k):
            for j in range(k):
                As[i] = As[i] + (A[i][j] @ ntt_s[j])
        
        ntt_t = [As[i] + ntt_e[i] for i in range(k)]

        # Encode
        pk, sk = bytes(0), bytes(0)
        for i in range(k):
            pk += ntt_t[i].encode()
            sk += ntt_s[i].encode()
        pk += rho

        ### pk and sk are int type
        return pk, sk

    def cca_keygen(self, z: bytearray, d: bytearray) -> bytearray:
        pk, sk_ = self.cpa_keygen(d)
        H = hash_H(pk)
        sk = sk_ + pk + H + z

        ### pk and sk are byte array
        return pk, sk
    
    def enc(self, pk: bytearray, m: bytearray, coin: bytearray) -> bytearray:
        t_ = [Rq.decode(pk[12*32*i:]) for i in range(k)]
        rho = pk[-32:]
        N = 0

        At = self.genA(rho, True)

        y = []
        for i in range(k):
            y.append(Rq.sample_cbd(coin, N, eta1))
            N += 1

        e1 = []
        for i in range(k):
            e1.append(Rq.sample_cbd(coin, N, eta2))
            N += 1

        e2 = Rq.sample_cbd(coin, N, eta2)

        ntt_y = [Rq.ntt(y[i]) for i in range(k)]

        # Matrix-vector multiplication
        Aty = [Rq() for x in range(k)]
        for i in range(k):
            for j in range(k):
                Aty[i] = Aty[i] + (At[i][j] @ ntt_y[j])

        u = [Rq.intt(Aty[i]) + e1[i] for i in range(k)]
        # print(hex(self.unpack(e1)))

        mu = Rq.msgdecode(m)

        # Vector-vector multiplication
        tty = Rq()
        for i in range(k):
                tty = tty + (t_[i] @ ntt_y[i])
        v = Rq.intt(tty) + e2 + mu

        c1 = Rq.polyvecCompEncode(u)    ### 640 bytes

        c2 = v.polyCompEncode()         ### 128 bytes

        return c1 + c2
    
    def cca_enc(self, pk: bytearray, m: bytearray) -> (bytearray, str):
        Kr = hash_G(m + hash_H(pk))
        K = Kr[:32]
        r = Kr[32:]
        c = self.enc(pk, m, r)

        return c, K

    def dec(self, dk: bytearray, c: bytearray) -> bytearray:
        u = Rq.polyvecDecodeDecomp(c)
        v = Rq.polyDecodeDecomp(c[32*du*k:])

        ntt_s = [Rq.decode(dk[384*i:]) for i in range(k)] 
        ntt_u = [Rq.ntt(u[i]) for i in range(k)]

        # Vector-vector multiplication
        stu = Rq()
        for i in range(k):
                stu = stu + (ntt_s[i] @ ntt_u[i])
        intt_stu = Rq.intt(stu)

        w = v - intt_stu
        m = Rq.msgencode(w)
        print(m)
        print(hex(m))
        print(w)
        input()

        return m.to_bytes(32, 'big')

    def cca_dec(self, c: bytearray, sk: bytearray):
        dk = sk[:384*k]
        pk = sk[384*k:768*k+32]
        h = sk[768*k+32:768*k+64]
        z = sk[768*k+64:768*k+96]
        mp = self.dec(dk, c)

        Kprp = hash_G(mp + h)
        Kp = Kprp[:32]
        rp = Kprp[32:]
        K_bar = hash_J(z + c)
        cp = self.enc(pk, mp, rp)
        
        if c != cp:
            Kp = K_bar

        return Kp

def KAT_ref(N):
    for idx in range(N):
        print(f'#{idx} test vector...')
        start = idx * 11
        with open(r'../kat/ml_kem_768_ref.kat', 'r') as f:
            tv = f.readlines()[start:start+10]

        tv_d = bytes.fromhex(tv[0].split(':')[1])
        tv_z = bytes.fromhex(tv[1].split(':')[1])
        tv_pk = bytes.fromhex(tv[2].split(':')[1])
        tv_sk = bytes.fromhex(tv[3].split(':')[1])
        tv_m = bytes.fromhex(tv[4].split(':')[1])
        tv_ct = bytes.fromhex(tv[5].split(':')[1])
        tv_ss = bytes.fromhex(tv[6].split(':')[1])
        tv_ss = bytes.fromhex(tv[7].split(':')[1])
        tv_ct_p = bytes.fromhex(tv[8].split(':')[1])
        tv_z_ss = bytes.fromhex(tv[9].split(':')[1])

        inst = my_ML_KEM()
        pk, sk = inst.cca_keygen(tv_z, tv_d)
        assert pk == tv_pk, f'{pk:x} != {tv_pk:x}'
        assert sk == tv_sk, f'{sk:x} != {tv_sk:x}'

        c, K = inst.cca_enc(pk, tv_m)
        assert c == tv_ct, f'{c:x} != {tv_ct:x}'
        assert K == tv_ss, f'{K:x} != {tv_ss:x}'
        print(f'Bob (enc) side shared secret K: {K}')

        K = inst.cca_dec(c, sk)
        assert K == tv_ss, f'{K:x} != {tv_ss:x}'
        print(f'Alice (dec) side shared secret K: {K}')

        K = inst.cca_dec(tv_ct_p, sk)
        assert K == tv_z_ss, f'{K:x} != {tv_z_ss:x}'
        print(f'Pseudorandom shared secret Kp: {K}')

def KAT_itzmeanjan(N):
    for idx in range(N):
        print(f'#{idx} test vector...')
        start = idx * 8
        with open(r'../kat/ml_kem_768.kat', 'r') as f:
            tv = f.readlines()[start:start+7]

        tv_d = bytes.fromhex(tv[0].split('=')[1])
        tv_z = bytes.fromhex(tv[1].split('=')[1])
        tv_pk = bytes.fromhex(tv[2].split('=')[1])
        tv_sk = bytes.fromhex(tv[3].split('=')[1])
        tv_m = bytes.fromhex(tv[4].split('=')[1])
        tv_ct = bytes.fromhex(tv[5].split('=')[1])
        tv_ss = bytes.fromhex(tv[6].split('=')[1])

        inst = my_ML_KEM()
        pk, sk = inst.cca_keygen(tv_z, tv_d)
        assert pk == tv_pk, f'{pk} != {tv_pk}'
        assert sk == tv_sk, f'{sk:x} != {tv_sk:x}'

        c, K = inst.cca_enc(pk, tv_m)
        assert c == tv_ct, f'{c:x} != {tv_ct:x}'
        assert K == tv_ss, f'{K:x} != {tv_ss:x}'
        print(f'Bob (enc) side shared secret K: {K}')

        K = inst.cca_dec(c, sk)
        assert K == tv_ss, f'{K:x} != {tv_ss:x}'
        print(f'Alice (dec) side shared secret K: {K}')


if __name__ == "__main__":
    N = 100
    # KAT_ref(N)
    KAT_itzmeanjan(N)