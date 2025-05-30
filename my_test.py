import my_ml_kem
import subprocess, sys

tv_d = bytes.fromhex('8085c873120ca7c2325f22418210332b9de6664656a6fa736864134165b489b6')
tv_z = bytes.fromhex('c91460e5aebd987585c9e7175155ecaa02008c8d9ae00aedfa11205ecc8b42c5')
tv_m = bytes.fromhex('530e0f36e99c8035e47b17656dc8e6c24e2510ab587b93c6a9a9547816026612')

def main():
    inst = my_ml_kem.my_ML_KEM()
    pk, sk = inst.cca_keygen(tv_z, tv_d)
    ct, K = inst.cca_enc(pk, tv_m)
    inst.cca_dec(ct, sk)

    cp = subprocess.run(['../bin/my_test_ml-kem-512_clean'], input = ct + pk,\
                        stdout=subprocess.PIPE)
    # print(cp)
    out = cp.stdout.decode('utf-8').split('/n')
    print(out[0])

if __name__ == "__main__":
    main()

