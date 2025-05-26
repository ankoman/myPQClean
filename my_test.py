import my_ml_kem
import subprocess, sys

tv_d = bytes.fromhex('8085c873120ca7c2325f22418210332b9de6664656a6fa736864134165b489b6')
tv_z = bytes.fromhex('c91460e5aebd987585c9e7175155ecaa02008c8d9ae00aedfa11205ecc8b42c5')
tv_m = bytes.fromhex('530e0f36e99c8035e47b17656dc8e6c24e2510ab587b93c6a9a9547816026612')

def main():
    inst = my_ml_kem.my_ML_KEM()
    cp = subprocess.run(['../bin/my_test_ml-kem-512_clean'], input = tv_d + tv_z,\
                        stdout=subprocess.PIPE)
    out = cp.stdout.split()
    print(out)
    
    pk = out[0]
    sk = out[1]
    ct = out[2]
    key_b = out[3]
    key_a = out[4]



    pk, sk = inst.cca_keygen(tv_z, tv_d)
    print(pk.hex())


if __name__ == "__main__":
    main()

