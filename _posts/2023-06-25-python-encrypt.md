---
title: python脚本加密和文本加密
tags: coding
---



有些时候我们要在其他服务器部署流程，需要对自己的程序进行加密来保护源码。



## 脚本加密
### Cython加密

先使用[cython](https://cython.org/)将源码转为C，然后再编译的加密方法。之后写一个运行脚本通过引用加密后的库来调用里面的函数，再使用[Py-Fuscate](https://github.com/Sl-Sanda-Ru/Py-Fuscate)来混淆这个运行脚本，加大逆向的难度。

cython的使用方式是先编写一个setup脚本如setup.py
```python
from distutils.core import setup
from Cython.Build import cythonize
module = cythonize(['script.py'], build_dir='directory', language_level='3')
setup(ext_modules = module)
```

然后执行
```bash
python setup.py build_ext --inplace
```

这样，便会把script.py编译生成对应的pyd或so格式（根据系统类型）。然后再编写一个运行脚本来import编译后的pyd或so中的函数，再混淆此脚本

```bash
python py_fuscate.py -i run.py -o output_path -c 50
```

### pyarmor加密

相比以上cython加密方式，使用[pyarmor](https://pyarmor.dashingsoft.com/index-zh.html)更加简便。当然pyarmor也是一种加密和混淆的方式，还是有对应的解包工具的，如[PyArmor-Unpacker](https://github.com/Svenskithesource/PyArmor-Unpacker)，不过逆向的难度还是不小的。

```bash
# pyarmor加密foo.py和同目录下所有的.py文件
pip install pyarmor
pyarmor gen foo.py

# 加密后的脚本会保存到dist文件夹下
python dist/foo.py
```



## 文本加密

有时还需要对用到的一些文本文件进行加密，这里可采用AES-256-GCM进行加密，并把密匙写到脚本中，再对脚本进行加密。

### AES-256bit

这里把最终的输出都编码为字符串了，仅适用于我自己的流程。

```python
import base64
import hashlib
from Crypto.Cipher import AES
from Crypto.Util.Padding import pad
from Crypto import Random

# 字符串加密
def encrypt_AES(raw, password):
    private_key = hashlib.sha256(password.encode("utf-8")).digest()
    raw = pad(raw.encode("utf-8"), 16, style="pkcs7")
    iv = Random.new().read(AES.block_size)
    cipher = AES.new(private_key, AES.MODE_GCM, iv)
    encrypt_data = base64.b64encode(iv + cipher.encrypt(raw))
    return encrypt_data.decode("utf-8")

# 字符串解密
def decrypt_AES(enc, password):
    enc = enc.encode("utf-8")
    unpad = lambda x: x[: -ord(x[len(x) - 1:])]
    private_key = hashlib.sha256(password.encode("utf-8")).digest()
    enc = base64.b64decode(enc)
    iv = enc[:16]
    cipher = AES.new(private_key, AES.MODE_GCM, iv)
    unppad = unpad(cipher.decrypt(enc[16:]))
    return unppad.decode("utf-8")

# 文本文件加密
def file_encrypt_AES(inputFile, outputFile, password):
    input = open(inputFile, "r", encoding="utf-8")
    output = open(outputFile, "w", encoding="utf-8")
    
    for line in input:
        encryptString = encrypt_AES(line.strip(), password)
        output.write(encryptString + "\n")
    
    input.close()
    output.close()

# 文本文件解密
def file_decrypt_AES(inputFile, outputFile, password):
    input = open(inputFile, "r", encoding="utf-8")
    output = open(outputFile, "w", encoding="utf-8")

    for line in input:
        outputString = decrypt_AES(line.strip(), password)
        output.write(outputString + "\n")

    output.close()
    input.close()
    
# 读取加密的文本文件
def file_decrypt_read_AES(inputFile, password):
    input = open(inputFile, "r", encoding="utf-8")
    lineList = []
    for line in input:
        outputString = decrypt_AES(line.strip(), password)
        lineList.append(outputString)
    input.close()
    return lineList
```
