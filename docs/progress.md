我们再进行一些修订

1. matrix.cppm 改名 array2d.cppm， matrix 改名 array2d
    1. 查找用到这个的代码，同步
    2. 增加 operator[][]，讨论是column major 还是 row major
2. linalg.cppm 中，eigenvectors 用 array2d 返回， array2d[i] 是一个eigenvector