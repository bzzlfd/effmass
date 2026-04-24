---
name: cpp23-std-module-manual
description: >
  不使用 CMake 的情况下，手动编译 libc++ 的 C++23 std 模块（std.pcm），
  用于命令行编译包含 `import std;` 的代码。
---

# 手动编译 C++23 std 模块（无 CMake）

## 1. 定位系统上的 std.cppm

```bash
clang++ -print-resource-dir
```
会输出类似 `/usr/lib/llvm-20/lib/clang/20` 的路径，`std.cppm` 就在这个路径下的 `share/libc++/v1/std.cppm`。 

## 2. 预编译 std.pcm

```bash
clang++ -std=c++23 -stdlib=libc++ -Wno-reserved-module-identifier \
    -c /usr/lib/llvm-20/share/libc++/v1/std.cppm \
    -o /path/to/std.pcm
```

参数说明：
- `-std=c++23`：`import std` 的最低要求
- `-stdlib=libc++`：目前只有 libc++ 提供 `std.cppm`
- `-Wno-reserved-module-identifier`：抑制 `std` 是保留模块名的警告
- `-c`：只编译不链接

## 3. 编译项目代码

### 简单程序（仅 import std）

```bash
clang++ -std=c++23 -stdlib=libc++ \
    -fmodule-file=std=/path/to/std.pcm \
    main.cpp -o main
```

### 含自定义 module

```bash
# 编译自定义模块接口
clang++ -std=c++23 -stdlib=libc++ \
    -fmodule-file=std=/path/to/std.pcm \
    -c src/io.cppm -Xclang -emit-module-interface \
    -o io.pcm

# 编译主程序
clang++ -std=c++23 -stdlib=libc++ \
    -fmodule-file=std=/path/to/std.pcm \
    -fmodule-file=io=io.pcm \
    -c main.cpp -o main.o

# 链接
clang++ -stdlib=libc++ main.o -o main
```
