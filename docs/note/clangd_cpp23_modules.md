# clangd 不需要任何额外配置就能自动找到 `std.pcm`

项目使用 C++23 模块（`import std`）和 CMake 3.31 的实验性 `CMAKE_CXX_MODULE_STD` 支持。起初认为需要在 `.clangd` 中硬编码 `std.pcm` 的路径，例如：

```yaml
# 不需要这样做
CompileFlags:
  Add:
    - -fprebuilt-module-path=build/CMakeFiles/__cmake_cxx23.dir
```

`build/CMakeFiles/__cmake_cxx23.dir` 是 CMake 内部临时目录，依赖它不是好的实践。实测发现 clangd **根本不需要这些配置**。

## 验证

```bash
mv .clangd .clangd.bak
clangd --check=src/io.cppm --compile-commands-dir=build
# All checks completed，无 import std 相关错误
```

## 为什么能自动找到

### 1. compile_commands.json 引用响应文件

CMake 生成的 `compile_commands.json` 中，编译命令使用了**响应文件**：

```json
{
    "directory": "/home/kai/projects/effmass/build",
    "command": "/usr/bin/clang++ -stdlib=libc++ -std=gnu++23 @CMakeFiles/io.dir/src/io.cppm.o.modmap -o ... -c src/io.cppm"
}
```

`@CMakeFiles/io.dir/src/io.cppm.o.modmap` 是编译器响应文件。clang/gcc/msvc 都支持 `@文件` 语法——编译器把文件内容展开为命令行参数。

### 2. .modmap 文件里写了 std.pcm 的位置

```
-x c++-module
-fmodule-output=CMakeFiles/io.dir/io.pcm
-fmodule-file=std=CMakeFiles/__cmake_cxx23.dir/std.pcm
```

`-fmodule-file=std=...` 显式绑定了 `std` 模块的 BMI（.pcm）路径。clangd 读取 `compile_commands.json` 时自动展开响应文件，获得这个参数，从而正确解析 `import std;`。

### 3. 路径是相对 directory 字段解析的

响应文件中的 `CMakeFiles/__cmake_cxx23.dir/std.pcm` 是相对于 `"directory": "/home/kai/projects/effmass/build"` 解析的，clangd 能正确处理。

## 结论

`.clangd` 只需保留最基本的 fallback：

```yaml
CompileFlags:
  Add: [-std=c++23, -stdlib=libc++]
```

甚至这两行也可以省略。绝对不要写 `-fprebuilt-module-path=build/CMakeFiles/...`。

## 前提

- `CMAKE_EXPORT_COMPILE_COMMANDS ON` 已启用
- **已执行过 `cmake --build build`**（生成了 `.modmap` 和 `.pcm`）
- clangd >= 21（实测 21.1.8）

### 必须先编译一遍

`.modmap` 响应文件和 `.pcm` 模块文件都是在**构建过程中**生成的，configure 阶段不会创建它们。如果只在 `build/` 里执行了 `cmake -B build -G Ninja` 而没有 `cmake --build build`，clangd 展开响应文件后会指向一个不存在的 `.pcm` 路径，`import std` 仍然会报错。

所以完整流程是：

```bash
cmake -B build -G Ninja
cmake --build build   # 必须有这一步，生成 .modmap 和 .pcm
```

然后 clangd 才能正确解析。

## VS Code 中的配置

在项目根目录创建 `.vscode/settings.json`：

```json
{
    "clangd.path": "/usr/bin/clangd",
    "clangd.arguments": [
        "--compile-commands-dir=build"
    ],
    "C_Cpp.intelliSenseEngine": "disabled"
}
```

`--compile-commands-dir=build` 告诉 clangd 去 `build/` 目录下找 `compile_commands.json`。clangd 默认也会搜索 `build/` 等常见子目录，显式指定更可靠。
