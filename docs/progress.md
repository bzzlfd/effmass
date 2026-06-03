# H_psi

H_psi 作为一个Functor（即实现了operator()的类）
构造是读取 GKK，ATOM，pseudo.ncpp文件们，Functor 对 WG 的数据进行变换。

Functor 对 WG 变换过程中，累积到自己内部的与 WG 同样尺寸的vector `H_psi`上，包含以下几项
- [ ] `H_psi` 的命名

```
H_psi[:] = 0
```

## 动能
验证 GKK中 kinetic 是否可靠，然后直接使用这个kinetic数据。进行数据累加
```
H_psi[g] += kinetic[g] * WG[g]
```

## 局域势
通过傅里叶变换快速实现：
```
WG --[FFT:G2R]--> WR --> WR * VR --[FFT:R2G]--> WVG
H_psi[g] += WVG[g]
```


## Pseudo Potential 非局域项

### 结构因子
- [ ] 命名
Pseudo Potential 和原子位置有关。对于平面波 exp(iqr), 和原子位置 tau，
相对于原子位于原点，pseudopotential对平面波的作用会额外多一个 exp(iq tau)的系数

具体实现上，
1. 我们 GKK IntegerView 时把其中的数据向量看成 -K=-(k+G)
2. exp函数内部的加法可以拆分成外面的乘法
所以我们不采用 n1*n2*n3 大小的数组记录，而是使用三个数组，大小分别是 n1 n2 n3 记录 exp(-i b1 tau)^iG exp(-i b2 tau)^jG  exp(-i b3 tau)^kG，然后再乘 exp(-i k tau)

### 非局域势计算

- 在 SphericalHarmonics 中加入 reset(l_max), 更改常驻 Y_lm(GKK) 的数量

- simpson积分 meshtype 从 ncpp 到 bessel_fourier 转换
```cpp
  auto mt = (ncpp.mesh.type == MeshType::Uniform)
            ? RadialMeshType::Uniform : RadialMeshType::General;
  simpson(f, rab, mt);
```


```cpp
for(ityp: ATOM.eachtype())
    for(iatm: ATOM.eachatom(ityp))
        V_psp = ncpp(ityp)
        Y_lm = SphericalHarmonics.reset(GKK, l_max)
        structure_factor = structure_factor(GKK)
        for(l=0; l<=l_max; ++l)
            V_NL_l = V_psp.projectorBlock(l)
            for(ib=0; ib<=nb; ++ib)
                beta_q_interpolaror = BetaQInterpolator(V_psp.mesh, V_NL_l.beta[ib])  // handle q_max 
                beta_q = (GKK => beta_q_interpolaror)
                for(m=-l; m<=l; ++m)
                    projector = (beta_q .* Y_lm .* structure_factor)
                    H_psi[:] +=
                        innerProduct(conj(projector), WG) * 
                        V_NL_l.B[ib, ib] * 
                        projector[:]
```




# 未来
我们会增加 $H^{\alpha} \psi$ 和 $H^{\alpha\beta} \psi$，alpha和beta指的是对于  以 k 为参数的哈密顿量 H(k)，对于 k_alpha 和 k_beta 分量的一阶偏导数和二阶偏导数
- [ ] 命名

# module 组织

#### 前置声明
`module_A` 和 `module_B`
	**A 类需要内置一个 B 类的对象**：为了知道 B 类对象的大小（以便为 A 分配内存空间），编译器在编译 A 时，**必须完整地看到 B 类的具体定义**。
	**B 类需要让 A 是自己的 friend**：在 C++ 模块中，要声明另一个模块的类为 `friend`，B 必须知道 A 的存在。
结果导致循环依赖

解决：
	利用“前置声明（Forward Declaration）”减弱某一侧的依赖，并把它们揉进同一个模块（可以分在不同的分区文件中）。
	就是说：让它们属于一个模块，而不是两个模块
```cpp
// MyModule-Types.cppm
export module MyModule:Types; // 属于 MyModule 的 Types 分区

// 1. 前置声明 A
export class A; 

// 2. 完整定义 B（因为 B 只需要知道 A 的名字就能设为友元）
export class B {
public:
    void do_something_with_A(A& a);
private:
    int secret_data = 42;
    
    friend class A; // 声明 A 是友元，正确！此时不需要 A 的完整定义
};
```

```cpp
// MyModule-A.cppm
export module MyModule:A; // 属于 MyModule 的 A 分区
import :Types;            // 导入 Types 分区，这样就能看到 B 的完整定义

export class A {
public:
    void ClassAMethod() {
        // 因为 A 是 B 的友元，且看到了 B 的定义，可以直接访问 B 的私有成员
        b_object.secret_data = 100; 
    }
private:
    B b_object; // 必须内置 B 的对象，因为有 import :Types，编译器知道 B 的大小，正确！
};
```

```cpp
// MyModule.cppm
export module MyModule; // 主模块

export import :Types;
export import :A;

// 顺便在这里实现一下 B 的成员函数（因为这里同时能看到 A 和 B 的完整定义了）
module :private; // 或是放在非导出的 .cpp 文件中
void B::do_something_with_A(A& a) {
    // 具体的实现...
}
```


> [!note]
> 类的定义不能拆开，因为要知道它的内存分布

#### 接口与实现分离
内部接口分区 (Internal Interface Partitions)

在 `.cppm` 文件中只保留声明（骨架），在非分区的 `.cpp` 模块实现文件中写实现。

1. 接口文件 (`MyModule.cppm`)
```cpp
export module MyModule;

export class A {
public:
    // 声明嵌套类
    class nest {
    public:
        void do_something(); // 只留声明，不写具体实现
        int get_value() const;
    private:
        int data_ = 0;
    };
    
    void foo();
};
```
2. 实现文件 (`MyModule.cpp`)**
	1. 没有 `export`，它是模块内部隐私
```cpp
module MyModule; // 声明我属于 MyModule 模块

// 使用 作用域解析运算符 (::) 来实现嵌套类的成员函数
void A::nest::do_something() {
    // 这里可以写几百行长代码...
}

int A::nest::get_value() const {
    return data_;
}

void A::foo() {
    // A 的其他成员函数实现
}
```

#### **Pimpl 模式** 或者 **前置声明+指针**

Pimpl 模式（Pointer to Implementation）

1. 接口文件 (`MyModule.cppm`)
```cpp
export module MyModule;
import std; // 用到 std::unique_ptr

export class WindowManager {
public:
    WindowManager();  // 构造函数
    ~WindowManager(); // 析构函数

    void init_system();
    void render_all();

private:
    // 核心：前置声明一个实现类，但不在此处定义它
    struct Impl; 
    
    // 类里唯一的私有成员就是一个指针。
    // 无论 Impl 里面有多庞大，指针在 64 位系统上永远只占 8 个字节。
    // 编译器不需要知道 Impl 的大小，就能确定 WindowManager 的大小！
    std::unique_ptr<Impl> pimpl; 
};
```
2. 内部实现分区或普通的实现文件 (`MyModule.cpp`)
```cpp
module MyModule;

// 揭开 Impl 的真面目！这个结构体只在这个 .cpp 文件里可见
struct WindowManager::Impl {
    // 把之前所有复杂的私有变量和复杂的嵌套类全部塞到这里
    struct WindowNode {
        int width = 0;
        int height = 0;
        void update_bounds(int w, int h) {}
    };

    WindowNode root_node;
    int system_status = 0;
};

// 构造函数：真正去创建这个隐藏的对象
WindowManager::WindowManager() 
    : pimpl(std::make_unique<Impl>()) {}

// 析构函数：必须要写在这里，因为这里编译器才知道 Impl 怎么销毁
WindowManager::~WindowManager() = default;

// 外部调用接口时，WindowManager 只是个“传话筒”，把工作转交给 pimpl
void WindowManager::init_system() {
    pimpl->system_status = 1; // 实际操作隐藏的数据
}

void WindowManager::render_all() {
    pimpl->root_node.update_bounds(1024, 768);
}
```

> [!note]
> 不需要把外围类当作 friend 哦

> [!note] cpp 这种实现写法，无法对这个函数进行inline展开
> **时下主流编译器都支持 LTO（Link-Time Optimization，链接期优化）**（在 MSVC 中被称为 LTCG，Link-Time Code Generation）。
> 
> - **CMake**: 只需要加一行 `set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)`
> - **GCC / Clang 命令行**: 编译和链接时都加上 `-flto`
> - **MSVC 命令行**: 加上 `/GL`（编译）和 `/LTCG`（链接）

