# io.EIGEN 已完成

## 完成

- 确定技术路线：球贝塞尔 $j_l(x)$ 自写递推 + `径向加权Riemann` 积分，实球谐 $Y_{lm}$ 自写连带勒让德递推，3D FFT 提供 PocketFFT 和自写两种实现。
- 实现了 `transform` 模块（`src/transform/`），包括三个子模块：
  - `transform.sph_bessel`：球贝塞尔函数 $j_l(x)$（向上递推 + 小量幂级数）+ 径向积分 $\int r^2 f(r) j_l(qr) dr$
  - `transform.sph_harmonics`：实球谐函数 $Y_{lm}(\theta, \phi)$
  - `transform.fft`：1D/3D FFT（PocketFFT 包装 + 自写混合基实现，Functor 风格）
- 添加 PocketFFT 为 git submodule（`vendor/pocketfft/`）
- 更新 CMake 构建，添加 `transform` library target
- 编写 `test/test_transform.cpp`，涵盖球贝塞尔、球谐函数、FFT 三个组件的验证测试
- 所有 4 个测试通过（`test_io_ncpp`, `test_ncpp`, `test_lattice`, `test_transform`）

## 技术说明

- `std::cyl_bessel_j` 在当前 libc++ / clang 20 中不可用，球贝塞尔改用自实现
- FFT 约定（FFTW 风格）：`forward=true` = 负指数无归一化，`forward=false` = 正指数除 $n$ 归一化；`fft(true) → fft(false)` 得回原序列
- 自写 FFT 覆盖：power-of-2 radix-2、小规模 direct DFT ($n \le 64$)、大非 power-of-2 Bluestein 算法

## 下一步

- 集成 `transform` 到伪势非局域项计算（`NCPP` 算符层）
- 集成 `fft` 到波函数 k→r 空间变换
