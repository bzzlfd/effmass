我们已经正确实现了 
1. $H\ket{\psi_{nk}}$： @src/H_psi/hamiltonian_callable.cpp 
2. $\frac{\partial H}{\partial k_\alpha} \ket{\psi_{nk}}$ ：@src/H_psi/hamiltonian_gradient.cpp

借由此，现在已经可以计算 $\frac{\partial \varepsilon_{nk}}{\partial k_\alpha}$: @test/H_psi/test_hgrad_eigen_fd.cpp 

接下来，我们来完成计算 $\frac{\partial \psi_{nk}}{\partial k_\alpha}$ 的代码。其中 $\psi_{nk}$ 是非简并态，我们在这里把它当作 WG 文件中的平面波系数向量。





