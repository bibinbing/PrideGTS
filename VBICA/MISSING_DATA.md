# NaN 数据直接提取

继续使用原来的调用，不增加参数，也不需要提前插值：

```matlab
VBICAResult = RunvbICA(X,NumReconComp,NumGuassModel,Source_type,SEE_PDFS,tol,max_steps,ARD);
```

`RunvbICA.m` 已恢复为修改前的原文件，参数、默认设置、分量数量、排序、
滤波方式以及五个输出字段均不变。`Components` 是直接估计的分量；
`OriginTimeSeries`、`FilteredData` 的缺失位置仍为 NaN。
没有补全输出字段，也不自动替换用户的数据。

仅在 `RunvbICA → vbica2` 训练路径保留 NaN 所必需的内部处理：

- 初始化时，仅使用通道均值为缺失位置构造有限初值。
- 训练始终使用原始缺失标记，源、混合矩阵、均值、噪声和自由能的
  观测项只累加实际有效的数据。
- 完整数据继续使用原计算分支。先验参数、收敛阈值、迭代次数及 ARD 等设置不变。

撤回了上一版新增的补全结果、vbICA1/unmix/似然函数修改，以及
`digamma` 的兼容性修复；未另行修复原项目的其他问题。

MATLAB 验证：`addpath(genpath('VBICA')); test_missing_data`。
本次 NaN 回归测试覆盖 `g` 源模型、随机/连续缺失、整个时刻缺失、
两种噪声设置与 ARD 开关。其他源模型保留原选项，未在此次缩小修改范围后
声明其所有组合已通过验证。

每个通道至少需要两个有效观测来初始化；整个时刻均缺失时，
该时刻没有数据约束，不能据此声称恢复了真实信号。
