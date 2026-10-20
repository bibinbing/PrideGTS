# PrideGTS

[English](README.md) | **简体中文**

![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-orange)
![Platform](https://img.shields.io/badge/platform-Windows%2010%2F11-blue)
![License](https://img.shields.io/badge/license-GPL--3.0-green)

PrideGTS（Pride GNSS Time Series Analysis Software）是基于MATLAB的开源GNSS坐标时间序列分析软件，在同一个图形界面中完成从多源坐标文件输入到参数估计的处理。

## 主要功能

- 读取PBO（`.pos`）、NGL（`.tenv3`）和CMONOC（`.txt`）坐标文件，转换为Hector兼容的`.mom`文件
- 预处理：异常值剔除、已知阶跃导入、三分量阶跃探测（MGLR）和缺失数据重构（KKF、PCA iteration、RegEM、VBPCA）
- 环境负荷改正（NTAL、NTOL、HYDL）和共模误差滤波（PCA、VBPCA、VBRPCA、FastICA、VBICA）
- 基于Hector的噪声模型选择和参数估计
- 所有模块均支持批量处理

## 运行环境

- MATLAB R2022b或更高版本；仅需MATLAB基础模块，无需其他工具箱
- Windows 10/11（`Hector/exe`中附带的Hector 2.1可执行程序为Windows程序）

## 安装

1. 下载或克隆本仓库到有写权限的本地目录。
2. 在MATLAB中将PrideGTS文件夹及其全部子文件夹添加到路径（**Home → Set Path → Add with Subfolders**）。
3. 双击`Install/PrideGTS.mlappinstall`。
4. 在**APPS**选项卡中启动PrideGTS。

## 使用说明

见[用户手册](docs/PrideGTS_Software_User_Manual_EN.pdf)（英文），其中包括用`sample_data/`样例数据完成的逐步示例。

## 第三方组件

PrideGTS包含以下第三方代码，各组件保留其原有的版权声明。

| 组件 | 参考文献 | 许可证 |
|---|---|---|
| Hector 2.1 | Bos et al. (2013) | GPL-3.0-or-later |
| FastICA 2.5 | Hyvärinen (1999) | GPL-2.0-or-later |
| VBRPCA | Babacan et al. (2012) | GPL-3.0-or-later |
| VBPCA | Ilin and Raiko (2010) | 见原始发布 |
| vbICA 1.0和Netlab | Choudrey (2003); Nabney (1996–2001) | 见原始发布；BSD-3-Clause |
| RegEM | Schneider (2001) | 见原始发布 |
| GMIS（KKF） | Liu et al. (2018) | 见原始发布 |

## 许可证

[GNU通用公共许可证第3版（GPL-3.0）](LICENSE)

## 联系方式

Xianpeng Huang，桂林理工大学（2640355140@qq.com）
