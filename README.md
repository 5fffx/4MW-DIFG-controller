# 4MW 双馈风电机组主控制器 C++ 源码

本项目提供 4MW 双馈风电机组背靠背 NPC 变流器主控制器的 C++ 实现与综合技术手册，涵盖网侧与机侧控制、参数与命令处理、逐步调度，以及 PI、PLL、坐标变换、滤波、调制、故障穿越和保护算法。

## 文件结构

```text
技术手册.pdf            A4 综合技术手册
src/                    7 个 C++ 实现文件和 1 个接口头文件
build_primitives.cmd    Windows x64 / MSVC 编译脚本
.gitignore              忽略编译产物
```

## 阅读指南

[技术手册（A4 打印版）](技术手册.pdf) 共 130 页、13 章，包含 187 条编号公式、25 组算法与源码对照、109 项面板参数、操作流程及历史验证记录。PDF 自带目录、书签和页码，无需安装 LaTeX。

| 想了解的内容 | 手册章节 |
| --- | --- |
| 双馈系统、功率流向与物理量基准 | 第 2 章 |
| 状态区、调度顺序与输入输出接口 | 第 4 章 |
| 参数换算、命令及完整参数表 | 第 5、13 章 |
| 公共算法的数学公式与源码对照 | 第 6 章 |
| 网侧 GSC 与机侧 RSC 控制 | 第 7、8 章 |
| 频率支撑、穿越恢复与保护 | 第 9、10 章 |
| 编译、接入与故障排查 | 第 3、12 章 |
| 历史工况、验证方法与结果适用范围 | 第 11 章 |

参数说明统一以手册第 5、13 章为准。手册保留原资料的覆盖索引与校订记录，便于追溯公式、参数和控制描述。

## 源码导航

| 文件 | 职责 |
| --- | --- |
| [recovered_primitives.h](src/recovered_primitives.h) | 状态结构与导出接口声明 |
| [recovered_primitives.cpp](src/recovered_primitives.cpp) | PI、PLL、变换、滤波、调制及输入输出算法 |
| [recovered_controller.cpp](src/recovered_controller.cpp) | 初始化入口、逐步调度、结束入口及快照 |
| [recovered_initialization.cpp](src/recovered_initialization.cpp) | 上电默认值与状态初始化 |
| [recovered_parameters.cpp](src/recovered_parameters.cpp) | 参数换算与保护定值装入 |
| [recovered_commands.cpp](src/recovered_commands.cpp) | 运行命令及功率给定处理 |
| [recovered_grid.cpp](src/recovered_grid.cpp) | 网侧控制、穿越、保护及频率支撑 |
| [recovered_motor.cpp](src/recovered_motor.cpp) | 机侧控制、穿越、谐波及撬棒相关处理 |

## 编译

使用 Windows x64 和 Visual Studio 的 MSVC C++ 工具链。在 **x64 Native Tools Command Prompt for VS** 中切换到本项目根目录，运行：

```bat
build_primitives.cmd
```

也可以将环境变量 `VCVARS64` 设为本机 `vcvars64.bat` 的完整路径，由脚本加载编译环境。脚本使用 `/O2 /fp:strict /EHsc /W4 /utf-8` 编译设置，输出 `build/recovered_primitives.dll`；中间文件也位于 `build/`，已加入忽略规则。

## 使用范围

本目录包含控制器源码、编译脚本和技术手册。接入 PSCAD 需要另行准备相应工程、Fortran 宿主包装及运行依赖；完成 DLL 编译后仍需核对接口、数组长度、状态布局与调用约定。

手册中的仿真结果和差分结果来自历史记录。修改控制器后，应在配套环境中重新验证。手册引用的历史目录、原始文件名和验证工具用于说明原工程，当前公开文件结构以本 README 为准。
