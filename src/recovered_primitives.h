// ============================================================================
// 文件：recovered_primitives.h —— 明码控制器的“接口说明书”
// ----------------------------------------------------------------------------
// 这个文件不包含具体算法，只声明两样东西：
//   1. 数据结构（struct）：控制器内部各模块使用的“变量打包盒”；
//   2. 函数清单（extern "C"）：这个 DLL 对外（对 PSCAD）提供哪些函数。
//
// 【背景：PSCAD 是怎么调用这个 DLL 的】
// PSCAD 的仿真内核是 Fortran 程序。它在每个仿真步都会调用一次本 DLL 里的函数。
// 调用方式不是传单个数值，而是传“一大块内存的首地址”（C++ 里叫指针）。
// 这块内存就是控制器的“状态区”（state）：控制器需要记住的所有中间量
// （积分器数值、滤波器历史值、锁相环角度、开关状态……）都存在这里。
// 每个变量在状态区里有一个固定的位置，用“字节偏移量 offset”表示。
// 比如偏移量 0xC0 处存的是网侧 A 相电压，0x550 处存的是 PLL 的比例系数。
// 每次调用函数时，函数就从这块内存里读出需要的量，算完再写回去。
//
// 【关于 0x…… 这种数字】
// 以 0x 开头的数字是十六进制（计算机常用的计数方式，0x10 = 十进制 16）。
// 在本工程里它们几乎都是“状态区内的字节偏移量”，可以理解成变量的“门牌号”：
// 0x270 = 状态区第 624 个字节处。偏移量之所以长这样，是因为它们是直接从
// 原 DLL 的二进制里逆向还原出来的，我们保留了原样，方便逐字节核对。
//
// 【extern "C" __declspec(dllexport) 是什么意思】
// - __declspec(dllexport)：把这个函数“导出”，让 DLL 外部的程序能调用它；
// - extern "C"：告诉编译器函数名不要被 C++ 改写（C++ 默认会把函数名和参数
//   类型混在一起编码），保证 PSCAD 侧能用原始函数名找到它。
//
// 【怎么阅读下面的结构体】
// struct 是“把若干变量打包成一个整体”的语法。比如 RegulatorState 里打包了
// reference、feedback、error……这就是一个 PI 调节器需要用到的全部变量。
// 结构体里变量的排列顺序 = 这些变量在状态区内存里的排列顺序，
// 这一点必须和原 DLL 完全一致（有 static_assert 在编译期检查，见 .cpp 文件）。
// ============================================================================
#pragma once
#include <cstdint>

// PI（比例-积分）调节器状态，另带抗饱和（antiwindup）和微分项。
// 控制里最基本的单元：给定值 reference 和反馈值 feedback 相减得到误差 error，
// 误差乘比例系数 kp 得比例项；误差不断累加、乘积分系数得积分项；两者相加输出。
// 积分和输出都有上下限（upper_limit / lower_limit，即“饱和”）。
// antiwindup_error 是“抗饱和量”：当输出被限幅时，把被削掉的部分反馈回积分器，
// 防止积分器越积越深（积分饱和）。ki_per_step 是折算到每步的积分增益。
struct RegulatorState {
    double reference;
    double feedback;
    double error;
    double kp;
    double proportional;
    double integral;
    double derivative;
    double unsaturated;
    double upper_limit;
    double lower_limit;
    double output;
    double antiwindup_error;
    double ki_per_step;
    double tracking_gain;
    double derivative_gain;
    double previous_error;
};

// 线电压 → αβ 变换（Clarke 变换）的状态。
// vab、vbc 是两路线电压采样值（第三路可由前两路算出），
// alpha、beta 是变换到两相静止坐标系的结果。
// 注意 reserved 是原 DLL 状态区里保留的字节，只占位、无意义。
struct LineVoltageState {
    double vab;
    double vbc;
    double reserved;
    double alpha;
    double beta;
};

// 三相电流 → αβ 变换（等幅值 Clarke 变换，系数 2/3）。
// ia/ib/ic 三相电流，alpha/beta 是变换结果。
struct CurrentTransformState {
    double ia;
    double ib;
    double ic;
    double alpha;
    double beta;
};

// Park 变换：把静止的 αβ 坐标系旋转到与电网电压同步的 dq 旋转坐标系。
// theta 是旋转角（通常来自锁相环 PLL），d、q 是变换结果。
// 在 dq 坐标系里，50Hz 的正弦量会变成直流量，PI 调节器才能无静差地控制它。
struct ParkState {
    double alpha;
    double beta;
    double theta;
    double d;
    double q;
};

// 两态滤波器（two-state filter）：用 x、z 两个内部状态实现的二阶滤波结构。
// step 是计算步长，scale/gain 决定滤波特性，x_coefficient/z_coefficient
// 是每步更新的系数。它常用于提取电压基波（对指定频率几乎无衰减）。
struct TwoStateFilter {
    double step;
    double scale;
    double input;
    double error;
    double x;
    double next_x;
    double z;
    double next_z;
    double x_coefficient;
    double z_coefficient;
    double gain;
};

// 一阶递归滤波器（IIR）：y[n] = a0*x[n] + a1*x[n-1] + b1*y[n-1]。
// 当前输出 = 当前输入×系数 + 上一步输入×系数 + 上一步输出×系数。
struct RecursiveFilterState {
    double input_coefficient;
    double previous_input_coefficient;
    double previous_output_coefficient;
    double input;
    double previous_input;
    double output;
    double previous_output;
};

// 残差滤波器：先对输入做低通滤波，再用“输入 - 低通结果”得到高频残差。
// 常用于提取谐波/纹波分量，然后反相注入以抑制谐波。
// 两通道（channel 0/1）分别对应 d 轴和 q 轴。
struct ResidualFilterState {
    double input[2];
    double lowpass[2];
    double residual[2];
};

// 锁相环 PLL（Phase Locked Loop）的完整状态。
// 功能：从三相电压里“锁定”电网的相位角 angle 和角速度 angular_speed。
// 结构：线电压 → Clarke 变换（stationary）→ 正负序解耦网络（positive_* / negative_*）
//       → 低通滤波 → q 轴分量送 PI 调节器（regulator）→ 角速度积分得角度。
// voltage_magnitude 是正序电压幅值，negative_magnitude 是负序电压幅值，
// feedforward 是给后续电流环使用的电压前馈量。
struct PllState {
    double line_voltage[3];
    double angular_speed;
    double angle;
    double positive_decoupled[2];
    double negative_decoupled[2];
    double positive_filtered[2];
    double negative_filtered[2];
    double positive_raw[2];
    double negative_raw[2];
    double voltage_magnitude;
    double negative_magnitude;
    double feedforward[2];
    double step;
    LineVoltageState stationary;
    ParkState rotating;
    double reserved[16];
    RegulatorState regulator;
};

// 调制模块（把控制环算出的电压参考值变成开关管的占空比）。
// alpha/beta：调制波（αβ 形式）；phase[3]：三相调制波；
// phase_max/phase_min/middle：三相的最大/最小/中间值，用于注入零序分量
// （相当于 SVPWM/三次谐波注入，提高直流电压利用率约 15%）；
// dc_plus/dc_minus：直流母线正、负电压（用于中点电位平衡）；
// correction：中点电位平衡的修正量；duty[3][2]：三相上、下管的占空比。
struct ModulationState {
    double alpha;
    double beta;
    double phase[3];
    double phase_max;
    double phase_min;
    double middle;
    double reserved_040[9];
    double modulation[3];
    double reserved_0A0[5];
    double dc_plus;
    double dc_minus;
    double delta;
    double reserved_0E0;
    double correction;
    double duty[3][2];
    double gain;
    double limit;
    double deadband;
    double reserved_138[14];
    std::int32_t polarity;
    std::uint32_t reserved_tail;
};

// ============================================================================
// 下面是 DLL 对外导出的全部函数清单（共 34 个）。
// 命名规律：
//   - 纯计算单元（regulator、pll、modulation……）只操作自己的结构体，
//     各自在状态区里占一小段；
//   - *_input / *_output 是 PSCAD 与状态区之间的“搬运工”（输入映射/输出映射）；
//   - *_defaults / *_initialize 是初始化；
//   - grid_* 是网侧变流器控制，motor_* / rotor_* 是机侧（转子侧）变流器控制。
// 每个函数的具体算法和公式，见同名的 .cpp 文件里的详细注释。
// ============================================================================
extern "C" __declspec(dllexport) void regulator(RegulatorState* state);
extern "C" __declspec(dllexport) void line_voltage_transform(LineVoltageState* state);
extern "C" __declspec(dllexport) void current_transform(CurrentTransformState* state);
extern "C" __declspec(dllexport) void park_transform(ParkState* state);
extern "C" __declspec(dllexport) void inverse_park_transform(ParkState* state);
extern "C" __declspec(dllexport) void two_state_filter(TwoStateFilter* state);
extern "C" __declspec(dllexport) void recursive_filter(RecursiveFilterState* state);
extern "C" __declspec(dllexport) void residual_filter(ResidualFilterState* state);
extern "C" __declspec(dllexport) void pll(double base_frequency, PllState* state);
extern "C" __declspec(dllexport) void modulation(ModulationState* state);
extern "C" __declspec(dllexport) void alternate_modulation(void* state);
extern "C" __declspec(dllexport) void motor_calculate(double base_frequency, void* state);
extern "C" __declspec(dllexport) void grid_input(const double* input, void* state);
extern "C" __declspec(dllexport) void rotor_input(const double* input, double base_frequency, void* state);
extern "C" __declspec(dllexport) void grid_output(double* output, const void* state);
extern "C" __declspec(dllexport) void rotor_output(double* output, const void* state);
extern "C" __declspec(dllexport) void motor_defaults(double base_frequency, void* state);
extern "C" __declspec(dllexport) void motor_post(double base_frequency, void* state);
extern "C" __declspec(dllexport) void grid_defaults(double base_frequency, void* state);
extern "C" __declspec(dllexport) void controller_initialize(double base_frequency, void* state);
extern "C" __declspec(dllexport) void average_converter_output(void* state);
extern "C" __declspec(dllexport) void view_output(double* output, const void* state);
extern "C" __declspec(dllexport) void rotor_nameplate(const double* parameters, void* state);
extern "C" __declspec(dllexport) void grid_post(double base_frequency, void* state);
extern "C" __declspec(dllexport) void command_input(const double* commands, void* state);
extern "C" __declspec(dllexport) void protection_input(const double* input, double base_frequency, void* state);
// Input is public90 + 10 (80 readable doubles); state spans 0x3748 bytes.
extern "C" __declspec(dllexport) void parameter_input(const double* input, double base_frequency, void* state);
// State spans 0x1AB8 bytes; ring lengths are 1..5, with indices in [0, length).
extern "C" __declspec(dllexport) void grid_frequency_support(double base_frequency, void* state);
extern "C" __declspec(dllexport) void grid_calculate(double base_frequency, void* state);
extern "C" __declspec(dllexport) void adaptive_oscillator(void* state, double base_frequency);
extern "C" __declspec(dllexport) void current_sequence(void* state);
extern "C" __declspec(dllexport) void rotor_current_sequence(void* state);
extern "C" __declspec(dllexport) void resonant_filter(void* state);
extern "C" __declspec(dllexport) void difference_filter(void* state);
