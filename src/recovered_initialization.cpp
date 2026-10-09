// ============================================================================
// 文件：recovered_initialization.cpp —— 控制器初始化（上电默认值）
// ----------------------------------------------------------------------------
// PSCAD 启动时（或复位时）把整个状态区填成默认值。
// 两个入口：
//   grid_defaults          ：网侧状态区（前 0x1AB8 字节）默认值；
//   controller_initialize  ：整块状态区清零后依次调用两侧默认值。
// ============================================================================
#include "recovered_primitives.h"

#include <cstddef>
#include <cstring>
#include <emmintrin.h>

namespace {

// ----------------------------------------------------------------------------
// 常数表：原 DLL 数据段里的浮点常数，按它们原来的内存地址命名
// （C61930 表示原地址 0x61930 处的常数）。
// 保留地址命名是为了可追溯：每个数都能对回原 DLL 的二进制。
// 下面挑几个重要的说明：
//   0.00025   = 控制周期 250µs（4kHz 控制频率）；
//   1.414     = √2；0.577 = 1/√3；2.3567 = 3π/4；
//   62.83185  = 2π×10；314.159 = 2π×50；471.2389 = 2π×75；
//   628.3185  = 2π×100；1884.9556 = 2π×300；6283.185 = 2π×1000；
//   930.6928  = 网侧线电压基准（幅值）；1611.96 ≈ √2×1140（机端相电压峰值）；
//   636.3     = 直流母线额定电压；989.8 = 机侧电流基准。
// ----------------------------------------------------------------------------
constexpr double C61930 = 1.2e-05;
constexpr double C61950 = 3.5e-05;
constexpr double C61970 = 0.00025;
constexpr double C61990 = 0.0012;
constexpr double C619D0 = 0.002827433388230811;
constexpr double C61A80 = 0.05;
constexpr double C61A90 = 0.06;
constexpr double C61AC0 = 0.088;
constexpr double C61AF0 = 0.12;
constexpr double C61B00 = 0.15;
constexpr double C61B60 = 0.2;
constexpr double C61CD0 = 0.5;
constexpr double C61D20 = 0.577;
constexpr double C61E00 = 0.8;
constexpr double C61E20 = 0.854;
constexpr double C61E70 = 0.9;
constexpr double C61E90 = 0.942;
constexpr double C61EF0 = 1.0;
constexpr double C61F00 = 1.02;
constexpr double C61F40 = 1.1;
constexpr double C61F50 = 1.123;
constexpr double C61FC0 = 1.2;
constexpr double C62000 = 1.25;
constexpr double C62060 = 1.414;
constexpr double C620E0 = 1.8;
constexpr double C62100 = 2.0;
constexpr double C62140 = 2.3566666666666665;
constexpr double C62200 = 4.0;
constexpr double C62210 = 5.0;
constexpr double C62250 = 7.069999999999999;
constexpr double C622A0 = 10.0;
constexpr double C622E0 = 14.139999999999999;
constexpr double C62320 = 30.0;
constexpr double C623A0 = 62.8318530717958;
constexpr double C623D0 = 80.0;
constexpr double C62480 = 114.0;
constexpr double C624E0 = 212.1;
constexpr double C62520 = 282.8;
constexpr double C62530 = 314.159265358979;
constexpr double C62540 = 318.15;
constexpr double C62560 = 471.2388980384685;
constexpr double C62590 = 565.6;
constexpr double C625A0 = 628.318530717958;
constexpr double C625B0 = 636.3;
constexpr double C625C0 = 658.1986143187067;
constexpr double C62600 = 930.6928406466512;
constexpr double C62640 = 1000.0;
constexpr double C62680 = 1030.0;
constexpr double C626A0 = 1140.0;
constexpr double C626B0 = 1170.0;
constexpr double C62710 = 1400.0;
constexpr double C62720 = 1596.0;
constexpr double C62730 = 1611.9599999999998;
constexpr double C62740 = 1692.558;
constexpr double C62760 = 1800.0;
constexpr double C62770 = 1884.955592153874;
constexpr double C62780 = 1900.0;
constexpr double C62790 = 1950.0;
constexpr double C627A0 = 2100.0;
constexpr double C627D0 = 2400.0;
constexpr double C62820 = 6283.18530717958;
constexpr double C62830 = 6363.0;
constexpr double C62840 = 8000.0;
constexpr double C62850 = 12019.0;
constexpr double C62860 = 12566.37061435916;
constexpr double C62880 = 141400.0;
constexpr double C628A0 = 4000000.0;
constexpr double C628B0 = 8388608.0;
constexpr double C628C0 = -0.05;
constexpr double C628D0 = -0.07853981633974474;
constexpr double C629F0 = -62.8318530717958;
constexpr double C62A20 = -282.8;
constexpr double C62A30 = -424.2;
constexpr double C62A40 = -471.2388980384685;

inline double load_double(const unsigned char* state, std::size_t offset)
{
    double value;
    std::memcpy(&value, state + offset, sizeof(value));
    return value;
}

inline void store_double(unsigned char* state, std::size_t offset, double value)
{
    std::memcpy(state + offset, &value, sizeof(value));
}

inline void store_int(unsigned char* state, std::size_t offset, std::int32_t value)
{
    std::memcpy(state + offset, &value, sizeof(value));
}

// period_count：把频率 f 换算成“每周期的采样步数”（1/f ÷ 250µs，向下取整）。
inline std::int32_t period_count(double f)
{
    return _mm_cvttsd_si32(_mm_set_sd((1.0 / f) / 0.00025));
}

// The R, H and F store sequences in the recovered specification.
// ----------------------------------------------------------------------------
// 三个“初始化序列”，被下面的默认值函数反复调用：
//   initialize_r：PI/谐振环节初值 —— 把 p 处值复制到 b，q 处值复制到 b+0x48，
//     算出比值，再放两个对称限幅（±0.577×v）；
//   initialize_h：谐振器初值 —— 1.414、314.159（2π×50）、0.00025、±62.83；
//   initialize_f：滤波器初值 —— 4.0、628.3185（2π×100）、0.00025。
// ----------------------------------------------------------------------------
inline void initialize_r(unsigned char* state, std::size_t b,
                         std::size_t p, std::size_t q, std::size_t v)
{
    store_double(state, b, load_double(state, p));
    store_double(state, b + 0x48, load_double(state, q));
    store_double(state, b + 0x50,
                 load_double(state, b + 0x48) / load_double(state, b));
    store_double(state, b + 0x28, C61D20 * load_double(state, v));
    store_double(state, b + 0x30, -load_double(state, b + 0x28));
}

inline void initialize_h(unsigned char* state, std::size_t b)
{
    store_double(state, b, C62060);
    store_double(state, b + 0x18, C62530);
    store_double(state, b + 0x08, C61970);
    store_double(state, b + 0x28, C623A0);
    store_double(state, b + 0x30, C629F0);
}

inline void initialize_f(unsigned char* state, std::size_t b)
{
    store_double(state, b + 0x08, C62200);
    store_double(state, b, C625A0);
    store_double(state, b + 0x50, C61970);
}

} // namespace

// Original RVA 0x363A0. Unspecified bytes must remain untouched.
// ============================================================================
// ----------------------------------------------------------------------------
// grid_defaults —— 网侧状态区默认值。
// 把网侧控制要用到的所有变量写一遍初值。按块说明如下。
// ----------------------------------------------------------------------------
// ============================================================================
extern "C" __declspec(dllexport) void grid_defaults(double f, void* state)
{
    auto* s = static_cast<unsigned char*>(state);

// 开关标志与计数器：0x000 启动计数、0x004 启动时长 400 步（0.1 秒）、
// 0x008 = 800 步（0.2 秒）、0x35C 停机时长 800 步；
// 0x350 运行状态、0x50C 输出使能、0x4B8/0x4D0/0x4E8 输出清零。
    store_int(s, 0x000, 0);
    store_int(s, 0x018, 0);
    store_int(s, 0x01C, 0);
    store_int(s, 0x004, 400);
    store_int(s, 0x008, 800);
    store_int(s, 0x35C, 800);
    store_int(s, 0x350, 0);
    store_double(s, 0x4B8, 0.0);
    store_double(s, 0x4D0, 0.0);
    store_double(s, 0x4E8, 0.0);
    store_int(s, 0x50C, 0);
    store_int(s, 0x500, 0);

// 直流电压基准与限幅：0x110 = 636.3V（直流母线额定电压），
// 0x120/0x138/0x140/0x148 = 636.3 的 0.5/1.2/1.2/1.2 倍（外环输出限幅），
// 0x118 = 5×636.3（积分限幅），0x080 = 8388608（大数，用作除法基准）。
    store_double(s, 0x110, C625B0);
    store_double(s, 0x120, C61CD0 * load_double(s, 0x110));
    store_double(s, 0x138, C62000 * load_double(s, 0x110));
    store_double(s, 0x140, C62000 * load_double(s, 0x110));
    store_double(s, 0x148, C62000 * load_double(s, 0x110));
    store_double(s, 0x118, C62210 * load_double(s, 0x110));
    store_int(s, 0x158, 1);
    store_double(s, 0x080, C628B0);
    store_double(s, 0x178, C61EF0);
    store_int(s, 0x160, 2);
    store_double(s, 0x180, 0.0);
    store_double(s, 0x188, 0.0);
    store_double(s, 0x190, 0.0);
    store_int(s, 0x198, 0);
    store_int(s, 0x19C, 40);
    store_double(s, 0x1A8, C62740);
    store_int(s, 0x1A0, 800);

    store_double(s, 0x1B8, 0.0);
    store_double(s, 0x1C0, 0.0);
    store_double(s, 0x1C8, 0.0);
    store_double(s, 0x1D0, C623D0);
    store_double(s, 0x168, C62760);
    store_double(s, 0x170, C627A0);
    store_double(s, 0x1E8, C626B0);
    store_double(s, 0x1F0, C626A0);
    store_double(s, 0x1D8, C62680);
    store_double(s, 0x1E0, C62640);
    store_double(s, 0x1F8, C62790);
    store_double(s, 0x200, C62780);

// 电流环 PI 参数组（0x1320 起）：kp = 0.15、ki = 0.0012、比值 0.008、
// 输出限幅 ±471.24（2π×75）。0x1230 = 0.00025 步长；
// 0x1120/0x1168 = 电压外环 kp 4.0 / ki 0.15；0x1148/0x1150 = ±1.2×636.3 限幅。
    store_double(s, 0x1320, C61AF0);
    store_double(s, 0x1368, C61990);
    store_double(s, 0x1370, load_double(s, 0x1368) / load_double(s, 0x1320));
    store_double(s, 0x1348, C62560);
    store_double(s, 0x1350, C62A40);
    store_double(s, 0x1230, C61970);
    store_double(s, 0x1120, C62200);
    store_double(s, 0x1168, C61B00);
    store_double(s, 0x1170, load_double(s, 0x1168) / load_double(s, 0x1120));
    store_double(s, 0x1148, load_double(s, 0x138));
    store_double(s, 0x1150, -load_double(s, 0x138));
    initialize_f(s, 0x14D0);

// 四个 R 组（PI/谐振环节）初始化：0xF20、0xFA0、0x1020、0x10A0，
// 限幅 1400/30（取自 0x038/0x040）。
// 三个 H 组（谐振器）：0x1890、0x1948、0x1A00；三个 F 组（滤波器）：0x1788、0x17E0、0x1838。
    store_double(s, 0x038, C62710);
    store_double(s, 0x040, C62320);
    initialize_r(s, 0xF20, 0x038, 0x040, 0x080);
    initialize_r(s, 0xFA0, 0x038, 0x040, 0x080);
    initialize_r(s, 0x1020, 0x038, 0x040, 0x080);
    initialize_r(s, 0x10A0, 0x038, 0x040, 0x080);
    initialize_h(s, 0x1890);
    initialize_h(s, 0x1948);
    initialize_h(s, 0x1A00);
    initialize_f(s, 0x1788);
    initialize_f(s, 0x17E0);
    initialize_f(s, 0x1838);

// 额定值与 LVRT/HVRT 门槛：
//   0x1780 = 930.69（网侧线电压基准幅值）；0x270 = 1611.96（机端相电压峰值基准）；
//   0x240/0x248 = LVRT 门槛 0.854/0.942（0.9 与 0.95 标幺的截断值）；
//   0x250/0x258 = HVRT 门槛 1.123/1.02（1.13/1.02 标幺）；
//   0x260/0x268 = 保护门槛 0.088/0.06；0x278/0x280 = LVRT/HVRT 无功系数初值。
    store_double(s, 0x1780, C62600);
    store_double(s, 0x220, C628D0);
    store_double(s, 0x218, 0.0);
    store_double(s, 0x270, C62730);
    store_double(s, 0x240, load_double(s, 0x1780) * C61E20);
    store_double(s, 0x248, load_double(s, 0x1780) * C61E90);
    store_double(s, 0x250, load_double(s, 0x1780) * C61F50);
    store_double(s, 0x258, load_double(s, 0x1780) * C61F00);
    store_double(s, 0x260, load_double(s, 0x270) * C61AC0);
    store_double(s, 0x268, load_double(s, 0x270) * C61A90);
    store_double(s, 0x278, C620E0);
    store_double(s, 0x280, C62100);
    store_double(s, 0x5E8, C61EF0);
    store_double(s, 0x288, C625B0);
    store_double(s, 0x290, C62A30);
    store_double(s, 0x298, C62A20);
    store_double(s, 0x2A0, C62540);
    store_double(s, 0x5F8, C624E0);
    store_double(s, 0x2A8, C62520);
    store_double(s, 0x2B0, 0.0);
    store_double(s, 0x5A0, C62590);
    store_int(s, 0x320, 40);
    store_int(s, 0x324, 40);
    store_int(s, 0x028, 20);

// 更多常数：0x068 = 3π/4、0x070 = 2π×2.25、0x150 = 2π×1.125（斜坡步长类），
// 0x088 = 35µs；0x3A0/0x464 = 每周期采样步数（由 f 算出）；
// 0x238 = 交叉解耦系数（按频率 f 折算）；0x550/0x560 = PLL 初值 0.15/0.0012。
    store_double(s, 0x068, C62140);
    store_double(s, 0x070, C622E0);
    store_double(s, 0x150, C62250);
    store_double(s, 0x088, C61950);
    store_int(s, 0x3A0, period_count(f));
    store_int(s, 0x464, period_count(f));
    store_double(s, 0x238, (C625C0 * (C619D0 * f)) * C62060);
    store_double(s, 0x550, C61AF0);
    store_double(s, 0x560, C61990);
    store_double(s, 0x558, C61AF0);
    store_double(s, 0x568, C61930);
    store_double(s, 0xE20, C62210);
    store_double(s, 0xE28, C61B60);
    store_int(s, 0xEA8, 1);
    store_double(s, 0xE30, C622A0);

// 两组谐振器初值（0x15A8 起、0x1660 起）：中心频率 2π×300 = 1884.96（300Hz），
// 带宽 ±62.83（±10Hz），步长 0.00025，输出限幅 ±0.05×8388608。
// 原 DLL 会重复写 0x15B0/0x15B8 两处（重复但不影响），原样保留。
    store_double(s, 0x15B0, C62770);
    store_double(s, 0x15B8, C623A0);
    store_double(s, 0x1628, 0.0);
    store_double(s, 0x15A8, C61970);
    store_double(s, 0x1650, C61A80 * load_double(s, 0x080));
    store_double(s, 0x1658, C628C0 * load_double(s, 0x080));

    // The original repeats these offsets; it does not initialize 0x1668/0x1670.
    store_double(s, 0x15B0, C62770);
    store_double(s, 0x15B8, C623A0);
    store_double(s, 0x16E0, 0.0);
    store_double(s, 0x1660, C61970);
    store_double(s, 0x1708, C61A80 * load_double(s, 0x080));
    store_double(s, 0x1710, C628C0 * load_double(s, 0x080));

// 保护门槛（有名值）：0x5B8 = 0.8×1611.96、0x5A8 = 1.2×1611.96；
// 0x1528/0x1568 两组差分滤波初值（2π×1000 等）；0x5C4/0x620/0x628 开关初值；
// 0x538 = 4,000,000（额定视在功率基准 4 MVA）。
    store_double(s, 0x5B8, C61E00 * load_double(s, 0x270));
    store_double(s, 0x5A8, C61FC0 * load_double(s, 0x270));
    store_double(s, 0x1528, C62820);
    store_double(s, 0x1530, C62860);
    store_double(s, 0x1538, C62840);
    store_double(s, 0x1540, C61EF0);
    store_double(s, 0x1568, C62820);
    store_double(s, 0x1570, C62860);
    store_double(s, 0x1578, C62840);
    store_double(s, 0x1580, C61EF0);
    store_int(s, 0x5C4, 1);
    store_int(s, 0x620, 0);
    store_int(s, 0x628, 1);
    store_double(s, 0x538, C628A0);
    store_int(s, 0x668, 0);
    store_int(s, 0x704, 0);
    store_int(s, 0x708, 5);
    store_int(s, 0x70C, 1);
    store_double(s, 0x6D8, 0.0);
    store_double(s, 0x6E0, 0.0);
    store_double(s, 0x6E8, 0.0);
    store_double(s, 0x6F0, 0.0);
    store_int(s, 0x6F8, 0);
    store_int(s, 0x6FC, 5);
    store_int(s, 0x700, 0);
    store_double(s, 0x6C0, 0.0);
    store_double(s, 0x6C8, 0.0);
// 频率支撑用的两个 5 点环形缓冲清零（0x670 起、0x698 起，各 5 个 double）。
    for (std::size_t j = 0; j < 5; ++j) {
        store_double(s, 0x670 + 8 * j, 0.0);
    }
    for (std::size_t j = 0; j < 5; ++j) {
        store_double(s, 0x698 + 8 * j, 0.0);
    }
}

// Original RVA 0x29DD0. Clear the complete controller before either initializer.
// ----------------------------------------------------------------------------
// controller_initialize —— 控制器整体初始化（PSCAD 启动时调用一次）。
// 步骤：
//   1) 把整块状态区（0x3748 字节）清零；
//   2) 网侧部分（前 0x1AB8 字节）写默认值，机侧部分（0x1AB8 起）写电机默认值；
//      —— 0x1AB8 = 6838 正是网侧状态区的长度，机侧状态区紧跟在网侧之后；
//   3) 各开关量初始为 1（使能）；
//   4) 写入额定电压 1596、低压门槛 114、限幅 1.1f/0.9f、12019、2400、1000 等。
// ----------------------------------------------------------------------------
extern "C" __declspec(dllexport) void controller_initialize(double f, void* state)
{
    auto* s = static_cast<unsigned char*>(state);
    std::memset(state, 0, 0x3748);
    grid_defaults(f, state);
    motor_defaults(f, s + 0x1AB8);

    store_int(s, 0x348, 1);
    store_int(s, 0x1FAC, 1);
    store_int(s, 0x518, 1);
    store_int(s, 0x1FEC, 1);
    store_int(s, 0x51C, 1);
    store_int(s, 0x1FF0, 1);
    store_int(s, 0x2088, 1);
    store_double(s, 0x370, C62720);
    store_double(s, 0x378, C62480);
    store_double(s, 0x380, C61F40 * f);
    store_double(s, 0x388, C61E70 * f);
    store_double(s, 0x368, C62850);
    store_double(s, 0x390, C627D0);
    store_double(s, 0x398, C62640);
    store_double(s, 0x1FC8, C62830);
    store_double(s, 0x1FD0, C62720);
    store_double(s, 0x1FD8, C62880);
}
