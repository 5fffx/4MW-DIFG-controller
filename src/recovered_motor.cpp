// ============================================================================
// 文件：recovered_motor.cpp —— 机侧变流器控制（MSC / RSC，转子侧）
// ----------------------------------------------------------------------------
// 双馈电机（DFIG）的转子绕组通过滑环接到机侧变流器，再经直流母线、
// 网侧变流器接入电网。本文件控制的就是转子电流：
//   - 控制转子电流的 d 轴分量 = 控制转矩/有功功率；
//   - 控制转子电流的 q 轴分量 = 控制无功/励磁；
//   - 转子电流的频率必须等于转差频率（由转速实时算出）。
// 本文件三个函数：
//   motor_defaults   ：机侧状态区默认值；
//   motor_post       ：状态机 + 保护（每控制周期最后执行）；
//   motor_calculate  ：主控制算法（PLL、序分解、功率/电流环、调制、保护）。
// 机侧状态区从总状态区的 0x1AB8 字节处开始。
// ============================================================================
#include "recovered_primitives.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <emmintrin.h>
#include <initializer_list>

namespace {

// MotorState：帮助类，与 GridState 同款（字节偏移量读写 + apply 结构体计算）。
class MotorState {
public:
    explicit MotorState(void* state) : bytes_(static_cast<unsigned char*>(state)) {}

    double d(std::size_t offset) const { return read<double>(offset); }
    std::int32_t i(std::size_t offset) const { return read<std::int32_t>(offset); }
    void d(std::size_t offset, double value) { write(offset, value); }
    void i(std::size_t offset, std::int32_t value) { write(offset, value); }

    void increment(std::size_t offset) {
        // Increment the bit pattern without signed-overflow undefined behavior.
        write(offset, read<std::uint32_t>(offset) + std::uint32_t{1});
    }

    void* at(std::size_t offset) { return bytes_ + offset; }

    template <typename T>
    void apply(std::size_t offset, void (*function)(T*)) {
        T value = read<T>(offset);
        function(&value);
        write(offset, value);
    }

private:
    template <typename T>
    T read(std::size_t offset) const {
        T value;
        std::memcpy(&value, bytes_ + offset, sizeof(value));
        return value;
    }

    template <typename T>
    void write(std::size_t offset, T value) {
        std::memcpy(bytes_ + offset, &value, sizeof(value));
    }

    unsigned char* bytes_;
};

std::int32_t cvtt_i32(double value) {
    // CVTT preserves the hardware INT_MIN result for invalid conversions.
    return _mm_cvttsd_si32(_mm_set_sd(value));
}

struct RmsRow {
    std::size_t source;
    std::size_t accumulator;
    std::size_t output;
};

// 两个 RMS 计算的小工具：accumulate 累加平方和，finish_rms 开方并清零累加器。
template <std::size_t N>
void accumulate(MotorState& s, const RmsRow (&rows)[N]) {
    for (const RmsRow& row : rows) {
        const double source = s.d(row.source);
        s.d(row.accumulator, s.d(row.accumulator) + source * source);
    }
}

template <std::size_t N>
void finish_rms(MotorState& s, const RmsRow (&rows)[N], std::int32_t count) {
    for (const RmsRow& row : rows) {
        s.d(row.output, std::sqrt(s.d(row.accumulator) / static_cast<double>(count)));
    }
    // All outputs precede every accumulator clear.
    for (const RmsRow& row : rows) {
        s.d(row.accumulator, 0.0);
    }
}

} // namespace

// ============================================================================
// ----------------------------------------------------------------------------
// motor_defaults —— 机侧状态区默认值。按块说明如下。
// ----------------------------------------------------------------------------
// ============================================================================
extern "C" __declspec(dllexport) void motor_defaults(double p, void* state) {
    MotorState s(state);
// 计数器与基准：0x008 = 400 步（0.1 秒启动时长）、0x00C = 800 步、
// 0x504 = 800 步（停机时长）；0x358 = 8388608（大数，除法基准）；
// 0x1F8/0x200 = 正序电流环 kp 22000 / ki 150；0x218/0x220 = 负序电流环 6000/60；
// 0x1E0 = 2（极对数）；0x48 = 989.8（机侧电流基准）；0x1C8 = 362.44。
    s.d(0x348, 0.0);
    s.d(0x350, 0.0);
    s.i(0x3A0, 0);
    s.i(0x8, 400);
    s.i(0xC, 800);
    s.i(0x504, 800);
    s.i(0x4FC, 0);
    s.d(0x358, 8388608.0);
    s.d(0x1F8, 22000.0);
    s.d(0x200, 150.0);
    s.d(0x208, 6000.0);
    s.d(0x210, 60.0);
    s.d(0x218, 6000.0);
    s.d(0x220, 60.0);
    s.d(0x1E0, 2.0);
    s.i(0x1C, 0);
// 更多基准与限幅：0x1F0 = 3.0665（转矩系数类）、0x50 = 0.95×989.8、
// 0x60 = 0.5×989.8、0x40 = 5×989.8（过流门槛）、0x28 = 3000；
// 0x1190/0x11D8 = 电流环 kp 0.12 / ki 0.0012；0x11B8/0x11C0 = ±2π×75 限幅。
    s.d(0x1F0, 3.0665434957540336);
    s.d(0x48, 989.8);
    s.d(0x1C8, 362.43647999999996);
    s.d(0x1C0, s.d(0x1C8));
    s.d(0x50, s.d(0x48) * .95);
    s.d(0x60, s.d(0x48) * .5);
    s.d(0x378, 2.828);
    s.d(0x370, 14.139999999999999);
    s.d(0x360, 2.828);
    s.d(0x368, .42419999999999997);
    s.d(0x40, 5.0 * s.d(0x48));
    s.d(0x28, 3000.0);
    s.i(0x24, 0);
    s.d(0x1190, .12);
    s.d(0x11D8, .0012);
    s.d(0x11E0, s.d(0x11D8) / s.d(0x1190));
    s.d(0x11B8, 471.2388980384685);
    s.d(0x11C0, -471.2388980384685);
    s.d(0x10A0, .00025);

// 四个 PI 组（0xE10、0xE90、0xF10、0xF90）：kp/ki 取自上面的 0x1F8/0x200 或
// 0x218/0x220，输出限幅 ±0.577×8388608。
    struct PiRow {
        std::size_t base;
        std::size_t kp;
        std::size_t ki;
    };
    const PiRow pi_rows[] = {
        {0xE10, 0x1F8, 0x200},
        {0xE90, 0x1F8, 0x200},
        {0xF10, 0x218, 0x220},
        {0xF90, 0x218, 0x220}
    };
    for (const PiRow& row : pi_rows) {
        const std::size_t b = row.base;
        s.d(b, s.d(row.kp));
        s.d(b + 0x48, s.d(row.ki));
        s.d(b + 0x50, s.d(b + 0x48) / s.d(b));
        s.d(b + 0x28, .577 * s.d(0x358));
        s.d(b + 0x30, -s.d(b + 0x28));
    }

// 三组谐振器（0x1790、0x1848、0x1900）：中心 2π×50 = 314.16，带宽 ±2π×10，步长 0.00025。
// 三组滤波器（0x1688、0x16E0、0x1738）：4.0、2π×100、0.00025。
    const std::size_t resonant_bases[] = {0x1790, 0x1848, 0x1900};
    for (std::size_t b : resonant_bases) {
        s.d(b, 1.414);
        s.d(b + 0x18, 314.159265358979);
        s.d(b + 0x8, .00025);
        s.d(b + 0x28, 62.8318530717958);
        s.d(b + 0x30, -62.8318530717958);
    }

    const std::size_t filter_bases[] = {0x1688, 0x16E0, 0x1738};
    for (std::size_t b : filter_bases) {
        s.d(b + 0x8, 4.0);
        s.d(b, 628.318530717958);
        s.d(b + 0x50, .00025);
    }

// 电压门槛（标幺）：0x3B0/0x3B8 = 0.854/0.942（LVRT），0x3C0/0x3C8 = 1.123/1.02（HVRT）；
// 0x3D0/0x3D8/0x3E8/0x3F0 = 0.088/0.06/0.08/0.02（保护值）；
// 0x418/0x420 = 1.8/1.5（LVRT/HVRT 无功系数）；0x410 = 1611.96（电压基准）。
    s.d(0x1A20, 930.6928406466512);
    s.d(0x410, 1611.9599999999998);
    s.d(0x3B0, s.d(0x1A20) * .854);
    s.d(0x3B8, s.d(0x1A20) * .942);
    s.d(0x3C0, s.d(0x1A20) * 1.123);
    s.d(0x3C8, s.d(0x1A20) * 1.02);
    s.d(0x3D0, s.d(0x410) * .088);
    s.d(0x3D8, s.d(0x410) * .06);
    s.d(0x3E8, s.d(0x410) * .08);
    s.d(0x3F0, s.d(0x410) * .02);
    s.d(0x418, 1.8);
// 保护定值段：0x428 = 989.8、0x430/0x438 = -424.2/-212.1、0x450 = 565.6、
// 0x640 = 494.9、0x440 = 282.8、0x458 = 707.0（各电流/电压门槛）；
// 0x30/0x4D0/0x4D4 = 40 步、0xCC = 80 步（确认计数）；
// 0xE4/0x1AC = 每周期采样步数（由频率 f 算出）；
// 0x548/0x558/0x560 = 转矩/励磁电流换算系数（含截断常数，原样保留）。
    s.d(0x420, 1.5);
    s.d(0x630, 1.0);
    s.d(0x428, 989.8);
    s.d(0x430, -424.2);
    s.d(0x438, -212.1);
    s.d(0x450, 565.6);
    s.d(0x640, 494.9);
    s.d(0x440, 282.8);
    s.d(0x448, 0.0);
    s.d(0x458, 707.0);
    s.d(0x1468, .2);
    s.d(0x14A0, 3.106140350877193);
    s.i(0x30, 40);
    s.i(0x4D0, 40);
    s.i(0x4D4, 40);
    s.i(0xCC, 80);
    s.i(0xE4, cvtt_i32((1.0 / p) / .00025));
    s.i(0x1AC, cvtt_i32((1.0 / p) / .00025));
    s.d(0x548, -823.4127048830165);
    s.d(0x558, 922.2222294689786);
    s.d(0x560, 922.2222294689786);
    s.d(0x550, 0.0);
    s.d(0x5A0, .12);
    s.d(0x5B0, .0012);
    s.d(0x5A8, .12);
    s.d(0x5B8, .000012);
    s.d(0x5F0, 3.106140350877193);
    s.d(0x600, 1.2 * s.d(0x410));
    s.d(0x5F8, 56.0 - p);
    s.i(0x6E8, 1);
    s.d(0xD10, 5.0);
    s.d(0xD18, .2);
    s.i(0xD98, 1);
    s.d(0xD20, 20.0);
    s.i(0x60C, 1);
    s.i(0x668, 0);
    s.i(0x66C, 1);
    s.i(0x690, 0);
    s.d(0x6B8, .3);
    s.d(0x6C0, .3);
    s.d(0x598, 4000000.0);
    s.d(0x14B0, 2199.114857512853);
    s.d(0x14B8, 62.8318530717958);
    s.d(0x1528, 0.0);
    s.d(0x14A8, .00025);
    s.d(0x1550, .05 * s.d(0x358));
    s.d(0x1558, -.05 * s.d(0x358));
    s.d(0x1568, 2199.114857512853);
    s.d(0x1570, 62.8318530717958);
    s.d(0x15E0, 0.0);
    s.d(0x1560, .00025);
    s.d(0x1608, .05 * s.d(0x358));
    s.d(0x1610, -.05 * s.d(0x358));
}

// ============================================================================
// ----------------------------------------------------------------------------
// motor_post —— 机侧状态机 + 保护（每个控制周期最后执行）
// ----------------------------------------------------------------------------
// 一、运行状态机（0x4FC）：0 = 停机，1 = 运行，2 = 停机过程中，3 = 跳闸。
//     跳闸时清输出和标志；复位命令（0x53C == 1）时整块复位重写默认值。
// 二、RMS 有效值：9 路网侧相关量 + 3 路转子电压，平方累加 → 窗口满后开方；
//     转子电压的窗口长度按转差频率自适应（频率太低时钳位到 ±0.1）。
// 三、故障锁存（0x934）：一类故障立即转 3；二类故障延时 4 秒转 3。
// ----------------------------------------------------------------------------
// ============================================================================
extern "C" __declspec(dllexport) void motor_post(double p, void* state) {
    MotorState s(state);
    if (s.i(0x4F8) == 1) {
        s.i(0x4FC, 3);
    } else if (s.i(0x4FC) == 1 && s.i(0x4F4) == 0) {
        s.i(0x4FC, 2);
    } else if (s.i(0x4FC) == 0 && s.i(0x4F4) == 1) {
        s.i(0x4FC, 1);
    }

    if (s.i(0x4FC) == 3) {
        s.i(0x10, 0);
        s.i(0x1C, 0);
        s.i(0x3A0, 0);
        s.i(0x3A4, 0);
        const std::size_t clear_offsets[] = {
            0x350, 0x348, 0x300, 0x318, 0x330,
            0x2F0, 0x2F8, 0xC78, 0xC80, 0xC88
        };
        for (std::size_t offset : clear_offsets) {
            s.d(offset, 0.0);
        }
    }

    if (s.i(0x53C) == 1) {
        std::memset(state, 0, 0x1C68);
        motor_defaults(p, state);
        return;
    }

    const RmsRow rows_a[] = {
        {0x258, 0x118, 0xE8},
        {0x260, 0x128, 0xF8},
        {0x268, 0x138, 0x108},
        {0x240, 0x190, 0x178},
        {0x248, 0x198, 0x180},
        {0x250, 0x1A0, 0x188},
        {0x19B8, 0x120, 0xF0},
        {0x19C0, 0x130, 0x100},
        {0x19C8, 0x140, 0x110}
    };
    s.increment(0x1A8);
    accumulate(s, rows_a);
    if (s.i(0x1A8) >= s.i(0x1AC)) {
        finish_rms(s, rows_a, s.i(0x1A8));
        s.i(0x1A8, 0);
        s.i(0x1B0, 1);
    }

    const RmsRow rows_b[] = {
        {0x228, 0x160, 0x148},
        {0x230, 0x168, 0x150},
        {0x238, 0x170, 0x158}
    };
    s.increment(0x1B4);
    accumulate(s, rows_b);
    double x = s.d(0x5F8);
    if (x < .1 && x > 0.0) {
        x = .1;
    } else if (x > -.1 && x < 0.0) {
        x = -.1;
    }
    s.i(0x1B8, cvtt_i32(std::fabs((1.0 / x) / .00025)));
    if (s.i(0x1B4) >= s.i(0x1B8)) {
        finish_rms(s, rows_b, s.i(0x1B4));
        s.i(0x1B4, 0);
        s.i(0x1BC, 1);
    }

// ---- 故障锁存：本拍先置状态，输出清零要等下一拍（原 DLL 的行为，原样保留） ----
    // Late faults latch now; output clearing waits until the next call.
    if ((s.i(0x934) & 1) != 0) {
        s.i(0x4FC, 3);
    } else if ((s.i(0x934) & 2) != 0) {
        if (s.i(0x20) == 1) {
            s.d(0x938, s.d(0x938) + .00025);
            if (s.d(0x938) >= 4.0) {
                s.i(0x4FC, 3);
                s.d(0x938, 4.0);
            }
        } else {
            s.i(0x4FC, 3);
        }
    } else {
        s.d(0x938, 0.0);
    }
}

// ============================================================================
// ----------------------------------------------------------------------------
// motor_calculate —— 机侧变流器主控制（每个控制周期执行一次）
// ----------------------------------------------------------------------------
// 整体信号流：
//   1) 功率/电流换算：0x1C8 功率折算；0x548/0x558/0x560 转矩/励磁电流换算系数
//      （含转速、额定功率、转差频率等）；
//   2) 锁相环 PLL（0xFF8 起，锁转子电压）→ 转子电压角度/频率；
//   3) 三相自适应振荡器 + 两态滤波器（0x1790/0x1688 起）：转子电压基波幅值，
//      及三相最大/最小（0x1A00 起，供 LVRT/HVRT 判定）；
//   4) 转子电流正负序分解（两套，0x11F8 与 0x12F0）；
//   5) 启动逻辑（0x1C 标志 + 计数）：先空跑 400 步，等电压建立后切换；
//   6) LVRT/HVRT 模式判定（0x10 状态机，20/4 步确认 + 滞回退出）：
//      - 模式 1（LVRT）：转子电压低 → 励磁电流按 LVRT 规则给，限幅与保持；
//      - 模式 2（HVRT）：转子电压高 → 按 HVRT 规则吸收无功；
//      - 模式 0（正常）：转速/有功外环 → 有功电流给定；
//   7) 电流内环 PI ×4（0xDF8/0xE78/0xEF8/0xF78）+ 300Hz 谐振滤波 + 谐波抑制；
//   8) 前馈 + 反 Park → 调制（0xBF0，普通/备用两版）→ 输出占空比 0x308 起；
//   9) 保护：过流（±0x40 门槛）→ 闭锁输出；crowbar（撬棒）逻辑（0x5F8 转速条件）；
//  10) 停机（0x4FC==2）持续 0x504 步后整块复位。
// ----------------------------------------------------------------------------
// ============================================================================
extern "C" __declspec(dllexport) void motor_calculate(double base_frequency, void* state) {
    MotorState s(state);
    const auto clip = [&s](std::size_t offset, double lower, double upper) {
        if (s.d(offset) > upper) s.d(offset, upper);
        else if (lower > s.d(offset)) s.d(offset, lower);
    };
    const auto wrap = [&s](std::size_t offset) {
        if (s.d(offset) > 6.28318530717958) s.d(offset, s.d(offset) - 6.28318530717958);
        else if (0.0 > s.d(offset)) s.d(offset, s.d(offset) + 6.28318530717958);
    };
    const auto slew = [&s](std::size_t destination, std::size_t target, std::size_t step) {
        if (s.d(target) > s.d(destination)) {
            s.d(destination, s.d(destination) + s.d(step));
            if (s.d(destination) > s.d(target)) s.d(destination, s.d(target));
        } else if (s.d(destination) > s.d(target)) {
            s.d(destination, s.d(destination) - s.d(step));
            if (s.d(target) > s.d(destination)) s.d(destination, s.d(target));
        }
    };
    const auto maximum = [](double a, double b, double c) {
        const double q = a > b ? a : b;
        return q > c ? q : c;
    };
    const auto minimum = [](double a, double b, double c) {
        const double q = a > b ? b : a;
        return q > c ? c : q;
    };

// ---- 第 1 步：功率与电流换算系数（0x1C8、0x548、0x558、0x560） ----
    s.d(0x1C8, (((s.d(0x6D0) / 1.732) / (s.d(0x6D8) / s.d(0x6D0))) / s.d(0x6C8)) * 1.414);
    s.d(0x548, ((((((-s.d(0x598)) * 3000.0) / s.d(0x6E0)) / s.d(0x1E0)) / s.d(0x6D8)) / 1.732) * 1.414 * (s.d(0xB28) / base_frequency));
    s.d(0x558, ((s.d(0x598) / s.d(0x6D8)) / 1.732) * 1.414);
    s.d(0x560, ((s.d(0x598) / s.d(0x6D8)) / 1.732) * 1.414);
    s.d(0x1190, s.d(0x5A0));
    s.d(0x11D8, s.d(0x5B0));
// ---- 电流环 PI 参数与限幅刷新（按模式 0x10 选 ±2π×1.25 或 ±2π×10） ----
    if (s.i(0x10) == 1) {
        s.d(0x11B8, 7.853981633974475);
        s.d(0x11C0, -7.853981633974475);
    } else {
        s.d(0x11B8, 62.8318530717958);
        s.d(0x11C0, -62.8318530717958);
    }
    if (s.d(0x2A0) >= 0.6) s.i(0x57C, 1);
    else if (0.6 > s.d(0x2A0)) s.i(0x57C, 0);
    if (s.i(0x6E8) == 1) {
        s.d(0x2A0, s.d(0x2A0) * 1.04);
        s.d(0x2A8, s.d(0x2A8) + 0.015);
    } else {
        s.d(0x2A0, s.d(0x2A0) * 1.045);
        s.d(0x2A8, s.d(0x2A8) + 0.015);
    }
// ---- 第 2 步：锁相环（0xFF8 起，锁转子电压 0x258/0x260/0x268） ----
    s.d(0xFF8, s.d(0x258));
    s.d(0x1000, s.d(0x260));
    s.d(0x1008, s.d(0x268));
    {
        PllState p;
        std::memcpy(&p, s.at(0xFF8), sizeof(p));
        pll(base_frequency, &p);
        std::memcpy(s.at(0xFF8), &p, sizeof(p));
    }
    s.d(0x1E8, s.d(0x1018));
    s.d(0x3A8, std::sqrt(s.d(0x1040) * s.d(0x1040) + s.d(0x1048) * s.d(0x1048)));
// ---- 第 3 步：三相自适应振荡器 + 两态滤波器（转子电压基波提取） ----
    for (int phase = 0; phase < 3; ++phase) {
        const std::size_t adaptive = 0x1790 + 0xB8 * phase;
        const std::size_t filter = 0x1688 + 0x58 * phase;
        s.d(adaptive + 0x40, s.d(0x19B8 + 8 * phase));
        adaptive_oscillator(s.at(adaptive), base_frequency);
        s.d(filter + 0x10, s.d(adaptive + 0x68));
        s.d(filter, s.d(0x1010) * 2.0);
        s.apply(filter, two_state_filter);
        s.d(0x19D0 + 8 * phase, s.d(filter + 0x18));
        s.d(0x19E8 + 8 * phase, s.d(0x19D0 + 8 * phase) / s.d(0x1A20));
    }
    s.d(0x1A00, maximum(s.d(0x19D0), s.d(0x19D8), s.d(0x19E0)));
    s.d(0x1A08, minimum(s.d(0x19D0), s.d(0x19D8), s.d(0x19E0)));
    s.d(0x1A10, maximum(s.d(0x19E8), s.d(0x19F0), s.d(0x19F8)));
    s.d(0x1A18, minimum(s.d(0x19E8), s.d(0x19F0), s.d(0x19F8)));
// ---- 直流母线电压平均（0xB8）与转子电压幅值（0xD0/0x1040） ----
    s.d(0x278, s.d(0x280) + s.d(0x288));
    s.increment(0xC8);
    s.d(0xC0, s.d(0xC0) + s.d(0x278));
    if (s.i(0xC8) >= s.i(0xCC)) {
        s.d(0xB8, s.d(0xC0) / static_cast<double>(s.i(0xCC)));
        s.i(0xC8, 0);
        s.d(0xC0, 0.0);
    }
    s.increment(0xE0);
    s.d(0xD8, s.d(0xD8) + s.d(0x1040));
    if (s.i(0xE0) >= s.i(0xE4)) {
        s.i(0xE0, 0);
        s.d(0xD8, 0.0);
    }
    s.d(0xD0, s.d(0x1040));
    wrap(0x290);
    s.d(0x2E0, (s.d(0x1E8) - s.d(0x290)) + s.d(0x1F0));
    s.d(0x2E8, ((-s.d(0x1E8)) - s.d(0x290)) + s.d(0x1F0));
    wrap(0x2E0);
    wrap(0x2E8);
// ---- 闸门：只有“运行/停机过程”状态才执行下面的控制 ----
    if (s.i(0x4FC) != 1 && s.i(0x4FC) != 2) goto final_status;

// ---- 第 4 步：转子电流正负序分解（两套：0x11F8 与 0x12F0） ----
    s.d(0x11F8, s.d(0x228));
    s.d(0x1200, s.d(0x230));
    s.d(0x1208, s.d(0x238));
    s.d(0x1210, s.d(0x2E0));
    s.d(0x1218, s.d(0x2E8));
    rotor_current_sequence(s.at(0x11F8));
    s.d(0x12F0, s.d(0x240));
    s.d(0x12F8, s.d(0x248));
    s.d(0x1300, s.d(0x250));
    s.d(0x1308, s.d(0x1E8));
    s.d(0x1310, -s.d(0x1E8));
    rotor_current_sequence(s.at(0x12F0));
    if (s.i(0x1C) == 0) {
        s.d(0x2C0, s.d(0x1260));
        s.d(0x2C8, s.d(0x1268));
        s.d(0x2D0, 0.0);
        s.d(0x2D8, 0.0);
    } else if (s.i(0x728) == 0) {
        s.d(0x2C0, s.d(0x1220));
        s.d(0x2C8, s.d(0x1228));
        s.d(0x2D0, s.d(0x1230));
        s.d(0x2D8, s.d(0x1238));
    } else {
        s.d(0x2C0, s.d(0x1280));
        s.d(0x2C8, s.d(0x1288));
        s.d(0x2D0, s.d(0x1290));
        s.d(0x2D8, s.d(0x1298));
    }
// ---- 第 5 步：启动逻辑。前 400 步直接复位返回（等电路建立）；
//      0x1C == 0 表示“预同步阶段”：等 0xA8（功率折算）与 0xA0 相等
//      满 200 步后置 0x1C = 1，并把四组 PI 参数按额定值刷新 ----
    if (s.i(0) < s.i(8)) {
        s.increment(0);
        goto exit_reset;
    }
    if (s.i(0x1C) == 0) {
        s.d(0xA0, s.d(0x1C8));
        if (s.d(0xA8) == s.d(0xA0)) {
            if (s.i(4) < 200) s.increment(4);
            else {
                s.i(0x1C, 1);
                s.i(0x3A0, 1);
                s.i(4, 0);
                s.d(0xE10, s.d(0x208));
                s.d(0xE58, s.d(0x210));
                s.d(0xE60, s.d(0xE58) / s.d(0xE10));
                s.d(0xE90, s.d(0x208));
                s.d(0xED8, s.d(0x210));
                s.d(0xEE0, s.d(0xED8) / s.d(0xE90));
                s.d(0xF10, s.d(0x218));
                s.d(0xF58, s.d(0x220));
                s.d(0xF60, s.d(0xF58) / s.d(0xF10));
                s.d(0xF90, s.d(0x218));
                s.d(0xFD8, s.d(0x220));
                s.d(0xFE0, s.d(0xFD8) / s.d(0xF90));
            }
        }
        // Startup still executes the common ramps, regulators, and modulation.
        goto common;
    }
// ---- 再等 100 步稳定，然后进入正常控制 ----
    if (s.i(4) < 100) {
        s.increment(4);
        goto common;
    }
// ---- 转矩/励磁电流给定：0x2B0 = 转矩电流（模式 1 时按系数 0x548），
//      0x2B8 = 励磁电流（按转速符号选 0x558/0x560）----
    if (s.i(0x4FC) == 1) {
        s.d(0x2B0, s.d(0x2A0) * s.d(0x548));
        if (s.d(0x2A8) >= 0.0) s.d(0x2B8, s.d(0x2A8) * s.d(0x558));
        else s.d(0x2B8, s.d(0x2A8) * s.d(0x560));
    } else {
        s.d(0x2B0, 0.0);
        s.d(0x2B8, 0.0);
    }
// ---- 过压/过流预判标志（0x3E0/0x3F8，20 步确认）与 LVRT/HVRT 模式判定 ----
    if (s.d(0x1088) > s.d(0x3D0)) {
        s.i(0x4BC, 0);
        if (s.i(0x4A4) < 20) s.increment(0x4A4);
        else s.i(0x3E0, 1);
    } else s.i(0x4A4, 0);
    if (s.i(0x3E0) == 1) {
        if (s.d(0x3D8) > s.d(0x1088)) {
            if (s.i(0x4BC) < 20) s.increment(0x4BC);
            else s.i(0x3E0, 0);
        } else s.i(0x4BC, 0);
    }
    if (s.d(0x1088) > s.d(0x3E8)) {
        s.i(0x404, 0);
        if (s.i(0x400) < 20) s.increment(0x400);
        else s.i(0x3F8, 1);
    } else s.i(0x400, 0);
    if (s.i(0x3F8) == 1) {
        if (s.d(0x3F0) > s.d(0x1088)) {
            if (s.i(0x404) < 20) s.increment(0x404);
            else s.i(0x3F8, 0);
        } else s.i(0x404, 0);
    }
    if (s.d(0x3B0) > s.d(0x1A08) && s.i(0x4C4) == 0 && s.i(0x534) == 1 && s.d(0x410) >= s.d(0x3A8)) {
        s.i(0x498, 0);
        s.i(0x4A8, 0);
        if (s.i(0x490) < 20) s.increment(0x490);
        else {
            if (s.i(0x10) != 1) s.increment(0x970);
            s.i(0x10, 1);
        }
    } else if (s.d(0x1A00) > s.d(0x3C0) && s.i(0x538) == 1 && s.d(0x3A8) > s.d(0x410) && s.d(0x580) == 0.0) {
        s.i(0x408, 0);
        s.i(0x490, 0);
        s.i(0x4B0, 0);
        if (s.i(0x498) < 4) s.increment(0x498);
        else {
            if (s.i(0x10) != 2) s.increment(0x974);
            s.i(0x10, 2);
        }
    } else {
        s.i(0x490, 0);
        s.i(0x498, 0);
        s.i(0x408, 0);
        if (s.i(0x10) == 1 && s.i(0x4C0) == 0) {
            if (s.d(0x1A08) > s.d(0x3B8)) {
                s.increment(0x4A8);
                if (s.i(0x4A8) > 4) {
                    s.i(0x4A8, 0);
                    s.i(0x10, 0);
                    s.d(0x480, 0.0);
                    s.d(0x488, 0.0);
                    s.d(0x468, 0.0);
                    s.i(0x408, 0);
                    s.d(0x88, 0.0);
                    s.i(0x570, 1);
                    s.d(0x1C10, 0.0);
                }
            } else s.i(0x4A8, 0);
        } else if (s.i(0x10) == 2 && s.i(0x4C4) == 0) {
            if (s.d(0x3C8) > s.d(0x1A00) && 56.559999999999995 > s.d(0x1088)) {
                s.increment(0x4B0);
                if (s.i(0x4B0) > 4) {
                    s.i(0x4B0, 0);
                    s.i(0x10, 0);
                    s.d(0x468, 0.0);
                }
            } else if (1.02 * s.d(0x1A20) > s.d(0x1A00)) {
                s.increment(0x4B0);
                if (s.i(0x4B0) > 4) {
                    s.i(0x4B0, 0);
                    s.i(0x10, 0);
                    s.d(0x468, 0.0);
                }
            } else s.i(0x4B0, 0);
        }
    }
// ---- 模式保持：穿越结束后维持 800 步（0.2 秒）再回正常 ----
    if (s.i(0x10) != 0) {
        s.i(0x14, 1);
        s.i(0x18, 0);
    } else if (s.i(0x14) == 1) {
        if (s.i(0x18) < 800) s.increment(0x18);
        else {
            s.i(0x18, 0);
            s.i(0x14, 0);
            s.d(0x478, 0.0);
        }
    }
// ---- 模式 1（LVRT）控制段：
//      0x4E0 = 1 - 有功/额定（功率偏差）；0x470 = 转速偏差经 PI（0x418 系数）
//      得有功电流给定 0x4D8；0x470/0x460 为 d/q 电流给定与限幅，
//      0x478 = √(1.21 - Id²)×基准（无功电流上限）----
    if (s.i(0x10) == 1) {
        s.increment(0x408);
        s.d(0x5C8, 0.0);
        s.i(0x4CC, 0);
        if (s.i(0x4C8) < s.i(0x4D4)) {
            s.increment(0x4C8);
            s.i(0x4C0, 1);
        } else s.i(0x4C0, 0);
        s.d(0x4E0, 1.0 - s.d(0x1040) / s.d(0x410));
        s.d(0x470, s.d(0x6F0) - s.d(0x3A8) / s.d(0x410));
        s.d(0x5C0, (30.0 * s.d(0x470) + 994.0 * s.d(0x5C0)) / 1024.0);
        if (s.i(0x60C) == 1) s.d(0x4D8, s.d(0x418) * s.d(0x5C0));
        else s.d(0x4D8, s.d(0x418) * s.d(0x470));
        if (s.d(0x4D8) > 0.9863) s.d(0x4D8, 0.9863);
        else if (0.0 > s.d(0x4D8)) s.d(0x4D8, s.d(0xB68));
        s.d(0x470, s.d(0x4D8) * s.d(0x48));
        if (s.i(0x60C) == 1 || s.d(0x470) > s.d(0x468)) {
            s.d(0x468, s.d(0x470));
            s.d(0x460, s.d(0x468));
            s.d(0x478, std::sqrt(1.21 - s.d(0x4D8) * s.d(0x4D8)) * s.d(0x48));
            if (s.d(0x478) > s.d(0x458)) s.d(0x478, s.d(0x458));
            if (212.1 > s.d(0x478)) s.d(0x478, 212.1);
        }
        {
            double u = s.d(0x1050) / s.d(0x410);
            if (u > 0.48828125) u = 0.48828125;
            else if (-0.48828125 > u) u = -0.48828125;
            s.d(0x610, s.d(0x630) * u);
        }
        clip(0x610, -0.48828125, 0.48828125);
        if (s.i(0x60C) == 1) s.d(0x618, (30.0 * s.d(0x610) + 994.0 * s.d(0x618)) / 1024.0);
        else s.d(0x618, s.d(0x610));
        if (s.i(0x3E0) == 1 && s.i(0x668) == 1) {
            s.d(0x610, s.d(0x618) * s.d(0x48));
            clip(0x610, -s.d(0x640), s.d(0x640));
        } else s.d(0x610, 0.0);
        {
            double v = s.d(0x1058) / s.d(0x410);
            if (v > 0.48828125) v = 0.48828125;
            else if (-0.48828125 > v) v = -0.48828125;
            s.d(0x620, (-s.d(0x630)) * v);
        }
        clip(0x620, -0.48828125, 0.48828125);
        if (s.i(0x60C) == 1) s.d(0x628, (30.0 * s.d(0x620) + 994.0 * s.d(0x628)) / 1024.0);
        else s.d(0x628, s.d(0x620));
        if (s.i(0x3E0) == 1 && s.i(0x668) == 1) {
            s.d(0x620, s.d(0x628) * s.d(0x48));
            clip(0x620, -s.d(0x640), s.d(0x640));
        } else s.d(0x620, 0.0);
        if (s.i(0x3E0) == 1 && s.d(0x460) > s.d(0x450)) {
            s.d(0x460, s.d(0x450));
            s.d(0x468, s.d(0x460));
        }
        if (s.d(0x1020) > s.d(0x410)) {
            if (s.d(0xA8) > s.d(0x440)) s.d(0xA8, s.d(0x440));
            if (s.d(0x460) > s.d(0x440)) {
                s.d(0x460, s.d(0x440));
                s.d(0x468, s.d(0x460));
            }
        }
        if (s.d(0x460) > s.d(0x428)) s.d(0x460, s.d(0x428));
// ---- 模式 2（HVRT）控制段：吸收无功（0x420 系数），限幅 0.3418，
//      0x478 按 √(I² - Id²) 重算 ----
    } else if (s.i(0x10) == 2) {
        if (s.i(0x7A8) == 1) s.i(0x570, 1);
        s.d(0x5C0, 0.0);
        s.d(0x480, 0.0);
        s.i(0x408, 0);
        s.i(0x4C8, 0);
        if (s.i(0x4CC) < s.i(0x4D0)) {
            s.increment(0x4CC);
            s.i(0x4C4, 1);
        } else s.i(0x4C4, 0);
        s.d(0x470, s.d(0x1040) / s.d(0x410) - s.d(0x6F8));
        s.d(0x5C8, (30.0 * s.d(0x470) + 994.0 * s.d(0x5C8)) / 1024.0);
        if (s.i(0x60C) == 1) s.d(0x470, s.d(0x420) * s.d(0x5C8));
        else s.d(0x470, s.d(0x420) * s.d(0x470));
        s.d(0x5E0, (20.0 * s.d(0x1088) + 1004.0 * s.d(0x5E0)) / 1024.0);
        s.d(0x5E8, (0.5 * s.d(0x5E0)) / s.d(0x410));
        if (0.0 > s.d(0x5E8)) s.d(0x5E8, 0.0);
        else if (s.d(0x5E8) > 0.1) s.d(0x5E8, 0.1);
        if (0.0 > s.d(0x470)) s.d(0x470, 0.0);
        else if (s.d(0x470) > 0.3418) s.d(0x470, 0.3418);
        s.d(0x470, (-s.d(0x470)) * s.d(0x48));
        if (s.i(0x60C) == 1 || s.d(0x468) > s.d(0x470)) {
            s.d(0x468, s.d(0x470));
            s.d(0x460, s.d(0x468));
        }
        if (s.d(0x460) > s.d(0x438)) s.d(0x460, s.d(0x438));
        if (s.d(0x430) > s.d(0x460)) s.d(0x460, s.d(0x430));
        s.d(0x478, s.d(0x48) * s.d(0x48) - s.d(0x460) * s.d(0x460));
        if (0.0 > s.d(0x478)) s.d(0x478, 0.0);
        s.d(0x478, std::sqrt(s.d(0x478)));
    }
// ---- 0x7A8 记住上一拍模式；0x1C0 是 80 步平均功率 ----
    s.i(0x7A8, s.i(0x10));
    s.d(0x1D8, s.d(0x1D8) + 1.0);
    s.d(0x1D0, s.d(0x1D0) + (s.d(0x1C8) * s.d(0x1040)) / s.d(0x410));
    if (s.d(0x1D8) >= 80.0) {
        s.d(0x1C0, s.d(0x1D0) / 80.0);
        s.d(0x1D0, 0.0);
        s.d(0x1D8, 0.0);
    }
// ---- 励磁电流给定 0x88（多分支：正常/惯量/恢复模式），并限幅 ----
    s.d(0x88, (s.d(0x2B0) * s.d(0x410)) / s.d(0xD0));
    if (s.d(0x88) > 0.0) s.d(0x88, (-s.d(0xB90)) * 1.414);
    else if (-s.d(0x50) > s.d(0x88)) s.d(0x88, -s.d(0x50));

// ---- 模式 1 的功率环：0xA0 = 功率电流给定合成 + 限幅；
//      下面的 0x5D4 状态机是“crowbar（撬棒）动作档位”编码：
//      121/481/861（并网运行）与 122/482/862（脱网）等档位按功率偏差选择 ----
    if (s.i(0x10) == 1) {
        if (s.i(0x568) < 0) s.increment(0x568);
        else {
            if (s.i(0x57C) == 1) {
                if (s.i(0x710) == 1) {
                    if (s.d(0x4E0) <= 0.2) {
                        if (s.d(0x1C10) == 1.0) s.d(0x88, s.d(0x1C30) * s.d(0x548));
                    } else {
                        s.d(0x88, s.d(0x1C30) * s.d(0x548));
                        s.d(0x1C10, 1.0);
                    }
                } else if (s.i(0x710) == 2) s.d(0x88, (s.d(0x2B0) * s.d(0x410)) / s.d(0xD0));
            } else {
                volatile double ratio = s.d(0x3A8) / s.d(0x410);
                (void)ratio;
                s.d(0x88, (s.d(0x2B0) * s.d(0x410)) / s.d(0xD0));
            }
        }
        s.d(0xA0, (s.d(0x1C0) + s.d(0x460)) + s.d(0x2B8) * s.d(0x700));
        if (s.d(0xA0) > 1.0 * s.d(0x48)) s.d(0xA0, 1.0 * s.d(0x48));
        if (s.i(0x60C) == 1) {
            s.d(0x468, s.d(0x470));
            s.d(0x460, s.d(0x468));
            s.d(0x478, (1.21 * s.d(0x48)) * s.d(0x48) - s.d(0xA0) * s.d(0xA0));
            if (0.0 > s.d(0x478)) s.d(0x478, 0.0);
            else s.d(0x478, std::sqrt(s.d(0x478)));
        }
        if (-s.d(0x478) > s.d(0x88)) s.d(0x88, -s.d(0x478));
        if (s.d(0x88) > (-s.d(0xB90)) * 1.414) s.d(0x88, (-s.d(0xB90)) * 1.414);
        s.d(0x648, s.d(0x620));
        s.d(0x650, s.d(0x610));
        if ((s.i(0x3E0) == 0 && 120.19 > s.d(0x1088)) || 0.55 * s.d(0x410) > s.d(0x3A8)) {
            if ((s.d(0x4E0) > 0.6 && s.i(0x57C) == 1) || s.i(0x5D4) == 121) {
                s.i(0x5D4, 121);
                s.i(0x774, 1);
            } else if ((s.d(0x4E0) > 0.35 && s.i(0x57C) == 1 && 0.45 > s.d(0x4E0)) || s.i(0x5D4) == 481) {
                s.i(0x5D4, 481);
                s.i(0x778, 1);
            } else if ((s.d(0x4E0) > 0.1 && s.i(0x57C) == 1 && 0.15 > s.d(0x4E0)) || s.i(0x5D4) == 861) {
                s.i(0x5D4, 861);
            } else if ((s.d(0x4E0) > 0.8 && s.i(0x57C) == 0) || s.i(0x5D4) == 122) {
                s.i(0x5D4, 122);
            } else if ((s.d(0x4E0) > 0.4 && s.i(0x57C) == 0) || s.i(0x5D4) == 482) {
                s.i(0x5D4, 482);
            } else if ((0.15 > s.d(0x4E0) && s.d(0x4E0) > 0.1 && s.i(0x57C) == 0) || s.i(0x5D4) == 862) {
                s.i(0x5D4, 862);
            }
        } else if (s.i(0x3E0) == 1 || s.d(0x1088) > 120.19 || 0.55 * s.d(0x410) > s.d(0x3A8)) {
            if ((s.d(0x4E0) > 0.41 && s.i(0x57C) == 1) || s.i(0x5D4) == 123) {
                s.i(0x5D4, 123);
                s.i(0x784, 1);
            } else if ((s.d(0x4E0) > 0.18 && s.i(0x57C) == 1) || s.i(0x5D4) == 483) {
                s.i(0x5D4, 483);
                s.i(0x788, 1);
            } else if ((s.d(0x4E0) > 0.12 && s.i(0x57C) == 1) || s.i(0x5D4) == 703) {
                s.i(0x5D4, 703);
                s.i(0x78C, 1);
            } else if ((s.d(0x4E0) > 0.1 && s.i(0x57C) == 1) || s.i(0x5D4) == 863) {
                s.i(0x5D4, 863);
                s.i(0x790, 1);
            } else if ((s.d(0x4E0) > 0.4 && s.i(0x57C) == 0) || s.i(0x5D4) == 124) {
                s.i(0x5D4, 124);
            } else if ((s.d(0x4E0) > 0.22 && s.i(0x57C) == 0) || s.i(0x5D4) == 484) {
                s.i(0x5D4, 484);
            } else if ((s.d(0x4E0) > 0.15 && s.i(0x57C) == 0) || s.i(0x5D4) == 704) {
                s.i(0x5D4, 704);
            } else if ((s.d(0x4E0) > 0.05 && s.i(0x57C) == 0) || s.i(0x5D4) == 864) {
                s.i(0x5D4, 864);
            }
        }
// ---- 模式 2 的保护档位状态机（0x5D8 编码 1361/1241/1131……） ----
    } else if (s.i(0x10) == 2) {
        s.d(0xA0, (s.d(0x1C0) + s.d(0x460)) + s.d(0x2B8) * s.d(0x708));
        if (-s.d(0x478) > s.d(0x88)) s.d(0x88, -s.d(0x478));
        s.d(0x648, 0.0);
        s.d(0x650, 0.0);
        if (s.i(0x3E0) == 0 && 42.419999999999995 > s.d(0x1088)) {
            if ((s.d(0x3A8) / s.d(0x410) > 1.2 && s.i(0x57C) == 1) || s.i(0x5D8) == 1361) {
                s.i(0x5D8, 1361);
            } else if ((s.d(0x3A8) / s.d(0x410) > 1.16 && s.i(0x57C) == 1) || s.i(0x5D8) == 1241) {
                s.i(0x5D8, 1241);
            } else if ((s.d(0x3A8) / s.d(0x410) > 1.05 && s.i(0x57C) == 1) || s.i(0x5D8) == 1131) {
                s.i(0x5D8, 1131);
            } else if ((s.d(0x3A8) / s.d(0x410) > 1.25 && s.i(0x57C) == 0) || s.i(0x5D8) == 1362) {
                s.i(0x5D8, 1362);
            } else if ((s.d(0x3A8) / s.d(0x410) > 1.2 && s.i(0x57C) == 0) || s.i(0x5D8) == 1242) {
                s.i(0x5D8, 1242);
            } else if ((s.d(0x3A8) / s.d(0x410) > 1.05 && s.i(0x57C) == 0) || s.i(0x5D8) == 1132) {
                s.i(0x5D8, 1132);
            }
        } else if (s.i(0x3E0) == 1 || s.d(0x1088) > 42.419999999999995) {
            if ((s.d(0x1088) > 219.17 && s.i(0x57C) == 1) || s.i(0x5D8) == 1363) {
                s.i(0x5D8, 1363);
            } else if ((s.d(0x1088) > 169.67999999999998 && s.i(0x57C) == 1) || s.i(0x5D8) == 1243) {
                s.i(0x5D8, 1243);
            } else if ((s.d(0x1088) > 98.97999999999999 && s.i(0x57C) == 1) || s.i(0x5D8) == 1133) {
                s.i(0x5D8, 1133);
            } else if ((s.d(0x1088) > 226.23999999999998 && s.i(0x57C) == 0) || s.i(0x5D8) == 1364) {
                s.i(0x5D8, 1364);
            } else if ((s.d(0x1088) > 169.67999999999998 && s.i(0x57C) == 0) || s.i(0x5D8) == 1244) {
                s.i(0x5D8, 1244);
            } else if ((s.d(0x1088) > 113.11999999999999 && s.i(0x57C) == 0) || s.i(0x5D8) == 1134) {
                s.i(0x5D8, 1134);
            }
        }
    } else {
        s.i(0x568, 0);
        s.d(0x550, 0.0);
        s.i(0x4CC, 0);
        s.i(0x4C8, 0);
        s.i(0x5D4, 0);
        s.i(0x5D8, 0);
        s.d(0x5C0, 0.0);
        s.d(0x5C8, 0.0);
        s.d(0x480, 0.0);
        s.d(0x488, 0.0);
        s.d(0x460, 0.0);
        s.d(0x4E0, 0.0);
        clip(0x2B8, -s.d(0x60), s.d(0x60));
        s.d(0x2B8, (s.d(0x2B8) * s.d(0x410)) / s.d(0xD0));
        s.d(0xA0, s.d(0x1C0) + s.d(0x2B8));
        s.d(0x648, 0.0);
        s.d(0x650, 0.0);
        s.d(0x7A0, s.d(0xD0));
    }

// ---- common：三种模式共用的后续处理（从这里起每步都执行） ----
common:
// ---- 斜坡步长选择：0x390/0x398 按运行状态给不同的爬升速度 ----
    if (s.i(0x14) == 0) {
        if (s.i(0x1C) == 1) {
            s.d(0x390, s.d(0x360));
            s.d(0x398, s.d(0x360));
        } else {
            s.d(0x390, s.d(0x368));
            s.d(0x398, s.d(0x368));
        }
    } else {
        s.d(0x390, s.d(0x378));
        if (s.i(0x408) < 180) s.d(0x398, s.d(0x370));
        else s.d(0x398, 2.0 * s.d(0x370));
    }
// ---- 惯量检测：电网电压相位突变检测（0x1BD0），用于惯量响应切换 ----
    if (s.i(0x57C) == 1 && s.i(0x10) == 1 && s.d(0x4E0) >= 0.15) {
        s.d(0xB88, s.d(0xB88) + 0.00025);
        if (s.d(0xB88) > 0.15) s.i(0x1BD0, 1);
        else s.i(0x1BD0, 0);
    } else if (s.i(0x57C) != 1) {
        s.i(0x1BD0, 0);
        s.d(0xB88, 0.0);
    }
    if (s.i(0x1BD0) == 1) {
        s.d(0x70, s.d(0x70) + 0.00025);
        if (s.d(0x70) > 0.02) s.d(0x90, 0.0);
        else if (s.d(0x70) > 0.01) s.d(0x90, 1.0);
    } else {
        s.d(0x98, 0.0);
        s.d(0x1C50, 0.0);
        s.d(0x1C58, 0.0);
        s.d(0x78, 0.0);
        s.d(0x70, 0.0);
    }
    if (s.d(0x90) == 1.0 && s.i(0x1BD0) == 1) {
        s.d(0x1C50, s.d(0x1C50) + s.d(0x1C40));
        s.d(0x1C58, s.d(0x1C58) + s.d(0xD0) / s.d(0x410));
        s.d(0x98, s.d(0x98) + s.d(0x88));
        s.d(0x78, s.d(0x78) + 1.0);
        s.d(0x80, s.d(0x98) / s.d(0x78));
        s.d(0x1C48, s.d(0x1C50) / s.d(0x78));
        s.d(0x1C60, s.d(0x1C58) / s.d(0x78));
    }
// ---- 特殊控制模式（0xB60 选 1/2/3）：把给定切换到平均值/额定值/记录值 ----
    if (s.i(0x10) != 1 && s.i(0x5D0) == 1 && s.i(0x570) == 1 && s.i(0x57C) == 1 && s.i(0x1BD0) == 1) {
        if (s.i(0x56C) < 80) {
            if (s.d(0xB60) == 1.0) {
                s.d(0x88, s.d(0x80));
                s.d(0xB0, s.d(0x88));
                s.increment(0x56C);
            } else if (s.d(0xB60) == 2.0) {
                s.d(0x88, s.d(0xB48) * s.d(0x548));
                s.d(0xB0, s.d(0x88));
                s.increment(0x56C);
            } else if (s.d(0xB60) == 3.0) {
                s.d(0x88, s.d(0x1C48) * s.d(0x548));
                s.d(0xB0, s.d(0x88));
                s.increment(0x56C);
            }
        } else {
            if (s.d(0xB60) == 1.0) {
                s.d(0x390, ((-((s.d(0x2A0) - s.d(0x80) / s.d(0x548)) / s.d(0xB58))) * s.d(0x548)) * 0.00025);
            } else if (s.d(0xB60) == 2.0) {
                s.d(0x390, ((-((s.d(0x2A0) - s.d(0xB48)) / s.d(0xB58))) * s.d(0x548)) * 0.00025);
            } else if (s.d(0xB60) == 3.0) {
                s.d(0x390, ((-((s.d(0x2A0) - s.d(0x1C48)) / s.d(0xB58))) * s.d(0x548)) * 0.00025);
            }
            s.d(0x80, 0.0);
            s.d(0x90, 0.0);
            s.d(0x1C48, 0.0);
            s.d(0x1C60, 0.0);
            if (s.d(0xB0) == s.d(0x88)) {
                for (std::size_t offset : {0x570, 0x56C, 0x1BD0, 0x774, 0x778, 0x77C, 0x780, 0x784, 0x788, 0x78C, 0x790}) s.i(offset, 0);
            }
            s.d(0x758, 0.0);
            for (std::size_t offset : {0x760, 0x764, 0x768, 0x76C, 0x770}) s.i(offset, 0);
        }
    } else s.i(0x570, 0);
    if (s.i(0x570) == 1) {
        s.d(0x580, s.d(0x580) + 0.00025);
        if (s.d(0x580) > s.d(0xB58)) {
            s.i(0x1BD0, 0);
            s.d(0x1C10, 0.0);
            s.d(0x580, 0.0);
            s.i(0x56C, 0);
        }
    }
    s.i(0x574, s.i(0x570));
// ---- 斜坡执行：0xA8（d 轴给定）、0xB0（q 轴给定）、0x658/0x660 缓缓逼近目标 ----
    slew(0xA8, 0xA0, 0x398);
    slew(0xB0, 0x88, 0x390);
    slew(0x658, 0x648, 0x390);
    slew(0x660, 0x650, 0x398);
// ---- 谐波提取：0x730/0x738 低通，0x1630/0x1668 残差，过递归滤波后限幅 ----
    s.d(0x730, 0.01 * s.d(0x1060) + 0.99 * s.d(0x730));
    s.d(0x1630, s.d(0x1060) - s.d(0x730));
    s.apply(0x1618, recursive_filter);
    s.d(0x1640, s.d(0x740) * s.d(0x1640));
    s.d(0x738, 0.01 * s.d(0x1068) + 0.99 * s.d(0x738));
    s.d(0x1668, s.d(0x1068) - s.d(0x738));
    s.apply(0x1650, recursive_filter);
    s.d(0x1678, s.d(0x748) * s.d(0x1678));
    if (s.d(0x1640) > std::fabs(s.d(0x88) * 0.05)) s.d(0x1640, std::fabs(s.d(0x88) * 0.05));
    else if (-std::fabs(s.d(0x88) * 0.05) > s.d(0x1640)) s.d(0x1640, -std::fabs(s.d(0x88) * 0.05));
    if (s.d(0x1678) > std::fabs(s.d(0x88) * 0.05)) s.d(0x1678, std::fabs(s.d(0x88) * 0.05));
    else if (-std::fabs(s.d(0x88) * 0.05) > s.d(0x1678)) s.d(0x1678, -std::fabs(s.d(0x88) * 0.05));
    if (s.i(0x750) == 0 || s.i(0x10) != 0) {
        s.d(0x1640, 0.0);
        s.d(0x1678, 0.0);
    }
// ---- 第 7 步：四个电流内环 PI：
//      0xDF8：q 轴（给定 0xB0+谐波补偿、反馈 0x2C0）；0xE78：d 轴（0xA8/0x2C8）；
//      0xEF8/0xF78：负序 d/q（0x658/0x660，过压标志 0x3F8 时清零）----
    s.d(0xDF8, s.d(0xB0) + s.d(0x1678));
    s.d(0xE00, s.d(0x2C0));
    s.apply(0xDF8, regulator);
    s.d(0xE78, s.d(0xA8) + s.d(0x1640));
    s.d(0xE80, s.d(0x2C8));
    s.apply(0xE78, regulator);
    if (s.i(0x3F8) == 0) {
        s.d(0xEF8, s.d(0x658));
        s.d(0xF00, s.d(0x2D0));
        s.apply(0xEF8, regulator);
        s.d(0xF78, s.d(0x660));
        s.d(0xF80, s.d(0x2D8));
        s.apply(0xF78, regulator);
    } else {
        s.d(0xF20, 0.0);
        s.d(0xFA0, 0.0);
        s.d(0xF48, 0.0);
        s.d(0xFC8, 0.0);
    }
// ---- 前馈与负序补偿（0x670 起，含 0x577 限幅），
//      以及残差滤波器的谐波抑制量（0x698/0x6A0/0x6A8/0x6B0）----
    if (s.i(0x66C) == 1) {
        s.d(0x670, (((s.d(0x1040) * s.d(0x5F0)) * (s.d(0x5F8) / 50.0)) / s.d(0xB8)) * 0.577);
        s.d(0x678, (((s.d(0x1048) * s.d(0x5F0)) * (s.d(0x5F8) / 50.0)) / s.d(0xB8)) * 0.577);
        if (s.d(0x670) > 0.577) s.d(0x670, 0.577);
        // The original lower-limit store targets 0x678 and leaves 0x670 below its bound.
        else if (-0.577 > s.d(0x670)) s.d(0x678, -0.577);
        clip(0x678, -0.577, 0.577);
        s.d(0x680, ((((-s.d(0x1050)) * s.d(0x5F0)) * (2.0 + s.d(0x5F8) / 50.0)) / s.d(0xB8)) * 0.577);
        s.d(0x688, ((((-s.d(0x1058)) * s.d(0x5F0)) * (2.0 + s.d(0x5F8) / 50.0)) / s.d(0xB8)) * 0.577);
        clip(0x680, -0.1, 0.1);
        clip(0x688, -0.1, 0.1);
    } else if (s.i(0x3F8) == 1) {
        s.d(0x670, 0.0);
        s.d(0x678, 0.0);
        s.d(0x680, ((((-s.d(0x1050)) * s.d(0x5F0)) * (2.0 + s.d(0x5F8) / 50.0)) / s.d(0xB8)) * 0.577);
        s.d(0x688, ((((-s.d(0x1058)) * s.d(0x5F0)) * (2.0 + s.d(0x5F8) / 50.0)) / s.d(0xB8)) * 0.577);
        clip(0x680, -0.1, 0.1);
        clip(0x688, -0.1, 0.1);
    } else {
        s.d(0x670, 0.0);
        s.d(0x678, 0.0);
        s.d(0x680, 0.0);
        s.d(0x688, 0.0);
    }
    if (s.i(0x690) == 0 || s.i(0x10) != 0) {
        s.d(0x698, 0.0);
        s.d(0x6A0, 0.0);
        s.d(0x6A8, 0.0);
        s.d(0x6B0, 0.0);
    } else if (s.i(0x690) == 1) {
        s.d(0x1470, s.d(0x1318));
        s.d(0x1478, s.d(0x1320));
        s.apply(0x1470, residual_filter);
        s.d(0x698, (((1.732 * s.d(0x1490)) * s.d(0x6B8)) / s.d(0xB8)) * 0.577);
        s.d(0x6A0, (((1.732 * s.d(0x1498)) * s.d(0x6B8)) / s.d(0xB8)) * 0.577);
        clip(0x698, -0.2, 0.2);
        clip(0x6A0, -0.2, 0.2);
        s.d(0x6A8, 0.0);
        s.d(0x6B0, 0.0);
    } else if (s.i(0x690) == 2) {
        s.d(0x1438, s.d(0x2C0));
        s.d(0x1440, s.d(0x2C8));
        s.apply(0x1438, residual_filter);
        s.d(0x6A8, (((1.732 * s.d(0x1458)) * s.d(0x6C0)) / s.d(0xB8)) * 0.577);
        s.d(0x6B0, (((1.732 * s.d(0x1460)) * s.d(0x6C0)) / s.d(0xB8)) * 0.577);
        clip(0x6A8, -0.2, 0.2);
        clip(0x6B0, -0.2, 0.2);
        s.d(0x698, 0.0);
        s.d(0x6A0, 0.0);
    }
// ---- 300Hz 谐振滤波（0x14A8/0x1560 两组） ----
    s.d(0x14C0, 0.0);
    s.d(0x14C8, s.d(0x1260));
    resonant_filter(s.at(0x14A8));
    s.d(0x1578, 0.0);
    s.d(0x1580, s.d(0x1268));
    resonant_filter(s.at(0x1560));
// ---- 调制波合成：正序/负序电压参考 + 各补偿量 → 反 Park（0x13E8/0x1410） ----
    s.d(0x1400, (((s.d(0x670) + s.d(0xE48) / s.d(0x358)) + s.d(0x6A0)) - s.d(0x6B0)) + s.d(0x1510) / s.d(0x358));
    s.d(0x1408, (((s.d(0x678) + s.d(0xEC8) / s.d(0x358)) - s.d(0x698)) + s.d(0x6A8)) + s.d(0x15C8) / s.d(0x358));
    s.d(0x13F8, s.d(0x2E0));
    s.apply(0x13E8, inverse_park_transform);
    s.d(0x1428, s.d(0x680) + s.d(0xF48) / s.d(0x358));
    s.d(0x1430, s.d(0x688) + s.d(0xFC8) / s.d(0x358));
    s.d(0x1420, s.d(0x2E8));
    s.apply(0x1410, inverse_park_transform);
// ---- crowbar（撬棒）逻辑：转速高于 0.2 或低于 -2.0（标幺）时
//      切换 0xD98/0xD10（撬棒投入状态与延时 5 秒） ----
    if (s.d(0x5F8) > 0.2) {
        s.i(0xD98, 1);
        s.d(0xD10, 5.0);
    } else if (-2.0 > s.d(0x5F8)) {
        s.i(0xD98, 0);
        s.d(0xD10, 5.0);
    } else s.d(0xD10, 0.0);
    s.d(0xDB0, s.d(0x228));
    s.d(0xDB8, s.d(0x230));
    s.d(0xDC0, s.d(0x238));
    s.d(0xCB8, s.d(0x280));
    s.d(0xCC0, s.d(0x288));
    s.d(0xBF0, s.d(0x13E8) + s.d(0x1410));
    s.d(0xBF8, s.d(0x13F0) + s.d(0x1418));
// ---- 第 8 步：调制（0x6E8==0/2 用普通版，==1 用备用版）→ 占空比输出 0x308 起 ----
    if (s.i(0x6E8) == 0 || s.i(0x6E8) == 2) s.apply(0xBF0, modulation);
    else alternate_modulation(s.at(0xBF0));
    s.d(0x308, s.d(0xCE0));
    s.d(0x310, s.d(0xCE8));
    s.d(0x320, s.d(0xCF0));
    s.d(0x328, s.d(0xCF8));
    s.d(0x338, s.d(0xD00));
    s.d(0x340, s.d(0xD08));
// ---- 过流保护：任一相转子电流超过 ±0x40 门槛 → 置 0x24 闭锁标志 ----
    if (std::fabs(s.d(0x228)) > s.d(0x40) || std::fabs(s.d(0x230)) > s.d(0x40) || std::fabs(s.d(0x238)) > s.d(0x40)) {
        s.i(0x24, 1);
        s.i(0x34, 0);
    }
    if (s.d(0x278) > s.d(0x28)) {
        s.i(0x24, 1);
        s.i(0x34, 0);
    }
// ---- 闭锁执行：0x24 == 1 时禁止输出（不擦除刚算出的占空比），
//      否则置运行标志 0x20、输出使能 0x350 = 1 ----
    if (s.i(0x24) == 1) {
        // Protection disables the output without erasing the duties just copied.
        s.d(0x350, 0.0);
        s.d(0xB0, 0.0);
        s.d(0xA8, 0.0);
        s.d(0xE20, 0.0);
        s.d(0xEA0, 0.0);
        s.d(0xF20, 0.0);
        s.d(0xFA0, 0.0);
        s.i(0x24, 0);
    } else {
        s.i(0x20, 1);
        s.d(0x350, 1.0);
        s.i(0x3A4, 1);
    }
// ---- exit_reset：停机过程（0x4FC==2）持续 0x504 步后整块复位 ----
exit_reset:
    if (s.i(0x4FC) == 2) {
        s.increment(0x508);
        if (s.i(0x508) > s.i(0x504)) {
            std::memset(state, 0, 0x1C68);
            motor_defaults(base_frequency, state);
        }
    }
// ---- final_status：0x528 = 运行标志（1/2 → 1.0），供输出与监控 ----
final_status:
    s.d(0x528, s.i(0x4FC) == 1 || s.i(0x4FC) == 2 ? 1.0 : 0.0);
}
