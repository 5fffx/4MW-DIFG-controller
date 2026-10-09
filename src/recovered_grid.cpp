// ============================================================================
// 文件：recovered_grid.cpp —— 网侧变流器控制（GSC，Grid Side Converter）
// ----------------------------------------------------------------------------
// 双馈风电机组有两台背靠背变流器：
//   - 网侧变流器（本文件）：接在电网和直流母线之间，任务是把转子回路
//     滑差功率送回电网，同时维持直流母线电压稳定；
//   - 机侧变流器（recovered_motor.cpp）：接在转子绕组和直流母线之间，
//     任务是根据转速和功率给定控制转子电流。
// 本文件三个函数：
//   grid_post              ：状态机 + 保护（每控制周期最后执行）；
//   grid_frequency_support ：惯量响应 + 一次调频（频率支撑）；
//   grid_calculate         ：主控制算法（PLL、序分解、外环、电流环、调制）。
// ============================================================================
#include "recovered_primitives.h"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <emmintrin.h>
#include <initializer_list>

namespace {

// ----------------------------------------------------------------------------
// GridState：帮助类，与前面文件里的一样，d(offset)/i(offset) 按字节偏移量
// 读写状态区；apply() 把状态区里某一段内存当成一个结构体（如滤波器）
// 传给基本算法单元计算；enabled(mask)/trip(bit) 操作保护使能位与跳闸位。
// ----------------------------------------------------------------------------
class GridState {
public:
    explicit GridState(void* state) : bytes_(static_cast<unsigned char*>(state)) {}

    double d(std::size_t offset) const { return read<double>(offset); }
    std::int32_t i(std::size_t offset) const { return read<std::int32_t>(offset); }
    std::uint32_t u(std::size_t offset) const { return read<std::uint32_t>(offset); }
    void d(std::size_t offset, double value) { write(offset, value); }
    void i(std::size_t offset, std::int32_t value) { write(offset, value); }
    void increment(std::size_t offset) { write(offset, u(offset) + std::uint32_t{1}); }
    void* at(std::size_t offset) { return bytes_ + offset; }

    template <typename T>
    void apply(std::size_t offset, void (*function)(T*)) {
        T value = read<T>(offset);
        function(&value);
        write(offset, value);
    }

    bool enabled(std::uint32_t mask) const {
        return (u(0x7E0) & mask) != 0 || (u(0x7D8) & mask) != 0;
    }

    void trip(unsigned bit) {
        write(0x7E8, u(0x7E8) | (std::uint32_t{1} << bit));
        write(0x7DC, u(0x7DC) | (std::uint32_t{1} << bit));
        write(0x7F0, u(0x7F0) | std::uint32_t{1});
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
    return _mm_cvttsd_si32(_mm_set_sd(value));
}

struct RmsRow {
    std::size_t source;
    std::size_t accumulator;
    std::size_t output;
};

struct VoltageStage {
    std::size_t limit;
    std::size_t timer_a;
    std::size_t timer_b;
    std::size_t duration;
};

struct FrequencyStage {
    std::size_t limit;
    std::size_t timer;
    std::size_t duration;
};

} // namespace

// ============================================================================
// ----------------------------------------------------------------------------
// grid_post —— 网侧状态机 + 保护（每个控制周期最后执行）
// ----------------------------------------------------------------------------
// 一、运行状态机（0x350）：0 = 停机，1 = 运行，2 = 停机过程中，3 = 跳闸。
//     迁移条件：跳闸命令 → 3；停机命令 → 2；启动命令 → 1。
//     跳闸时：清所有输出（PWM 使能、调制波），只保留直流电压极性标志。
//     复位命令（0x520 == 1）：整块网侧状态区清零并重写默认值 = 重启。
//
// 二、测量：
//   - 频率测量：累加角频率 0x11A0，每满一个周期算一次 f（0x468），
//     并限幅在 f±10Hz 内（0x47C 置“本轮已更新”标志）；
//   - RMS 计算：10 路信号（三相电压、三相电流、折算后的量）平方累加，
//     累加满一个窗口（0x464 步）后开方，得到各量有效值。
//
// 三、保护（都在 0x7E4/0x7D8 的保护使能位与 0x2C 运行标志满足时动作）：
//   - 过压快速段：电压越限就累积时间，超过定值 0x918 跳闸 trip(10)；
//   - 过压三段式：0x920/0x940/0x960 三组“阈值+延时”，超时 trip(10)；
//   - 欠压四段式：0x980/0x9A0/0x9C0/0x9E0 四组，超时 trip(11)；
//   - 过频三段：0xA78/0xA90/0xAA8 起，超时 trip(14)；
//   - 欠频三段：0xAC0/0xAD8/0xAF0 起，超时 trip(15)。
//   trip() 把跳闸位写进 0x7E8/0x7DC，并在 0x7F0 置锁存位。
//
// 四、故障锁存：一类故障立即把状态机打回 3；二类故障延时 4 秒后打回 3。
// ----------------------------------------------------------------------------
// ============================================================================
extern "C" __declspec(dllexport) void grid_post(double f, void* state) {
    GridState s(state);
    if (s.i(0x34C) == 1) {
        s.i(0x350, 3);
    } else if (s.i(0x350) == 1 && s.i(0x348) == 0) {
        s.i(0x350, 2);
    } else if (s.i(0x350) == 0 && s.i(0x348) == 1) {
        s.i(0x350, 1);
    }

    if (s.i(0x350) == 3) {
        for (std::size_t offset : {0x18, 0x50C, 0x2C, 0x510, 0x514}) {
            s.i(offset, 0);
        }
        for (std::size_t offset : {0x4B8, 0x4D0, 0x4E8, 0x4A8, 0x4B0, 0xD88, 0xD90, 0xD98}) {
            s.d(offset, 0.0);
        }
        s.i(0x504, s.d(0x100) > 0.0 ? 1 : 0);
        s.i(0x32C, s.i(0x504));
        s.i(0x508, s.d(0x108) > 0.0 ? 1 : 0);
        s.i(0x330, s.i(0x508));
    }

    if (s.i(0x520) == 1) {
        std::memset(state, 0, 0x1AB8);
        grid_defaults(f, state);
        return;
    }

// ---- 频率测量：每步累加 PLL 给出的角频率，满一个周期算一次平均值 ----
    s.d(0x470, s.d(0x470) + s.d(0x11A0));
    s.increment(0x478);
    if (s.i(0x478) >= s.i(0x3A0)) {
        s.d(0x468, ((s.d(0x470) / static_cast<double>(s.i(0x478))) / 2.0) / 3.14159265358979);
        s.i(0x464, cvtt_i32(4000.0 / s.d(0x468)));
        s.i(0x478, 0);
        s.d(0x470, 0.0);
        s.i(0x47C, 1);
        if (s.i(0x464) < cvtt_i32(4000.0 / (f + 10.0))) {
            s.i(0x464, cvtt_i32(4000.0 / (f + 10.0)));
        } else if (s.i(0x464) > cvtt_i32(4000.0 / std::fabs(f - 10.0))) {
            s.i(0x464, cvtt_i32(4000.0 / std::fabs(f - 10.0)));
        }
    }

    s.d(0x3B0, s.d(0x3B0) + s.d(0x11D0));
    s.increment(0x3B8);
    if (s.i(0x3B8) >= s.i(0x464)) {
        s.d(0x3A8, (s.d(0x3B0) / static_cast<double>(s.i(0x3B8))) / 1.414);
        s.i(0x3B8, 0);
        s.d(0x3B0, 0.0);
    }

// ---- RMS 有效值计算：10 路信号的平方累加 → 窗口满后开方 ----
// 每路信号一个“源地址 → 累加器 → 结果”三元组，避免手写 10 遍同样的代码。
    const RmsRow rms_rows[] = {
        {0xC0, 0x3F0, 0x3C0}, {0xC8, 0x400, 0x3D0}, {0xD0, 0x410, 0x3E0},
        {0xD8, 0x438, 0x420}, {0xE0, 0x440, 0x428}, {0xE8, 0x448, 0x430},
        {0xF0, 0x458, 0x450}, {0x1718, 0x3F8, 0x3C8},
        {0x1720, 0x408, 0x3D8}, {0x1728, 0x418, 0x3E8}
    };
    s.increment(0x460);
    for (const RmsRow& row : rms_rows) {
        const double source = s.d(row.source);
        s.d(row.accumulator, s.d(row.accumulator) + source * source);
    }
    if (s.i(0x460) >= s.i(0x464)) {
        for (const RmsRow& row : rms_rows) {
            s.d(row.output, std::sqrt(s.d(row.accumulator) / static_cast<double>(s.i(0x460))));
        }
        for (const RmsRow& row : rms_rows) {
            s.d(row.accumulator, 0.0);
        }
        s.i(0x460, 0);
        s.i(0x480, 1);
    }

// ---- 保护段：以下四组分别对应过压快速段、过压三段、欠压四段、过频/欠频三段 ----
// 0x3C0/0x3D0 是 A/B 相电压有效值；0x900、0x920…… 是定值（来自保护输入），
// 0x928/0x930 等是各段的“累积计时器”和“允许时长”。
    if ((s.u(0x7E4) & 0x8000) != 0 || (s.u(0x7D8) & 0x40000000) != 0) {
        if (s.enabled(0x400)) {
            s.d(0x908, s.d(0x3C0) > s.d(0x900) ? s.d(0x908) + 0.00025 : 0.0);
            s.d(0x910, s.d(0x3D0) > s.d(0x900) ? s.d(0x910) + 0.00025 : 0.0);
            if (s.d(0x908) > s.d(0x918)) {
                s.trip(10);
            }
            if (s.d(0x910) > s.d(0x918)) {
                s.trip(10);
            }
        }

        if (s.enabled(0x400) && s.i(0x2C) == 1 && s.i(0x480) == 1) {
            const VoltageStage stages[] = {
                {0x920, 0x928, 0x930, 0x938},
                {0x940, 0x948, 0x950, 0x958},
                {0x960, 0x968, 0x970, 0x978}
            };
            for (const VoltageStage& stage : stages) {
                s.d(stage.timer_a, s.d(0x3C0) > s.d(stage.limit) ? s.d(stage.timer_a) + 0.00025 : 0.0);
                s.d(stage.timer_b, s.d(0x3D0) > s.d(stage.limit) ? s.d(stage.timer_b) + 0.00025 : 0.0);
            }
            for (const VoltageStage& stage : stages) {
                if (s.d(stage.timer_a) >= s.d(stage.duration) || s.d(stage.timer_b) >= s.d(stage.duration)) {
                    s.d(stage.timer_a, 0.0);
                    s.d(stage.timer_b, 0.0);
                    s.trip(10);
                }
            }
        }

        if (s.enabled(0x800) && s.i(0x2C) == 1 && s.i(0x480) == 1) {
            const VoltageStage stages[] = {
                {0x980, 0x988, 0x990, 0x998},
                {0x9A0, 0x9A8, 0x9B0, 0x9B8},
                {0x9C0, 0x9C8, 0x9D0, 0x9D8},
                {0x9E0, 0x9E8, 0x9F0, 0x9F8}
            };
            for (const VoltageStage& stage : stages) {
                s.d(stage.timer_a, s.d(stage.limit) > s.d(0x3C0) ? s.d(stage.timer_a) + 0.00025 : 0.0);
                s.d(stage.timer_b, s.d(stage.limit) > s.d(0x3D0) ? s.d(stage.timer_b) + 0.00025 : 0.0);
            }
            for (const VoltageStage& stage : stages) {
                if (s.d(stage.timer_a) > s.d(stage.duration) || s.d(stage.timer_b) > s.d(stage.duration)) {
                    s.d(stage.timer_a, 0.0);
                    s.d(stage.timer_b, 0.0);
                    s.trip(11);
                }
            }
        }

        if (s.i(0x2C) == 1 && s.i(0x47C) == 1 && s.enabled(0x4000)) {
            const FrequencyStage stages[] = {
                {0xA78, 0xA80, 0xA88}, {0xA90, 0xA98, 0xAA0}, {0xAA8, 0xAB0, 0xAB8}
            };
            for (const FrequencyStage& stage : stages) {
                if (s.d(0x468) > s.d(stage.limit)) {
                    s.d(stage.timer, s.d(stage.timer) + 0.00025);
                    if (s.d(stage.timer) > s.d(stage.duration)) {
                        s.d(stage.timer, 0.0);
                        s.trip(14);
                    }
                } else {
                    s.d(stage.timer, 0.0);
                }
            }
        }

        if (s.i(0x2C) == 1 && s.i(0x47C) == 1 && s.enabled(0x8000)) {
            const FrequencyStage stages[] = {
                {0xAC0, 0xAC8, 0xAD0}, {0xAD8, 0xAE0, 0xAE8}, {0xAF0, 0xAF8, 0xB00}
            };
            for (const FrequencyStage& stage : stages) {
                if (s.d(stage.limit) > s.d(0x468)) {
                    s.d(stage.timer, s.d(stage.timer) + 0.00025);
                    if (s.d(stage.timer) > s.d(stage.duration)) {
                        s.d(stage.timer, 0.0);
                        s.trip(15);
                    }
                } else {
                    s.d(stage.timer, 0.0);
                }
            }
        }
    }

// ---- 故障锁存：置位后本拍不立刻清输出，下一拍状态机转 3 才清 ----
    // A fault latched here clears outputs on the next call, not this one.
    if ((s.u(0x7F0) & 1) != 0) {
        s.i(0x350, 3);
    } else if ((s.u(0x7F0) & 2) != 0) {
        if (s.i(0x2C) == 1) {
            s.d(0x7F8, s.d(0x7F8) + 0.00025);
            if (s.d(0x7F8) >= 4.0) {
                s.i(0x350, 3);
                s.d(0x7F8, 4.0);
            }
        } else {
            s.i(0x350, 3);
        }
    } else {
        s.d(0x7F8, 0.0);
    }
}

// ============================================================================
// ----------------------------------------------------------------------------
// grid_frequency_support —— 惯量响应 + 一次调频（频率支撑）
// ----------------------------------------------------------------------------
// 先做测量：
//   5 点滑动平均（0x670 起的环形缓冲）得到平稳的频率 0x6E0；
//   用相邻两次平均值的差分近似算 df/dt（0x6D0）。
//
// 模式 1（惯量响应，commands[10]==1）：模拟同步机的转动惯量——
//   df/dt 或频偏超过死区后，输出功率指令 0x750 = ±幅值（方向由频偏符号定），
//   维持一段时间（0xBC8 秒）后退出。
//
// 模式 2（一次调频，commands[10]==2）：功频下垂——
//   频偏超过死区后，0x758 = -(f - fn) × 下垂系数，带功率上/下限
//   （0xC40/0xC48，标幺值换算），并带最小功率门槛（0.1 倍基准）。
//
// 最后：0x768 = 惯量功率 + 调频功率，限幅（不超过 1 - 有功给定），
// 由总调度叠加到有功给定上（0x1D58 += 0x768）。
// 0xC60 是给监控输出用的“当前支撑模式”标志。
// ----------------------------------------------------------------------------
// ============================================================================
extern "C" __declspec(dllexport) void grid_frequency_support(double f, void* state) {
    GridState s(state);
    if (s.i(0x478) == 1) {
        s.d(0x670 + 8 * static_cast<std::ptrdiff_t>(s.i(0x704)), s.d(0x468));
        s.increment(0x704);
        if (s.i(0x704) >= s.i(0x708)) {
            s.i(0x704, 0);
            if (s.i(0x70C) == 0) {
                s.i(0x70C, 1);
            }
        }
        if (s.i(0x70C) == 1) {
            s.d(0x6D8, 0.0);
            for (std::int32_t j = 0; j < s.i(0x708); ++j) {
                s.d(0x6D8, s.d(0x6D8) + s.d(0x670 + 8 * static_cast<std::ptrdiff_t>(j)));
            }
            s.d(0x6E0, s.d(0x6D8) / static_cast<double>(s.i(0x708)));
            s.d(0x6F0, ((s.d(0x6E0) - s.d(0x6E8)) / 0.00025) / static_cast<double>(s.i(0x464)));
            s.d(0x6E8, s.d(0x6E0));
            s.d(0x698 + 8 * static_cast<std::ptrdiff_t>(s.i(0x6F8)), s.d(0x6F0));
            s.increment(0x6F8);
            if (s.i(0x6F8) >= s.i(0x6FC)) {
                s.i(0x6F8, 0);
                if (s.i(0x700) == 0) {
                    s.i(0x700, 1);
                }
            }
            if (s.i(0x700) == 1) {
                s.d(0x6C0, 0.0);
                for (std::int32_t j = 0; j < s.i(0x6FC); ++j) {
                    s.d(0x6C0, s.d(0x6C0) + s.d(0x698 + 8 * static_cast<std::ptrdiff_t>(j)));
                }
                s.d(0x6C8, s.d(0x6C0) / static_cast<double>(s.i(0x6FC)));
                s.d(0x6D0, s.d(0x6C8));
            }
        }
    }

    const double measured = s.d(0x468) * static_cast<double>(s.i(0x70C) == 0)
                          + s.d(0x6E0) * static_cast<double>(s.i(0x70C) == 1);
    const double threshold_a = s.d(0xC38);
    const double threshold_b = s.d(0xC38);
    const double amplitude = s.d(0x720);
    const double derivative = s.d(0x6D0);
    const double product = (measured - f) * derivative;
    const double amplitude_scaled = amplitude * 1.0;
    if (s.i(0x740) == 1) {
        if (s.i(0xBD0) == 1) {
            s.d(0xBE0, s.d(0xBE0) + 0.00025);
            if (s.d(0xBE0) >= 2.0) {
                if (s.i(0x788) == 1) {
                    s.d(0xBC0, s.d(0xBC0) + 0.00025);
                    if (s.d(0xBC0) > s.d(0xBC8)) {
                        s.d(0xBC0, 0.0);
                        s.i(0x788, 0);
                        s.d(0x790, 0.0);
                    } else {
                        s.d(0x750, s.d(0x790));
                    }
                } else if (s.i(0x780) == 1) {
                    s.d(0x750, s.d(0x790));
                } else {
                    s.d(0x750, 0.0);
                    s.d(0x790, 0.0);
                }
                if ((std::fabs(derivative) > threshold_a && product > 0.0)
                    || std::fabs(measured - f) > s.d(0xC38)) {
                    if (0.005 > s.d(0x748)) {
                        s.d(0x748, s.d(0x748) + 0.00025);
                    } else {
                        s.i(0x780, 1);
                        if (measured - f > s.d(0xC38)) {
                            s.d(0x790, -amplitude_scaled);
                        } else if (-s.d(0xC38) > measured - f) {
                            s.d(0x790, amplitude_scaled);
                        }
                    }
                } else if (threshold_b > std::fabs(derivative)) {
                    s.d(0x748, 0.0);
                    s.i(0x780, 0);
                }
                if (s.i(0x780) == 1 && s.i(0x784) == 0) {
                    s.i(0x788, 1);
                }
                s.i(0x784, s.i(0x780));
                if (s.i(0x788) == 0) {
                    s.d(0x790, 0.0);
                    s.d(0x750, 0.0);
                }
            }
        }
    } else {
        s.d(0x750, 0.0);
        s.d(0x748, 0.0);
    }

    const double upper_frequency = (f + s.d(0x738)) + 0.005;
    const double lower_frequency = (f - s.d(0x738)) - 0.005;
    const double positive_max = (s.d(0xC40) * 1.05) * s.d(0x7C0);
    const double negative_max = ((-s.d(0xC40)) * 1.05) * s.d(0x7C0);
    const double positive_min = (s.d(0xC48) * 1.1) * s.d(0x7C0);
    const double negative_min = ((-s.d(0xC48)) * 1.1) * s.d(0x7C0);
    if (s.i(0x740) == 2) {
        if (s.i(0xBD0) == 1) {
            s.d(0xBE0, s.d(0xBE0) + 0.00025);
            if (s.d(0xBE0) >= 2.0) {
                if (lower_frequency > measured) {
                    if (0.05 > s.d(0xBF0)) {
                        s.d(0xBF0, s.d(0xBF0) + 0.00025);
                    } else {
                        s.i(0x798, 1);
                    }
                } else if (measured > upper_frequency) {
                    if (0.05 > s.d(0xBF8)) {
                        s.d(0xBF8, s.d(0xBF8) + 0.00025);
                    } else {
                        s.i(0x798, 2);
                    }
                } else {
                    s.i(0x798, 0);
                    s.i(0x7A0, 0);
                    s.d(0xBF0, 0.0);
                    s.d(0xBF8, 0.0);
                    s.d(0x758, 0.0);
                }
            }
            if (s.i(0x7A0) == 1) {
                s.d(0x7A8, s.d(0x7A8) + 0.00025);
                if (s.d(0x7A8) > s.d(0x7B0)) {
                    s.d(0x7A8, 0.0);
                    s.i(0x7A0, 0);
                    s.d(0x758, 0.0);
                } else {
                    s.d(0x758, (-(measured - f)) * (((1.0 / f) * s.d(0x7C0)) / s.d(0xC30)));
                }
            } else if (s.i(0x7A0) == 2) {
                s.d(0x7B8, s.d(0x7B8) + 0.00025);
                if (s.d(0x7B8) > s.d(0x7B0)) {
                    s.d(0x7B8, 0.0);
                    s.i(0x7A0, 0);
                    s.d(0x758, 0.0);
                } else {
                    s.d(0x758, (-(measured - f)) * (((1.0 / f) * s.d(0x7C0)) / s.d(0xC30)));
                }
            }
            if (s.i(0x798) == 1 && s.i(0x79C) == 0) {
                s.i(0x7A0, 1);
            }
            if (s.i(0x798) == 2 && s.i(0x79C) == 0) {
                s.i(0x7A0, 2);
            }
            s.i(0x79C, s.i(0x798));
            if (s.d(0x758) > positive_max && std::fabs(s.d(0x758)) >= 0.1 * s.d(0x7C0) && s.i(0x7A0) == 1) {
                s.d(0x758, positive_max);
            } else if (positive_min >= s.d(0x758) && s.d(0x758) >= 0.0 && s.i(0x7A0) == 1) {
                s.d(0x758, positive_min);
            } else if (s.d(0x758) >= negative_min && 0.0 >= s.d(0x758) && s.i(0x7A0) == 2) {
                s.d(0x758, negative_min);
            } else if (positive_max > s.d(0x758) && std::fabs(s.d(0x758)) >= 0.1 * s.d(0x7C0) && s.i(0x7A0) == 2) {
                s.d(0x758, negative_max);
            } else {
                s.d(0x758, 0.0);
            }
        }
    } else {
        s.d(0x758, 0.0);
        s.d(0xBD8, 0.0);
    }

    if (s.d(0x750) != 0.0 && s.i(0x740) == 1) {
        s.i(0xC60, 1);
    } else if (s.d(0x758) != 0.0 && s.i(0x740) == 2) {
        s.i(0xC60, 2);
    } else {
        s.i(0xC60, 0);
    }
    s.d(0x760, s.d(0x750) + s.d(0x758));
    if (1.0 >= s.d(0x760) + s.d(0x7C0)) {
        s.d(0x768, s.d(0x760));
    } else if (s.d(0x760) + s.d(0x7C0) > 1.0) {
        s.d(0x768, 1.0 - s.d(0x7C0));
    }
}

// ============================================================================
// ----------------------------------------------------------------------------
// grid_calculate —— 网侧变流器主控制（每个控制周期执行一次）
// ----------------------------------------------------------------------------
// 整体信号流（从上到下）：
//   1) 刷新 PI 参数（参数区 0x550/0x560 → 0x1320/0x1368），并按模式
//      （0x18：0=正常 1=LVRT 2=HVRT）选择电流环积分限幅 7.85 / 62.83 / 15.71；
//   2) 锁相环 PLL（状态区 0x1188 起）→ 网侧角度 0x11A8、频率 0x11A0；
//   3) 三相自适应振荡器 + 两态滤波器（每相一组，0x1890/0x1788 起）：
//      从三相电压里提取基波幅值，得到三相电压标幺 0x1730 起，
//      及其最大/最小值 0x1760/0x1768/0x1770/0x1778（供 LVRT/HVRT 判定）；
//   4) 电流正负序分解 current_sequence：正序 dq（0x90/0x98）、负序 dq（0xA0/0xA8）；
//   5) 直流母线电压平均（0x1B8，每 0x1D0 步更新一次）；
//   6) 启动逻辑：前 400 步走软启动斜坡（0x14D0 两态滤波 + 电压外环 PI 0x1108）；
//   7) LVRT/HVRT 模式判定（0x18 状态机，带 20 步确认计数和退出滞回）：
//      - 低压穿越：电压低于 0.9 标幺 → 模式 1，无功优先（按系数 3.0 发无功），
//        有功电流受 √(Imax² - Iq²) 限制；
//      - 高压穿越：电压高于 1.13 标幺 → 模式 2，吸收无功（系数 2.2841）；
//   8) 外环/内环：
//      - 模式 1/2：无功电流给定直接算（0x2D0），再合成 dq 电流给定；
//      - 模式 0：直流电压外环（0x1108 的 PI）→ 有功电流给定 + 交叉解耦；
//   9) 四个电流内环 PI（0xF08 正序 d、0xF88 正序 q、0x1008 负序 d、0x1088 负序 q），
//      加上 300Hz 谐振滤波器（0x15A8/0x1660）和差分滤波修正；
//  10) 电压前馈（0x1480/0x14A8 反 Park）叠加负序补偿 → 调制（0xD00）→
//      输出占空比到 0x4C0 起；0x31C 软启动斜坡控制 PWM 使能；
//  11) 停机时（0x350==2）持续 0x35C 步后整块复位。
// ----------------------------------------------------------------------------
// ============================================================================
extern "C" __declspec(dllexport) void grid_calculate(double f, void* state) {
    GridState s(state);
// 四个小工具 lambda：ramp = 斜坡（按固定步长逼近目标值，上升/下降步长可不同）；
// clip = 限幅；maximum/minimum = 三数取最大/最小。
    const auto ramp = [&s](std::size_t value, std::size_t target, std::size_t step, double divisor) {
        if (s.d(target) > s.d(value)) {
            s.d(value, s.d(value) + s.d(step));
            if (s.d(value) > s.d(target)) s.d(value, s.d(target));
        } else if (s.d(value) > s.d(target)) {
            s.d(value, s.d(value) - s.d(step) / divisor);
            if (s.d(target) > s.d(value)) s.d(value, s.d(target));
        }
    };
    const auto clip = [](double value, double lower, double upper) {
        if (value > upper) return upper;
        if (lower > value) return lower;
        return value;
    };
    const auto maximum = [](double a, double b, double c) {
        const double t = a > b ? a : b;
        return t > c ? (a > b ? a : b) : c;
    };
    const auto minimum = [](double a, double b, double c) {
        const double t = a > b ? b : a;
        return t > c ? c : (a > b ? b : a);
    };
// ---- 第 1 步：PI 参数刷新 + 按运行模式选积分限幅 ----
    s.d(0x1320, s.d(0x550));
    s.d(0x1368, s.d(0x560));
    if (s.i(0x18) == 1) {
        s.d(0x1348, 7.853981633974475);
        s.d(0x1350, -7.853981633974475);
    } else if (s.i(0x18) == 2) {
        s.d(0x1348, 62.8318530717958);
        s.d(0x1350, -62.8318530717958);
    } else {
        s.d(0x1348, 15.70796326794895);
        s.d(0x1350, -15.70796326794895);
    }
    s.d(0x1188, s.d(0xC0));
    s.d(0x1190, s.d(0xC8));
    s.d(0x1198, s.d(0xD0));
// ---- 第 2 步：锁相环（把 0x1188 起的一段状态区当 PllState 使用） ----
    {
        PllState p;
        std::memcpy(&p, s.at(0x1188), sizeof(p));
        pll(f, &p);
        std::memcpy(s.at(0x1188), &p, sizeof(p));
    }
    s.d(0x218, s.d(0x11A8));
    s.d(0x230, std::sqrt(s.d(0x11D0) * s.d(0x11D0) + s.d(0x11D8) * s.d(0x11D8)));
// ---- 第 3 步：三相自适应振荡器 + 两态滤波器，提取每相基波幅值 ----
// 0x1718 起是三相电压，每相先过自适应振荡器（0x1890/0x1948/0x1A00），
// 再过两态滤波器（0x1788/0x17E0/0x1838），结果除以线电压基准 0x1780 得标幺。
    for (int phase = 0; phase < 3; ++phase) {
        const std::size_t adaptive = 0x1890 + 0xB8 * phase;
        const std::size_t filter = 0x1788 + 0x58 * phase;
        s.d(adaptive + 0x40, s.d(0x1718 + 8 * phase));
        adaptive_oscillator(s.at(adaptive), f);
        s.d(filter + 0x10, s.d(adaptive + 0x68));
        s.d(filter, s.d(0x11A0) * 2.0);
        s.apply(filter, two_state_filter);
        s.d(0x1730 + 8 * phase, s.d(filter + 0x18));
        s.d(0x1748 + 8 * phase, s.d(0x1730 + 8 * phase) / s.d(0x1780));
    }
    s.d(0x1760, maximum(s.d(0x1730), s.d(0x1738), s.d(0x1740)));
    s.d(0x1768, minimum(s.d(0x1730), s.d(0x1738), s.d(0x1740)));
    s.d(0x1770, maximum(s.d(0x1748), s.d(0x1750), s.d(0x1758)));
    s.d(0x1778, minimum(s.d(0x1748), s.d(0x1750), s.d(0x1758)));
    if (s.i(0x228) == 1) s.d(0x218, s.d(0x11A8) + s.d(0x220));
    else s.d(0x218, s.d(0x11A8));
    if (s.d(0x218) > 6.28318530717958) s.d(0x218, s.d(0x218) - 6.28318530717958);
    else if (0.0 > s.d(0x218)) s.d(0x218, s.d(0x218) + 6.28318530717958);
// ---- 第 4 步：电流正负序分解（0x1388 起的状态区，θ 用 0x218 及其负值） ----
// 0x668 选择用哪一组滤波结果：0 = 短时间常数组（0x13B0 起），1 = 另一组（0x1410 起）。
    s.d(0x1388, s.d(0xD8));
    s.d(0x1390, s.d(0xE0));
    s.d(0x1398, s.d(0xE8));
    s.d(0x13A0, s.d(0x218));
    s.d(0x13A8, -s.d(0x218));
    current_sequence(s.at(0x1388));
    if (s.i(0x668) == 0) {
        s.d(0x90, s.d(0x13B0));
        s.d(0x98, s.d(0x13B8));
        s.d(0xA0, s.d(0x13C0));
        s.d(0xA8, s.d(0x13C8));
    } else {
        s.d(0x90, s.d(0x1410));
        s.d(0x98, s.d(0x1418));
        s.d(0xA0, s.d(0x1420));
        s.d(0xA8, s.d(0x1428));
    }
    s.d(0xF8, s.d(0x100) + s.d(0x108));
// ---- 第 5 步：直流母线电压平均（0xF8 = 正负母线电压和，每 0x1D0 步平均） ----
    s.d(0x1C8, s.d(0x1C8) + 1.0);
    s.d(0x1C0, s.d(0x1C0) + s.d(0xF8));
    if (s.d(0x1C8) >= s.d(0x1D0)) {
        s.d(0x1B8, s.d(0x1C0) / s.d(0x1D0));
        s.d(0x1C8, 0.0);
        s.d(0x1C0, 0.0);
    }
// ---- 闸门：只有“运行/停机过程”状态才执行下面的控制，否则直接跳到最后 ----
    if (s.i(0x350) != 1 && s.i(0x350) != 2) goto final_flags;

    if (s.d(0xC98) == 1.0) {
        if (0.2 > s.d(0x7C8)) s.d(0x7C8, s.d(0x7C8) + 0.00025);
        else {
            s.d(0x7D0, 1.0);
            s.d(0x7C8, 0.0);
        }
    }
// ---- 启动计数：0x34 累计运行步数；前 0x008 步（800 步 = 0.2 秒）内 0x228=1，
//      锁相角加 0x220 修正（启动阶段的角度补偿）；0x190 累计锁相频率误差 ----
    if (s.i(0x2C) == 1) {
        if (s.i(0x34) < 4000) s.increment(0x34);
        if (s.i(0x328) < s.i(0x8)) {
            s.increment(0x328);
            s.i(0x228, 1);
        } else {
            s.increment(0x198);
            s.d(0x190, s.d(0x190) + (s.d(0x11B0) - s.d(0x1A8)) * 1.5);
            if (s.i(0x198) >= s.i(0x19C)) {
                s.d(0x180, s.d(0x190) / static_cast<double>(s.i(0x19C)));
                if (0.0 > s.d(0x180)) s.d(0x180, 0.0);
                s.i(0x198, 0);
                s.d(0x190, 0.0);
            }
        }
    }
    if (s.d(0x180) > s.d(0x188)) {
        s.d(0x188, s.d(0x180));
        s.i(0x1B0, 0);
    } else {
        s.increment(0x1B0);
        if (s.i(0x1B0) >= s.i(0x1A0)) {
            s.i(0x1B0, 0);
            s.d(0x188, s.d(0x180));
        }
    }
    if (s.i(0x18) != 0 && s.d(0xF8) > s.d(0x130)) {
        if (s.d(0xF8) > s.d(0x170)) s.d(0x130, s.d(0x170));
        else s.d(0x130, s.d(0xF8));
        if (s.d(0x130) - s.d(0x168) > s.d(0x188)) s.d(0x188, s.d(0x130) - s.d(0x168));
    }
// ---- 第 6 步：软启动。前 400 步：直流电压给定从当前值缓慢爬升
//      （0x14E0 = 当前直流电压 → 两态滤波 → 电压外环 PI 0x1108），
//      并让输出调制波为零，避免冲击。400 步后转入正常控制 ----
    if (s.i(0) < s.i(4)) {
        s.increment(0);
        s.i(0x510, 1);
        s.i(0x514, 0);
        s.d(0x130, s.d(0xF8));
        s.d(0x14E0, s.d(0xF8));
        s.apply(0x14D0, two_state_filter);
        goto exit_tail;
    }
    s.i(0x510, 0);
    s.i(0x514, 1);
    if (static_cast<double>(s.i(4)) * 1.5 > static_cast<double>(s.i(0))) s.increment(0);
    else s.i(0x30, 1);
// ---- 直流电压水平判定：根据当前直流电压落在哪个区间，
//      选择母线过压/欠压的限幅基准（0x208/0x210）与标志 0x504/0x508 ----
    if (s.d(0x128) > s.d(0x1F8)) {
        s.d(0x208, s.d(0x1E8));
        s.d(0x210, s.d(0x1F0));
    } else if (s.d(0x200) > s.d(0x130)) {
        s.d(0x208, s.d(0x1D8));
        s.d(0x210, s.d(0x1E0));
    }
    if (s.d(0x100) > s.d(0x208)) {
        s.i(0x504, 1);
        s.i(0x32C, 1);
    } else if (s.d(0x210) > s.d(0x100)) {
        s.i(0x504, 0);
        s.i(0x32C, 0);
    }
    if (s.d(0x108) > s.d(0x208)) {
        s.i(0x508, 1);
        s.i(0x330, 1);
    } else if (s.d(0x210) > s.d(0x108)) {
        s.i(0x508, 0);
        s.i(0x330, 0);
    }
// ---- 电压外环：0x138 是有功电流给定，向目标 0x140 斜坡（步长 0x150，
//      上升限速、下降快 64 倍）；0x1148 是限幅（母线正常时 1.25×基准） ----
    if (s.i(0x158) == 1) {
        ramp(0x138, 0x140, 0x150, 1.0);
        s.d(0x1150, -s.d(0x138));
        if (s.i(0x32C) == 1 || s.i(0x330) == 1) s.d(0x1148, 0.0);
        else s.d(0x1148, 1.25 * s.d(0x110));
    }
    s.d(0x128, s.d(0x168) + s.d(0x188));
    if (s.d(0x128) > s.d(0x170)) s.d(0x128, s.d(0x170));
    ramp(0x130, 0x128, 0x178, 64.0);
// ---- 电压外环 PI（0x1108）：输入 0x130（直流电压给定）与 0x14E8（滤波后实际值） ----
    s.increment(0x15C);
    s.d(0x14E0, s.d(0xF8));
    s.apply(0x14D0, two_state_filter);
    if (s.i(0x15C) >= s.i(0x160)) {
        s.i(0x15C, 0);
        s.d(0x1108, s.d(0x130));
        s.d(0x1110, s.d(0x14E8));
        s.apply(0x1108, regulator);
    }
// ---- 电流给定合成：0x58 = 外环输出（d 轴给定），
//      0x588 对外环输出做平滑（8/1024 的低通）----
    s.d(0x58, s.d(0x1158));
    s.d(0x588, (8.0 * s.d(0x58) + 1016.0 * s.d(0x588)) / 1024.0);
    s.d(0x78, s.i(0x1C) == 0 ? s.d(0x68) : s.d(0x70));
    ramp(0x60, 0x50, 0x78, 1.0);

// ------------------------------------------------------------------
// ---- 第 7 步：LVRT/HVRT 模式判定（0x18 状态机，带滞回与确认计数） ----
//   低压穿越（0x240 门槛）：三相最小电压低于门槛且条件满足时，连续 20 步
//   确认后进入模式 1，有功电流给定切到 0x5A0；恢复条件满足并维持 20 步才退出；
//   高压穿越（0x250 门槛）：对称逻辑进入模式 2；
//   0x1C/0x20：穿越结束后保持 800 步（0.2 秒）再回到正常模式。
// ------------------------------------------------------------------
    if (s.d(0x1218) > s.d(0x260)) {
        s.i(0x304, 0);
        if (s.i(0x2EC) < 20) s.increment(0x2EC);
        else s.i(0x308, 1);
    } else s.i(0x2EC, 0);
    if (s.i(0x308) == 1) {
        if (s.d(0x268) > s.d(0x1218)) {
            if (s.i(0x304) < 20) s.increment(0x304);
            else s.i(0x308, 0);
        } else s.i(0x304, 0);
    }
    if (s.d(0x240) > s.d(0x1768) && s.i(0x310) == 0 && s.i(0x518) == 1 && s.d(0x270) >= s.d(0x230)) {
        s.i(0x2E0, 0);
        s.i(0x2F0, 0);
        s.i(0x2F4, 0);
        s.increment(0x334);
        if (s.i(0x2D8) < 20) s.increment(0x2D8);
        else {
            if (s.i(0x18) != 1) s.increment(0xB38);
            s.i(0x18, 1);
            s.d(0x140, s.d(0x5A0));
        }
    } else if (s.d(0x1760) > s.d(0x250) && s.i(0x51C) == 1 && s.d(0x230) > s.d(0x270) && s.d(0xCA0) == 0.0) {
        s.i(0x334, 0);
        s.i(0x2D8, 0);
        s.i(0x2DC, 0);
        s.i(0x2F8, 0);
        s.i(0x2FC, 0);
        if (s.i(0x2E0) < 4) s.increment(0x2E0);
        else {
            if (s.i(0x18) != 2) s.increment(0xB3C);
            s.i(0x18, 2);
            s.d(0x140, 1.25 * s.d(0x110));
        }
    } else {
        s.i(0x2DC, 0);
        s.i(0x2E4, 0);
        s.i(0x2D8, 0);
        s.i(0x2E0, 0);
        s.i(0x334, 0);
        if (s.i(0x18) == 1 && s.i(0x30C) == 0) {
            if (s.d(0x1768) > s.d(0x248)) {
                s.increment(0x2F0);
                if (s.i(0x2F0) > 20) {
                    s.i(0x2F0, 0);
                    s.i(0x18, 0);
                    s.d(0x2C8, 0.0);
                    s.d(0x140, 1.25 * s.d(0x110));
                }
            } else s.i(0x2F0, 0);
        } else if (s.i(0x18) == 2 && s.i(0x310) == 0) {
            if (s.d(0x258) > s.d(0x1760) && 56.559999999999995 > s.d(0x1218)) {
                s.increment(0x2F8);
                if (s.i(0x2F8) > 4) {
                    s.i(0x2F8, 0);
                    s.i(0x18, 0);
                    s.d(0x2C8, 0.0);
                    s.d(0x140, 1.25 * s.d(0x110));
                }
            } else if (1.02 * s.d(0x1780) > s.d(0x1760)) {
                s.increment(0x2F8);
                if (s.i(0x2F8) > 4) {
                    s.i(0x2F8, 0);
                    s.i(0x18, 0);
                    s.d(0x2C8, 0.0);
                    s.d(0x140, 1.25 * s.d(0x110));
                }
            } else s.i(0x2F8, 0);
        }
    }
    if (s.i(0x18) != 0) {
        s.i(0x1C, 1);
        s.i(0x20, 0);
    } else if (s.i(0x1C) == 1) {
        if (s.i(0x20) < 800) s.increment(0x20);
        else {
            s.i(0x20, 0);
            s.i(0x1C, 0);
        }
    }
// ---- 模式 1（LVRT）：无功优先。
//   0x2D0 = (额定无功系数 0x638 - 电压标幺) 的 PI/低通处理 × LVRT 无功系数 0x278，
//   限幅后乘电压基准得无功电流给定；0x5C8/0x5D0 为 d/q 轴电流给定（限幅 ±0.488）；
//   0x2C8 取历史最大值（无功保持），0x2C0 = √(Imax² - Iq²) 有功电流上限 ----
    if (s.i(0x18) == 1) {
        s.d(0x590, 0.0);
        s.d(0x578, 0.0);
        s.i(0x318, 0);
        if (s.i(0x314) < s.i(0x324)) {
            s.increment(0x314);
            s.i(0x30C, 1);
        } else s.i(0x30C, 0);
        s.d(0x2D0, s.d(0x638) - s.d(0x230) / s.d(0x270));
        s.d(0x570, (30.0 * s.d(0x2D0) + 994.0 * s.d(0x570)) / 1024.0);
        if (s.i(0x5C4) == 1) s.d(0x2D0, s.d(0x278) * s.d(0x570));
        else s.d(0x2D0, s.d(0x278) * s.d(0x2D0));
        if (0.0 > s.d(0x2D0)) s.d(0x2D0, s.d(0xC58));
        else if (s.d(0x2D0) > 1.1719) s.d(0x2D0, 1.1719);
        s.d(0x2D0, s.d(0x2D0) * s.d(0x110));
        for (int axis = 0; axis < 2; ++axis) {
            const std::size_t output = 0x5C8 + 0x10 * axis;
            const std::size_t filtered = output + 8;
            const double t = clip(s.d(0x11E0 + 8 * axis) / s.d(0x270), -0.48828125, 0.48828125);
            const double multiplier = axis == 0 ? s.d(0x5E8) : -s.d(0x5E8);
            s.d(output, multiplier * t);
            s.d(output, clip(s.d(output), -0.48828125, 0.48828125));
            if (s.i(0x5C4) == 1) s.d(filtered, (30.0 * s.d(output) + 994.0 * s.d(filtered)) / 1024.0);
            else s.d(filtered, s.d(output));
            if (s.i(0x308) == 1 && s.i(0x620) == 1) {
                s.d(output, s.d(filtered) * s.d(0x110));
                s.d(output, clip(s.d(output), -s.d(0x5F8), s.d(0x5F8)));
            } else s.d(output, 0.0);
        }
        if (s.i(0x5C4) == 1) s.d(0x2C8, s.d(0x2D0));
        else if (s.d(0x2D0) > s.d(0x2C8)) s.d(0x2C8, s.d(0x2D0));
        s.d(0x2C0, (2.0 * s.d(0x110)) * s.d(0x110) - s.d(0x58) * s.d(0x58));
        if (0.0 > s.d(0x2C0)) s.d(0x2C0, 0.0);
        s.d(0x2C0, std::sqrt(s.d(0x2C0)));
        if (s.d(0x2C0) > s.d(0x2C8)) s.d(0x2C0, s.d(0x2C8));
        if (s.i(0x308) == 1 && s.d(0x2C0) > s.d(0x2A0)) s.d(0x2C0, s.d(0x2A0));
        if (s.d(0x11B0) > s.d(0x270) * 0.95) {
            if (s.d(0x60) > s.d(0x2A8)) s.d(0x60, s.d(0x2A8));
            if (s.d(0x2C0) > s.d(0x2A8)) s.d(0x2C0, s.d(0x2A8));
        }
        if (s.d(0x2C0) > s.d(0x288)) s.d(0x2C0, s.d(0x288));
// ---- 模式 2（HVRT）：吸收无功。
//   0x2D0 = (电压标幺 - 0x640) × HVRT 无功系数 0x280，限幅 0.4883 后
//   取负乘电压基准（吸收无功即负的无功电流）；0x590/0x598 是电压平滑项 ----
    } else if (s.i(0x18) == 2) {
        s.d(0x570, 0.0);
        s.i(0x314, 0);
        if (s.i(0x318) < s.i(0x320)) {
            s.increment(0x318);
            s.i(0x310, 1);
        } else s.i(0x310, 0);
        s.d(0x2D0, s.d(0x11D0) / s.d(0x270) - s.d(0x640));
        s.d(0x578, (30.0 * s.d(0x2D0) + 994.0 * s.d(0x578)) / 1024.0);
        if (s.i(0x5C4) == 1) s.d(0x2D0, s.d(0x280) * s.d(0x578));
        else s.d(0x2D0, s.d(0x280) * s.d(0x2D0));
        s.d(0x590, (20.0 * s.d(0x1218) + 1004.0 * s.d(0x590)) / 1024.0);
        s.d(0x598, (0.5 * s.d(0x590)) / s.d(0x270));
        if (0.0 > s.d(0x598)) s.d(0x598, 0.0);
        else if (s.d(0x598) > 0.1) s.d(0x598, 0.1);
        s.d(0x2D0, s.d(0x2D0) + s.d(0x598));
        if (0.0 > s.d(0x2D0)) s.d(0x2D0, 0.0);
        else if (s.d(0x2D0) > 0.4883) s.d(0x2D0, 0.4883);
        s.d(0x2D0, (-s.d(0x2D0)) * s.d(0x110));
        if (s.d(0x2D0) > s.d(0x298)) s.d(0x2D0, s.d(0x298));
        if (s.d(0x60) > s.d(0x298)) s.d(0x60, s.d(0x298));
        if (s.i(0x5C4) == 1) s.d(0x2C8, s.d(0x2D0));
        else if (s.d(0x2C8) > s.d(0x2D0)) s.d(0x2C8, s.d(0x2D0));
        if (s.d(0x290) > s.d(0x2C8)) s.d(0x2C8, s.d(0x290));
        if (s.d(0x290) > s.d(0x60)) s.d(0x60, s.d(0x290));
        s.d(0x2C0, s.d(0x2C8));
        s.d(0x148, (2.0 * s.d(0x110)) * s.d(0x110) - s.d(0x2C0) * s.d(0x2C0));
        if (0.0 > s.d(0x148)) s.d(0x148, 0.0);
        s.d(0x148, std::sqrt(s.d(0x148)));
        if (s.d(0x140) > s.d(0x148)) s.d(0x140, s.d(0x148));
    }
// ---- 三种模式下的 d/q 电流给定合成（0x50/0x600/0x608），
//      含交叉解耦项（0x238×电压/基准）与限幅（±0x120）；
//      0x610/0x618 是给负序电流环的给定（正常模式为 0）----
    if (s.i(0x18) == 1) {
        s.d(0x50, (s.d(0x2C0) + s.d(0x48) * s.d(0x648)) - (s.d(0x238) * s.d(0x230)) / s.d(0x270));
        s.d(0x600, s.d(0x5D8));
        s.d(0x608, s.d(0x5C8));
    } else if (s.i(0x18) == 2) {
        s.d(0x50, (s.d(0x2C0) + s.d(0x48) * s.d(0x650)) - (s.d(0x238) * s.d(0x230)) / s.d(0x270));
        s.d(0x600, 0.0);
        s.d(0x608, 0.0);
    } else {
        s.d(0x600, 0.0);
        s.d(0x608, 0.0);
        s.i(0x318, 0);
        s.i(0x314, 0);
        s.d(0x570, 0.0);
        s.d(0x578, 0.0);
        s.d(0x2C0, 0.0);
        s.d(0x590, 0.0);
        s.d(0x630, (s.d(0x238) * s.d(0x230)) / s.d(0x270));
        s.d(0x50, s.d(0x48) - s.d(0x630));
        s.d(0x50, clip(s.d(0x50), -s.d(0x120), s.d(0x120)));
    }
    ramp(0x610, 0x600, 0x78, 1.0);
    ramp(0x618, 0x608, 0x78, 1.0);
// ---- 300Hz 谐振滤波（0x15A8/0x1660 两组），用于抑制特定次谐波 ----
    s.d(0x15C0, 0.0);
    s.d(0x15C8, s.d(0x13F0));
    resonant_filter(s.at(0x15A8));
    s.d(0x1678, 0.0);
    s.d(0x1680, s.d(0x13F8));
    resonant_filter(s.at(0x1660));
// ---- 第 9 步：四个电流内环 PI：
//      0xF08：正序 d 轴（给定 0x58、反馈 0x90）；0xF88：正序 q 轴（0x60/0x98）；
//      0x1008：负序 d 轴（0x610/0xA0）；0x1088：负序 q 轴（0x618/0xA8） ----
    s.d(0xF08, s.d(0x58));
    s.d(0xF10, s.d(0x90));
    s.apply(0xF08, regulator);
    s.d(0xF88, s.d(0x60));
    s.d(0xF90, s.d(0x98));
    s.apply(0xF88, regulator);
    s.d(0x1008, s.d(0x610));
    s.d(0x1010, s.d(0xA0));
    s.apply(0x1008, regulator);
    s.d(0x1088, s.d(0x618));
    s.d(0x1090, s.d(0xA8));
    s.apply(0x1088, regulator);
// ---- 谐波抑制与负序补偿：0xB0/0xB8 为谐波补偿量（限幅 0.61），
//      0x540/0x548 为差分滤波修正（限幅 ±0.1），
//      0x1480/0x14A8 把正序、负序电压参考分别反 Park 回 αβ ----
    s.d(0xB0, (s.d(0x1220) / s.d(0x1B8)) * 0.577);
    s.d(0xB8, (s.d(0x1228) / s.d(0x1B8)) * 0.577);
    if (s.d(0xB0) > 0.61) s.d(0xB0, 0.61);
    if (s.d(0xB8) > 0.61) s.d(0xB8, 0.61);
    s.d(0x1548, s.d(0x1260));
    s.d(0x1588, s.d(0x1268));
    difference_filter(s.at(0x1528));
    difference_filter(s.at(0x1568));
    s.d(0x540, (0.577 * s.d(0x1558)) / s.d(0x1B8));
    s.d(0x548, (0.577 * s.d(0x1598)) / s.d(0x1B8));
    s.d(0x540, clip(s.d(0x540), -0.1, 0.1));
    s.d(0x548, clip(s.d(0x548), -0.1, 0.1));
    s.d(0x1498, s.d(0xB0) - s.d(0xF58) / s.d(0x80));
    s.d(0x14A0, s.d(0xB8) - s.d(0xFD8) / s.d(0x80));
    s.d(0x1490, s.d(0x218));
    s.apply(0x1480, inverse_park_transform);
    s.d(0x14C0, (-s.d(0x1058)) / s.d(0x80));
    s.d(0x14C8, (-s.d(0x10D8)) / s.d(0x80));
    s.d(0x14B8, -s.d(0x218));
    s.apply(0x14A8, inverse_park_transform);
// ---- 第 10 步：调制。正序 + 负序 + 修正 → 0xD00 调制模块 →
//      输出占空比 0xDF0 起 → 拷贝到输出区 0x4C0 起 ----
    s.d(0xEC0, s.d(0xD8));
    s.d(0xEC8, s.d(0xE0));
    s.d(0xED0, s.d(0xE8));
    s.d(0xDC8, s.d(0x100));
    s.d(0xDD0, s.d(0x108));
    s.d(0xD00, (s.d(0x1480) + s.d(0x14A8)) + s.d(0x540));
    s.d(0xD08, (s.d(0x1488) + s.d(0x14B0)) + s.d(0x548));
    s.apply(0xD00, modulation);
    s.d(0x4C0, s.d(0xDF0));
    s.d(0x4C8, s.d(0xDF8));
    s.d(0x4D8, s.d(0xE00));
    s.d(0x4E0, s.d(0xE08));
    s.d(0x4F0, s.d(0xE10));
    s.d(0x4F8, s.d(0xE18));
// ---- 软启动使能（0x31C）：启动阶段 PWM 使能 1，超过 0x28 步后转正常 ----
    if (s.i(0x31C) == 1) {
        s.i(0x50C, 1);
        if (s.i(0x28) > 1) {
            for (std::size_t offset : {0x60, 0xF30, 0xFB0, 0x1030, 0x10B0, 0x1130, 0x138}) s.d(offset, 0.0);
            s.i(0x158, 0);
        }
        s.increment(0x24);
        if (s.i(0x24) > s.i(0x28)) {
            s.i(0x24, 0);
            s.i(0x31C, 0);
            s.i(0x158, 1);
        }
    } else {
        s.i(0x50C, 1);
        s.i(0x2C, 1);
    }
// ---- exit_tail：停机过程（0x350==2）持续 0x35C 步（0.2 秒）后整块复位 ----
exit_tail:
    if (s.i(0x350) == 2) {
        s.increment(0x358);
        if (s.i(0x358) > s.i(0x35C)) {
            std::memset(state, 0, 0x1AB8);
            grid_defaults(f, state);
        }
    }
// ---- final_flags：0x360 = 运行标志（1/2 → 1），供输出与监控使用 ----
final_flags:
    s.i(0x360, s.i(0x350) == 1 || s.i(0x350) == 2 ? 1 : 0);
}
