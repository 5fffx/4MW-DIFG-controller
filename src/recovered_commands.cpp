// ============================================================================
// 文件：recovered_commands.cpp —— 命令与给定值的入口
// ----------------------------------------------------------------------------
// PSCAD 每步给控制器送 11 路命令 commands[0..10]（启动、模式、复位等），
// 本函数把它们写进状态区相应位置，并完成“无功指令三种模式”的换算。
// 这是“操作员/上位机 → 控制器”的唯一入口。
// ============================================================================
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <emmintrin.h>

namespace {

// CommandState：帮助类。把 void* 状态指针包装成带偏移量读写的对象：
// d(offset) 读/写 double，i(offset) 读/写 32 位整数。
class CommandState {
public:
    explicit CommandState(void* state) : bytes_(static_cast<unsigned char*>(state)) {}

    double d(std::size_t offset) const { return read<double>(offset); }
    std::int32_t i(std::size_t offset) const { return read<std::int32_t>(offset); }
    void d(std::size_t offset, double value) { write(offset, value); }
    void i(std::size_t offset, std::int32_t value) { write(offset, value); }

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

// cvtt_i32：双精度 → 32 位整数的“向零取整”转换。
// 用的是 x86 的 CVTT 指令（通过 <emmintrin.h> 内建函数），
// 保证和原 DLL 的截断行为逐位一致（包括异常输入时产生 INT_MIN 的细节）。
std::int32_t cvtt_i32(double value) {
    // CVTT preserves the hardware INT_MIN result for invalid conversions.
    return _mm_cvttsd_si32(_mm_set_sd(value));
}

} // namespace

// ----------------------------------------------------------------------------
// command_input —— 命令与给定值的入口
// ----------------------------------------------------------------------------
// commands[] 含义（下标从 0 开始，对应 PSCAD 里的英文标签）：
//   [0] 启动命令 Start_Cmd        [1] 有功给定 P*       [2] 无功相关给定 Q*
//   [3] LVRT 使能                  [4] HVRT 使能         [5] 复位命令 Reset_Cmd
//   [6] LVRT 恢复标志              [7] 其他标志          [8] 无功控制模式 Q_ctrol_mode
//   [9] 输出模式 Flag_AVG_IGBT    [10] 频率支撑模式 Flag_Interia_Freq
//
// 无功给定的三种模式（0x21A4 存的是 commands[8]，注意用的是“上一拍”的
// 模式值，本拍的新值存到最后才更新——这是原 DLL 的执行顺序，原样保留）：
//   模式 0：恒功率因数。给定有功 P 和功率因数 pf，先算视在功率 S = P/pf，
//           再反推无功 Q = √(S² - P²)。pf 被限制在 |pf| ≥ 0.95。
//   模式 2：直接无功给定 + LVRT 无功增发。当电压偏差（1.0 - 机端电压标幺）
//           超过死区（0x36A0）时，按系数（0x3698）比例增发无功。
//   模式 1：直接无功给定，限幅 ±0.33。
// ----------------------------------------------------------------------------
extern "C" __declspec(dllexport) void command_input(const double* commands, void* state) {
    CommandState s(state);
    s.i(0x348, cvtt_i32(commands[0]));
    s.i(0x1FAC, cvtt_i32(commands[0]));

    double p;
    double q;
    if (s.d(0xC98) == 1.0) {
        if (s.d(0x7D0) == 1.0) {
            s.d(0x7D0, 0.0);
            s.d(0x36D8, commands[1]);
            s.d(0x36E0, commands[2]);
        }
        p = s.d(0x36D8);
        q = s.d(0x36E0);
    } else {
        p = commands[1];
        q = commands[2];
    }
    s.d(0x1D58, p);

    // The previous mode applies now; commands[8] is stored only afterward.
    switch (s.i(0x21A4)) {
    case 0: {
        if (q > 0.0 && q < .95) {
            q = .95;
        } else if (q < 0.0 && q > -.95) {
            q = -.95;
        }
        // Zero is intentionally not guarded: retain IEEE division behavior.
        const double apparent = p / q;
        s.d(0x1D60, apparent * apparent - p * p);
        if (s.d(0x1D60) < 0.0) {
            s.d(0x1D60, 0.0);
        }
        s.d(0x1D60, std::sqrt(s.d(0x1D60)));
        if (q < 0.0) {
            s.d(0x1D60, -s.d(0x1D60));
        }
        break;
    }
    case 2:
        s.d(0x1D60, q);
        if (s.i(0x1AD8) == 1) {
            const double deviation = s.d(0x36A8) - s.d(0x1E60) / s.d(0x1EC8);
            if (deviation >= s.d(0x36A0)) {
                s.d(0x3690, s.d(0x3698) * deviation);
                s.d(0x1D60, s.d(0x1D60) + s.d(0x3690));
            } else if (deviation <= -s.d(0x36A0)) {
                s.d(0x3690, s.d(0x3698) * deviation);
                s.d(0x1D60, s.d(0x1D60) + s.d(0x3690));
            } else {
                s.d(0x3690, 0.0);
                s.d(0x1D60, q);
            }
        }
        break;
    case 1:
        s.d(0x1D60, q);
        if (s.d(0x1D60) > .33) {
            s.d(0x1D60, .33);
        } else if (s.d(0x1D60) < -.33) {
            s.d(0x1D60, -.33);
        }
        break;
    default:
        break;
    }

    s.i(0x518, cvtt_i32(commands[3]));
    s.i(0x1FEC, cvtt_i32(commands[3]));
    s.i(0x51C, cvtt_i32(commands[4]));
    s.i(0x1FF0, cvtt_i32(commands[4]));

    if (commands[5] == 1.0 && s.i(0x3720) == 0) {
        s.i(0x520, cvtt_i32(commands[5]));
        s.i(0x1FF4, cvtt_i32(commands[5]));
    }
    s.i(0x3720, cvtt_i32(commands[5]));

    s.i(0x2088, cvtt_i32(commands[6]));
    s.i(0x7D8, cvtt_i32(commands[7]));
    s.i(0x23D0, cvtt_i32(commands[7]));
    s.i(0x21A4, cvtt_i32(commands[8]));
    s.i(0x628, cvtt_i32(commands[9]));
    s.i(0x21A0, cvtt_i32(commands[9]));
    s.i(0x740, cvtt_i32(commands[10]));
    s.d(0xCA0, s.d(0x2038));
}
