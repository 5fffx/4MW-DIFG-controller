// ============================================================================
// 文件：recovered_parameters.cpp —— 参数与保护定值的输入换算
// ----------------------------------------------------------------------------
// 两个入口：
//   parameter_input ：PSCAD 的 90 个参数（PI 系数、LVRT/HVRT 定值、限幅……）
//                     从标幺值换算成内部控制用的有名值，写进状态区；
//   protection_input：28 路保护定值（过压/欠压/过频/欠频各段阈值和时间）。
// 参数名与位置的对照表见《参数说明.md》。
// ============================================================================
#include "recovered_primitives.h"

#include <cstddef>
#include <cstring>
#include <emmintrin.h>

namespace {

// ParameterState：帮助类，和 CommandState 一样是“状态区读写器”。
// 构造时记住状态区首地址，之后 d(offset)/i(offset) 按字节偏移量读写。
class ParameterState {
public:
    explicit ParameterState(void* state) : bytes_(static_cast<unsigned char*>(state)) {}

    double d(std::size_t offset) const {
        double value;
        std::memcpy(&value, bytes_ + offset, sizeof(value));
        return value;
    }

    void d(std::size_t offset, double value) {
        std::memcpy(bytes_ + offset, &value, sizeof(value));
    }

    std::int32_t i(std::size_t offset) const {
        std::int32_t value;
        std::memcpy(&value, bytes_ + offset, sizeof(value));
        return value;
    }

    void i(std::size_t offset, std::int32_t value) {
        std::memcpy(bytes_ + offset, &value, sizeof(value));
    }

private:
    unsigned char* bytes_;
};

} // namespace

// ----------------------------------------------------------------------------
// protection_input —— 保护定值入口。
// 把 28 路保护相关输入写进状态区：前 16 个是电压类定值（0x900 起），
// 后 12 个是频率类定值（0xA78 起）。
// 电压偶数下标项乘电压基准 voltage_base（= 0x270/1.414，有名值化），
// 奇数下标项（时间类）直接存入。
// 开头的 settings 表是这些定值的“默认值”：
// 过压四段（1.33/1.28/1.23/1.13 标幺，对应时间 0.15/0.6/1.1/10.4 秒）、
// 欠压四段（0.12/0.4/0.6/0.83，对应 0.3/1.4/2.2/3.0 秒）、
// 过频/欠频各三段（54/53/52 Hz 与 46/47/47.5 Hz）等。
// 这些默认值随后会被 PSCAD 传来的实际定值覆盖（见参数说明 Protect_*）。
// ----------------------------------------------------------------------------
extern "C" __declspec(dllexport) void protection_input(const double* input, double, void* state) {
    ParameterState s(state);
    const double voltage_base = s.d(0x270) / 1.414;
    struct Setting {
        std::size_t offset;
        double value;
    };
    const Setting settings[] = {
        {0x800, 1350.0}, {0x808, 0.5529}, {0x820, 0.004},
        {0x828, 1300.0}, {0x830, 1280.0}, {0x838, 500.0},
        {0x850, 0.01}, {0x858, 0.0025}, {0x860, 500.0},
        {0x870, 0.005}, {0x888, 0.0316}, {0x890, 0.002},
        {0x8C8, 0.4}, {0x8D0, 0.7}, {0x8E8, 0.0395}, {0x8F0, 12.0},
        {0xA18, 0.7503}, {0xA20, 0.004}, {0xA28, 0.4836},
        {0xA30, 1.9344}, {0xA38, 2.5792}, {0xA40, 3.8689},
        {0xA60, 0.2}, {0xA68, 0.1}, {0xA70, 0.00025},
        {0xB20, 3.0}, {0xB28, 12.0}, {0xB58, 15.0}, {0xB48, 2.0},
        {0xB78, 0.002}, {0xB80, 0.6}, {0xB88, 0.15},
        {0xBA0, 15.0}, {0xBA8, 0.15}, {0xBB8, 4.0},
        {0x2410, 3.0}, {0x2418, 12.0}, {0x2448, 15.0}, {0x2438, 2.0},
        {0x2468, 0.002}, {0x2470, 0.6}, {0x2478, 0.15},
        {0x2490, 15.0}, {0x2498, 0.15}, {0x24A8, 4.0},
        {0x24B0, 1.5796}, {0x24C0, 0.004}, {0x24D8, 0.1974}, {0x24E0, 0.005},
        {0x24F8, 0.158}, {0x2500, 0.6}, {0x2528, 1.3821}, {0x2530, 0.00025},
        {0x2578, 1.1847}, {0x2570, 0.1}, {0x2548, 0.1947}, {0x2550, 0.00125},
        {0x2598, 0.3949}, {0x25A0, 0.1}, {0x25B0, 3.949}, {0x25C0, 1.0},
        {0x25C8, 3.949}, {0x25D8, 1.0}
    };
    for (const Setting& setting : settings) {
        s.d(setting.offset, setting.value);
    }
    const std::size_t voltage_offsets[] = {
        0x900, 0x918, 0x920, 0x938, 0x940, 0x958, 0x980, 0x998,
        0x960, 0x978, 0x9A0, 0x9B8, 0x9C0, 0x9D8, 0x9E0, 0x9F8
    };
    for (std::size_t j = 0; j < 16; ++j) {
        s.d(voltage_offsets[j], j % 2 == 0 ? input[j] * voltage_base : input[j]);
    }
    const std::size_t frequency_offsets[] = {
        0xA78, 0xA88, 0xA90, 0xAA0, 0xAA8, 0xAB8,
        0xAC0, 0xAD0, 0xAD8, 0xAE8, 0xAF0, 0xB00
    };
    for (std::size_t j = 0; j < 12; ++j) {
        s.d(frequency_offsets[j], input[j + 16]);
    }
}

// The caller supplies the full parameter array starting at zero-based index 10.
// ----------------------------------------------------------------------------
// parameter_input —— 参数换算入口：PSCAD 参数 → 内部物理量
// ----------------------------------------------------------------------------
// p[] 是从 PSCAD 传来的参数数组（从 parameters+10 开始，即跳过前 10 个命令位）。
// 大部分参数在 PSCAD 里是标幺值，这里按基准换算成有名值。
// ----------------------------------------------------------------------------
extern "C" __declspec(dllexport) void parameter_input(const double* p, double, void* state) {
    ParameterState s(state);
// 三个换算公式（lambda 是 C++ 里的“小函数”，h/f/r 是函数名）：
//   h(n) = p[n] * (0x538/0x270)/√3 * √2   电流类：标幺 → 有名值电流
//   f(n) = h(n) * 1.414                     再乘 √2（峰值化）
//   r(n) = p[n]*(0x538/1140)/√3 / (a+b)     转子侧换算（extra 时再 ×0.75）
// 其中 0x538 = 额定视在功率基准（4 MVA），0x270 = 网侧电压基准，
// 1140 = 转子额定电压，a/b 是机侧电流基准的组合。
    const auto h = [&s, p](std::size_t n) {
        double q = s.d(0x538) / s.d(0x270);
        q = q / 1.732;
        q = q * 1.414;
        return p[n] * q;
    };
    const auto f = [&h](std::size_t n) {
        return h(n) * 1.414;
    };
    const auto r = [&s, p](std::size_t n, bool extra) {
        double q = s.d(0x538) / 1140.0;
        q = q / 1.732;
        q = p[n] * q;
        double a = s.d(0x110) / 1.414;
        double b = s.d(0x1B00) / 1.414;
        b = b * s.d(0x2F58);
        if (extra) {
            b = b * 0.75;
        }
        a = a + b;
        return q / a;
    };
    const auto cvtt = [](double value) {
        return _mm_cvttsd_si32(_mm_set_sd(value));
    };

// p[0..11]：各控制环的 PI 参数（kp 和 ki 成对出现）。
//   p[0]/p[1] = PLL_kp / PLL_ki（锁相环）；
//   p[2]/p[3] = 网侧电压外环 kp / ki；
//   p[4]/p[5] = 网侧正序电流环 kp / ki（存两处：0xF20 与 0xFA0）；
//   p[6]/p[7] = 网侧负序电流环 kp / ki；
//   p[8]/p[9] = 机侧正序电流环 kp / ki；p[10]/p[11] = 机侧负序电流环。
// 0x1170、0xF70…… 这些“相除”是把 ki/kp 折算成调节器内部的每步积分增益。
    s.d(0x550, p[0]);
    s.d(0x560, p[1]);
    s.d(0x1120, p[2]);
    s.d(0x1168, p[3]);
    s.d(0x1170, s.d(0x1168) / s.d(0x1120));
    s.d(0xF20, p[4]);
    s.d(0xF68, p[5]);
    s.d(0xF70, s.d(0xF68) / s.d(0xF20));
    s.d(0xFA0, p[4]);
    s.d(0xFE8, p[5]);
    s.d(0xFF0, s.d(0xFE8) / s.d(0xFA0));
    s.d(0x1020, p[6]);
    s.d(0x1068, p[7]);
    s.d(0x1070, s.d(0x1068) / s.d(0x1020));
    s.d(0x10A0, p[6]);
    s.d(0x10E8, p[7]);
    s.d(0x10F0, s.d(0x10E8) / s.d(0x10A0));
    s.d(0x2058, p[0]);
    s.d(0x2068, p[1]);
    s.d(0x1CC0, p[8]);
    s.d(0x1CC8, p[9]);
    if (s.i(0x1AC8) == 1) {
        s.d(0x560, p[1] / 10.0);
        s.d(0x2068, p[1] / 10.0);
    }
    s.d(0x1CD0, p[10]);
    s.d(0x1CD8, p[11]);
// p[12..19]：LVRT/HVRT 门槛与无功系数。
//   p[12]/p[13] = 低压穿越门槛 0.9 / 滞回 0.05 → 0x240/0x248（有名值）；
//   p[14] = LVRT 无功系数 3.0 → 0x278；
//   p[16]/p[17] = 高压穿越门槛 1.13 / 滞回 0.05 → 0x250/0x258；
//   p[18] = HVRT 无功系数 2.2841 → 0x280。
// 后面 0x1E68 起的一串是给“机侧”用的同一组定值。
    s.d(0x240, (s.d(0x270) * p[12]) / 1.732);
    s.d(0x248, (s.d(0x270) * (p[12] + p[13])) / 1.732);
    s.d(0x278, r(14, false));
    s.d(0x250, (s.d(0x270) * p[16]) / 1.732);
    s.d(0x258, (s.d(0x270) * (p[16] - p[17])) / 1.732);
    s.d(0x280, r(18, true));
    s.d(0x1E68, (s.d(0x1EC8) * p[12]) / 1.732);
    s.d(0x1E70, (s.d(0x1EC8) * (p[12] + p[13])) / 1.732);
    s.d(0x1ED0, r(14, false));
    s.d(0x1E78, (s.d(0x1EC8) * p[16]) / 1.732);
    s.d(0x1E80, (s.d(0x1EC8) * (p[16] - p[17])) / 1.732);
    s.d(0x1ED8, (r(19, true) * 3.0) / 4.0);
// p[20..25]：电流限幅（标幺 → 有名值）。
//   网侧有功电流限幅 0.2777 / 无功电流限幅 0.1111（0x140/0x120）；
//   机侧有功电流限幅 0.4 / 0.0，无功电流限幅 ±0.2895（0x1B08 起）。
// 注意 0x140 连续写两次——原 DLL 就是这样（后写覆盖前写），原样保留。
    s.d(0x140, f(20));
    s.d(0x140, f(20));
    s.d(0x120, f(21));
    s.d(0x1B08, f(22));
    s.d(0x1B10, f(23));
    s.d(0x1B18, f(24));
    s.d(0x1B20, f(25));
// p[26..45]：LVRT/HVRT 控制细节参数。
//   p[26]/p[27]/p[28] = LVRT 无功电流计算电压/初始系数/负序系数；
//   p[29..32] = 网侧/机侧 LVRT 无功电流最大/最小值；
//   p[33..39] = HVRT 对应的一组参数；
//   p[40]/p[41]/p[42] = LVRT 有功电流控制模式/系数、HVRT 有功功率控制系数；
//   p[43] = 转子 d 轴电流限幅上限；p[44] = LVRT 最小 d 轴电流。
// （p[45..55] 在原 DLL 里没有被使用，参数表里也不存在，跳过。）
    s.d(0x638, p[26]);
    s.d(0x648, p[27]);
    s.d(0x5E8, p[28]);
    s.d(0x21A8, p[26]);
    s.d(0x21B8, p[27]);
    s.d(0x20E8, p[28]);
    s.d(0x288, f(29));
    s.d(0x2A0, f(30));
    s.d(0x1EE0, f(31));
    s.d(0x1F08, f(32));
    s.d(0x640, p[33]);
    s.d(0x650, p[34]);
    s.d(0x5F0, p[35]);
    s.d(0x21B0, p[33]);
    s.d(0x21C0, p[34]);
    s.d(0x20F0, p[35]);
    s.d(0x298, f(36));
    s.d(0x290, f(37));
    s.d(0x1EF0, f(38));
    s.d(0x1EE8, f(39));
    s.i(0x21C8, cvtt(p[40]));
    s.d(0x21D0, p[41]);
    s.d(0x21D8, p[42]);
    s.d(0x1F10, f(43));
    s.d(0x2648, h(44));
    s.d(0x370, 1.4);
    s.i(0x658, 100);
    s.d(0x378, 0.1);
    s.i(0x65C, 600);
    s.d(0x380, 1.1);
    s.i(0x660, 100);
    s.d(0x388, 0.9);
    s.i(0x664, 100);
    s.d(0x368, 2.6889);
    s.d(0x1FC8, 2.6889);
    s.d(0x390, 2500.0);
    s.d(0x398, 1200.0);
// p[56..67]：频率支撑（惯量响应 + 一次调频）参数。
//   p[56] = 惯量响应时间 10s；p[57] = 惯量死区 0.3；p[58] = 惯量时间常数 0.1；
//   p[59] = 一次调频响应时间 5s；p[60] = 功频下垂系数 0.06；
//   p[61] = 一次调频死区 0.03；p[62]/p[63] = 一次调频功率上/下限；
//   p[64]/p[65] = 恢复斜坡起点/时间；p[66]/p[67] = 网侧/机侧 LVRT 无功最小值。
    s.d(0xBC8, p[56]);
    s.d(0xC38, p[57]);
    s.d(0x720, p[58]);
    s.d(0x7B0, p[59]);
    s.d(0xC30, p[60]);
    s.d(0x738, p[61]);
    s.d(0xC40, p[62]);
    s.d(0xC48, p[63]);
    s.d(0x2600, p[64]);
    s.d(0x2610, p[65]);
    s.d(0xC58, p[66]);
    s.d(0x2620, p[67]);
// p[68..74]：LVRT 无功下垂与恢复参数。
//   p[68] = 电压偏差死区（0x36A0）；p[69] = 无功增发比例系数（0x3698）；
//   p[70] = LVRT 无功控制参考电压（0x36A8）；p[71] = PWM 200ms 断路器使能；
//   p[72] = LVRT d 轴电流模式（2 = 保持故障前值，覆盖 p[40]）；
//   p[73] = LVRT d 轴电流限幅 0.7；p[74] = 恢复斜坡模式 3；
//   p[79] = 另一个恢复相关定值（0x36F8）。
    s.d(0x3698, p[69]);
    s.d(0x36A0, p[68]);
    s.d(0x36A8, p[70]);
    s.d(0xC98, p[71]);
    s.i(0x21C8, cvtt(p[72]));
    s.d(0x36E8, p[73]);
    s.d(0x2618, p[74]);
    s.d(0x36F8, p[79]);
// 尾部固定常数：0x778 = 1.0（使能标志），
// 0x1FD0/0x1FD8 = 机侧低压/高压限幅初值（0x370 与 0x368+0x1FC8）。
    s.d(0x778, 1.0);
    s.d(0x1FD0, s.d(0x370));
    s.d(0x1FD8, s.d(0x368) + s.d(0x1FC8));
}
