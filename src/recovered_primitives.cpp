// ============================================================================
// 文件：recovered_primitives.cpp —— 控制器的“数学积木箱”
// ----------------------------------------------------------------------------
// 本文件只包含最基本的算法单元，每个函数对应原 DLL 里的一个计算模块。
// 规律：每个函数拿到自己的状态区指针，从状态区读入变量 → 计算 → 写回状态区。
// 它们自己不决定“什么时候被调用”，调用顺序由 recovered_controller.cpp 调度。
// 阅读顺序建议：clip → regulator（PI）→ 坐标变换 → 滤波器 → pll → modulation，
// 其余是输入输出映射和一些高级算法单元。
// ============================================================================
#include "recovered_primitives.h"
#include <cmath>
#include <cstddef>
#include <cstring>

// ----------------------------------------------------------------------------
// 编译期检查：确认结构体字段在内存里的位置和原 DLL 状态区完全一致。
// offsetof(结构体, 字段) = 该字段距结构体开头的字节数；sizeof = 结构体总字节数。
// 例如 PllState 里 rotating 必须从第 0xD8 个字节开始（0xD8 = 十进制 216），
// 整个 PLL 状态区恰好 0x200 = 512 字节。任何一项对不上，编译直接报错——
// 这是防止我们“自己改错内存布局”的第一道保险。
// ----------------------------------------------------------------------------
static_assert(offsetof(PllState, stationary) == 0xB0, "PLL stationary offset");
static_assert(offsetof(PllState, rotating) == 0xD8, "PLL rotating offset");
static_assert(offsetof(PllState, regulator) == 0x180, "PLL regulator offset");
static_assert(sizeof(PllState) == 0x200, "PLL state extent");
static_assert(offsetof(ModulationState, modulation) == 0x88, "Modulation output offset");
static_assert(offsetof(ModulationState, dc_plus) == 0xC8, "Modulation DC offset");
static_assert(offsetof(ModulationState, correction) == 0xE8, "Modulation correction offset");
static_assert(offsetof(ModulationState, duty) == 0xF0, "Modulation duty offset");
static_assert(offsetof(ModulationState, gain) == 0x120, "Modulation gain offset");
static_assert(offsetof(ModulationState, polarity) == 0x1A8, "Modulation polarity offset");
static_assert(sizeof(ModulationState) == 0x1B0, "Modulation state extent");

// clip：限幅（饱和）函数。把 value 限制在 [lower, upper] 区间内。
// 控制系统中到处都要限幅：防止输出过大、防止积分器发散等。
static double clip(double value, double lower, double upper) {
    if (value > upper) return upper;
    if (value < lower) return lower;
    return value;
}

// ----------------------------------------------------------------------------
// regulator —— PI 调节器（带抗饱和、带可选的微分项），控制器里用得最多的单元
// ----------------------------------------------------------------------------
// 每步执行一次，公式：
//   error        = reference - feedback                      （误差）
//   proportional = kp * error                                （比例项）
//   integral     = integral + ki_per_step*error + tracking_gain*antiwindup_error
//                                                           （积分项，逐步累加）
//   derivative   = derivative_gain * (previous_error - error)（微分项，差分近似）
//   unsaturated  = P + I + D                                 （未限幅的总输出）
//   output       = clip(unsaturated, lower, upper)           （限幅后的最终输出）
//   antiwindup_error = output - unsaturated    （被限幅削掉的量，反馈给积分器）
// 控制直觉：
//   - 比例项立刻响应误差；积分项消除稳态误差（只要有误差就不断累加，直到误差为零）；
//   - 若输出撞到限幅（比如调制波已经到顶），被削掉的部分通过 tracking_gain
//     反向灌回积分器，积分值就不会无意义地继续涨——这就是“抗饱和”；
//   - ki_per_step 是“每步的积分增益”（原 DLL 已经把 ki 折算成每步增量）。
// ----------------------------------------------------------------------------
void regulator(RegulatorState* s) {
    const double error = s->reference - s->feedback;
    const double proportional = s->kp * error;
    const double integral = clip(s->integral + s->ki_per_step * error
        + s->tracking_gain * s->antiwindup_error, s->lower_limit, s->upper_limit);
    const double derivative = s->derivative_gain * (s->previous_error - error);
    const double unsaturated = proportional + integral + derivative;
    const double output = clip(unsaturated, s->lower_limit, s->upper_limit);
    s->error = error;
    s->proportional = proportional;
    s->integral = integral;
    s->derivative = derivative;
    s->unsaturated = unsaturated;
    s->output = output;
    s->antiwindup_error = output - unsaturated;
    s->previous_error = error;
}

// ----------------------------------------------------------------------------
// 线电压 → αβ（Clarke 变换）。只用 vab、vbc 两路：
//   alpha = 1.1547*vab + 0.577*vbc；beta = vbc。
// 注意：1.1547 ≈ 2/√3、0.577 ≈ 1/√3 是原 DLL 里被“截断”的常数——
// 原程序用的是截断后的数值而不是精确的 √3。为了和原 DLL 输出逐字节一致，
// 我们原样保留截断值，不做“更精确”的修正。
// ----------------------------------------------------------------------------
void line_voltage_transform(LineVoltageState* s) {
    // Preserve the DLL's truncated coefficients rather than ideal square roots.
    s->alpha = 1.1547 * s->vab + 0.577 * s->vbc;
    s->beta = s->vbc;
}

// ----------------------------------------------------------------------------
// 三相电流 → αβ（等幅值 Clarke 变换，系数 2/3）：
//   alpha = (2/3)*(ia - 0.5*ib - 0.5*ic)；beta = (2/3)*(0.866*ib - 0.866*ic)。
// 0.666667、0.866 同样是原 DLL 的截断常数，原样保留。
// ----------------------------------------------------------------------------
void current_transform(CurrentTransformState* s) {
    s->alpha = (s->ia - 0.5 * s->ib - 0.5 * s->ic) * 0.666667;
    s->beta = (0.866 * s->ib - 0.866 * s->ic) * 0.666667;
}

// Park 变换：αβ → dq。
//   d = alpha*cosθ + beta*sinθ；q = beta*cosθ - alpha*sinθ。
// θ 取锁相环角度时，三相正序基波量在 dq 下变成直流量，便于 PI 无静差控制。
void park_transform(ParkState* s) {
    s->d = s->alpha * std::cos(s->theta) + s->beta * std::sin(s->theta);
    s->q = s->beta * std::cos(s->theta) - s->alpha * std::sin(s->theta);
}

// 反 Park 变换：dq → αβ。把控制环在旋转坐标系里算出的电压参考
// 变回静止坐标系，交给调制模块转成占空比。
void inverse_park_transform(ParkState* s) {
    s->alpha = s->d * std::cos(s->theta) - s->q * std::sin(s->theta);
    s->beta = s->q * std::cos(s->theta) + s->d * std::sin(s->theta);
}

// ----------------------------------------------------------------------------
// 两态滤波器（two-state filter）。
// 用 x、z 两个内部状态实现二阶滤波：先算系数 a = (gain/scale)*step、
// b = (gain*scale)*step；x 追踪“输入误差减 z”（error = input - x，
// next_x = x + a*(error - z)），z 再追踪 x（next_z = z + b*next_x）。
// 两级串联形成谐振特性：对某频率（由 gain 决定）的信号几乎无衰减、
// 对其他频率衰减——所以它常被当作“基波提取器/陷波器”使用。
// ----------------------------------------------------------------------------
void two_state_filter(TwoStateFilter* s) {
    const double a = (s->gain / s->scale) * s->step;
    const double b = (s->gain * s->scale) * s->step;
    s->x_coefficient = a;
    s->z_coefficient = b;
    const double error = s->input - s->x;
    const double next_x = s->x + a * (error - s->z);
    const double next_z = s->z + b * next_x;
    s->error = error;
    s->next_x = next_x;
    s->next_z = next_z;
    s->x = next_x;
    s->z = next_z;
}

// 一阶递归滤波器：y[n] = a0*x[n] + a1*x[n-1] + b1*y[n-1]。
// input_coefficient = a0；previous_input_coefficient = a1；
// -previous_output_coefficient = b1。每步算完把输出和输入存成“上一步”。
void recursive_filter(RecursiveFilterState* s) {
    const double output = ((-s->previous_output_coefficient) * s->previous_output
        + s->input_coefficient * s->input)
        + s->previous_input_coefficient * s->previous_input;
    s->output = output;
    s->previous_output = output;
    s->previous_input = s->input;
}

// ----------------------------------------------------------------------------
// 残差滤波器（两通道，d/q 轴各一路）：
//   1) lowpass = 0.2042*input + 0.7958*lowpass   —— 一阶低通滤波；
//   2) residual = input - lowpass                —— 原始信号减低频分量 = 高频残差；
//   3) 限幅 ±1131.2。
// 0.2042 + 0.7958 = 1（保证直流增益为 1），这两个数同样来自原 DLL。
// 残差分量在电机控制里用于提取谐波，反相注入以抑制谐波。
// ----------------------------------------------------------------------------
void residual_filter(ResidualFilterState* s) {
    for (int channel = 0; channel < 2; ++channel) {
        s->lowpass[channel] = s->input[channel] * 0.20420352248333634
            + s->lowpass[channel] * 0.7957964775166637;
    }
    for (int channel = 0; channel < 2; ++channel) {
        s->residual[channel] = clip(s->input[channel] - s->lowpass[channel], -1131.2, 1131.2);
    }
}

// ----------------------------------------------------------------------------
// pll —— 锁相环：从三相电压里提取电网相位角、角速度和电压幅值
// ----------------------------------------------------------------------------
// 这是网侧和机侧控制共用的“基准生成器”：后面所有 Park 变换用的 θ 都来自它。
// 算法结构（正负序解耦 PLL，业界也常叫 DSOGI / 解耦网络类 PLL）：
//   1) 线电压 → Clarke 变换（stationary），并算正序电压幅值 voltage_magnitude；
//   2) 用当前角度 θ 做一次 Park 得正序分量（positive_raw），
//      用 -θ 再做一次得负序分量（negative_raw）；
//   3) 解耦：正序分量里减去负序滤波值旋转 2θ 后的投影（反之亦然），
//      把电压不平衡/谐波造成的 2 倍频纹波互相抵消；
//   4) 解耦结果过一阶低通（系数 0.03665，即每步保留 96.3% 的历史值）；
//   5) 正序 q 轴分量作为 PI 调节器的输入（reference = q，feedback = 0），
//      目标是把 q 轴分量压到 0 —— q = 0 意味着锁相角度与电压矢量完全重合；
//   6) 角速度 = PI 输出 + 2π*f_base（额定角速度作为前馈，PI 只修偏差），
//      角度积分后归一化到 [0, 2π)；
//   7) 最后对 dq 电压做前馈滤波（系数 0.10472），供电流环解耦用。
// ----------------------------------------------------------------------------
void pll(double base_frequency, PllState* s) {
    const double tau = 6.28318530717958;
    s->stationary.vab = s->line_voltage[0];
    s->stationary.vbc = s->line_voltage[1];
    s->stationary.reserved = s->line_voltage[2];
    line_voltage_transform(&s->stationary);
    s->rotating.alpha = s->stationary.alpha;
    s->rotating.beta = s->stationary.beta;
    s->voltage_magnitude = std::sqrt(s->rotating.alpha * s->rotating.alpha
        + s->rotating.beta * s->rotating.beta);

    s->rotating.theta = s->angle;
    park_transform(&s->rotating);
    s->positive_raw[0] = s->rotating.d;
    s->positive_raw[1] = s->rotating.q;
    s->rotating.theta = -s->angle;
    park_transform(&s->rotating);
    s->negative_raw[0] = s->rotating.d;
    s->negative_raw[1] = s->rotating.q;

    const double cosine = std::cos(2.0 * s->angle);
    const double sine = std::sin(2.0 * s->angle);
    s->positive_decoupled[0] = (s->positive_raw[0] - s->negative_filtered[0] * cosine)
        - s->negative_filtered[1] * sine;
    s->positive_decoupled[1] = (s->positive_raw[1] - s->negative_filtered[1] * cosine)
        + s->negative_filtered[0] * sine;
    s->negative_decoupled[0] = (s->negative_raw[0] - s->positive_filtered[0] * cosine)
        + s->positive_filtered[1] * sine;
    s->negative_decoupled[1] = (s->negative_raw[1] - s->positive_filtered[1] * cosine)
        - s->positive_filtered[0] * sine;

    // All four decouplings consume the previous filter histories.
    for (int axis = 0; axis < 2; ++axis) {
        s->positive_filtered[axis] = s->positive_decoupled[axis] * 0.03665191429188089
            + s->positive_filtered[axis] * 0.9633480857081191;
    }
    for (int axis = 0; axis < 2; ++axis) {
        s->negative_filtered[axis] = s->negative_decoupled[axis] * 0.03665191429188089
            + s->negative_filtered[axis] * 0.9633480857081191;
    }
    s->negative_magnitude = std::sqrt(s->negative_filtered[0] * s->negative_filtered[0]
        + s->negative_filtered[1] * s->negative_filtered[1]);
    s->regulator.reference = s->positive_decoupled[1];
    s->regulator.feedback = 0.0;
    regulator(&s->regulator);
    s->angular_speed = s->regulator.output + tau * base_frequency;
    s->angle = s->angle + s->angular_speed * s->step;
    if (s->angle > tau) s->angle = s->angle - tau;
    else if (s->angle < 0.0) s->angle = s->angle + tau;

    s->rotating.theta = s->angle;
    park_transform(&s->rotating);
    s->feedforward[0] = s->rotating.d * 0.10471975511965965
        + s->feedforward[0] * 0.8952802448803403;
    s->feedforward[1] = s->rotating.q * 0.10471975511965965
        + s->feedforward[1] * 0.8952802448803403;
}

// ----------------------------------------------------------------------------
// modulation —— 调制模块：电压参考值（αβ）→ 三相开关管占空比
// ----------------------------------------------------------------------------
// 输入是控制环给出的 αβ 电压参考（标幺），输出是三相上/下管占空比 duty[3][2]。
// 步骤：
//   1) 由 αβ 反推三相调制波 phase[0..2]（系数 2 是标幺基准的换算）；
//   2) 三相调制波限幅 ±1.1547（= 2/√3，正弦调制波线性区上限）；
//   3) 求三相最大 phase_max / 最小 phase_min / 中间值 middle；
//   4) 零序注入：每相都加 middle/2 —— 等价于注入三次谐波（零序分量），
//      把线性调制范围从 ±1 扩到 ±1.1547，直流电压利用率提高约 15%；
//   5) 中点电位平衡：delta = 正母线电压 - 负母线电压，超过死区 deadband 时
//      按 gain 比例产生 correction（极性 polarity 可反转），把中点电位拉回平衡；
//   6) 调制波 = 零序注入后的值 + correction，再限幅到 ±1；
//   7) 转占空比：调制波 ≥ 0 时上管占空比 = 调制波、下管 = 1；
//      调制波 < 0 时上管 = 0、下管 = 1 + 调制波。
// ----------------------------------------------------------------------------
void modulation(ModulationState* s) {
    s->phase[0] = 2.0 * s->alpha;
    s->phase[1] = 2.0 * ((-0.5 * s->alpha) + (0.866 * s->beta));
    s->phase[2] = 2.0 * ((-0.5 * s->alpha) - (0.866 * s->beta));
    for (int phase = 0; phase < 3; ++phase) {
        if (s->phase[phase] > 1.1547) s->phase[phase] = 1.1547;
        else if (s->phase[phase] <= -1.1547) s->phase[phase] = -1.1547;
    }
    s->phase_max = s->phase[0];
    if (s->phase[1] > s->phase_max) s->phase_max = s->phase[1];
    if (s->phase[2] > s->phase_max) s->phase_max = s->phase[2];
    s->phase_min = s->phase[0];
    if (s->phase_min > s->phase[1]) s->phase_min = s->phase[1];
    if (s->phase_min > s->phase[2]) s->phase_min = s->phase[2];
    s->middle = (((s->phase[0] + s->phase[1]) + s->phase[2]) - s->phase_max) - s->phase_min;
    s->delta = s->dc_plus - s->dc_minus;
    if (std::fabs(s->delta) > s->deadband) {
        const double signed_gain = s->polarity == 1 ? s->gain : -s->gain;
        s->correction = (signed_gain * s->delta) / (s->dc_plus + s->dc_minus);
        s->correction = clip(s->correction, -s->limit, s->limit);
    } else {
        s->correction = 0.0;
    }
    for (int phase = 0; phase < 3; ++phase) {
        s->modulation[phase] = (s->phase[phase] + s->middle / 2.0) + s->correction;
    }
    for (int phase = 0; phase < 3; ++phase) {
        s->modulation[phase] = clip(s->modulation[phase], -1.0, 1.0);
    }
    for (int phase = 0; phase < 3; ++phase) {
        if (s->modulation[phase] >= 0.0) {
            s->duty[phase][0] = s->modulation[phase];
            s->duty[phase][1] = 1.0;
        } else {
            s->duty[phase][0] = 0.0;
            s->duty[phase][1] = 1.0 + s->modulation[phase];
        }
    }
}

// ----------------------------------------------------------------------------
// 下面两个小工具：从状态区任意偏移量读出/写入一个 double。
// memcpy 把内存里的 8 个字节原样解释成浮点数（避免类型转换引入精度变化）。
// 后面的函数大量使用 d(0x...) / put(0x..., v) 这种写法：
// d(0x550) 表示“读状态区第 0x550 字节处的浮点数”。
// ----------------------------------------------------------------------------
template<typename T>
static T state_value(const void* state, std::size_t offset) {
    T value;
    std::memcpy(&value, static_cast<const unsigned char*>(state) + offset, sizeof(value));
    return value;
}

static void set_state_double(void* state, std::size_t offset, double value) {
    std::memcpy(static_cast<unsigned char*>(state) + offset, &value, sizeof(value));
}

// ----------------------------------------------------------------------------
// alternate_modulation —— 另一个调制实现（alternate = 备用的）。
// 与上面的 modulation 做的是同一件事，区别在于：
//   1) 完全用字节偏移量写（0x10、0x18……对应 ModulationState 的 alpha/beta/phase 等字段）；
//   2) 占空比分配上处理了“某相恰好等于最大值/最小值”的边界情形。
// 原 DLL 在机侧某些模式（0x6E8 == 1）下调用这个版本。
// ----------------------------------------------------------------------------
void alternate_modulation(void* state) {
    const auto d = [state](std::size_t offset) { return state_value<double>(state, offset); };
    const auto put = [state](std::size_t offset, double value) { set_state_double(state, offset, value); };
    put(0x10, 2.0 * d(0));
    put(0x18, 2.0 * ((-0.5 * d(0)) + (0.866 * d(8))));
    put(0x20, 2.0 * ((-0.5 * d(0)) - (0.866 * d(8))));
    for (std::size_t offset = 0x10; offset <= 0x20; offset += 8) {
        if (d(offset) > 1.1547) put(offset, 1.1547);
        else if (-1.1547 >= d(offset)) put(offset, -1.1547);
    }
    put(0x28, d(0x10));
    if (d(0x18) > d(0x28)) put(0x28, d(0x18));
    if (d(0x20) > d(0x28)) put(0x28, d(0x20));
    put(0x30, d(0x10));
    if (d(0x30) > d(0x18)) put(0x30, d(0x18));
    if (d(0x30) > d(0x20)) put(0x30, d(0x20));
    put(0x38, ((((d(0x10) + d(0x18)) + d(0x20)) - d(0x28)) - d(0x30)) / 2.0);
    put(0x88, d(0x10) + d(0x38));
    put(0x90, d(0x18) + d(0x38));
    put(0x98, d(0x20) + d(0x38));
    put(0x138, d(0x88));
    if (d(0x90) > d(0x138)) put(0x138, d(0x90));
    if (d(0x98) > d(0x138)) put(0x138, d(0x98));
    put(0x140, d(0x88));
    if (d(0x140) > d(0x90)) put(0x140, d(0x90));
    if (d(0x140) > d(0x98)) put(0x140, d(0x98));
    put(0x158, (((d(0x88) + d(0x90)) + d(0x98)) - d(0x138)) - d(0x140));
    put(0x148, (d(0x158) + d(0x138)) / 2.0);
    put(0x150, (d(0x158) - d(0x138)) / 2.0);
    for (std::size_t phase = 0; phase < 3; ++phase) {
        const double value = d(0x88 + phase * 8);
        const std::size_t duty = 0xF0 + phase * 16;
        if (value == d(0x138)) {
            put(duty, d(0x138));
            put(duty + 8, 1.0);
        } else if (value == d(0x140)) {
            put(duty, 0.0);
            put(duty + 8, 1.0 + d(0x140));
        } else {
            put(duty, d(0x148));
            put(duty + 8, d(0x150) + 1.0);
        }
    }
}

// status_times_100：把状态区里某个 32 位状态字（bit 位组合）
// 解释成“无符号数 × 100”再转成带符号整数输出给 PSCAD——
// 原 DLL 把多个开关状态打包成一个整数×100 送出来画波形，这里原样复现。
static double status_times_100(const void* state, std::size_t offset) {
    const std::uint32_t product = state_value<std::uint32_t>(state, offset) * 100u;
    std::int32_t signed_product;
    std::memcpy(&signed_product, &product, sizeof(product));
    return static_cast<double>(signed_product);
}

// ----------------------------------------------------------------------------
// grid_input —— 网侧输入映射。
// 把 PSCAD 每步送来的 8 个输入量（电压、电流、直流电压等）按原 DLL 的位置
// 写进状态区（0xC0/0xC8/... 是各自的门牌号）。
// 最后三行把两路线电压 vab/vbc 反推成三相电压存入 0x1718 起，供保护使用。
// ----------------------------------------------------------------------------
void grid_input(const double* input, void* state) {
    const std::size_t offsets[] = {0xC0, 0xC8, 0xD0, 0x100, 0x108, 0xD8, 0xE0, 0xE8};
    for (int index = 0; index < 8; ++index) set_state_double(state, offsets[index], input[index]);
    const double vab = state_value<double>(state, 0xC0);
    const double vbc = state_value<double>(state, 0xC8);
    set_state_double(state, 0x1718, (2.0 * vab + vbc) / 3.0);
    set_state_double(state, 0x1720, (-vab + vbc) / 3.0);
    set_state_double(state, 0x1728, (-vab - 2.0 * vbc) / 3.0);
}

// ----------------------------------------------------------------------------
// rotor_input —— 机侧（转子侧）输入映射 + 转差频率计算：
//   slip = speed * pole_pairs / 60 - base_frequency
// speed 是机械转速（rpm），pole_pairs 是极对数，除 60 换算到 Hz。
// 转差频率是双馈电机控制的关键量：转子电流的电气频率 = 转差频率。
// 末尾三行用“上一拍”的电网电压反推三相电压（原 DLL 的执行顺序，原样保留）。
// ----------------------------------------------------------------------------
void rotor_input(const double* input, double base_frequency, void* state) {
    const std::size_t offsets[] = {0x228, 0x230, 0x238, 0x240, 0x248, 0x250, 0x290, 0x6E0};
    for (int index = 0; index < 8; ++index) set_state_double(state, offsets[index], input[index]);
    const double speed = state_value<double>(state, 0x6E0);
    const double pole_pairs = state_value<double>(state, 0x1E0);
    set_state_double(state, 0x5F8, (speed * pole_pairs) / 60.0 - base_frequency);
    // The main routine copies this step's grid voltages only after this helper returns.
    const double previous_vab = state_value<double>(state, 0x258);
    const double previous_vbc = state_value<double>(state, 0x260);
    set_state_double(state, 0x19B8, (2.0 * previous_vab + previous_vbc) / 3.0);
    set_state_double(state, 0x19C0, (-previous_vab + previous_vbc) / 3.0);
    set_state_double(state, 0x19C8, (-previous_vab - 2.0 * previous_vbc) / 3.0);
}

// ----------------------------------------------------------------------------
// grid_output —— 网侧输出映射：把状态区里的结果打包成 PSCAD 期望的 20 个输出。
// case 0/1/2（由 0x628 处的模式字决定）对应原 DLL 的三种输出模式：
// 调制波输出、占空比输出、平均值输出。
// output[0] 是“使能”标志，其余信号乘以使能——停机时输出自动清零。
// ----------------------------------------------------------------------------
void grid_output(double* output, const void* state) {
    const auto d = [state](std::size_t offset) { return state_value<double>(state, offset); };
    const auto i = [state](std::size_t offset) { return state_value<std::int32_t>(state, offset); };
    const double enable = static_cast<double>(i(0x50C));
    output[0] = enable;
    switch (i(0x628)) {
    case 0:
        output[1] = d(0x4B8) * enable;
        output[2] = d(0x4D0) * enable;
        output[3] = d(0x4E8) * enable;
        output[4] = d(0x4A8) * enable;
        output[5] = d(0x4B0) * enable;
        output[6] = 0.0;
        break;
    case 1: {
        const std::size_t offsets[] = {0x4C0, 0x4C8, 0x4D8, 0x4E0, 0x4F0, 0x4F8};
        for (int index = 0; index < 6; ++index) output[index + 1] = d(offsets[index]) * enable;
        break;
    }
    case 2:
        output[1] = d(0xD88) * enable;
        output[2] = d(0xD90) * enable;
        output[3] = d(0xD98) * enable;
        output[4] = 0.0;
        output[5] = 0.0;
        output[6] = 0.0;
        break;
    }
    output[7] = static_cast<double>(i(0x504));
    output[8] = static_cast<double>(i(0x508));
    output[9] = static_cast<double>(i(0x510));
    output[10] = static_cast<double>(i(0x514));
    output[11] = status_times_100(state, 0x18);
    output[12] = d(0xF08) / 1.414;
    output[13] = d(0xF10) / 1.414;
    output[14] = d(0xF88) / 1.414;
    output[15] = d(0xF90) / 1.414;
    output[16] = status_times_100(state, 0xC60);
    output[17] = d(0x1218);
    output[18] = d(0x270);
    output[19] = d(0x278);
}

// ----------------------------------------------------------------------------
// rotor_output —— 机侧输出映射，与网侧对称。
// 注意 output[19] 这一行：原 DLL 的两个分支“都会被计算”（即使权重为 0），
// 这在浮点上可被观察到（可能产生 NaN/Inf），所以明码也原样保留这个行为。
// ----------------------------------------------------------------------------
void rotor_output(double* output, const void* state) {
    const auto d = [state](std::size_t offset) { return state_value<double>(state, offset); };
    const auto i = [state](std::size_t offset) { return state_value<std::int32_t>(state, offset); };
    const double enable = d(0x350);
    output[0] = enable;
    switch (i(0x6E8)) {
    case 0:
        output[1] = d(0x300) * enable;
        output[2] = d(0x318) * enable;
        output[3] = d(0x330) * enable;
        output[4] = d(0x2F0) * enable;
        output[5] = d(0x2F8) * enable;
        output[6] = 0.0;
        break;
    case 1: {
        const std::size_t offsets[] = {0x308, 0x310, 0x320, 0x328, 0x338, 0x340};
        for (int index = 0; index < 6; ++index) output[index + 1] = d(offsets[index]) * enable;
        break;
    }
    case 2:
        output[1] = d(0xC78) * enable;
        output[2] = d(0xC80) * enable;
        output[3] = d(0xC88) * enable;
        output[4] = 0.0;
        output[5] = 0.0;
        output[6] = 0.0;
        break;
    }
    output[7] = 0.0;
    output[8] = static_cast<double>(i(0x3A0));
    output[9] = static_cast<double>(i(0x3A4));
    output[10] = status_times_100(state, 0x10);
    output[11] = d(0xDF8) / 1.414;
    output[12] = d(0xE00) / 1.414;
    output[13] = d(0xE78) / 1.414;
    output[14] = d(0xE80) / 1.414;
    output[15] = static_cast<double>(i(0x3E0));
    output[16] = static_cast<double>(i(0x5D4) | i(0x5D8));
    output[17] = d(0x6C8);
    output[18] = (((d(0x6D0) / 1.732) / (d(0x6D8) / d(0x6D0))) / d(0x6C8)) * 1.414;
    // Both terms are evaluated even when their weight is zero.
    output[19] = d(0x4E0) * static_cast<double>(i(0x10) == 1)
        + (d(0x3A8) / d(0x410)) * static_cast<double>(i(0x10) != 1);
}

// ----------------------------------------------------------------------------
// average_converter_output —— 直流母线电流折算。
// 把两侧变流器的交流侧功率按当前调制波折算成直流母线电流，
// 供 PSCAD 里的“平均模型”（AVG 模式，不用开关管细节）使用。
// 原理：瞬时功率 = 相电压×相电流，按调制波 m 的符号把功率分到正/负母线。
// ----------------------------------------------------------------------------
void average_converter_output(void* state) {
    const auto d = [state](std::size_t offset) { return state_value<double>(state, offset); };
    const double half_dc = (d(0x100) + d(0x108)) / 2.0;
    for (int side = 0; side < 2; ++side) {
        const std::size_t modulation_base = side == 0 ? 0xD88 : 0x2730;
        const std::size_t current_base = side == 0 ? 0xD8 : 0x1CE0;
        const std::size_t phase_base = side == 0 ? 0x4B8 : 0x1DB8;
        const std::size_t dc_current_base = side == 0 ? 0x4A8 : 0x1DA8;
        double nonpositive[3];
        double positive[3];
        for (int phase = 0; phase < 3; ++phase) {
            const double m = d(modulation_base + 8 * phase);
            const double current = d(current_base + 8 * phase);
            if (m > 0.0) {
                nonpositive[phase] = 0.0;
                positive[phase] = m * current;
            } else {
                nonpositive[phase] = m * current;
                positive[phase] = 0.0;
            }
        }
        for (int phase = 0; phase < 3; ++phase) {
            set_state_double(state, phase_base + 24 * phase, half_dc * d(modulation_base + 8 * phase));
        }
        const double nonpositive_sum = (nonpositive[0] + nonpositive[1]) + nonpositive[2];
        const double positive_sum = (positive[0] + positive[1]) + positive[2];
        set_state_double(state, dc_current_base, side == 0 ? -nonpositive_sum : nonpositive_sum);
        set_state_double(state, dc_current_base + 8, side == 0 ? -positive_sum : positive_sum);
    }
}

// view_output —— 监控输出：把 41 个内部变量打包给 PSCAD 画波形/监视用。
// 第 29 个（index 28）是整数状态量，其余都是浮点数。
void view_output(double* output, const void* state) {
    const std::size_t offsets[] = {
        0xB18, 0x468, 0x908, 0x910, 0x918, 0x928, 0x930, 0x938,
        0x948, 0x950, 0x958, 0x988, 0x990, 0x998, 0xA80, 0xA88,
        0xA98, 0xAA0, 0xAB0, 0xAB8, 0xAC8, 0xAD0, 0xAE0, 0xAE8,
        0xAF8, 0xB00, 0xB10, 0xB20, 0x7DC, 0x968, 0x970, 0x978,
        0x9A8, 0x9B0, 0x9B8, 0x9C8, 0x9D0, 0x9D8, 0x9E8, 0x9F0, 0x9F8
    };
    for (int index = 0; index < 41; ++index) {
        if (index == 28) output[index] = static_cast<double>(state_value<std::int32_t>(state, offsets[index]));
        else output[index] = state_value<double>(state, offsets[index]);
    }
}

// rotor_nameplate —— 电机铭牌参数换算。
// 把 PSCAD 传入的额定功率、额定电压、额定电流、极对数等
// 换算成内部使用的基准值（如 0x6D0/0x6D8 乘 1000 是 kW→W 的换算）。
void rotor_nameplate(const double* parameters, void* state) {
    set_state_double(state, 0x598, parameters[0]);
    set_state_double(state, 0x6D0, parameters[1] * 1000.0);
    set_state_double(state, 0x6D8, parameters[2] * 1000.0);
    set_state_double(state, 0x6C8, parameters[3]);
    set_state_double(state, 0x1E0, parameters[4]);
    set_state_double(state, 0x5F0, state_value<double>(state, 0x6D8) / state_value<double>(state, 0x6D0));
}

// ----------------------------------------------------------------------------
// adaptive_oscillator —— 自适应振荡器：一个“跟随输入信号”的振荡器，
// 用于从每相电压里提取基波。
// 核心：内部振荡器（幅值 0x68、角速度 0x18、相位 0x38）不断与输入比对：
//   - 幅值按输入幅值自适应（下限 5.0）；
//   - 角速度 = 内部频率修正 0x20 + 2π*f_base（额定角速度打底）；
//   - 相位积分并归一化到 [0, 2π)。
// 输出 0x48/0x50 是与电网同步的基波两相正交分量。
// ----------------------------------------------------------------------------
void adaptive_oscillator(void* state, double base_frequency) {
    const auto d = [state](std::size_t offset) { return state_value<double>(state, offset); };
    const auto put = [state](std::size_t offset, double value) { set_state_double(state, offset, value); };
    const double tau = 6.28318530717958;
    put(0x68, std::sqrt(d(0x48) * d(0x48) + d(0x50) * d(0x50)));
    if (5.0 > d(0x68)) put(0x68, 5.0);
    double q = (-d(0)) * 20.0;
    q = q / d(0x68);
    q = q * (d(0x40) - d(0x48));
    q = q * d(0x50);
    q = q * d(8);
    put(0x20, d(0x20) + q);
    put(0x20, clip(d(0x20), d(0x30), d(0x28)));
    put(0x18, d(0x20) + tau * base_frequency);
    put(0x10, d(0x40) - d(0x48));
    put(0x10, d(0) * d(0x10) - d(0x50));
    put(0x48, d(0x48) + (d(0x18) * d(8)) * d(0x10));
    put(0x50, d(0x50) + (d(0x18) * d(8)) * d(0x48));
    put(0x38, d(0x38) + d(0x18) * d(8));
    if (d(0x38) > tau) put(0x38, d(0x38) - tau);
    else if (0.0 > d(0x38)) put(0x38, d(0x38) + tau);
}

// ----------------------------------------------------------------------------
// current_sequence —— 电流正负序分解（网侧版）。
// 流程：三相电流 → Clarke → 用 +θ 做一次 Park 得正序 dq（0x68/0x70），
// 用 -θ 做一次 Park 得负序 dq（0x78/0x80）→ 解耦（消除 2 倍频纹波）→ 低通。
// 输出：0x28..0x40 为解耦后的正/负序 dq；0x48..0x60 为低通滤波后的值
// （系数 0.0314，即时间常数较短的滤波）；0x88 起为另一组滤波（系数 0.5115）。
// 正序 dq 给电流环做反馈，负序 dq 用于负序电流抑制。
// ----------------------------------------------------------------------------
void current_sequence(void* state) {
    const auto d = [state](std::size_t offset) { return state_value<double>(state, offset); };
    const auto put = [state](std::size_t offset, double value) { set_state_double(state, offset, value); };
    put(0xA8, d(0));
    put(0xB0, d(8));
    put(0xB8, d(0x10));
    CurrentTransformState stationary = state_value<CurrentTransformState>(state, 0xA8);
    current_transform(&stationary);
    put(0xC0, stationary.alpha);
    put(0xC8, stationary.beta);
    put(0xD0, d(0xC0));
    put(0xD8, d(0xC8));
    put(0xE0, d(0x18));
    ParkState rotating = state_value<ParkState>(state, 0xD0);
    park_transform(&rotating);
    put(0xE8, rotating.d);
    put(0xF0, rotating.q);
    put(0x68, d(0xE8));
    put(0x70, d(0xF0));
    put(0xE0, d(0x20));
    rotating.theta = d(0xE0);
    park_transform(&rotating);
    put(0xE8, rotating.d);
    put(0xF0, rotating.q);
    put(0x78, d(0xE8));
    put(0x80, d(0xF0));

    double t = d(0x68) - d(0x58) * std::cos(d(0x18) - d(0x20));
    put(0x28, t - d(0x60) * std::sin(d(0x18) - d(0x20)));
    t = d(0x70) - d(0x60) * std::cos(d(0x18) - d(0x20));
    put(0x30, t + d(0x58) * std::sin(d(0x18) - d(0x20)));
    t = d(0x78) - d(0x48) * std::cos(d(0x18) - d(0x20));
    put(0x38, t + d(0x50) * std::sin(d(0x18) - d(0x20)));
    t = d(0x80) - d(0x50) * std::cos(d(0x18) - d(0x20));
    put(0x40, t - d(0x48) * std::sin(d(0x18) - d(0x20)));
    const std::size_t raw[] = {0x28, 0x30, 0x38, 0x40};
    const std::size_t filtered[] = {0x48, 0x50, 0x58, 0x60};
    for (int axis = 0; axis < 4; ++axis) {
        put(filtered[axis], d(raw[axis]) * 0.0314159265358979 + d(filtered[axis]) * 0.9685840734641021);
    }
    for (int axis = 0; axis < 4; ++axis) {
        const std::size_t offset = 0x88 + 8 * axis;
        put(offset, d(raw[axis]) * 0.5115273563044781 + d(offset) * 0.4884726436955219);
    }
}

// ----------------------------------------------------------------------------
// rotor_current_sequence —— 电流正负序分解（机侧版）。
// 算法与 current_sequence 完全一样，只是滤波系数不同：
// 0x48 起的低通系数 0.0785（时间常数更短），0x88 起的滤波系数 0.5237。
// ----------------------------------------------------------------------------
void rotor_current_sequence(void* state) {
    const auto d = [state](std::size_t offset) { return state_value<double>(state, offset); };
    const auto put = [state](std::size_t offset, double value) { set_state_double(state, offset, value); };
    put(0xA8, d(0));
    put(0xB0, d(8));
    put(0xB8, d(0x10));
    CurrentTransformState stationary = state_value<CurrentTransformState>(state, 0xA8);
    current_transform(&stationary);
    put(0xC0, stationary.alpha);
    put(0xC8, stationary.beta);
    put(0xD0, d(0xC0));
    put(0xD8, d(0xC8));
    put(0xE0, d(0x18));
    ParkState rotating = state_value<ParkState>(state, 0xD0);
    park_transform(&rotating);
    put(0xE8, rotating.d);
    put(0xF0, rotating.q);
    put(0x68, d(0xE8));
    put(0x70, d(0xF0));
    put(0xE0, d(0x20));
    rotating.theta = d(0xE0);
    park_transform(&rotating);
    put(0xE8, rotating.d);
    put(0xF0, rotating.q);
    put(0x78, d(0xE8));
    put(0x80, d(0xF0));

    double t = d(0x68) - d(0x58) * std::cos(d(0x18) - d(0x20));
    put(0x28, t - d(0x60) * std::sin(d(0x18) - d(0x20)));
    t = d(0x70) - d(0x60) * std::cos(d(0x18) - d(0x20));
    put(0x30, t + d(0x58) * std::sin(d(0x18) - d(0x20)));
    t = d(0x78) - d(0x48) * std::cos(d(0x18) - d(0x20));
    put(0x38, t + d(0x50) * std::sin(d(0x18) - d(0x20)));
    t = d(0x80) - d(0x50) * std::cos(d(0x18) - d(0x20));
    put(0x40, t - d(0x48) * std::sin(d(0x18) - d(0x20)));
    const std::size_t raw[] = {0x28, 0x30, 0x38, 0x40};
    const std::size_t filtered[] = {0x48, 0x50, 0x58, 0x60};
    for (int axis = 0; axis < 4; ++axis) {
        put(filtered[axis], d(raw[axis]) * 0.07853981633974474 + d(filtered[axis]) * 0.9214601836602553);
    }
    for (int axis = 0; axis < 4; ++axis) {
        const std::size_t offset = 0x88 + 8 * axis;
        put(offset, d(raw[axis]) * 0.523709146364097 + d(offset) * 0.47629085363590296);
    }
}

// ----------------------------------------------------------------------------
// resonant_filter —— 谐振滤波器（二阶数字滤波器，用于谐振控制器）。
// 每步按当前频率 d(0) 在线重算分母 0x90 和三个系数 0x88/0x98/0xA0，
// 然后更新滤波器状态（0x28 输入误差、0x40 输出、0x48/0x50 历史值）。
// 这种“系数随频率在线更新”的结构用于在指定频率处产生高增益，
// 作为谐振控制器抑制特定次谐波（本工程里谐振点由初始化给出）。
// ----------------------------------------------------------------------------
void resonant_filter(void* state) {
    const auto d = [state](std::size_t offset) { return state_value<double>(state, offset); };
    const auto put = [state](std::size_t offset, double value) { set_state_double(state, offset, value); };
    put(0x90, (4.0 + (4.0 * d(0x10)) * d(0)) + ((d(8) * d(8)) * d(0)) * d(0));
    put(0x88, (((4.0 * d(0x80)) * d(0x10)) * d(0)) / d(0x90));
    put(0x98, (((((2.0 * d(8)) * d(8)) * d(0)) * d(0)) - 8.0) / d(0x90));
    put(0xA0, ((4.0 - (4.0 * d(0x10)) * d(0)) + ((d(8) * d(8)) * d(0)) * d(0)) / d(0x90));
    put(0x28, d(0x18) - d(0x20));
    put(0x40, (d(0x88) * (d(0x28) - d(0x38)) - d(0x98) * d(0x48)) - d(0xA0) * d(0x50));
    put(0x38, d(0x30));
    put(0x30, d(0x28));
    put(0x50, d(0x48));
    put(0x48, d(0x40));
    put(0x58, d(0x70) * d(0x28));
    put(0x60, d(0x60) + (d(0x78) * d(0)) * d(0x28));
    put(0x68, -d(0x40));
    put(0x68, clip(d(0x68), d(0xB0), d(0xA8)));
}

// difference_filter —— 差分滤波器：提取信号的变化量（近似微分），
// 用于电压突变检测与动态补偿。
void difference_filter(void* state) {
    const auto d = [state](std::size_t offset) { return state_value<double>(state, offset); };
    double u = (((d(8) * d(0x10)) * d(0x18)) * d(0x20)) / (d(0) * (d(8) + d(0x10)));
    double v = (((d(8) * d(0x10)) * d(0x18)) * d(0x28)) / (d(0) * (d(8) + d(0x10)));
    u = u - v;
    v = ((d(8) - d(0x10)) * d(0x38)) / (d(8) + d(0x10));
    set_state_double(state, 0x30, u - v);
    set_state_double(state, 0x38, d(0x30));
    set_state_double(state, 0x28, d(0x20));
}
