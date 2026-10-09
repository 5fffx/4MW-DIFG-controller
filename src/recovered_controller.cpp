// ============================================================================
// 文件：recovered_controller.cpp —— 控制器的“总调度”（顶层入口）
// ----------------------------------------------------------------------------
// PSCAD 每步只调用下面三个函数：
//   dfig_controller_sy_wg4000kdf_npc_init ：启动时初始化一次；
//   dfig_controller_sy_wg4000kdf_npc     ：每个仿真步调用一次（核心）；
//   dfig_controller_sy_wg4000kdf_npc_end ：仿真结束时清理。
// 总调度负责：把输入分发给两侧、按采样周期决定是否执行控制、
// 把两侧输出打包返回给 PSCAD。
// ============================================================================
#include "recovered_primitives.h"
#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <emmintrin.h>
#include <map>

namespace {

// ----------------------------------------------------------------------------
// ControllerState：管理一整块状态区（0x3748 字节）。
// 和前面文件里的帮助类一样提供 d(offset)/i(offset) 读写；
// 另外 snapshot() 把整块状态区存成文件/从文件恢复——
// 这是原 DLL 自带的“断点快照”机制（PSCAD 的 snapshot 功能）。
// ----------------------------------------------------------------------------
class ControllerState {
public:
    void* at(std::size_t offset = 0) { return bytes_.data() + offset; }
    double d(std::size_t offset) const { return read<double>(offset); }
    std::int32_t i(std::size_t offset) const { return read<std::int32_t>(offset); }
    void d(std::size_t offset, double value) { write(offset, value); }
    void i(std::size_t offset, std::int32_t value) { write(offset, value); }
    void increment(std::size_t offset) { write(offset, read<std::uint32_t>(offset) + std::uint32_t{1}); }

    void clear() { bytes_.fill(0); }

    void snapshot(bool restore) {
        std::FILE* file = nullptr;
        if (fopen_s(&file, "snapshot.txt", restore ? "rb" : "wb") != 0 || file == nullptr) {
            std::abort();
        }
        if (restore) {
            std::fread(bytes_.data(), bytes_.size(), 1, file);
        } else {
            std::fwrite(bytes_.data(), bytes_.size(), 1, file);
        }
        std::fclose(file);
    }

private:
    template <typename T>
    T read(std::size_t offset) const {
        T value;
        std::memcpy(&value, bytes_.data() + offset, sizeof(value));
        return value;
    }

    template <typename T>
    void write(std::size_t offset, T value) {
        std::memcpy(bytes_.data() + offset, &value, sizeof(value));
    }

    std::array<unsigned char, 0x3748> bytes_{};
};

// The original array capacity is unknown; storage remains stable per instance.
// instances：按“实例号”管理多份状态。
// 一台仿真里可能放多个控制器（多个机组），每个实例号对应一块独立状态区。
// 原 DLL 内部同样是按实例号分开存储的，这里用 map 等价实现。
std::map<std::int32_t, ControllerState> instances;

std::int32_t run_number(double value) {
    return _mm_cvttsd_si32(_mm_set_sd(value));
}

} // namespace

// ----------------------------------------------------------------------------
// init：初始化入口。PSCAD 传参数数组（第 66 个是额定频率）和实例号，
// 我们取出该实例的状态区，调用 controller_initialize 填默认值。
// ----------------------------------------------------------------------------
extern "C" __declspec(dllexport) void dfig_controller_sy_wg4000kdf_npc_init(
    const double* parameters, const std::int32_t* instance) {
    ControllerState& s = instances[*instance];
    controller_initialize(parameters[65], s.at());
}

extern "C" __declspec(dllexport) void dfig_controller_sy_wg4000kdf_npc_end(
    const std::int32_t* instance) {
    instances[*instance].clear();
}

// ----------------------------------------------------------------------------
// 主入口：每个仿真步（10µs）被 PSCAD 调用一次。
// 参数：grid/motor = 网侧/机侧输入数组；parameters = 参数数组；
// protection/commands = 保护定值与命令；grid_result/motor_result = 两侧输出；
// view_result = 监控输出；delta = 本步时间增量；sample_period = 控制周期。
// ----------------------------------------------------------------------------
extern "C" __declspec(dllexport) void dfig_controller_sy_wg4000kdf_npc(
    const double* grid, const double* motor, const double* parameters,
    const double* protection, const double* commands,
    double* grid_result, double* motor_result, double* view_result,
    const double* delta, const double* sample_period, const std::int32_t* instance) {
    ControllerState& s = instances[*instance];
// 运行号检测：PSCAD 每重新运行一次仿真，grid[9] 里的运行号就变一次，
// 检测到变化就重新初始化控制器（除非机侧在复位过程中）。
    const double frequency = parameters[65];
    s.i(0x3730, run_number(grid[9]));
    if (s.i(0x3734) != s.i(0x3730)) {
        s.i(0x3740, 0);
        if (motor[9] != 1000.0) {
            controller_initialize(frequency, s.at());
            s.i(0x3730, run_number(grid[9]));
        }
    }
    s.i(0x3734, s.i(0x3730));
    if (s.i(0x3740) == 0 && motor[9] == 1000.0) {
        s.snapshot(true);
        s.i(0x3740, 1);
        s.i(0x3730, run_number(grid[9]));
        s.i(0x3734, s.i(0x3730));
    }

// ----------------------------------------------------------------------------
// 采样周期闸门：仿真步长 10µs，但控制算法每 250µs 才执行一次。
// 0x3728 是“距上次控制的累计时间”，只有它归零的那一步才执行下面的控制。
// ----------------------------------------------------------------------------
    if (s.d(0x3728) == 0.0) {
// 一个控制周期内的完整调用序列（顺序即原 DLL 的顺序）：
//   1) grid_input / rotor_input      —— 把 PSCAD 输入搬进状态区；
//   2) rotor_nameplate / parameter_input / protection_input / command_input
//                                    —— 铭牌、参数、保护定值、命令；
//   3) grid_calculate / motor_calculate —— 两侧主控制算法；
//   4) grid_post / motor_post        —— 两侧状态机与保护；
//   5) average_converter_output      —— 直流母线电流折算；
//   6) grid_output / rotor_output / view_output —— 打包输出给 PSCAD。
// 中间还有频率支撑（0x1E08 使能时）和 0x1D58 += 0x768（调频功率叠加到有功给定）。
        void* rotor = s.at(0x1AB8);
        grid_input(grid, s.at());
        rotor_input(motor, frequency, rotor);
        rotor_nameplate(parameters, rotor);
        parameter_input(parameters + 10, frequency, s.at());
        protection_input(protection, frequency, s.at());
        command_input(commands, s.at());
        s.d(0x1D10, s.d(0xC0));
        s.d(0x1D18, s.d(0xC8));
        s.d(0x1D20, s.d(0xD0));
        s.d(0x1D38, s.d(0x100));
        s.d(0x1D40, s.d(0x108));
        s.d(0x538, s.d(0x2050));
        s.i(0xBD0, s.i(0x1AD4));
        s.d(0x7C0, s.d(0x1D58));
// 保护计数器：0x3738 使能时，若两侧有保护标志，计数 8000 步（2 秒）
// 后置 0x520/0x1FF4 = 1，触发两侧复位。
        if (s.i(0x3738) == 1 && (s.i(0x1FB8) != 0 || s.i(0x354) != 0)) {
            s.increment(0x373C);
            if (s.i(0x373C) >= 8000) {
                s.i(0x520, 1);
                s.i(0x1FF4, 1);
            }
        } else {
            s.i(0x373C, 0);
        }
        grid_calculate(frequency, s.at());
        if (s.d(0x1E08) == 1.0) {
            grid_frequency_support(frequency, s.at());
        }
        s.d(0x1D58, s.d(0x1D58) + s.d(0x768));
        s.d(0x25E0, s.d(0x468));
        motor_calculate(frequency, rotor);
        grid_post(frequency, s.at());
        motor_post(frequency, rotor);
        s.i(0x34C, s.i(0x1FB4) == 3 ? 1 : 0);
        s.i(0x1FB0, s.i(0x350) == 3 ? 1 : 0);
        average_converter_output(s.at());
        grid_output(grid_result, s.at());
        rotor_output(motor_result, rotor);
        view_output(view_result, s.at());
    }

// 采样计时器推进：累计本步时间 delta，到满一个采样周期就清零，
// 下一个 250µs 窗口再执行一次控制。
// 末尾：当两侧“运行号”相等且机侧不处于复位过程时，把状态区存成快照。
    s.d(0x3728, s.d(0x3728) + *delta);
    if (s.d(0x3728) >= *sample_period) {
        s.d(0x3728, 0.0);
    }
    if (motor[8] == grid[8] && motor[9] != 1000.0) {
        s.i(0x3740, 0);
        s.snapshot(false);
    }
}
