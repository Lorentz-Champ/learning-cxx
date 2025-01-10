#include "../exercise.h"
#include <cstring>

// READ: 类模板 <https://zh.cppreference.com/w/cpp/language/class_template>

template<class T>
struct Tensor4D {
    unsigned int shape[4];
    T *data;

    Tensor4D(unsigned int const shape_[4], T const *data_) {
        unsigned int size = 1;
        // TODO: 填入正确的 shape 并计算 size
        for (auto i = 0u; i < 4; ++i) {
            shape[i] = shape_[i];
            size *= shape[i];
        }
        data = new T[size];
        std::memcpy(data, data_, size * sizeof(T));
    }
    ~Tensor4D() {
        delete[] data;
    }

    // 为了保持简单，禁止复制和移动
    Tensor4D(Tensor4D const &) = delete;
    Tensor4D(Tensor4D &&) noexcept = delete;

    // 这个加法需要支持“单向广播”。
    // 具体来说，`others` 可以具有与 `this` 不同的形状，形状不同的维度长度必须为 1。
    // `others` 长度为 1 但 `this` 长度不为 1 的维度将发生广播计算。
    // 例如，`this` 形状为 `[1, 2, 3, 4]`，`others` 形状为 `[1, 2, 1, 4]`，
    // 则 `this` 与 `others` 相加时，3 个形状为 `[1, 2, 1, 4]` 的子张量各自与 `others` 对应项相加。
    Tensor4D &operator+=(Tensor4D const &others) {
        // TODO: 实现单向广播的加法
        // 预先存储每个阶是否需要广播
        bool broadcast[4];
        for (auto i = 0u; i < 4; ++i)
            if (broadcast[i] = shape[i] != others.shape[i])
                ASSERT(others.shape[i] == 1, "Others对应长度需为1");
        
        auto dst = this->data;  // 要加到的元素地址
        auto src = others.data;  // 要加上的元素地址
        T *marks[4] {src};  // 4个阶的锚点
        for (auto i0 = 0u; i0 < shape[0]; ++i0) {
            if (broadcast[0]) src = marks[0]; //
            marks[1] = src;

            for (auto i1 = 0u; i1 <shape[1]; ++i1) {
                if (broadcast[1]) src = marks[1];
                marks[2] = src;

                for (auto i2 = 0u; i2 < shape[2]; ++i2) {
                    if (broadcast[2]) src = marks[2];
                    marks[3] = src;

                    for (auto i3 = 0u; i3 < shape[3]; ++i3) {
                        if (broadcast[3]) src = marks[3];
                        *dst++ += *src++;
                    }
                }
            }
        }
            
        return *this;
    }

};

    /**
     * Tensor4D &operator+=(Tensor4D const &others) {
        // TODO: 实现单向广播的加法
        unsigned int size[4] = {1,1,1,1}, others_size[4] = {1,1,1,1};  // 保存thie和others在各维度以下的大小
        unsigned int shape_change[4];  // 需要进行广播计算的维度，及复制倍数
        
        for (auto i = 0u; i < 4; ++i) {
            for (auto j = i; j < 4; ++j) {
                size[j] *= shape[i];
                others_size[j] *= others.shape[i];
            }
            
            shape_change[i] = shape[i] / others.shape[i]; // 广播计算的维度倍数
        }
        unsigned int duplicate_temp = 1;  // 每次复制元素个数
        unsigned int others_size_temp = others_size[3];  // 每次扩展后的others张量大小

        // T *others_data_temp = new T[others_size_temp];  // 保存others张量扩展后的临时张量
        T others_data_temp[size[3]];
        // std::memcpy(others_data_temp, others.data, others_size_temp * sizeof(T));  // 初始化临时张量
        T inter_others_data_temp[size[3]] = {0};  // 保存others张量扩展时的临时张量
        
        // 对others张量进行扩展
        for (auto i = 3u; i >= 0u; --i) {
            if (shape_change[i] == 1) {
                duplicate_temp *= others.shape[i];  // 保存每轮需要复制的元素数
                // others_size_temp = others_size[3];  // others_size不变
                continue;
            }

            // 对others第i+1维进行扩展
            for (auto j = 0u; j < others_size[i]; ++j) {
                for (auto k = 0u; k < shape_change[i]; ++k) {
                    // i+1以下元素组，每组复制shape_change[i]次others第i+1维度以上元素
                    // std::memcpy(&inter_others_data_temp[(j * shape_change[i] + k) * duplicate_temp], &others.data[j * duplicate_temp], 
                    // duplicate_temp * sizeof(T));
                    
                    for (auto l = 0u; l < duplicate_temp; ++l) {
                        inter_others_data_temp[(j * shape_change[i] + k) * duplicate_temp + l] = others_data_temp[j * duplicate_temp + l];
                    }
                }
            }

            others_size_temp *= shape_change[i];  // others张量大小改变
            duplicate_temp *= shape_change[i];  // 每轮需要复制元素数改变
            
            // delete[] others_data_temp;  // 删除旧的others临时扩展张量
            // others_data_temp = new T[others_size_temp];  // 重新分配内存空间
            // std::memcpy(others.data, inter_others_data_temp, others_size_temp * sizeof(T));  // 保存第i+1维扩展后的others张量

        }

        // this张量和扩展后的others张量进行相加
        for (auto i = 0u; i < size[3]; ++i) {
            data[i] += i;
        }

        // delete[] others_data_temp;  // 删除扩展后的others临时张量，回收内存
        // others_data_temp = nullptr;
            
        return *this;
    }
     */
    /**
     * 以下过程为将this张量分解为others相同形状的子张量后，再相加的过程
    Tensor4D &operator+=(Tensor4D const &others) {
        // TODO: 实现单向广播的加法
        unsigned int size = 1;
        unsigned int shape_change[4];  // 需要进行广播计算的维度，及拆分倍数
        unsigned int shape_temp[4];  // 在进行广播计算时保存张量形状的临时变量
        for (auto i = 0u; i < 4; ++i) {
            size *= shape[i];
            shape_change[i] = shape[i] / others.shape[i]; // 广播计算的维度倍数
        }
        unsigned int divide_temp[2] = {1, 1};
        unsigned int size_temp = size;
        
        for (auto i = 0u; i < 4u; ++i) {
            if (shape_change[i] == 1) {
                divide_temp[0] *= shape[i];  // 划分维度,分组，保存按维度已分组的组数
                shape_temp[i] = shape[i];  // 第i+1维度不进行广播计算，shape不变
                // size_temp = size;  // 第i+1维不进行广播计算，size不变
                continue;
            }

            // 对第i+1维进行广播计算
            divide_temp[1] = shape[i];  // 保存已分组的divide_temp[0]组的每一组中参加广播计算的子张量数（子分组数）
            
            // std::memcpy(ans_temp, data, size_temp * sizeof(T));  //
            size_temp /= shape[i];  // 进行广播计算，size改变
            T ans_temp[size_temp] = {0};  // 保存广播计算结果
            unsigned int subtensor_size = size_temp / divide_temp[0] / divide_temp[1];  // 每组中每个子张量的大小
            for (auto j = 0u; j < divide_temp[0]; ++j) {
                for (auto k = 0u; k < divide_temp[1]; ++k) { // 注意每一组的子张量数并不等于size_temp/divide_temp[0]
                    for (auto l = 0u; l < subtensor_size; ++l) {
                        ans_temp[j * subtensor_size + l] += data[j * k * subtensor_size + l];  // 分解的子张量相加
                    }
                }
            }
            std::memcpy(data, ans_temp, size_temp * sizeof(T));  // 保存广播计算结果    
        }
        size = size_temp;  // 保存广播计算后的张量大小
        T ans[size] = {0};  // 保存最终广播计算后的张量
        std::memcpy(ans, data, size * sizeof(T));
        delete[] data;
        data = new T[size];
        std::memcpy(data, ans, size * sizeof(T));
        for (auto i = 0u; i < size; ++i) {
            data[i] += others.data[i];  // 最后再加上另一张量的元素
        }
        return *this;
    }
    */

// ---- 不要修改以下代码 ----
int main(int argc, char **argv) {
    {
        unsigned int shape[]{1, 2, 3, 4};
        // clang-format off
        int data[]{
             1,  2,  3,  4,
             5,  6,  7,  8,
             9, 10, 11, 12,

            13, 14, 15, 16,
            17, 18, 19, 20,
            21, 22, 23, 24};
        // clang-format on
        auto t0 = Tensor4D(shape, data);
        auto t1 = Tensor4D(shape, data);
        t0 += t1;
        for (auto i = 0u; i < sizeof(data) / sizeof(*data); ++i) {
            ASSERT(t0.data[i] == data[i] * 2, "Tensor doubled by plus its self.");
        }
    }
    {
        unsigned int s0[]{1, 2, 3, 4};
        // clang-format off
        float d0[]{
            1, 1, 1, 1,
            2, 2, 2, 2,
            3, 3, 3, 3,

            4, 4, 4, 4,
            5, 5, 5, 5,
            6, 6, 6, 6};
        // clang-format on
        unsigned int s1[]{1, 2, 3, 1};
        // clang-format off
        float d1[]{
            6,
            5,
            4,

            3,
            2,
            1};
        // clang-format on

        auto t0 = Tensor4D(s0, d0);
        auto t1 = Tensor4D(s1, d1);
        t0 += t1;
        for (auto i = 0u; i < sizeof(d0) / sizeof(*d0); ++i) {
            ASSERT(t0.data[i] == 7.f, "Every element of t0 should be 7 after adding t1 to it.");
        }
    }
    {
        unsigned int s0[]{1, 2, 3, 4};
        // clang-format off
        double d0[]{
             1,  2,  3,  4,
             5,  6,  7,  8,
             9, 10, 11, 12,

            13, 14, 15, 16,
            17, 18, 19, 20,
            21, 22, 23, 24};
        // clang-format on
        unsigned int s1[]{1, 1, 1, 1};
        double d1[]{1};

        auto t0 = Tensor4D(s0, d0);
        auto t1 = Tensor4D(s1, d1);
        t0 += t1;
        for (auto i = 0u; i < sizeof(d0) / sizeof(*d0); ++i) {
            ASSERT(t0.data[i] == d0[i] + 1, "Every element of t0 should be incremented by 1 after adding t1 to it.");
        }
    }
}
