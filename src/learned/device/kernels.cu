typedef unsigned long long u64;

#define THREADS 256
#define MOST_MATES 32
#define MASKED -3.0e38f

__device__ __forceinline__ u64 mix(u64 value) {
    u64 mixed = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    mixed = (mixed ^ (mixed >> 27)) * 0x94d049bb133111ebULL;
    return mixed ^ (mixed >> 31);
}

__device__ __forceinline__ int flat_index() {
    return blockIdx.x * THREADS + threadIdx.x;
}

extern "C" __global__ void block_sums(const float* x, const float* y, const float* mean, int rows,
                                      int cols, int block, int mode, double* partial) {
    int col = flat_index();
    if (col >= cols) return;
    int start = blockIdx.y * block;
    int end = min(rows, start + block);
    double sum = 0.0;
    for (int row = start; row < end; ++row) {
        long long at = (long long)row * cols + col;
        float value = x[at];
        if (mode == 0) {
            sum += value;
        } else if (mode == 1) {
            double gap = (double)value - mean[col];
            sum += gap * gap;
        } else {
            sum += value * y[at];
        }
    }
    partial[(long long)blockIdx.y * cols + col] = sum;
}

extern "C" __global__ void fold(const double* partial, int blocks, int cols, float divisor,
                                float* out) {
    int col = flat_index();
    if (col >= cols) return;
    double total = 0.0;
    for (int at = 0; at < blocks; ++at) total += partial[(long long)at * cols + col];
    out[col] = (float)(total / divisor);
}

extern "C" __global__ void finish(const double* partial, int blocks, int cols, int rows,
                                  const float* mean, float momentum, float epsilon,
                                  float* running_mean, float* running_variance, float* inverse) {
    int col = flat_index();
    if (col >= cols) return;
    double total = 0.0;
    for (int at = 0; at < blocks; ++at) total += partial[(long long)at * cols + col];
    double variance = total / rows;
    double unbiased = (double)rows / (max(rows, 2) - 1);
    running_mean[col] = (1.0f - momentum) * running_mean[col] + momentum * mean[col];
    running_variance[col] =
        (1.0f - momentum) * running_variance[col] + momentum * (float)(variance * unbiased);
    // The float argument serves Metal. CUDA adds the CPU's own double epsilon, set at compile.
    inverse[col] = (float)(1.0 / sqrt(variance + NORM_EPSILON));
}

extern "C" __global__ void normalise(float* table, float* active, const float* mean,
                                     const float* inverse, const float* scale, const float* shift,
                                     int rows, int cols, int dropping, u64 bar, u64 key,
                                     float kept) {
    long long at = (long long)blockIdx.x * THREADS + threadIdx.x;
    if (at >= (long long)rows * cols) return;
    int row = at / cols;
    int col = at % cols;
    float value = (table[at] - mean[col]) * inverse[col];
    table[at] = value;
    float relu = fmaxf(value * scale[col] + shift[col], 0.0f);
    if (dropping) {
        u64 seed = mix(key ^ (u64)row);
        active[at] = mix(seed ^ (u64)col) < bar ? 0.0f : relu * kept;
    } else {
        active[at] = relu;
    }
}

extern "C" __global__ void mask(float* gradient, const float* active, int count, float kept) {
    int at = flat_index();
    if (at >= count) return;
    gradient[at] = active[at] > 0.0f ? gradient[at] * kept : 0.0f;
}

extern "C" __global__ void rewrite(float* gradient, const float* normed, const float* scale,
                                   const float* inverse, const float* shift, const float* scaled,
                                   int rows, int cols) {
    int at = flat_index();
    if (at >= rows * cols) return;
    int col = at % cols;
    float count = (float)rows;
    float factor = scale[col] * inverse[col] / count;
    gradient[at] = factor * (count * gradient[at] - shift[col] - normed[at] * scaled[col]);
}

extern "C" __global__ void broadcast(float* table, const float* row, int rows, int cols) {
    int at = flat_index();
    if (at >= rows * cols) return;
    table[at] = row[at % cols];
}

extern "C" __global__ void unit(float* table, float* lengths, int rows, int cols) {
    int row = flat_index();
    if (row >= rows) return;
    float* values = table + (long long)row * cols;
    float sum = 0.0f;
    for (int col = 0; col < cols; ++col) sum += values[col] * values[col];
    float length = fmaxf(sqrtf(sum), 1e-12f);
    for (int col = 0; col < cols; ++col) values[col] /= length;
    lengths[row] = length;
}

extern "C" __global__ void unit_backward(float* gradient, const float* output,
                                         const float* lengths, int rows, int cols) {
    int row = flat_index();
    if (row >= rows) return;
    float* values = gradient + (long long)row * cols;
    const float* unit = output + (long long)row * cols;
    float along = 0.0f;
    for (int col = 0; col < cols; ++col) along += values[col] * unit[col];
    for (int col = 0; col < cols; ++col) values[col] = (values[col] - unit[col] * along) / lengths[row];
}

extern "C" __global__ void join(const float* features, int width, int split, const float* depth,
                                int depth_cols, float lift, float* joined, int rows) {
    int cols = split + depth_cols;
    int at = flat_index();
    if (at >= rows * cols) return;
    int row = at / cols;
    int col = at % cols;
    joined[at] = col < split ? features[row * width + col]
                             : depth[row * depth_cols + col - split] * lift;
}

extern "C" __global__ void scaled_copy(const float* from, float* to, int count, float by) {
    int at = flat_index();
    if (at >= count) return;
    to[at] = from[at] * by;
}

extern "C" __global__ void adamw(float* params, const float* grads, float* first, float* second,
                                 int count, float shrink, float step, float momentum, float keep,
                                 float correction, float epsilon) {
    int at = flat_index();
    if (at >= count) return;
    float grad = grads[at];
    float moment = momentum * first[at] + (1.0f - momentum) * grad;
    float square = keep * second[at] + (1.0f - keep) * grad * grad;
    first[at] = moment;
    second[at] = square;
    params[at] = params[at] * shrink - step * moment / (sqrtf(square) / correction + epsilon);
}

__device__ double block_sum(double* scratch, double value) {
    int lane = threadIdx.x;
    scratch[lane] = value;
    __syncthreads();
    for (int half = THREADS / 2; half > 0; half >>= 1) {
        if (lane < half) scratch[lane] += scratch[lane + half];
        __syncthreads();
    }
    double total = scratch[0];
    __syncthreads();
    return total;
}

// Rows are unit length, so no weight passes the temperature's inverse and it serves as every row's
// peak.
__device__ __forceinline__ float lifted(float value, float peak) {
    return value == MASKED ? 0.0f : expf(value - peak);
}

// One block per row. Each positive pair is scored against the negatives alone, and the row
// leaves scaled so the gradient products need no factor of their own.
extern "C" __global__ void per_positive(float* weights, int n, const int* group_start,
                                        const int* group_size, const int* members, float peak,
                                        float scale, double* row_loss) {
    __shared__ double scratch[THREADS];
    __shared__ float own[MOST_MATES];
    __shared__ int own_at[MOST_MATES];
    __shared__ float shares[MOST_MATES];
    __shared__ float spread_shared;
    __shared__ int own_count;
    int row = blockIdx.x;
    int lane = threadIdx.x;
    float* values = weights + (long long)row * n;
    const int* mates = members + group_start[row];
    int count = group_size[row];
    if (lane == 0) {
        int kept = 0;
        for (int at = 0; at < count; ++at) {
            int mate = mates[at];
            if (mate == row) continue;
            own_at[kept] = mate;
            own[kept] = values[mate];
            ++kept;
        }
        own_count = kept;
    }
    __syncthreads();
    for (int at = lane; at < count; at += THREADS) values[mates[at]] = MASKED;
    __syncthreads();
    double sum = 0.0;
    for (int col = lane; col < n; col += THREADS) sum += lifted(values[col], peak);
    double total = block_sum(scratch, sum);
    if (lane == 0) {
        double log_negatives = peak + log(total);
        double loss = 0.0;
        double spread = 0.0;
        for (int at = 0; at < own_count; ++at) {
            double mine = own[at];
            double high = fmax(mine, log_negatives);
            double log_total = high + log(exp(mine - high) + exp(log_negatives - high));
            loss += log_total - mine;
            spread += exp(peak - log_total);
            shares[at] = (float)((exp(mine - log_total) - 1.0) * scale);
        }
        row_loss[row] = loss;
        spread_shared = (float)(spread * scale);
    }
    __syncthreads();
    float spread = spread_shared;
    for (int col = lane; col < n; col += THREADS) values[col] = lifted(values[col], peak) * spread;
    __syncthreads();
    for (int at = lane; at < own_count; at += THREADS) values[own_at[at]] = shares[at];
}

// One block per row. Every other row with the label is a positive (Khosla et al. 2020), and the
// row leaves scaled so the gradient products need no factor of their own.
extern "C" __global__ void contrastive(float* weights, int n, const int* group_start,
                                       const int* group_size, const int* members, float peak,
                                       float scale, double* row_loss) {
    __shared__ double scratch[THREADS];
    int row = blockIdx.x;
    int lane = threadIdx.x;
    float* values = weights + (long long)row * n;
    const int* mates = members + group_start[row];
    int count = group_size[row];
    if (count < 2) {
        for (int col = lane; col < n; col += THREADS) values[col] = 0.0f;
        if (lane == 0) row_loss[row] = 0.0;
        return;
    }
    double mated = 0.0;
    for (int at = lane; at < count; at += THREADS) {
        int mate = mates[at];
        if (mate != row) mated += values[mate];
    }
    mated = block_sum(scratch, mated);
    if (lane == 0) values[row] = MASKED;
    __syncthreads();
    double sum = 0.0;
    for (int col = lane; col < n; col += THREADS) sum += lifted(values[col], peak);
    double total = block_sum(scratch, sum);
    float share = (float)(scale / total);
    for (int col = lane; col < n; col += THREADS) values[col] = lifted(values[col], peak) * share;
    __syncthreads();
    float positives = (float)(count - 1);
    float pull = scale / positives;
    for (int at = lane; at < count; at += THREADS) {
        int mate = mates[at];
        if (mate != row) values[mate] -= pull;
    }
    if (lane == 0) row_loss[row] = peak + log(total) - mated / (count - 1);
}

#define FULL 0xffffffffu
#define NEAREST_SLOTS 4

__device__ __forceinline__ bool ahead(float d, unsigned int at, float e, unsigned int bt) {
    return d < e || (d == e && at < bt);
}

// A warp per query holds its k nearest in registers, slot s in lane s % 32, so an insert is a few
// shuffles and not a walk through memory. Ties go to the lower index as the host lists do.
extern "C" __global__ void nearest(const float* dots, int rows, int cols, int first, const unsigned int* selves,
                                   const float* query_variance, const float* base_variance, int k,
                                   float* distances, unsigned int* indices) {
    int row = flat_index() / 32;
    int lane = threadIdx.x & 31;
    if (row >= rows) return;
    float* best = distances + (long long)row * k;
    unsigned int* held = indices + (long long)row * k;
    float near[NEAREST_SLOTS];
    unsigned int near_at[NEAREST_SLOTS];
#pragma unroll
    for (int r = 0; r < NEAREST_SLOTS; ++r) {
        int slot = r * 32 + lane;
        near[r] = slot < k ? best[slot] : __int_as_float(0x7f800000);
        near_at[r] = slot < k ? held[slot] : 0xffffffffu;
    }
    int last_reg = (k - 1) / 32, last_lane = (k - 1) % 32;
    float worst = 0.0f;
    unsigned int worst_at = 0;
#pragma unroll
    for (int r = 0; r < NEAREST_SLOTS; ++r) {
        float d = __shfl_sync(FULL, near[r], last_lane);
        unsigned int a = __shfl_sync(FULL, near_at[r], last_lane);
        if (r == last_reg) { worst = d; worst_at = a; }
    }
    float mine = query_variance[row];
    unsigned int own = selves[row];
    const float* line = dots + (long long)row * cols;
    for (int start = 0; start < cols; start += 32) {
        int col = start + lane;
        unsigned int at = (unsigned int)(first + col);
        bool live = col < cols && at != own;
        float distance = __int_as_float(0x7f800000);
        if (live) {
            float total = mine + base_variance[at];
            distance = 0.0f;
            if (total != 0.0f) distance = (total - 2.0f * line[col]) / total;
            if (isnan(distance)) distance = 2.0f;
        }
        unsigned int wanted = __ballot_sync(FULL, live && ahead(distance, at, worst, worst_at));
        while (wanted) {
            int from = __ffs(wanted) - 1;
            wanted &= wanted - 1;
            float d = __shfl_sync(FULL, distance, from);
            unsigned int a = __shfl_sync(FULL, at, from);
            if (!ahead(d, a, worst, worst_at)) continue;
            int place = 0;
#pragma unroll
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                bool before = r * 32 + lane < k && ahead(near[r], near_at[r], d, a);
                place += __popc(__ballot_sync(FULL, before));
            }
            float up[NEAREST_SLOTS];
            unsigned int up_at[NEAREST_SLOTS];
#pragma unroll
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                up[r] = __shfl_up_sync(FULL, near[r], 1);
                up_at[r] = __shfl_up_sync(FULL, near_at[r], 1);
                float wrapped = __shfl_sync(FULL, r > 0 ? near[r - 1] : 0.0f, 31);
                unsigned int wrapped_at = __shfl_sync(FULL, r > 0 ? near_at[r - 1] : 0u, 31);
                if (lane == 0) { up[r] = wrapped; up_at[r] = wrapped_at; }
            }
#pragma unroll
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                int slot = r * 32 + lane;
                if (slot > place) { near[r] = up[r]; near_at[r] = up_at[r]; }
                else if (slot == place) { near[r] = d; near_at[r] = a; }
            }
#pragma unroll
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                float e = __shfl_sync(FULL, near[r], last_lane);
                unsigned int b = __shfl_sync(FULL, near_at[r], last_lane);
                if (r == last_reg) { worst = e; worst_at = b; }
            }
        }
    }
#pragma unroll
    for (int r = 0; r < NEAREST_SLOTS; ++r) {
        int slot = r * 32 + lane;
        if (slot < k) { best[slot] = near[r]; held[slot] = near_at[r]; }
    }
}
