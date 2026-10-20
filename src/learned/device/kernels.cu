typedef unsigned long long u64;

#define MAX_LANES 32
#define MAX_MATES 32
#define NEGATIVE_INFINITY __int_as_float(0xff800000)

__device__ __forceinline__ u64 mix(u64 value) {
    u64 mixed = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    mixed = (mixed ^ (mixed >> 27)) * 0x94d049bb133111ebULL;
    return mixed ^ (mixed >> 31);
}

// The host's polynomial exp step for step, so a row of the loss reads the same on both.
__device__ __forceinline__ float polynomial_exp(float x) {
    x = fminf(fmaxf(x, -87.0f), 88.0f);
    float k = roundf(x * 1.44269504088896340736f);
    float r = x - k * 0.6933594f + k * 2.1219444e-4f;
    float p = 1.9875691e-4f;
    p = p * r + 1.3981999e-3f;
    p = p * r + 8.333452e-3f;
    p = p * r + 4.1665796e-2f;
    p = p * r + 1.6666665e-1f;
    p = p * r + 0.5f;
    float y = p * r * r + r + 1.0f;
    return y * __int_as_float(((int)k + 127) << 23);
}

__device__ __forceinline__ long long flat_index() {
    return (long long)blockIdx.x * blockDim.x + threadIdx.x;
}

// Each block of rows is summed in row order and the blocks are folded in order, as on the host.
extern "C" __global__ void block_sums(const float* x, const float* y, const double* mean,
                                      int rows, int cols, int block, int mode, double* partial) {
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    if (col >= cols) return;
    int start = blockIdx.y * block;
    int end = min(rows, start + block);
    double sum = 0.0;
    for (int row = start; row < end; ++row) {
        long long at = (long long)row * cols + col;
        float value = x[at];
        if (mode == 0) {
            sum += (double)value;
        } else if (mode == 1) {
            double gap = (double)value - mean[col];
            sum += gap * gap;
        } else {
            sum += (double)(value * y[at]);
        }
    }
    partial[(long long)blockIdx.y * cols + col] = sum;
}

extern "C" __global__ void fold(const double* partial, int blocks, int cols, double divisor,
                                double* out64, float* out32) {
    int col = flat_index();
    if (col >= cols) return;
    double total = 0.0;
    for (int at = 0; at < blocks; ++at) total += partial[(long long)at * cols + col];
    if (out64) out64[col] = total / divisor;
    if (out32) out32[col] = (float)total;
}

extern "C" __global__ void finish(const double* partial, int blocks, int cols, int rows,
                                  const double* mean, float momentum, double epsilon,
                                  float* running_mean, float* running_variance, float* inverse,
                                  float* mean32) {
    int col = flat_index();
    if (col >= cols) return;
    double total = 0.0;
    for (int at = 0; at < blocks; ++at) total += partial[(long long)at * cols + col];
    double variance = total / (double)rows;
    double unbiased = (double)rows / (double)(max(rows, 2) - 1);
    running_mean[col] = (1.0f - momentum) * running_mean[col] + momentum * (float)mean[col];
    running_variance[col] =
        (1.0f - momentum) * running_variance[col] + momentum * (float)(variance * unbiased);
    inverse[col] = (float)(1.0 / sqrt(variance + epsilon));
    mean32[col] = (float)mean[col];
}

extern "C" __global__ void normalise(float* table, float* active, const float* mean,
                                     const float* inverse, const float* scale, const float* shift,
                                     int rows, int cols, int dropping, u64 bar, u64 key,
                                     float kept) {
    long long at = flat_index();
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

extern "C" __global__ void mask(float* gradient, const float* active, long long count,
                                float kept) {
    long long at = flat_index();
    if (at >= count) return;
    gradient[at] = active[at] > 0.0f ? gradient[at] * kept : 0.0f;
}

extern "C" __global__ void rewrite(float* gradient, const float* normed, const float* scale,
                                   const float* inverse, const float* shift, const float* scaled,
                                   int rows, int cols) {
    long long at = flat_index();
    if (at >= (long long)rows * cols) return;
    int col = at % cols;
    float count = (float)rows;
    float factor = scale[col] * inverse[col] / count;
    gradient[at] = factor * (count * gradient[at] - shift[col] - normed[at] * scaled[col]);
}

extern "C" __global__ void broadcast(float* table, const float* row, int rows, int cols) {
    long long at = flat_index();
    if (at >= (long long)rows * cols) return;
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
    long long at = flat_index();
    if (at >= (long long)rows * cols) return;
    int row = at / cols;
    int col = at % cols;
    joined[at] = col < split ? features[(long long)row * width + col]
                             : depth[(long long)row * depth_cols + col - split] * lift;
}

extern "C" __global__ void scaled_copy(const float* from, float* to, long long count, float by) {
    long long at = flat_index();
    if (at >= count) return;
    to[at] = from[at] * by;
}

extern "C" __global__ void adamw(float* params, const float* grads, float* first, float* second,
                                 long long count, float shrink, float step, float momentum,
                                 float keep, float correction, float epsilon) {
    long long at = flat_index();
    if (at >= count) return;
    float grad = grads[at];
    float param = params[at] * shrink;
    float moment = momentum * first[at] + (1.0f - momentum) * grad;
    float square = keep * second[at] + (1.0f - keep) * grad * grad;
    first[at] = moment;
    second[at] = square;
    params[at] = param - step * moment / (sqrtf(square) / correction + epsilon);
}

// One block per row. The negatives are summed in the host's lanes so the totals agree.
extern "C" __global__ void per_positive(float* weights, int n, const int* group_start,
                                        const int* group_size, const int* members, int lanes,
                                        double* row_loss) {
    __shared__ float peaks[256];
    __shared__ float lane_sums[MAX_LANES];
    __shared__ double own[MAX_MATES];
    __shared__ int own_at[MAX_MATES];
    __shared__ float shares[MAX_MATES];
    __shared__ float spread_shared;
    __shared__ int own_count;
    int row = blockIdx.x;
    int thread = threadIdx.x;
    float* values = weights + (long long)row * n;
    const int* mates = members + group_start[row];
    int count = group_size[row];
    if (thread == 0) values[row] = NEGATIVE_INFINITY;
    __syncthreads();
    float peak = NEGATIVE_INFINITY;
    for (int col = thread; col < n; col += blockDim.x) peak = fmaxf(peak, values[col]);
    peaks[thread] = peak;
    __syncthreads();
    for (int half = blockDim.x / 2; half > 0; half >>= 1) {
        if (thread < half) peaks[thread] = fmaxf(peaks[thread], peaks[thread + half]);
        __syncthreads();
    }
    peak = peaks[0];
    if (thread == 0) {
        int kept = 0;
        for (int at = 0; at < count; ++at) {
            int mate = mates[at];
            if (mate == row) continue;
            own_at[kept] = mate;
            own[kept] = (double)values[mate];
            ++kept;
        }
        own_count = kept;
    }
    __syncthreads();
    for (int col = thread; col < n; col += blockDim.x) values[col] = polynomial_exp(values[col] - peak);
    __syncthreads();
    for (int at = thread; at < count; at += blockDim.x) values[mates[at]] = 0.0f;
    __syncthreads();
    int full = n / lanes * lanes;
    if (thread < lanes) {
        float sum = 0.0f;
        for (int col = thread; col < full; col += lanes) sum += values[col];
        lane_sums[thread] = sum;
    }
    __syncthreads();
    if (thread == 0) {
        double tail = 0.0;
        for (int col = full; col < n; ++col) tail += (double)values[col];
        double total = 0.0;
        for (int lane = 0; lane < lanes; ++lane) total += (double)lane_sums[lane];
        double log_negatives = (double)peak + log(total + tail);
        double loss = 0.0;
        double spread = 0.0;
        for (int at = 0; at < own_count; ++at) {
            double mine = own[at];
            double high = fmax(mine, log_negatives);
            double log_total = high + log(exp(mine - high) + exp(log_negatives - high));
            loss += log_total - mine;
            spread += exp((double)peak - log_total);
            shares[at] = (float)(exp(mine - log_total) - 1.0);
        }
        row_loss[row] = loss;
        spread_shared = (float)spread;
    }
    __syncthreads();
    float spread = spread_shared;
    for (int col = thread; col < n; col += blockDim.x) values[col] *= spread;
    __syncthreads();
    for (int at = thread; at < own_count; at += blockDim.x) values[own_at[at]] = shares[at];
}
