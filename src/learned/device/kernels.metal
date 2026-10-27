#include <metal_stdlib>
using namespace metal;

#define THREADS 256
#define MOST_MATES 32
#define MASKED -3.0e38f

static inline ulong mix(ulong value) {
    ulong mixed = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9UL;
    mixed = (mixed ^ (mixed >> 27)) * 0x94d049bb133111ebUL;
    return mixed ^ (mixed >> 31);
}

// Metal has no doubles, so sums the host keeps in f64 carry their rounding in a second float.
static inline void add(thread float2& sum, float value) {
    float total = sum.x + value;
    sum.y += fabs(sum.x) >= fabs(value) ? (sum.x - total) + value : (value - total) + sum.x;
    sum.x = total;
}

static inline void add(thread float2& sum, float2 other) {
    add(sum, other.x);
    sum.y += other.y;
}

#define INDEX uint2 group [[threadgroup_position_in_grid]], uint2 local [[thread_position_in_threadgroup]]
#define FLAT (int)(group.x * THREADS + local.x)

kernel void block_sums(device const float* x [[buffer(0)]], device const float* y [[buffer(1)]],
                       device const float* mean [[buffer(2)]], constant int& rows [[buffer(3)]],
                       constant int& cols [[buffer(4)]], constant int& block [[buffer(5)]],
                       constant int& mode [[buffer(6)]], device float2* partial [[buffer(7)]],
                       INDEX) {
    int col = FLAT;
    if (col >= cols) return;
    int start = group.y * block;
    int end = min(rows, start + block);
    float2 sum = 0.0f;
    for (int row = start; row < end; ++row) {
        long at = (long)row * cols + col;
        float value = x[at];
        if (mode == 0) {
            add(sum, value);
        } else if (mode == 1) {
            float gap = value - mean[col];
            add(sum, gap * gap);
        } else {
            add(sum, value * y[at]);
        }
    }
    partial[(long)group.y * cols + col] = sum;
}

kernel void fold(device const float2* partial [[buffer(0)]], constant int& blocks [[buffer(1)]],
                 constant int& cols [[buffer(2)]], constant float& divisor [[buffer(3)]],
                 device float* out [[buffer(4)]], INDEX) {
    int col = FLAT;
    if (col >= cols) return;
    float2 total = 0.0f;
    for (int at = 0; at < blocks; ++at) add(total, partial[(long)at * cols + col]);
    out[col] = (total.x + total.y) / divisor;
}

kernel void finish(device const float2* partial [[buffer(0)]], constant int& blocks [[buffer(1)]],
                   constant int& cols [[buffer(2)]], constant int& rows [[buffer(3)]],
                   device const float* mean [[buffer(4)]], constant float& momentum [[buffer(5)]],
                   constant float& epsilon [[buffer(6)]], device float* running_mean [[buffer(7)]],
                   device float* running_variance [[buffer(8)]], device float* inverse [[buffer(9)]],
                   INDEX) {
    int col = FLAT;
    if (col >= cols) return;
    float2 total = 0.0f;
    for (int at = 0; at < blocks; ++at) add(total, partial[(long)at * cols + col]);
    float variance = (total.x + total.y) / rows;
    float unbiased = (float)rows / (float)(max(rows, 2) - 1);
    running_mean[col] = (1.0f - momentum) * running_mean[col] + momentum * mean[col];
    running_variance[col] =
        (1.0f - momentum) * running_variance[col] + momentum * variance * unbiased;
    inverse[col] = rsqrt(variance + epsilon);
}

kernel void normalise(device float* table [[buffer(0)]], device float* active [[buffer(1)]],
                      device const float* mean [[buffer(2)]], device const float* inverse [[buffer(3)]],
                      device const float* scale [[buffer(4)]], device const float* shift [[buffer(5)]],
                      constant int& rows [[buffer(6)]], constant int& cols [[buffer(7)]],
                      constant int& dropping [[buffer(8)]], constant ulong& bar [[buffer(9)]],
                      constant ulong& key [[buffer(10)]], constant float& kept [[buffer(11)]],
                      INDEX) {
    int at = FLAT;
    if (at >= rows * cols) return;
    int row = at / cols;
    int col = at % cols;
    float value = (table[at] - mean[col]) * inverse[col];
    table[at] = value;
    float relu = max(value * scale[col] + shift[col], 0.0f);
    if (dropping) {
        ulong seed = mix(key ^ (ulong)row);
        active[at] = mix(seed ^ (ulong)col) < bar ? 0.0f : relu * kept;
    } else {
        active[at] = relu;
    }
}

kernel void mask(device float* gradient [[buffer(0)]], device const float* active [[buffer(1)]],
                 constant int& count [[buffer(2)]], constant float& kept [[buffer(3)]], INDEX) {
    int at = FLAT;
    if (at >= count) return;
    gradient[at] = active[at] > 0.0f ? gradient[at] * kept : 0.0f;
}

kernel void rewrite(device float* gradient [[buffer(0)]], device const float* normed [[buffer(1)]],
                    device const float* scale [[buffer(2)]], device const float* inverse [[buffer(3)]],
                    device const float* shift [[buffer(4)]], device const float* scaled [[buffer(5)]],
                    constant int& rows [[buffer(6)]], constant int& cols [[buffer(7)]], INDEX) {
    int at = FLAT;
    if (at >= rows * cols) return;
    int col = at % cols;
    float count = (float)rows;
    float factor = scale[col] * inverse[col] / count;
    gradient[at] = factor * (count * gradient[at] - shift[col] - normed[at] * scaled[col]);
}

kernel void broadcast(device float* table [[buffer(0)]], device const float* row [[buffer(1)]],
                      constant int& rows [[buffer(2)]], constant int& cols [[buffer(3)]], INDEX) {
    int at = FLAT;
    if (at >= rows * cols) return;
    table[at] = row[at % cols];
}

kernel void unit(device float* table [[buffer(0)]], device float* lengths [[buffer(1)]],
                 constant int& rows [[buffer(2)]], constant int& cols [[buffer(3)]], INDEX) {
    int row = FLAT;
    if (row >= rows) return;
    device float* values = table + (long)row * cols;
    float sum = 0.0f;
    for (int col = 0; col < cols; ++col) sum += values[col] * values[col];
    float length = max(sqrt(sum), 1e-12f);
    for (int col = 0; col < cols; ++col) values[col] /= length;
    lengths[row] = length;
}

kernel void unit_backward(device float* gradient [[buffer(0)]], device const float* output [[buffer(1)]],
                          device const float* lengths [[buffer(2)]], constant int& rows [[buffer(3)]],
                          constant int& cols [[buffer(4)]], INDEX) {
    int row = FLAT;
    if (row >= rows) return;
    device float* values = gradient + (long)row * cols;
    device const float* unit = output + (long)row * cols;
    float along = 0.0f;
    for (int col = 0; col < cols; ++col) along += values[col] * unit[col];
    for (int col = 0; col < cols; ++col) values[col] = (values[col] - unit[col] * along) / lengths[row];
}

kernel void join(device const float* features [[buffer(0)]], constant int& width [[buffer(1)]],
                 constant int& split [[buffer(2)]], device const float* depth [[buffer(3)]],
                 constant int& depth_cols [[buffer(4)]], constant float& lift [[buffer(5)]],
                 device float* joined [[buffer(6)]], constant int& rows [[buffer(7)]], INDEX) {
    int cols = split + depth_cols;
    int at = FLAT;
    if (at >= rows * cols) return;
    int row = at / cols;
    int col = at % cols;
    joined[at] = col < split ? features[row * width + col]
                             : depth[row * depth_cols + col - split] * lift;
}

kernel void scaled_copy(device const float* from [[buffer(0)]], device float* to [[buffer(1)]],
                        constant int& count [[buffer(2)]], constant float& by [[buffer(3)]], INDEX) {
    int at = FLAT;
    if (at >= count) return;
    to[at] = from[at] * by;
}

kernel void adamw(device float* params [[buffer(0)]], device const float* grads [[buffer(1)]],
                  device float* first [[buffer(2)]], device float* second [[buffer(3)]],
                  constant int& count [[buffer(4)]], constant float& shrink [[buffer(5)]],
                  constant float& step [[buffer(6)]], constant float& momentum [[buffer(7)]],
                  constant float& keep [[buffer(8)]], constant float& correction [[buffer(9)]],
                  constant float& epsilon [[buffer(10)]], INDEX) {
    int at = FLAT;
    if (at >= count) return;
    float grad = grads[at];
    float moment = momentum * first[at] + (1.0f - momentum) * grad;
    float square = keep * second[at] + (1.0f - keep) * grad * grad;
    first[at] = moment;
    second[at] = square;
    params[at] = params[at] * shrink - step * moment / (sqrt(square) / correction + epsilon);
}

static float block_sum(threadgroup float* scratch, float value, uint lane) {
    scratch[lane] = value;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint span = THREADS / 2; span > 0; span >>= 1) {
        if (lane < span) scratch[lane] += scratch[lane + span];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    float total = scratch[0];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    return total;
}

// Rows are unit length, so no weight passes the temperature's inverse and it serves as every row's
// peak.
static inline float lifted(float value, float peak) {
    return value == MASKED ? 0.0f : exp(value - peak);
}

// One threadgroup per row. Each positive pair is scored against the negatives alone, and the row
// leaves scaled so the gradient products need no factor of their own.
kernel void per_positive(device float* weights [[buffer(0)]], constant int& n [[buffer(1)]],
                         device const int* group_start [[buffer(2)]],
                         device const int* group_size [[buffer(3)]],
                         device const int* members [[buffer(4)]], constant float& peak [[buffer(5)]],
                         constant float& scale [[buffer(6)]], device float* row_loss [[buffer(7)]],
                         INDEX) {
    threadgroup float scratch[THREADS];
    threadgroup float own[MOST_MATES];
    threadgroup int own_at[MOST_MATES];
    threadgroup float shares[MOST_MATES];
    threadgroup float spread_shared;
    threadgroup int own_count;
    int row = group.x;
    uint lane = local.x;
    device float* values = weights + (long)row * n;
    device const int* mates = members + group_start[row];
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
    threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
    for (int at = lane; at < count; at += THREADS) values[mates[at]] = MASKED;
    threadgroup_barrier(mem_flags::mem_device);
    float sum = 0.0f;
    for (int col = lane; col < n; col += THREADS) sum += lifted(values[col], peak);
    float total = block_sum(scratch, sum, lane);
    if (lane == 0) {
        float log_negatives = peak + log(total);
        float loss = 0.0f;
        float spread = 0.0f;
        for (int at = 0; at < own_count; ++at) {
            float mine = own[at];
            float high = max(mine, log_negatives);
            float log_total = high + log(exp(mine - high) + exp(log_negatives - high));
            loss += log_total - mine;
            spread += exp(peak - log_total);
            shares[at] = (exp(mine - log_total) - 1.0f) * scale;
        }
        row_loss[row] = loss;
        spread_shared = spread * scale;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float spread = spread_shared;
    for (int col = lane; col < n; col += THREADS) values[col] = lifted(values[col], peak) * spread;
    threadgroup_barrier(mem_flags::mem_device);
    for (int at = lane; at < own_count; at += THREADS) values[own_at[at]] = shares[at];
}

// One threadgroup per row. Every other row with the label is a positive (Khosla et al. 2020), and
// the row leaves scaled so the gradient products need no factor of their own.
kernel void contrastive(device float* weights [[buffer(0)]], constant int& n [[buffer(1)]],
                        device const int* group_start [[buffer(2)]],
                        device const int* group_size [[buffer(3)]],
                        device const int* members [[buffer(4)]], constant float& peak [[buffer(5)]],
                        constant float& scale [[buffer(6)]], device float* row_loss [[buffer(7)]],
                        INDEX) {
    threadgroup float scratch[THREADS];
    int row = group.x;
    uint lane = local.x;
    device float* values = weights + (long)row * n;
    device const int* mates = members + group_start[row];
    int count = group_size[row];
    if (count < 2) {
        for (int col = lane; col < n; col += THREADS) values[col] = 0.0f;
        if (lane == 0) row_loss[row] = 0.0f;
        return;
    }
    float mated = 0.0f;
    for (int at = lane; at < count; at += THREADS) {
        int mate = mates[at];
        if (mate != row) mated += values[mate];
    }
    mated = block_sum(scratch, mated, lane);
    if (lane == 0) values[row] = MASKED;
    threadgroup_barrier(mem_flags::mem_device);
    float sum = 0.0f;
    for (int col = lane; col < n; col += THREADS) sum += lifted(values[col], peak);
    float total = block_sum(scratch, sum, lane);
    float share = scale / total;
    for (int col = lane; col < n; col += THREADS) values[col] = lifted(values[col], peak) * share;
    threadgroup_barrier(mem_flags::mem_device);
    float positives = (float)(count - 1);
    float pull = scale / positives;
    for (int at = lane; at < count; at += THREADS) {
        int mate = mates[at];
        if (mate != row) values[mate] -= pull;
    }
    if (lane == 0) row_loss[row] = peak + log(total) - mated / positives;
}

// MPS ran the depth tower's two-column products under 1 GFLOP/s, so each thread sums one cell.
kernel void dot(device const float* a [[buffer(0)]], constant int& a_rows [[buffer(1)]],
                constant int& a_cols [[buffer(2)]], device const float* b [[buffer(3)]],
                constant int& b_rows [[buffer(4)]], constant int& b_cols [[buffer(5)]],
                device float* c [[buffer(6)]], constant int& c_rows [[buffer(7)]],
                constant int& m [[buffer(8)]], constant int& n [[buffer(9)]],
                constant int& k [[buffer(10)]], constant float& alpha [[buffer(11)]],
                constant float& beta [[buffer(12)]], INDEX) {
    int at = FLAT;
    if (at >= m * n) return;
    int i = at / n;
    int j = at % n;
    float2 total = 0.0f;
    for (int p = 0; p < k; ++p) add(total, a[i * a_rows + p * a_cols] * b[p * b_rows + j * b_cols]);
    float sum = total.x + total.y;
    long out = (long)i * c_rows + j;
    c[out] = beta == 0.0f ? alpha * sum : alpha * sum + beta * c[out];
}

#define NEAREST_SLOTS 4

static inline bool ahead(float d, uint at, float e, uint bt) {
    return d < e || (d == e && at < bt);
}

static inline uint ballot(bool flag) {
    return (uint)(simd_vote::vote_t)simd_ballot(flag);
}

// A SIMD group per query holds its k nearest in registers, slot s in lane s % 32, so an insert is
// a few shuffles and not a walk through memory. Ties go to the lower index as the host lists do.
kernel void nearest(device const float* dots [[buffer(0)]], constant int& rows [[buffer(1)]],
                    constant int& cols [[buffer(2)]], constant int& first [[buffer(3)]],
                    device const uint* selves [[buffer(4)]], device const float* query_variance [[buffer(5)]],
                    device const float* base_variance [[buffer(6)]], constant int& k [[buffer(7)]],
                    device float* distances [[buffer(8)]], device uint* indices [[buffer(9)]], INDEX) {
    int row = FLAT / 32;
    uint lane = local.x & 31;
    if (row >= rows) return;
    device float* best = distances + (long)row * k;
    device uint* held = indices + (long)row * k;
    float near[NEAREST_SLOTS];
    uint near_at[NEAREST_SLOTS];
    for (int r = 0; r < NEAREST_SLOTS; ++r) {
        int slot = r * 32 + lane;
        near[r] = slot < k ? best[slot] : INFINITY;
        near_at[r] = slot < k ? held[slot] : 0xffffffffu;
    }
    int last_reg = (k - 1) / 32;
    ushort last_lane = (k - 1) % 32;
    float worst = 0.0f;
    uint worst_at = 0;
    for (int r = 0; r < NEAREST_SLOTS; ++r) {
        float d = simd_shuffle(near[r], last_lane);
        uint a = simd_shuffle(near_at[r], last_lane);
        if (r == last_reg) { worst = d; worst_at = a; }
    }
    float mine = query_variance[row];
    uint own = selves[row];
    device const float* line = dots + (long)row * cols;
    for (int start = 0; start < cols; start += 32) {
        int col = start + lane;
        uint at = (uint)(first + col);
        bool live = col < cols && at != own;
        float distance = INFINITY;
        if (live) {
            float total = mine + base_variance[at];
            distance = 0.0f;
            if (total != 0.0f) distance = (total - 2.0f * line[col]) / total;
            if (isnan(distance)) distance = 2.0f;
        }
        uint wanted = ballot(live && ahead(distance, at, worst, worst_at));
        while (wanted) {
            ushort from = ctz(wanted);
            wanted &= wanted - 1;
            float d = simd_shuffle(distance, from);
            uint a = simd_shuffle(at, from);
            if (!ahead(d, a, worst, worst_at)) continue;
            int place = 0;
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                bool before = r * 32 + (int)lane < k && ahead(near[r], near_at[r], d, a);
                place += popcount(ballot(before));
            }
            float up[NEAREST_SLOTS];
            uint up_at[NEAREST_SLOTS];
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                up[r] = simd_shuffle_up(near[r], 1);
                up_at[r] = simd_shuffle_up(near_at[r], 1);
                float wrapped = simd_shuffle(r > 0 ? near[r - 1] : 0.0f, 31);
                uint wrapped_at = simd_shuffle(r > 0 ? near_at[r - 1] : 0u, 31);
                if (lane == 0) { up[r] = wrapped; up_at[r] = wrapped_at; }
            }
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                int slot = r * 32 + lane;
                if (slot > place) { near[r] = up[r]; near_at[r] = up_at[r]; }
                else if (slot == place) { near[r] = d; near_at[r] = a; }
            }
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                float e = simd_shuffle(near[r], last_lane);
                uint b = simd_shuffle(near_at[r], last_lane);
                if (r == last_reg) { worst = e; worst_at = b; }
            }
        }
    }
    for (int r = 0; r < NEAREST_SLOTS; ++r) {
        int slot = r * 32 + lane;
        if (slot < k) { best[slot] = near[r]; held[slot] = near_at[r]; }
    }
}
