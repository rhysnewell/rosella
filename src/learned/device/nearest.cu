#define FULL 0xffffffffu
#ifndef COVERAGE
#define COVERAGE float
#endif
typedef COVERAGE real;

// The host's coverage constants, with its erfc table after them.
struct Terms {
    real weight, heaviest, epsilon, min_var_epsilon, low, high, steps, width;
};

struct Moment {
    real mean, variance, deviation, log_variance;
};

__device__ __forceinline__ real erfc_of(real z, const Terms* terms) {
    if (isnan(z)) return z;
    z = z < terms->low ? terms->low : (z > terms->high ? terms->high : z);
    real at = round((z - terms->low) * terms->steps);
    real offset = z - (terms->low + at / terms->steps);
    int width = (int)terms->width;
    const real* row = (const real*)(terms + 1) + (long long)at * width;
    real sum = row[width - 1];
    for (int i = width - 2; i >= 0; --i) sum = fma(sum, offset, row[i]);
    return sum;
}

__device__ __forceinline__ real normal_cdf(Moment m, real x, const Terms* terms) {
    real z = -(x - m.mean) / (m.deviation * (real)1.4142135623730951);
    return fmin((real)0.5 * erfc_of(z, terms), (real)1.0);
}

// The host measures one pair from either end, and both have to round alike.
__device__ real overlap(Moment a, Moment b, const Terms* terms) {
    if (b.variance < a.variance || (b.variance == a.variance && b.mean < a.mean)) {
        Moment swap = a;
        a = b;
        b = swap;
    }
    real k1, k2;
    if (fabs(a.variance - b.variance) < terms->min_var_epsilon) {
        k1 = k2 = (a.mean + b.mean) / 2;
    } else {
        real gap = (a.mean - b.mean) * (a.mean - b.mean)
            - (a.variance - b.variance) * (b.log_variance - a.log_variance);
        real tmp = sqrt(a.variance * b.variance * gap);
        k1 = (tmp - a.mean * b.variance + b.mean * a.variance) / (a.variance - b.variance);
        k2 = (tmp + a.mean * b.variance - b.mean * a.variance) / (b.variance - a.variance);
    }
    if (k1 > k2) {
        real swap = k1;
        k1 = k2;
        k2 = swap;
    }
    Moment narrow = a.variance > b.variance ? b : a;
    Moment wide = a.variance > b.variance ? a : b;
    if (k1 == k2) return fabs(normal_cdf(narrow, k1, terms) - normal_cdf(wide, k1, terms));
    return fabs(normal_cdf(narrow, k2, terms) - normal_cdf(narrow, k1, terms)
                + normal_cdf(wide, k1, terms) - normal_cdf(wide, k2, terms));
}

__device__ __forceinline__ float composition_only(float dot, float mine, float theirs) {
    float total = mine + theirs;
    float distance = 0.0f;
    if (total != 0.0f) distance = (total - 2.0f * dot) / total;
    return isnan(distance) ? 2.0f : distance;
}

// The composition floor turns a pair away before its coverage is read, as the host's bound does.
__device__ float mixed(float dot, float mine, float theirs, unsigned int a, unsigned int b, float worst,
                       const real* moments, const real* presence, int stride, int samples,
                       const Terms* terms) {
    real total = (real)mine + (real)theirs;
    real composition = total == 0 ? (real)0 : (total - 2 * (real)dot) / total;
    if (isnan(composition)) composition = 2;
    if (composition > 0) {
        float floor = (float)((1 - terms->heaviest) * composition);
        if (floor > worst) return floor;
    }
    real epsilon = terms->epsilon, sum = 0, mine_bar = presence[a], their_bar = presence[b];
    int scored = 0;
    for (int s = 0; s < samples; ++s) {
        const real* sample = moments + (long long)s * 4 * stride;
        real x_mean = sample[a], y_mean = sample[b];
        if (x_mean - epsilon <= mine_bar && y_mean - epsilon <= their_bar) continue;
        Moment x = {x_mean, sample[stride + a], sample[2 * stride + a], sample[3 * stride + a]};
        Moment y = {y_mean, sample[stride + b], sample[2 * stride + b], sample[3 * stride + b]};
        real shared = overlap(x, y, terms);
        sum += shared < epsilon ? epsilon : (shared > 1 - epsilon ? 1 - epsilon : shared);
        scored += 1;
    }
    real coverage = epsilon;
    if (scored > 0) {
        coverage = sum / scored;
        if (isnan(coverage)) coverage = 1;
    }
    real weight = terms->weight < 0 ? (real)scored / ((real)scored + 1) : terms->weight;
    real distance = weight * coverage + (1 - weight) * composition;
    return isnan(distance) ? 1.0f : (float)distance;
}

#define NEAREST_SLOTS 4

__device__ __forceinline__ bool ahead(float d, unsigned int at, float e, unsigned int bt) {
    return d < e || (d == e && at < bt);
}

// A warp per query holds its k nearest in registers, slot s in lane s % 32, so an insert is a few
// shuffles and not a walk through memory. Ties go to the lower index as the host lists do.
extern "C" __global__ void nearest(const float* dots, int rows, int cols, int first, const unsigned int* selves,
                                   const float* query_variance, const float* base_variance, int k,
                                   float* distances, unsigned int* indices, const unsigned int* query_rows,
                                   const real* moments, const real* presence, int stride, int samples,
                                   const real* terms) {
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
    unsigned int query = samples > 0 ? query_rows[row] : 0u;
    const float* line = dots + (long long)row * cols;
    for (int start = 0; start < cols; start += 32) {
        int col = start + lane;
        unsigned int at = (unsigned int)(first + col);
        bool live = col < cols && at != own;
        float distance = __int_as_float(0x7f800000);
        if (live) {
            distance = samples > 0
                ? mixed(line[col], mine, base_variance[at], query, at, worst, moments, presence, stride,
                        samples, (const Terms*)terms)
                : composition_only(line[col], mine, base_variance[at]);
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
