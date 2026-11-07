#define FULL 0xffffffffu
#ifndef COVERAGE
#define COVERAGE float
#endif
typedef COVERAGE real;

struct Terms {
    real weight, heaviest, epsilon, min_var_epsilon;
};

// `scale` is 1 / (deviation * sqrt 2), so a cdf costs a product and not a division.
struct Moment {
    real mean, variance, scale, log_variance;
};

__device__ __forceinline__ real normal_cdf(Moment m, real x) {
    return fmin((real)0.5 * erfc((m.mean - x) * m.scale), (real)1.0);
}

// The host measures one pair from either end, and both have to round alike.
__device__ real overlap(Moment a, Moment b, const Terms* __restrict__ terms) {
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
        real across = 1 / (a.variance - b.variance);
        k1 = (tmp - a.mean * b.variance + b.mean * a.variance) * across;
        k2 = -(tmp + a.mean * b.variance - b.mean * a.variance) * across;
    }
    if (k1 > k2) {
        real swap = k1;
        k1 = k2;
        k2 = swap;
    }
    Moment narrow = a.variance > b.variance ? b : a;
    Moment wide = a.variance > b.variance ? a : b;
    if (k1 == k2) return fabs(normal_cdf(narrow, k1) - normal_cdf(wide, k1));
    return fabs(normal_cdf(narrow, k2) - normal_cdf(narrow, k1)
                + normal_cdf(wide, k1) - normal_cdf(wide, k2));
}

__device__ __forceinline__ float composition_only(float dot, float mine, float theirs) {
    float total = mine + theirs;
    float distance = 0.0f;
    if (total != 0.0f) distance = (total - 2.0f * dot) / total;
    return isnan(distance) ? 2.0f : distance;
}

__device__ __forceinline__ bool absent(real x_mean, real y_mean, real mine_bar, real their_bar, real epsilon) {
    return x_mean - epsilon <= mine_bar && y_mean - epsilon <= their_bar;
}

// Each overlap only adds to the sum, so once the sum alone puts the pair past `worst` the rest
// cannot bring it back. The composition floor is the same bound before any sample is read.
__device__ float mixed(float dot, float mine, float theirs, unsigned int a, unsigned int b, float worst,
                       const real* __restrict__ moments, const real* __restrict__ presence, int stride,
                       int samples, const Terms* __restrict__ terms) {
    const float past = __int_as_float(0x7f800000);
    real total = (real)mine + (real)theirs;
    real composition = total == 0 ? (real)0 : (total - 2 * (real)dot) / total;
    if (isnan(composition)) composition = 2;
    if (composition > 0 && (float)((1 - terms->heaviest) * composition) > worst) return past;
    real epsilon = terms->epsilon, mine_bar = presence[a], their_bar = presence[b];
    int scored = 0;
    for (int s = 0; s < samples; ++s) {
        const real* __restrict__ sample = moments + (long long)s * 4 * stride;
        scored += !absent(sample[a], sample[b], mine_bar, their_bar, epsilon);
    }
    real weight = terms->weight < 0 ? (real)scored / ((real)scored + 1) : terms->weight;
    real sum = 0, enough = past;
    if (scored > 0 && weight > 0 && isfinite(worst)) {
        real slack = sizeof(real) == 8 ? (real)1e-12 : (real)1e-5;
        enough = ((real)worst - (1 - weight) * composition) / weight * scored * (1 + slack) + slack;
    }
    for (int s = 0; s < samples && scored > 0; ++s) {
        const real* __restrict__ sample = moments + (long long)s * 4 * stride;
        real x_mean = sample[a], y_mean = sample[b];
        if (absent(x_mean, y_mean, mine_bar, their_bar, epsilon)) continue;
        Moment x = {x_mean, sample[stride + a], sample[2 * stride + a], sample[3 * stride + a]};
        Moment y = {y_mean, sample[stride + b], sample[2 * stride + b], sample[3 * stride + b]};
        real shared = overlap(x, y, terms);
        sum += shared < epsilon ? epsilon : (shared > 1 - epsilon ? 1 - epsilon : shared);
        if (sum > enough) return past;
    }
    real coverage = epsilon;
    if (scored > 0) {
        coverage = sum / scored;
        if (isnan(coverage)) coverage = 1;
    }
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
                                   const real* __restrict__ moments, const real* __restrict__ presence,
                                   int stride, int samples, const real* __restrict__ terms) {
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
