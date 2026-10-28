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
